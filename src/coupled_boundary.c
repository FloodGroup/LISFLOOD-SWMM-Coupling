/* Conservative finite-volume COUPLED receiving cells.
 * Interior pipe endpoint areas represent half lengths. At an outfall the
 * boundary is a face, not a storage node: the whole terminal pipe storage
 * belongs to its interior endpoint. The receiving cell is a finite volume.
 * a prescribed bed elevation cannot create water in an initially dry pipe.
 * Only components containing COUPLED outfalls use the volume-based continuity
 * update; conduit momentum is still solved by SWMM's dynamic-wave equations.
 */
#include "headers.h"
#include "coupled_boundary.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "hydraulic_restart_io.h"

extern TRoutingTotals FlowTotals, StepFlowTotals;
typedef struct {
    int first, marked, group, nextMember;
    double oldVolume, oldNet, oldBoundaryQ, oldEnd, exchange, standardExchange;
    double oldLat, oldOverflow, oldLoss;
    double evapRate, seepRate;
} CNode;
typedef struct { int link, next, side; } CPort;
typedef struct {
    int id, configured, first, last;
    double bed, area, water, oldWater, oldTotal, oldNet, head;
} CCell;
static CNode* CN;
static CPort* CP;
static CCell* CC;
static int NC, Active;
static double* InletPrevious;
static double* LinkEvapPrevious;
static double* LinkSeepPrevious;
static int SIInlets;
static double RemainingInterval;
/* Receiving-cell lookup: open addressing, 0 = empty, otherwise CC index + 1. */
static int* CellHash;
static int CellHashMask;

static int parent(int* p, int i)
{
    while (p[i] != i) { p[i] = p[p[i]]; i = p[i]; }
    return i;
}

int coupled_restart(void* stream, int reading)
{
    FILE* f=(FILE*)stream;
    int active=Active, n=NC, g;
    HR_VALUE(active);
    if(active!=Active) return 1;
    HR_VALUE(SIInlets); HR_VALUE(RemainingInterval);
    HR_ARRAY(InletPrevious,Nobjects[NODE]);
    if(!Active) return 0;
    HR_VALUE(n);
    if(n<0 || n>Nobjects[NODE]) return 1;
    NC=n;
    HR_ARRAY(CN,Nobjects[NODE]); HR_ARRAY(CC,NC);
    HR_ARRAY(LinkEvapPrevious,Nobjects[LINK]); HR_ARRAY(LinkSeepPrevious,Nobjects[LINK]);
    if(reading)
    {
        memset(CellHash,0,((size_t)CellHashMask+1)*sizeof(int));
        for(g=0;g<NC;++g)
        {
            unsigned int h=((unsigned int)CC[g].id*2654435761u)&(unsigned int)CellHashMask;
            while(CellHash[h]) h=(h+1)&(unsigned int)CellHashMask;
            CellHash[h]=g+1;
        }
    }
    return 0;
}

static double end_area(int link, double head, int side)
{
    TXsect* xs = &Link[link].xsect;
    int n = side ? Link[link].node2 : Link[link].node1;
    double offset = side ? Link[link].offset2 : Link[link].offset1;
    double y = MAX(0.0, head - Node[n].invertElev - offset);
    double a = xsect_getAofY(xs, MIN(y, xs->yFull));
    if (SurchargeMethod == SLOT && !xsect_isOpen(xs->type) && y > xs->yFull)
    {
        /* Integrate SWMM's slot width to obtain a monotone storage curve.
         * Above 1.78 diameters SWMM uses the constant 1% width. */
        double top = MIN(y, 1.78 * xs->yFull);
        double dy = (top - xs->yFull) / 16.0;
        int j;
        for (j = 0; j < 16; ++j)
        {
            double yn = (xs->yFull + (j + 0.5) * dy) / xs->yFull;
            a += dy * xs->wMax * 0.5423 * exp(-pow(yn, 2.4));
        }
        a += MAX(0.0, y - top) * 0.01 * xs->wMax;
    }
    return a;
}

static double own_volume(int n, double head)
{
    double d = MAX(0.0, head - Node[n].invertElev);
    double v;
    if (Node[n].type == OUTFALL) return 0.0;
    if (d <= Node[n].fullDepth) return node_getVolume(n, d);
    v = node_getVolume(n, Node[n].fullDepth);
    if (AllowPonding && Node[n].pondedArea > 0.0)
        return v + (d - Node[n].fullDepth) * Node[n].pondedArea;
    return v + (d - Node[n].fullDepth) * MAX(MinSurfArea,
        node_getSurfArea(n, Node[n].fullDepth));
}

static double end_volume(int n, double head)
{
    int p;
    double v = 0.0;
    for (p = CN[n].first; p >= 0; p = CP[p].next)
    {
        int l = CP[p].link;
        if (Link[l].type == CONDUIT && Link[l].xsect.type != DUMMY)
        {
            int other = CP[p].side ? Link[l].node1 : Link[l].node2;
            double weight = Node[n].type == OUTFALL ? 0.0 :
                (Node[other].type == OUTFALL ? 1.0 : 0.5);
            v += weight * link_getLength(l) * Conduit[Link[l].subIndex].barrels *
                end_area(l, head, CP[p].side);
        }
    }
    return v;
}

static double node_volume(int n, double head)
{
    return own_volume(n, head) + end_volume(n, head);
}

static int cell_group(int cell)
{
    unsigned int h = ((unsigned int)cell * 2654435761u) & (unsigned int)CellHashMask;
    while (CellHash[h])
    {
        if (CC[CellHash[h]-1].id == cell) return CellHash[h]-1;
        h = (h + 1) & (unsigned int)CellHashMask;
    }
    if (NC >= Nobjects[NODE]) return -1;
    CC[NC].id = cell; CC[NC].first = CC[NC].last = -1;
    CellHash[h] = NC + 1;
    return NC++;
}

static double cell_volume(int g, double head)
{
    return CC[g].area * MAX(0.0, head - CC[g].bed);
}

static double inverse_volume(int n, int g, double v, double lo, double hi)
{
    /* Bisect stalled Newton steps and retain the iterate with least residual. */
    int i;
    double x, best_x, best_r = HUGE_VAL, prev_r = HUGE_VAL;
    if (v <= 0.0) return lo;
    for (i = 0; i < 40; ++i)
    {
        double w = g < 0 ? node_volume(n, hi) : cell_volume(g, hi);
        if (w >= v) break;
        hi = lo + 2.0 * (hi - lo);
    }
    x = g < 0 ? MAX(lo,MIN(hi,Node[n].invertElev+Node[n].newDepth)) : 0.5*(lo+hi);
    best_x = x;
    for (i = 0; i < 200; ++i)
    {
        double w = g < 0 ? node_volume(n,x) : cell_volume(g,x);
        double r = fabs(w-v), eps=MAX(1.e-6,fabs(x)*1.e-7), slope, next;
        if (r < best_r) { best_r = r; best_x = x; }
        if (r <= 1.e-11+1.e-14*v) return x;
        if (w > v) hi=x; else lo=x;
        if (hi - lo <= 4.0*DBL_EPSILON*MAX(1.0,fabs(x))) break;
        slope=((g < 0 ? node_volume(n,x+eps) : cell_volume(g,x+eps))-w)/eps;
        next = slope > 0.0 ? x-(w-v)/slope : 0.5*(lo+hi);
        if (!(next > lo && next < hi) || r > 0.5*prev_r) next = 0.5*(lo+hi);
        prev_r = r;
        x = next;
    }
    return best_x;
}

int coupled_open(void)
{
    int i, l, k, *parents, *marks;
    coupled_close();
    for (i = 0; i < Nobjects[NODE]; ++i)
        if (Node[i].type == OUTFALL && Outfall[Node[i].subIndex].type == COUPLED_OUTFALL)
            Active = 1;
    InletPrevious = (double*)calloc(Nobjects[NODE], sizeof(double));
    if (Nobjects[NODE] && !InletPrevious) { ErrorCode=ERR_MEMORY; return ErrorCode; }
    if (!Active) return 0;
    if (RouteModel != DW || IgnoreRouting || NumEvents > 0 ||
        (Nobjects[POLLUT] > 0 && !IgnoreQuality))
    {
        fprintf(stderr, "COUPLED requires continuous DYNWAVE and hydraulic-only routing (no routing EVENTS).\n");
        ErrorCode = ERR_SYSTEM; return ErrorCode;
    }
    CN = (CNode*)calloc(Nobjects[NODE], sizeof(CNode));
    CP = (CPort*)calloc(MAX(1, 2*Nobjects[LINK]), sizeof(CPort));
    CC = (CCell*)calloc(Nobjects[NODE], sizeof(CCell));
    for (CellHashMask = 1; CellHashMask < 2*Nobjects[NODE]+1; CellHashMask <<= 1);
    CellHash = (int*)calloc(CellHashMask, sizeof(int));
    CellHashMask -= 1;
    LinkEvapPrevious=(double*)calloc(MAX(1,Nobjects[LINK]),sizeof(double));
    LinkSeepPrevious=(double*)calloc(MAX(1,Nobjects[LINK]),sizeof(double));
    parents = (int*)malloc(Nobjects[NODE]*sizeof(int));
    marks = (int*)calloc(Nobjects[NODE], sizeof(int));
    if (!CN || !CP || !CC || !CellHash || !InletPrevious || !parents || !marks || !LinkEvapPrevious || !LinkSeepPrevious)
    { FREE(parents); FREE(marks); coupled_close(); ErrorCode = ERR_MEMORY; return ErrorCode; }
    for (i = 0; i < Nobjects[NODE]; ++i)
    { parents[i] = i; CN[i].first = -1; CN[i].group = -1; }
    for (l = 0; l < Nobjects[LINK]; ++l)
    {
        int n1 = Link[l].node1, n2 = Link[l].node2;
        parents[parent(parents,n1)] = parent(parents,n2);
        CP[2*l].link = CP[2*l+1].link = l;
        CP[2*l].side = 0; CP[2*l+1].side = 1;
        CP[2*l].next = CN[n1].first; CN[n1].first = 2*l;
        CP[2*l+1].next = CN[n2].first; CN[n2].first = 2*l+1;
    }
    for (i = 0; i < Nobjects[NODE]; ++i)
        if (Node[i].type == OUTFALL && Outfall[Node[i].subIndex].type == COUPLED_OUTFALL)
            marks[parent(parents,i)] = 1;
    for (i = 0; i < Nobjects[NODE]; ++i)
    {
        CN[i].marked = marks[parent(parents,i)];
        if (!CN[i].marked) continue;
        /* Account explicitly for the minimum nodal area already used by
         * SWMM's depth equation, instead of unreported numerical storage. */
        if (Node[i].type == JUNCTION && Node[i].fullVolume == 0.0)
            Node[i].fullVolume = MinSurfArea * Node[i].fullDepth;
        Node[i].newVolume = Node[i].oldVolume = own_volume(i,
            Node[i].invertElev + Node[i].newDepth);
    }
    for (l = 0; l < Nobjects[LINK]; ++l)
    {
        i = Link[l].node1; k = Link[l].node2;
        if (!CN[i].marked || Link[l].type != CONDUIT || Link[l].xsect.type == DUMMY) continue;
        Link[l].newVolume = Link[l].oldVolume = link_getLength(l) *
            Conduit[Link[l].subIndex].barrels *
            ((Node[i].type == OUTFALL ? 0.0 : (Node[k].type == OUTFALL ? 1.0 : 0.5))*end_area(l, Node[i].invertElev+Node[i].newDepth,0) +
             (Node[k].type == OUTFALL ? 0.0 : (Node[i].type == OUTFALL ? 1.0 : 0.5))*end_area(l, Node[k].invertElev+Node[k].newDepth,1));
        /* Preserve midpoint flow-area history used by dwflow's inertial term. */
    }
    FREE(parents); FREE(marks);
    return 0;
}

void coupled_close(void)
{ FREE(CN); FREE(CP); FREE(CC); FREE(CellHash); FREE(InletPrevious); FREE(LinkEvapPrevious); FREE(LinkSeepPrevious);
  NC = 0; CellHashMask = 0; Active = 0; SIInlets = 0; RemainingInterval=0.0; }
void coupled_si_inlets(int enabled) { SIInlets=enabled; }
void coupled_remaining_interval(double seconds) { RemainingInterval=seconds; }
double coupled_inlet_factor(void) { return SIInlets ? 0.028316846592 : UCF(FLOW); }
int coupled_ready(void)
{
    int n;
    if (!Active) return 1;
    for (n=0;n<Nobjects[NODE];++n)
        if (coupled_outfall(n) && (CN[n].group < 0 || !CC[CN[n].group].configured)) return 0;
    return 1;
}
double coupled_available(int n)
{
    if (n < 0 || n >= Nobjects[NODE]) return 0.0;
    if (coupled_node(n)) return node_volume(n,Node[n].invertElev+Node[n].newDepth)*0.028316846592;
    return MAX(0.0,Node[n].newVolume)*0.028316846592;
}
double coupled_old_volume(int n) { return coupled_node(n) ? CN[n].oldVolume : Node[n].oldVolume; }
int coupled_active(void) { return Active; }
int coupled_node(int n) { return Active && CN && CN[n].marked; }
int coupled_outfall(int n)
{ return coupled_node(n) && Node[n].type == OUTFALL && Outfall[Node[n].subIndex].type == COUPLED_OUTFALL; }

int coupled_set_cell(int n, int cell, double bed, double depth, double area)
{
    int g, i;
    if (!coupled_outfall(n) || cell < 0 || !isfinite(bed) || !isfinite(depth) ||
        !isfinite(area) || depth < -1.e-9 || area <= 0.0) return 1;
    g = cell_group(cell);
    if (g < 0) return 1;
    if (CN[n].group >= 0 && CN[n].group != g) return 1;
    if (CN[n].group < 0)
    {
        CN[n].group=g; CN[n].nextMember=-1;
        if (CC[g].last < 0) CC[g].first=n;
        else CN[CC[g].last].nextMember=n;
        CC[g].last=n;
    }
    if (CC[g].configured && (fabs(CC[g].bed - bed/0.3048) > 1.e-7 ||
        fabs(CC[g].area - area/(0.3048*0.3048)) > 1.e-7 ||
        fabs(CC[g].water - MAX(0.0,depth)*area/0.028316846592) >
            1.e-7+1.e-12*MAX(1.0,CC[g].water))) return 1;
    CC[g].bed = bed/0.3048; CC[g].area = area/(0.3048*0.3048);
    CC[g].water = MAX(0.0,depth)*area/0.028316846592;
    CC[g].configured = 1;
    /* Bed is the dry receiving-face elevation, with ZERO supply. It never
     * contributes a fictitious downstream half-pipe storage volume. */
    for (i = CC[g].first; i >= 0; i=CN[i].nextMember)
        {
            Outfall[Node[i].subIndex].externalStageActive = 1;
            Outfall[Node[i].subIndex].externalStage = CC[g].bed+CC[g].water/CC[g].area;
            Node[i].newDepth = MAX(0.0,Outfall[Node[i].subIndex].externalStage-Node[i].invertElev);
        }
    return 0;
}

void coupled_begin_interval(void)
{
    int n, g;
    if (!Active) return;
    for (n = 0; n < Nobjects[NODE]; ++n)
        CN[n].exchange = CN[n].standardExchange = 0.0;
    for (g = 0; g < NC; ++g) CC[g].configured = 0;
}

void coupled_sync_inlets(void)
{
    int n;
    if (!InletPrevious || ADDLINK != 1 || !Node_Inflow) return;
    for (n = 0; n < Nobjects[NODE]; ++n)
    {
        double now = Node_Inflow[n]/coupled_inlet_factor(), delta = now - InletPrevious[n];
        /* The interface supplies an interval-average rate, constant on
         * both ends of the SWMM trapezoid, including the first substep. */
        Node[n].newLatFlow += delta;
        Node[n].inflow += MAX(now,0.0)-MAX(InletPrevious[n],0.0);
        Node[n].outflow += MAX(-now,0.0)-MAX(-InletPrevious[n],0.0);
        StepFlowTotals.exInflow += MAX(now,0.0)-MAX(InletPrevious[n],0.0);
        StepFlowTotals.outflow += MAX(-now,0.0)-MAX(-InletPrevious[n],0.0);
        InletPrevious[n] = now;
    }
}

void coupled_begin_step(void)
{
    int n, g;
    if (!Active) return;
    for (n = 0; n < Nobjects[NODE]; ++n)
    {
        if (!CN[n].marked) continue;
        CN[n].oldVolume = node_volume(n, Node[n].invertElev+Node[n].newDepth);
        CN[n].oldEnd = end_volume(n, Node[n].invertElev+Node[n].newDepth);
        CN[n].oldLat = Node[n].newLatFlow;
        CN[n].oldOverflow = Node[n].overflow;
        CN[n].oldLoss = Node[n].losses;
        /* Preserve signed pipe flow independently of the reporting mirror
         * installed by node_getSystemOutflow for reverse-flow outfalls. */
        if (Node[n].type != OUTFALL) CN[n].oldNet = Node[n].inflow-Node[n].outflow;
    }
    for (g = 0; g < NC; ++g)
    {
        CC[g].oldWater = CC[g].water;
        CC[g].oldTotal = CC[g].water; CC[g].oldNet = 0.0;
    }
}

double coupled_solve_node(int n, double dt)
{
    double target, lo, hi, head, cap;
    if (!coupled_node(n) || Node[n].type == OUTFALL) return Node[n].newDepth;
    target = CN[n].oldVolume + dt*(Node[n].inflow-Node[n].outflow);
    lo = Node[n].invertElev; hi = lo + MAX(1.0,Node[n].fullDepth+Node[n].surDepth);
    head = inverse_volume(n,-1,MAX(0.0,target),lo,hi);
    Node[n].overflow = 0.0;
    if (!(AllowPonding && Node[n].pondedArea > 0.0))
    {
        cap = lo + Node[n].fullDepth + Node[n].surDepth;
        if (head > cap)
        {
            head = cap;
            Node[n].overflow = MAX(0.0,(target-node_volume(n,head))/dt);
        }
    }
    Node[n].newVolume = own_volume(n,head);
    return head-lo;
}

double coupled_solve_cells(double dt)
{
    int n, g;
    double maxChange = 0.0;
    if (!Active) return 0.0;
    for (g = 0; g < NC; ++g)
    {
        double net = 0.0, lo = CC[g].bed, hi, target, head;
        for (n = CC[g].first; n >= 0; n=CN[n].nextMember)
            { lo = MIN(lo,Node[n].invertElev); net += Node[n].inflow-Node[n].outflow; }
        target = CC[g].oldTotal + dt*net;
        hi = MAX(CC[g].bed+MAX(1.0,CC[g].oldWater/CC[g].area), lo+1.0);
        head = CC[g].bed+MAX(0.0,target)/CC[g].area;
        CC[g].head = head;
        CC[g].water = CC[g].area*MAX(0.0,head-CC[g].bed);
        for (n = CC[g].first; n >= 0; n=CN[n].nextMember)
            {
                double y = MAX(0.0,head-Node[n].invertElev);
                maxChange = MAX(maxChange,fabs(y-Node[n].newDepth));
                Node[n].newDepth = y;
                Outfall[Node[n].subIndex].externalStage = head;
            }
    }
    return maxChange;
}

/* Reduce donor fluxes before continuity is solved. Seepage/evaporation are
 * also withdrawals, so a dry donor cannot pay them from an empty volume. */
static void scale_donor(int n, double factor)
{
    int p;
    for (p=CN[n].first;p>=0;p=CP[p].next)
    {
        int l=CP[p].link, n1=Link[l].node1, n2=Link[l].node2;
        double old=Link[l].newFlow, q=old*factor;
        int donor=old>=0.0?n1:n2;
        if (donor != n) continue;
        Node[n1].outflow += MAX(q,0.0)-MAX(old,0.0);
        Node[n2].inflow += MAX(q,0.0)-MAX(old,0.0);
        Node[n1].inflow += MAX(-q,0.0)-MAX(-old,0.0);
        Node[n2].outflow += MAX(-q,0.0)-MAX(-old,0.0);
        Link[l].newFlow=q;
        if (Link[l].type == CONDUIT)
        {
            int k=Link[l].subIndex;
            double loss=(Conduit[k].evapLossRate+Conduit[k].seepLossRate)*Conduit[k].barrels;
            Node[n].outflow -= loss*(1.0-factor);
            Conduit[k].evapLossRate *= factor;
            Conduit[k].seepLossRate *= factor;
            Conduit[k].q1=Conduit[k].q2=q/Conduit[k].barrels;
        }
    }
    Node[n].outflow -= Node[n].losses*(1.0-factor);
    Node[n].losses *= factor;
    if (Node[n].type == STORAGE)
    {
        Storage[Node[n].subIndex].evapLoss *= factor;
        Storage[Node[n].subIndex].exfilLoss *= factor;
    }
}

void coupled_limit_flows(double dt)
{
    int g, n, pass;
    if (!Active) return;
    for (pass=0;pass<8;++pass)
    {
    int changed=0;
    for (g = 0; g < NC; ++g)
    {
        double incoming = 0.0, withdrawing = 0.0, factor;
        for (n = CC[g].first; n >= 0; n=CN[n].nextMember)
        { incoming += Node[n].inflow; withdrawing += Node[n].outflow; }
        if (withdrawing <= incoming+CC[g].oldWater/dt+1.e-12) continue;
        changed=1;
        factor = (incoming+CC[g].oldWater/dt)/withdrawing;
        for (n = CC[g].first; n >= 0; n=CN[n].nextMember) scale_donor(n,factor);
    }
    /* Several pipes can withdraw from one junction. Limit their SUM, not
     * each pipe independently. Reserve water already committed to the 2D
     * manhole over the remaining exchange interval. */
    for (n=0;n<Nobjects[NODE];++n)
    {
        double total, available, reserve=0.0, factor, fixed;
        if (!CN[n].marked || Node[n].type == OUTFALL) continue;
        fixed=MAX(0.0,-Node[n].newLatFlow);
        total=MAX(0.0,Node[n].outflow-fixed);
        if (Node_Inflow && ADDLINK == 1)
            reserve=MAX(0.0,-Node_Inflow[n])/coupled_inlet_factor()*MAX(0.0,RemainingInterval-dt);
        available=MAX(0.0,(CN[n].oldVolume-reserve)/dt+Node[n].inflow-fixed);
        if (total <= available+1.e-12) continue;
        changed=1; factor=available/total;
        scale_donor(n,factor);
    }
    if (!changed) break;
    }
}

int coupled_commit(double dt)
{
    int n, l, g, bad;
    double tolerance = 1.e-7; /* ft3 absolute + relative below */
    if (!Active) return 0;
    for (g = 0; g < NC; ++g)
    {
        double net = 0.0, actual, target;
        for (n = CC[g].first; n >= 0; n=CN[n].nextMember) net += Node[n].inflow-Node[n].outflow;
        target = CC[g].oldTotal+dt*net;
        actual = cell_volume(g,CC[g].head);
        if (!isfinite(actual) || !isfinite(target) || fabs(actual-target) > tolerance + 1.e-8*MAX(1.0,fabs(target)))
        {
            fprintf(stderr,"finite receiving cell %d did not converge (volume residual %.9g m3). Reduce ROUTING_STEP / coupleint.\n",CC[g].id,(actual-target)*0.028316846592);
            ErrorCode = ERR_TIMESTEP; return ErrorCode;
        }
    }
    /* Each thread owns its node/link state; report errors in node order. */
    bad = Nobjects[NODE];
#pragma omp parallel for schedule(dynamic, 64) num_threads(NumThreads)
    for (n = 0; n < Nobjects[NODE]; ++n)
    {
        if (!CN[n].marked) continue;
        if (Node[n].type != OUTFALL)
        {
            double target=CN[n].oldVolume+dt*(Node[n].inflow-Node[n].outflow-Node[n].overflow);
            double actual=node_volume(n,Node[n].invertElev+Node[n].newDepth);
            if (!isfinite(actual) || !isfinite(target) || fabs(actual-target)>1.e-6+1.e-8*MAX(1.0,fabs(target)))
            {
#pragma omp critical (coupled_commit_bad)
                { if (n < bad) bad = n; }
            }
        }
        if (Node[n].type == OUTFALL)
        {
            double q = Node[n].inflow-Node[n].outflow;
            double integral = q*dt;
            double end = end_volume(n,Node[n].invertElev+Node[n].newDepth);
            double dv = integral-(end-CN[n].oldEnd);
            CN[n].exchange += dv;
            CN[n].standardExchange += dv;
            CN[n].oldBoundaryQ = CN[n].oldNet;
            CN[n].oldNet = q;
        }
    }
    if (bad < Nobjects[NODE])
    {
        double target=CN[bad].oldVolume+dt*(Node[bad].inflow-Node[bad].outflow-Node[bad].overflow);
        double actual=node_volume(bad,Node[bad].invertElev+Node[bad].newDepth);
        fprintf(stderr,"node %s volume equation failed (%.9g m3); reduce ROUTING_STEP.\n",Node[bad].ID,(actual-target)*0.028316846592);
        ErrorCode=ERR_TIMESTEP; return ErrorCode;
    }
#pragma omp parallel for schedule(dynamic, 64) num_threads(NumThreads)
    for (l = 0; l < Nobjects[LINK]; ++l)
    {
        int n1 = Link[l].node1, n2 = Link[l].node2;
        if (!CN[n1].marked || Link[l].type != CONDUIT || Link[l].xsect.type == DUMMY) continue;
        Link[l].newVolume = link_getLength(l)*Conduit[Link[l].subIndex].barrels*
            ((Node[n1].type == OUTFALL ? 0.0 : (Node[n2].type == OUTFALL ? 1.0 : 0.5))*end_area(l,Node[n1].invertElev+Node[n1].newDepth,0)+
             (Node[n2].type == OUTFALL ? 0.0 : (Node[n1].type == OUTFALL ? 1.0 : 0.5))*end_area(l,Node[n2].invertElev+Node[n2].newDepth,1));
        /* Keep the midpoint flow area maintained by dwflow. */
    }
    return 0;
}

void coupled_account(double dt)
{
    int n, l;
    if (!Active) return;
    for (n = 0; n < Nobjects[NODE]; ++n)
    {
        double oldq, q, rawOut, rawIn, end, dv;
        if (!CN[n].marked) continue;
        /* Volume continuity in coupled components uses backward Euler.
         * Replace only their integrated lateral/overflow terms; unrelated
         * components retain SWMM's native trapezoidal continuity. */
        FlowTotals.exInflow += 0.5*dt*(MAX(Node[n].newLatFlow,0.0)-MAX(CN[n].oldLat,0.0));
        FlowTotals.outflow += 0.5*dt*(MAX(-Node[n].newLatFlow,0.0)-MAX(-CN[n].oldLat,0.0));
        if (Node[n].type != OUTFALL)
        {
            FlowTotals.flooding += 0.5*dt*(Node[n].overflow-CN[n].oldOverflow);
            if (Node[n].type == STORAGE)
            {
                double evap=Storage[Node[n].subIndex].evapLoss/dt;
                double seep=Storage[Node[n].subIndex].exfilLoss/dt;
                FlowTotals.evapLoss += 0.5*dt*(evap-CN[n].evapRate);
                FlowTotals.seepLoss += 0.5*dt*(seep-CN[n].seepRate);
                CN[n].evapRate=evap; CN[n].seepRate=seep;
            }
            continue;
        }
        oldq = CN[n].oldBoundaryQ; q = CN[n].oldNet;
        rawOut = 0.5*dt*(MAX(oldq,0.0)+MAX(q,0.0));
        rawIn = 0.5*dt*(MAX(-oldq,0.0)+MAX(-q,0.0));
        end = end_volume(n,Node[n].invertElev+Node[n].newDepth);
        dv = q*dt;
        /* Replace infinite-reservoir pipe-end flux by the boundary CV flux.
         * This is a storage term, NOT a correction using system residuals. */
        FlowTotals.outflow += MAX(dv,0.0)-rawOut;
        FlowTotals.exInflow += MAX(-dv,0.0)-rawIn;
    }
    for (l=0;l<Nobjects[LINK];++l)
    {
        int k;
        double evap, seep;
        if (!coupled_node(Link[l].node1) || Link[l].type != CONDUIT) continue;
        k=Link[l].subIndex;
        evap=Conduit[k].evapLossRate*Conduit[k].barrels;
        seep=Conduit[k].seepLossRate*Conduit[k].barrels;
        FlowTotals.evapLoss += 0.5*dt*(evap-LinkEvapPrevious[l]);
        FlowTotals.seepLoss += 0.5*dt*(seep-LinkSeepPrevious[l]);
        LinkEvapPrevious[l]=evap; LinkSeepPrevious[l]=seep;
    }
}

double coupled_exchange(int n) { return coupled_outfall(n) ? CN[n].exchange*0.028316846592 : 0.0; }
double coupled_standard_exchange(int n) { return coupled_node(n) ? CN[n].standardExchange*0.028316846592 : 0.0; }
