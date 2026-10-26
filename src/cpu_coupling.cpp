#include "cpu_coupling.h"
#include "swmm5.h"
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <vector>
#include <unordered_map>

extern "C" {
#include "headers.h"
}

namespace lis { namespace cpu_coupling {

void map_internal_outfalls(Fnames* filenames, Vconnect* coupling,
    Pars* pars, Arrays* arrays, int verbose)
{
    coupling->outfall_count = 0;
    coupling->outfall_dem_count = 0;
    coupling->outfall_node_Vindex = nullptr;
    coupling->outfall_dem_Vindex = nullptr;
    coupling->outfall_dem_Vindex_uni = nullptr;
    coupling->outfall_dynamic_stage = nullptr;
    coupling->dynamic_outfall_count = 0;
    coupling->outfall_slot = nullptr;

    int total = 0;
    for (int n=0; n<Nobjects[NODE]; ++n) if (Node[n].type == OUTFALL) ++total;
    if (total == 0) return;

    coupling->outfall_node_Vindex = new int[total]();
    coupling->outfall_dem_Vindex = new int[total]();
    coupling->outfall_dem_Vindex_uni = new int[total]();
    coupling->outfall_dynamic_stage = new unsigned char[total]();
    coupling->outfall_slot = new int[total]();

    FILE* log = fopen(filenames->flooderrorfile, "a");
    if (log) fprintf(log, "SWMM OUTFALL TO 2D MAPPING START\n");
    const double xmax=pars->blx+pars->xsz*pars->dx;
    const double ymax=pars->bly+pars->ysz*pars->dy;
    for (int n=0; n<Nobjects[NODE]; ++n)
    {
        if (Node[n].type != OUTFALL) continue;
        const double x=NodeXY[n][0], y=NodeXY[n][1];
        if (x < pars->blx-1.e-6 || x > xmax+1.e-6 ||
            y < pars->bly-1.e-6 || y > ymax+1.e-6) continue;
        int i=(int)((x-pars->blx)/pars->dx);
        int j=(int)((ymax-y)/pars->dy);
        if(i>=pars->xsz) i=pars->xsz-1;
        if(j>=pars->ysz) j=pars->ysz-1;
        if(i<0 || j<0) continue;
        const int cell=i+j*pars->xsz;
        if(FABS(arrays->DEM[cell]-pars->nodata_elevation)<1.e-6) continue;

        const int k=coupling->outfall_count++;
        coupling->outfall_node_Vindex[k]=n;
        coupling->outfall_dem_Vindex[k]=cell;
        int slot=-1;
        for(int u=0;u<coupling->outfall_dem_count;++u)
            if(coupling->outfall_dem_Vindex_uni[u]==cell){slot=u;break;}
        if(slot<0){slot=coupling->outfall_dem_count;coupling->outfall_dem_Vindex_uni[coupling->outfall_dem_count++]=cell;}
        coupling->outfall_slot[k]=slot;

        const int oi=Node[n].subIndex;
        if(Outfall[oi].type==COUPLED_OUTFALL)
        {
            if(RouteModel!=DW){fprintf(stderr,"COUPLED outfall requires DYNWAVE\n");exit(EXIT_FAILURE);}
            coupling->outfall_dynamic_stage[k]=1;
            ++coupling->dynamic_outfall_count;
        }
        if(log) fprintf(log,"outfall:%-20s cell=%d coupled_stage=%d\n",
            Node[n].ID,cell,(int)coupling->outfall_dynamic_stage[k]);
    }
    for(int n=0;n<Nobjects[NODE];++n)
    {
        if(Node[n].type!=OUTFALL || Outfall[Node[n].subIndex].type!=COUPLED_OUTFALL) continue;
        bool mapped=false;
        for(int k=0;k<coupling->outfall_count;++k)
            if(coupling->outfall_node_Vindex[k]==n){mapped=true;break;}
        if(!mapped){fprintf(stderr,"COUPLED outfall %s has no valid DEM cell\n",Node[n].ID);exit(EXIT_FAILURE);}
    }
    if(log){fprintf(log,"mapped=%d unique_cells=%d coupled_stages=%d\n",
        coupling->outfall_count,coupling->outfall_dem_count,
        coupling->dynamic_outfall_count);fclose(log);}
    if(verbose) printf("coupled outfalls: mapped=%d cells=%d dynamic=%d\n",
        coupling->outfall_count,coupling->outfall_dem_count,
        coupling->dynamic_outfall_count);
}

void free_outfalls(Vconnect* coupling)
{
    delete[] coupling->outfall_node_Vindex;
    delete[] coupling->outfall_dem_Vindex;
    delete[] coupling->outfall_dem_Vindex_uni;
    delete[] coupling->outfall_dynamic_stage;
    delete[] coupling->outfall_slot;
    coupling->outfall_slot=nullptr;
    coupling->outfall_node_Vindex=nullptr;
    coupling->outfall_dem_Vindex=nullptr;
    coupling->outfall_dem_Vindex_uni=nullptr;
    coupling->outfall_dynamic_stage=nullptr;
}

void update_outfall_stages(Vconnect* coupling, Pars* pars, Arrays* arrays)
{
    for(int k=0;k<coupling->outfall_count;++k)
    {
        if(!coupling->outfall_dynamic_stage[k]) continue;
        const int cell=coupling->outfall_dem_Vindex[k];
        const double bed=arrays->DEM[cell], depth=arrays->H[cell];
        if(!std::isfinite(bed)||!std::isfinite(depth)||depth < -1.e-5)
        {fprintf(stderr,"invalid receiving-cell stage at %s\n",Node[coupling->outfall_node_Vindex[k]].ID);exit(EXIT_FAILURE);}
        if(swmm_setOutfallReceivingCell(coupling->outfall_node_Vindex[k],cell,
            bed,depth>0?depth:0,pars->dA)!=0)
        {fprintf(stderr,"failed to set outfall stage at %s\n",Node[coupling->outfall_node_Vindex[k]].ID);exit(EXIT_FAILURE);}
    }
}

void collect_return_sources(Vconnect* coupling, BoundCs* boundaries,
    Pars* pars, Arrays* arrays, int surface_source_start,
    int outfall_source_start, double exchange_seconds,
    double application_seconds)
{
    if(exchange_seconds<=0 || application_seconds<=0) return;
    const double scale=exchange_seconds/application_seconds;
    for(int u=0;u<coupling->count_dem_index;++u)
        boundaries->PS_Val[surface_source_start+u]=C(0.0);
    for(int u=0;u<coupling->outfall_dem_count;++u)
        boundaries->PS_Val[outfall_source_start+u]=C(0.0);

    // Map receiving cells to source slots.
    std::unordered_map<int,int> slot_of_cell;
    slot_of_cell.reserve((size_t)coupling->count_dem_index*2+1);
    for(int u=0;u<coupling->count_dem_index;++u) slot_of_cell.emplace(coupling->dem_Vindex_uni[u],u);
    std::vector<int> slots(Nobjects[NODE],-1);
    for(int i=0;i<coupling->count_index;++i)
    {
        std::unordered_map<int,int>::const_iterator f=slot_of_cell.find(coupling->dem_Vindex[i]);
        slots[coupling->node_Vindex[i]]= f==slot_of_cell.end() ? -1 : f->second;
    }
    for(int n=0;n<Nobjects[NODE];++n)
    {
        const double volume=swmm_getNodeFloodExchangeVolume(n);
        if(volume<=0) continue;
        if(slots[n]<0){fprintf(stderr,"flooded node %s has no mapped 2D receiving cell\n",Node[n].ID);exit(EXIT_FAILURE);}
        boundaries->PS_Val[surface_source_start+slots[n]] +=
            (NUMERIC_TYPE)(volume/(exchange_seconds*pars->dy)*scale);
    }

    for(int k=0;k<coupling->outfall_count;++k)
    {
        const double volume=swmm_getOutfallExchangeVolume(coupling->outfall_node_Vindex[k]);
        if(volume<=0 && !coupling->outfall_dynamic_stage[k]) continue;
        // Use the precomputed receiving-cell slot.
        boundaries->PS_Val[outfall_source_start+coupling->outfall_slot[k]] +=
            (NUMERIC_TYPE)(volume/(exchange_seconds*pars->dy)*scale);
    }
    for(int u=0;u<coupling->outfall_dem_count;++u)
    {
        const int cell=coupling->outfall_dem_Vindex_uni[u];
        const double q=boundaries->PS_Val[outfall_source_start+u]*pars->dy;
        const double available=std::max(0.0,(double)arrays->H[cell])*pars->dA;
        const double roundoff=(sizeof(NUMERIC_TYPE)==sizeof(float)?8.*FLT_EPSILON:8.*DBL_EPSILON)*std::max(1.0,available);
        if(q<0 && -q*application_seconds>available+1.e-7+roundoff)
        {fprintf(stderr,"reverse outfall flow exceeds receiving-cell water at cell %d\n",cell);exit(EXIT_FAILURE);}
    }
}

void apply_coupled_sources(BoundCs* boundaries, Pars* pars, Solver* solver,
    Arrays* arrays, int source_start, int source_end)
{
    for(int p=source_start;p<source_end;++p)
    {
        if(boundaries->PS_Ident[p]!=QFIX4) continue;
        const int cell=boundaries->xpi[p]+boundaries->ypi[p]*pars->xsz;
        const NUMERIC_TYPE old_h=arrays->H[cell];
        NUMERIC_TYPE new_h=old_h+boundaries->PS_Val[p]*pars->dx*solver->Tstep/pars->dA;
        if(new_h<C(0.0)) new_h=C(0.0);
        arrays->H[cell]=new_h;
        const NUMERIC_TYPE actual_q=(new_h-old_h)*pars->dA/solver->Tstep;
        if(actual_q>C(0.0)) boundaries->Qpoint_pos+=actual_q;
        else boundaries->Qpoint_neg+=actual_q;
    }
}

} }
