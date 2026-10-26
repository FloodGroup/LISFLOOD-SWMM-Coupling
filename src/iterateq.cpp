/*
*****************************************************************************
ITERATEQ and UPDATEH
---------------------


*****************************************************************************
*/

#include "lisflood.h"
#include "rain/rain.h"
#include "sgc.h"
#include "geometry.h"
#include "cpu_coupling.h"
#include <float.h>
#include <chrono>
extern "C" {
#include "headers.h"
};
#include "swmm_couple.h"
#include "coupling_restart.h"
#include "swmm5.h"


namespace flood {
	struct verticalcon
	{
		static void load_floodinfo(Fnames* Fnameptr, Vconnect* verlink, Pars* Parptr, Arrays* Arrptr, int verbose);
		static void vcon_init(Vconnect* verlink, Pars* Parptr);
		static void read_vconinfo(Vconnect* verlink, FILE* flood_file, int verbose);
		static void cal_dem_index(Vconnect* verlink, FILE* flood_error, Pars* Parptr, Arrays* Arrptr);
		static void demcountuni_cal(Vconnect* verlink);
		static void check_nodedepth(Vconnect* verlink, FILE* flood_error, Arrays* Arrptr);
		static void cal_exchange_flow(Vconnect* verlink, Arrays* Arrptr, Pars* Parptr,
			Solver* Solverptr, const double* pending_volumes);
		static NUMERIC_TYPE orificeflow_cal(NUMERIC_TYPE c0, NUMERIC_TYPE area, NUMERIC_TYPE H1D, NUMERIC_TYPE H2D);
		static NUMERIC_TYPE weirflow_cal(NUMERIC_TYPE cw, NUMERIC_TYPE w, NUMERIC_TYPE h2d);
		static void free_vcon(Vconnect* verlink);
		static void nodeinflow_cal(Vconnect* verlink);
	};
}

void flood::verticalcon::load_floodinfo
(
	Fnames* Fnameptr,
	Vconnect* verlink,
	Pars* Parptr,
	Arrays* Arrptr,
	int verbose
)
{
	ADDLINK = 1;
	FILE* flood_file = fopen_or_die(Fnameptr->floodfilename, "rb", "Loading FloodConfig\n", verbose);
	FILE* flood_error = fopen_or_die(Fnameptr->flooderrorfile, "w", "Loading FloodError\n", verbose);
	fscanf(flood_file, "%*s %d", &(verlink->count_index));
	vcon_init(verlink, Parptr);
	fscanf(flood_file, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
	read_vconinfo(verlink, flood_file, verbose);
	fclose(flood_file);
	cal_dem_index(verlink, flood_error, Parptr, Arrptr);
	check_nodedepth(verlink, flood_error, Arrptr);
	demcountuni_cal(verlink);
	fclose(flood_error);

	return;
}

void flood::verticalcon::vcon_init
(
	Vconnect* verlink,
	Pars* Parptr
)
{

	verlink->node_Vindex = new int[verlink->count_index]();
	verlink->dem_Vindex = new int[verlink->count_index]();
	verlink->dem_count = new int[Parptr->xsz * Parptr->ysz]();
	verlink->node_surfelev = new double[verlink->count_index]();
	verlink->node_c0 = new double[verlink->count_index]();
	verlink->node_area = new double[verlink->count_index]();
	verlink->node_qoutmax = new double[verlink->count_index]();
	verlink->node_cw = new double[verlink->count_index]();
	verlink->node_w = new double[verlink->count_index]();
	verlink->node_qinmax = new double[verlink->count_index]();
	verlink->node_c1 = new double[verlink->count_index]();
	verlink->node_area1 = new double[verlink->count_index]();

	verlink->count_dem_index = verlink->count_index;
	verlink->dem_Vindex_uni = new int[verlink->count_index]();


	verlink->dem_xi = new int[verlink->count_index]();
	verlink->dem_yi = new int[verlink->count_index]();
	//----------------------------------------------------------------------

	verlink->flow_1D_node = new double[verlink->count_index]();
	verlink->flow_2D_node = new double[verlink->count_index]();
	verlink->flow_2D_dem = new double[Parptr->xsz * Parptr->ysz]();
	verlink->flow_neg_count = new int[Parptr->xsz * Parptr->ysz]();
	//----------------------------------------------------------------------
	return;
}

void flood::verticalcon::read_vconinfo(Vconnect* verlink, FILE* flood_file, int verbose)
{
	for (int i = 0; i < verlink->count_index; ++i) {
		fscanf(flood_file,
			"%*s %d %lf %lf %lf %lf %lf %lf %lf %lf",
			&(verlink->node_Vindex[i]),
			&(verlink->node_c0[i]),
			&(verlink->node_area[i]),
			&(verlink->node_qoutmax[i]),
			&(verlink->node_cw[i]),
			&(verlink->node_w[i]),
			&(verlink->node_qinmax[i]),
			&(verlink->node_c1[i]),
			&(verlink->node_area1[i]));

	}
	if (verbose)
	{
		printf("read flood config done!\n");
	}
	return;
}

void flood::verticalcon::cal_dem_index
(
	Vconnect* verlink,
	FILE* flood_error,
	Pars* Parptr,
	Arrays* Arrptr
)
{
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "THE CORRESPONDING INDEX CALCULATION START:\n");
	fprintf(flood_error, "====================================================================\n");

	int x_index, y_index, dem_index;
	double x_min = Parptr->blx;
	double x_max = Parptr->blx + Parptr->xsz * Parptr->dx;
	double y_min = Parptr->bly;
	double y_max = Parptr->bly + Parptr->ysz * Parptr->dy;

	for (int j = 0; j < verlink->count_index; ++j)
	{
		verlink->node_surfelev[j] = (Node[verlink->node_Vindex[j]].invertElev + Node[verlink->node_Vindex[j]].fullDepth) * UCF(LENGTH);

		if ((x_min - NodeXY[verlink->node_Vindex[j]][0]) > 1e-6 || (NodeXY[verlink->node_Vindex[j]][0] - x_max) > 1e-6 ||
			(y_min - NodeXY[verlink->node_Vindex[j]][1]) > 1e-6 || (NodeXY[verlink->node_Vindex[j]][1] - y_max) > 1e-6)
		{
			verlink->dem_Vindex[j] = -9999;
			verlink->dem_xi[j] = -9999;
			verlink->dem_yi[j] = -9999;
			fprintf(flood_error, "node: %-15s is out of DEM range\n", Node[verlink->node_Vindex[j]].ID);
		}
		else if (Node[verlink->node_Vindex[j]].type == OUTFALL  || Node[verlink->node_Vindex[j]].type == DIVIDER || Node[verlink->node_Vindex[j]].type == STORAGE)
		{
			verlink->dem_Vindex[j] = -9999;
			verlink->dem_xi[j] = -9999;
			verlink->dem_yi[j] = -9999;
			fprintf(flood_error, "node: %-15s is the wrong node type\n", Node[verlink->node_Vindex[j]].ID);
		}
		else
		{
			x_index = (int)((NodeXY[verlink->node_Vindex[j]][0] - x_min) / Parptr->dx);
			y_index = (int)((y_max - NodeXY[verlink->node_Vindex[j]][1]) / Parptr->dy);

			if (x_index >= Parptr->xsz) x_index -= 1;
			if (y_index >= Parptr->ysz) y_index -= 1;

			dem_index = x_index + y_index * Parptr->xsz;

			if (FABS(Arrptr->DEM[dem_index] - Parptr->nodata_elevation) < 1e-6)
			{
				verlink->dem_Vindex[j] = -9999;
				verlink->dem_xi[j] = -9999;
				verlink->dem_yi[j] = -9999;
				fprintf(flood_error, "node: %-15s is out of DEM range\n", Node[verlink->node_Vindex[j]].ID);
			}
			else
			{
				verlink->dem_Vindex[j] = dem_index;
				verlink->dem_xi[j] = x_index;
				verlink->dem_yi[j] = y_index;
			}
		}
	}
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "THE CORRESPONDING INDEX CALCULATION END:\n");
	fprintf(flood_error, "====================================================================\n");
	return;
}


void flood::verticalcon::demcountuni_cal
(
	Vconnect* verlink
)
{
	int count_init = verlink->count_index;
	int k_temp = 0;
	memcpy(verlink->dem_Vindex_uni, verlink->dem_Vindex, verlink->count_index * sizeof(int));

	for (size_t m = 0; m < verlink->count_index; m++)
	{
		if (verlink->dem_Vindex_uni[m] == -9999)
		{

			for (size_t n = m; n < verlink->count_index - 1; n++)
			{
				verlink->dem_Vindex_uni[n] = verlink->dem_Vindex_uni[n + 1];

				verlink->dem_xi[n] = verlink->dem_xi[n + 1];
				verlink->dem_yi[n] = verlink->dem_yi[n + 1];
			}
			verlink->dem_Vindex_uni[verlink->count_index - 1] = 0;
			verlink->dem_xi[verlink->count_index - 1] = 0;
			verlink->dem_yi[verlink->count_index - 1] = 0;
			m--;
			count_init--;
		}
	}


	for (int k = 0; k < verlink->count_index; ++k)
	{
		if (verlink->dem_Vindex[k] != -9999)
			verlink->dem_count[verlink->dem_Vindex[k]] += 1;
	}


	for (size_t i = 0; i < count_init; i++)
	{
		for (size_t j = i + 1; j < count_init; j++)
		{
			if (verlink->dem_Vindex_uni[i] == verlink->dem_Vindex_uni[j])
			{
				for (size_t temp = j; temp < count_init - 1; temp++)
				{
					verlink->dem_Vindex_uni[temp] = verlink->dem_Vindex_uni[temp + 1];
					verlink->dem_xi[temp] = verlink->dem_xi[temp + 1];
					verlink->dem_yi[temp] = verlink->dem_yi[temp + 1];
				}
				j--;
				count_init--;
			}
		}
	}
	verlink->count_dem_index = count_init;
	return;
}


void flood::verticalcon::check_nodedepth
(
	Vconnect* verlink,
	FILE* flood_error,
	Arrays* Arrptr
)
{
	int count = 0;
	int count1 = 0;
	fprintf(flood_error, "\n");
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "FIX NODE MAX DEPTH START:\n");
	fprintf(flood_error, "====================================================================\n");

	for (size_t i = 0; i < verlink->count_index; i++)
	{
		if (verlink->dem_Vindex[i] != -9999)
		{
			++count1;
			if (FABS(Arrptr->DEM[verlink->dem_Vindex[i]] - verlink->node_surfelev[i]) > 0.1)
			{
				fprintf(flood_error, "The difference between node: %-15s invertElev:%lf + fullDepth:%lf = %lf and DEM:%lf is over %lf\n",
					Node[verlink->node_Vindex[i]].ID,
					Node[verlink->node_Vindex[i]].invertElev * UCF(LENGTH),
					Node[verlink->node_Vindex[i]].fullDepth * UCF(LENGTH),
					(Node[verlink->node_Vindex[i]].invertElev + Node[verlink->node_Vindex[i]].fullDepth) * UCF(LENGTH),
					Arrptr->DEM[verlink->dem_Vindex[i]],
					0.1);
				++count;
			}
		}
	}
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "FIX NODE MAX DEPTH END:\n");
	fprintf(flood_error, "====================================================================\n");
	return;
}

void flood::verticalcon::cal_exchange_flow
(
	Vconnect* verlink,
	Arrays* Arrptr,
	Pars* Parptr,
	Solver* Solverptr,
	const double* pending_volumes
)
{
	double H1D, H2D, Z, h2D;
	double flow_exchange;

	for (size_t i = 0; i < verlink->count_index; i++)
	{
		if (verlink->dem_Vindex[i] == -9999)
		{
			flow_exchange = 0.0;
			verlink->flow_1D_node[i] = 0.0;
			verlink->flow_2D_node[i] = 0.0;
			continue;
		}

		H1D = (Node[verlink->node_Vindex[i]].newDepth + Node[verlink->node_Vindex[i]].invertElev) * UCF(LENGTH);
		H2D = Arrptr->DEM[verlink->dem_Vindex[i]] + Arrptr->H[verlink->dem_Vindex[i]];
		Z = getmax(verlink->node_surfelev[i], Arrptr->DEM[verlink->dem_Vindex[i]]);

		if ((H1D - Z > COUPLING_HEAD_TOL) && (H1D - H2D > COUPLING_HEAD_TOL))
		{
				flow_exchange = orificeflow_cal(verlink->node_c0[i],
					verlink->node_area[i],
					H1D,
					H2D);
				flow_exchange = getmin(flow_exchange, verlink->node_qoutmax[i]);
				const double stored = swmm_getNodeAvailableCoupledVolume(
					verlink->node_Vindex[i]);
				if (stored >= 0.0)
				{
					const double pending = pending_volumes ? pending_volumes[i] : 0.0;
					flow_exchange = getmin(flow_exchange,
						getmax(0.0, stored + pending) / Solverptr->Tstep);
				}
				verlink->flow_1D_node[i] = -1 * flow_exchange;
				verlink->flow_2D_node[i] = flow_exchange;
		}
		else if ((H1D <= Z + COUPLING_HEAD_TOL) && (H2D - Z > COUPLING_HEAD_TOL))
		{
			// Free inlet: the surface head is above the crest while the sewer
			// water level remains below it, so use the free-weir relation.
			flow_exchange = weirflow_cal(verlink->node_cw[i],
				verlink->node_w[i],
				(H2D - Z));
			verlink->flow_1D_node[i] = flow_exchange;
			verlink->flow_2D_node[i] = flow_exchange * -1;
		}
		else if ((H1D >= Z - COUPLING_HEAD_TOL) && (H2D - H1D > COUPLING_HEAD_TOL))
		{
				flow_exchange = orificeflow_cal(verlink->node_c1[i],
					verlink->node_area1[i],
					H2D,
					H1D);
				verlink->flow_1D_node[i] = flow_exchange;
				verlink->flow_2D_node[i] = flow_exchange * -1;
		}
		else
		{
			flow_exchange = 0.0;
			verlink->flow_1D_node[i] = 0.0;
			verlink->flow_2D_node[i] = 0.0;
		}
	}

	for (size_t k = 0; k < verlink->count_index; k++)
	{
		if (verlink->dem_Vindex[k] == -9999)
		{
			continue;
		}
		if (verlink->flow_2D_node[k] < -DBL_EPSILON)
		{
			verlink->flow_neg_count[verlink->dem_Vindex[k]] += 1;
		}
	}

	for (size_t j = 0; j < verlink->count_index; j++)
	{
		if (verlink->dem_Vindex[j] == -9999) continue;
		if (verlink->flow_2D_node[j] > DBL_EPSILON)
		{
			verlink->flow_2D_dem[verlink->dem_Vindex[j]] += verlink->flow_2D_node[j] / Parptr->dy;
		}
		else if (verlink->flow_2D_node[j] < -DBL_EPSILON)
		{
			if (verlink->flow_neg_count[verlink->dem_Vindex[j]] == 0 || Arrptr->H[verlink->dem_Vindex[j]] < -DBL_EPSILON)
			{
				printf("wrong-------------------------------------------------------------\n");
				printf("%lf", (double)verlink->flow_2D_node[j]);
				continue;
			}
			else
			{
				flow_exchange = getmin((Parptr->dx * Parptr->dy * Arrptr->H[verlink->dem_Vindex[j]]) / ((verlink->flow_neg_count[verlink->dem_Vindex[j]]) * Solverptr->Tstep),
						verlink->flow_1D_node[j]);
				flow_exchange = getmin(flow_exchange, verlink->node_qinmax[j]);
				verlink->flow_1D_node[j] = flow_exchange;
				verlink->flow_2D_node[j] = flow_exchange * -1;

				verlink->flow_2D_dem[verlink->dem_Vindex[j]] += verlink->flow_2D_node[j] / Parptr->dy;
			}
		}
		else
		{
			verlink->flow_2D_dem[verlink->dem_Vindex[j]] += 0;
		}
	}
	return;
}

NUMERIC_TYPE flood::verticalcon::orificeflow_cal
(
	NUMERIC_TYPE c0,
	NUMERIC_TYPE area,
	NUMERIC_TYPE H1D,
	NUMERIC_TYPE H2D
)
{
	return c0 * area * (sqrt(2.0 * SI_GRAVITY * (H1D - H2D)));
}

NUMERIC_TYPE flood::verticalcon::weirflow_cal
(
	NUMERIC_TYPE cw,
	NUMERIC_TYPE w,
	NUMERIC_TYPE h2d
)
{
	return cw * w * h2d * (sqrt(2.0 * SI_GRAVITY * h2d));
}

void flood::verticalcon::free_vcon
(
	Vconnect* verlink
)
{
	delete[] verlink->node_Vindex;
	delete[] verlink->dem_Vindex;
	delete[] verlink->dem_count;
	delete[] verlink->node_surfelev;
	delete[] verlink->node_c0;
	delete[] verlink->node_area;
	delete[] verlink->node_qoutmax;
	delete[] verlink->node_cw;
	delete[] verlink->node_w;
	delete[] verlink->node_qinmax;
	//delete[] verlink->flow_type;
	delete[] verlink->node_c1;
	delete[] verlink->node_area1;

	delete[] verlink->flow_1D_node;
	delete[] verlink->flow_2D_node;
	delete[] verlink->flow_2D_dem;

	delete[] verlink->flow_neg_count;
	delete[] verlink->dem_Vindex_uni;
	delete[] verlink->dem_xi;
	delete[] verlink->dem_yi;
	return;
}

void flood::verticalcon::nodeinflow_cal
(
	Vconnect* verlink
)
{
	for (size_t i = 0; i < verlink->count_index; i++)
	{
		Node_Inflow[verlink->node_Vindex[i]] = verlink->flow_1D_node[i];
	}
	return;
}


//-----------------------------------------------------------------------------
// ITERATE THROUGH TIME STEPS
void IterateQ(Fnames *Fnameptr, Files *Fptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, Stage *Locptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, SGCprams *SGCptr, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, const int verbose)
{
	int i, j, ptr;
	int tstep_counter = -1;   // start at -1 so that in first run through we calculate river
	int steadyCount = 0;
	NUMERIC_TYPE tstep_channel = 0; // channel timestep
	NUMERIC_TYPE Previous_t;      // previous time channel was calculated
	NUMERIC_TYPE FloodArea, dA;
	NUMERIC_TYPE discharge = C(0.0); // value of discharge for virtual gauge output
	NUMERIC_TYPE loss; //temp variable to keep track of losses since last mass interval
	NUMERIC_TYPE Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin;
    DynamicRain<> dynamic_rain(Fnameptr->dynamicrainfilename, verbose);


	//--------------------------------------------------------------------------------------------

	double elapsedTime = 0.0;
	ErrorCode = 0;
	const bool swmm_coupled = swmm_coupling_enabled(Fnameptr);

	Vconnect vconnection = {};
	Vconnect* verlink = &vconnection;

	FILE* timefile = nullptr;
	if (swmm_coupled)
	{
		if (swmm_open(Fnameptr->swmminputfilename,
				Fnameptr->swmmreportfilename,
				Fnameptr->swmmoutputfilename) != 0 || swmm_start(1) != 0)
		{
			fprintf(stderr, "ERROR: unable to initialise SWMM coupling.\n");
			swmm_close();
			throw std::runtime_error("Unable to initialise SWMM coupling; see the SWMM report");
		}

		char passtime[256] = { 0 };
		const size_t input_name_length = strlen(Finp.name);
		if (input_name_length >= 3)
		{
			strncat(passtime, Finp.name, input_name_length - 3);
			strcat(passtime, "time");
			timefile = fopen(passtime, "w");
		}

		flood::verticalcon::load_floodinfo(Fnameptr, verlink, Parptr, Arrptr,
				verbose);
		lis::cpu_coupling::map_internal_outfalls(Fnameptr, verlink, Parptr,
			Arrptr, verbose);
	}

	//--------------------------------------------------------------------------------------------

	Solverptr->vol1 = DomainVol(Statesptr, Parptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr);

	if (verbose == ON)
	{
		printf("\nStarting time steps: ");
		fflush(stdout);
	}
	Solverptr->itrn_time_now = Solverptr->itrn_time;

	// Populating Tstep variables prior to start of simulation
	if (Statesptr->adaptive_ts == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("adaptive mode\n\n");
		fflush(stdout);
	}
	else if (Statesptr->acceleration == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("acceleration mode\n\n");
		fflush(stdout);
	}
	else if (Statesptr->Roe == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("Roe mode\n\n");
		fflush(stdout);
	}
	else
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("non-adaptive mode\n\n");
		fflush(stdout);
	}
	if (Statesptr->SGC == ON)
	{
		// because the SGC model calculates the time step in UpdateH rather than during calcFPflow it needs to initalise
		// SGCtmpTstep, which would usually be calculated in UpdateH
		Solverptr->Tstep = Solverptr->InitTstep;
		CalcT(Parptr, Solverptr, Arrptr);
		Solverptr->SGCtmpTstep = Solverptr->Tstep;
		if (verbose == ON) printf("SGC mode\n\n");
		fflush(stdout);
	}

	// set previous time to one timestep backwards so that first river calcs uses Solverptr->Tstep for 1st iteration
	// this is because of the way the timestep is calculated as the time difference from the last time the river was run
	//Previous_t = -Solverptr->Tstep;
	Previous_t = Solverptr->t - Solverptr->Tstep;

	//--------------------------------------------------------------------
	const int temp_numPS = BCptr->numPS > 0 ? BCptr->numPS : 0;
	const int surface_source_start = temp_numPS;
	const int outfall_source_start = surface_source_start +
		(swmm_coupled ? verlink->count_dem_index : 0);
	const int coupled_source_end = outfall_source_start +
		(swmm_coupled ? verlink->outfall_dem_count : 0);
	if (swmm_coupled)
	{
		int* temp_xpi = new int[coupled_source_end]();
		int* temp_ypi = new int[coupled_source_end]();
		ESourceType* temp_PS_Ident = new ESourceType[coupled_source_end];
		NUMERIC_TYPE* temp_PS_Val = new NUMERIC_TYPE[coupled_source_end]();
		if (temp_numPS > 0)
		{
			memcpy(temp_xpi, BCptr->xpi, sizeof(int) * temp_numPS);
			memcpy(temp_ypi, BCptr->ypi, sizeof(int) * temp_numPS);
			memcpy(temp_PS_Ident, BCptr->PS_Ident,
				sizeof(ESourceType) * temp_numPS);
			memcpy(temp_PS_Val, BCptr->PS_Val,
				sizeof(NUMERIC_TYPE) * temp_numPS);
		}
		for (int source = 0; source < verlink->count_dem_index; ++source)
		{
			const int p = surface_source_start + source;
			temp_xpi[p] = verlink->dem_xi[source];
			temp_ypi[p] = verlink->dem_yi[source];
			temp_PS_Ident[p] = QFIX4;
		}
		for (int source = 0; source < verlink->outfall_dem_count; ++source)
		{
			const int p = outfall_source_start + source;
			const int cell = verlink->outfall_dem_Vindex_uni[source];
			temp_xpi[p] = cell % Parptr->xsz;
			temp_ypi[p] = cell / Parptr->xsz;
			temp_PS_Ident[p] = QFIX4;
		}

		delete[] BCptr->xpi;
		delete[] BCptr->ypi;
		delete[] BCptr->PS_Ident;
		memory_free_legacy(&BCptr->PS_Val);
		BCptr->xpi = temp_xpi;
		BCptr->ypi = temp_ypi;
		BCptr->PS_Ident = temp_PS_Ident;
		BCptr->PS_Val = temp_PS_Val;
		BCptr->numPS = coupled_source_end;
	}

	if (swmm_coupled)
	{
		Node_Inflow = (double*)calloc(Nobjects[NODE], sizeof(double));
	}
	SwmmCoupleClock couple_clock;
	swmm_couple_clock_init(&couple_clock,
		!swmm_coupled || Parptr->swmm_couple_every_step,
		(double)Parptr->swmm_coupleint,
		swmm_coupled ? verlink->count_index : 0);
	// Queue positive returns for the next fixed exchange interval.
	SwmmReturnSpread return_spread;
	swmm_return_init(&return_spread, swmm_coupled && !couple_clock.every_step
		? coupled_source_end - surface_source_start : 0);
	double t_exact = (double)Solverptr->t;
	double swmm_wall = 0.0;
	if (swmm_coupled)
		printf("SWMM exchange %s, SWMM threads %d\n",
			couple_clock.every_step ? "every 2D step" : "fixed interval (equal sub-steps, returns spread over the next interval)",
			swmm_getNumThreads());
	//--------------------------------------------------------------------


    int countnum=0;
    if(swmm_coupled && (Fnameptr->hsffloodtemp[0] || Fnameptr->hsffloodend[0] || coupling_restart_file(Fnameptr->hsffloodread)) &&
        !swmm_hydraulicRestartSupported())
        throw std::runtime_error("Hydraulic restart is unavailable with SWMM runoff, RDII, storage exfiltration, routing-interface input, events or active quality routing");
    const auto restart_fields=cpu_restart_fields(*Parptr,*Arrptr,false);
    if(coupling_restart_file(Fnameptr->hsffloodread)) {
        coupling_restart_io(Fnameptr->hsffloodread,true,1,*Parptr,*Solverptr,verlink,
            swmm_coupled,couple_clock,return_spread,t_exact,restart_fields);
        Solverptr->vol1=DomainVol(Statesptr,Parptr,nullptr,Arrptr,nullptr);
    }
    auto save_restart=[&](const char* path) {
        coupling_restart_io(path,false,1,*Parptr,*Solverptr,verlink,
            swmm_coupled,couple_clock,return_spread,t_exact,restart_fields);
    };
    if(swmm_coupled) printf("SWMM threads %d; pending returns included in hydraulic restart\n",swmm_getNumThreads());

	time_t loop_start;
	time(&loop_start);
    const auto precise_loop_start=std::chrono::steady_clock::now();

	// main iteration loop
	while (t_exact < (double)Solverptr->Sim_Time-1.e-9
			&& (!swmm_coupled
					|| (NewRoutingTime < TotalDuration && !ErrorCode)))
	{


		if (timefile != nullptr)
		{
			fseek(timefile, 0, SEEK_SET);
			fprintf(timefile, "%-10.2f", NewRoutingTime / TotalDuration * 100.0);
		}

		// ChannelQ routine is run using Tstep calculated for previous iteration
		// Use default kinematic solver unless diffusive flag set
		if (Statesptr->ChannelPresent == ON)
		{
			if (tstep_counter == -1 || tstep_counter == Solverptr->ts_multiple)
				// do river calc if start of run (-1) or if timestep counter is equal to ts_multiple
			{
				tstep_channel = Solverptr->t - Previous_t; // calculate river timestep (deltaT since last river calc)
				if (Statesptr->diffusive == ON) ChannelQ_Diff(tstep_channel, Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, RiversIndexVecPtr, RiversIndexPtr);
				else ChannelQ(tstep_channel, Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, RiversIndexVecPtr, RiversIndexPtr);
				tstep_counter = 0;  // set timestep counter to zero
				Previous_t = Solverptr->t; // record previous timestep
			}
			tstep_counter++; // increment timestep counter
		}
		// Chose between a sub-grid or conventional floodplain model
		if (Statesptr->SGC == ON)
		{
			printf("SGC not processed here\n");
			//// SGC uses time step calculated in update H and stored in Solverptr->SGCtmpTstep
			//Solverptr->Tstep = Solverptr->SGCtmpTstep;
			//// sub grid floodplain models
			//SGC_FloodplainQ(Statesptr, Parptr, Solverptr, Arrptr, SGCptr);

			//SGC_BCs(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, SGCptr);
			//// Infiltration, evaporation and rainfall routines after time step update (TJF)11
			//if (Statesptr->calc_evap == ON) SGC_Evaporation(Parptr, Solverptr, Arrptr, SGCptr);
			//if (Statesptr->rainfall == ON && Statesptr->routing == OFF) SGC_Rainfall(Parptr, Solverptr, Arrptr); // CCS rainfall with routing scheme disabled
			//if (Statesptr->routing == ON) SGC_Routing(Statesptr, Parptr, Solverptr, Arrptr);	// CCS Routing scheme (controls rainfall if enabled; can be used without rainfall)
			//if (Statesptr->hazard == ON) UpdateV(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);
			//NUMERIC_TYPE Hmax = SGC_UpdateH(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, SGCptr);

			//// now calulate the next time step to be used if it is smaller than inital time step
			////Solverptr->SGCtmpTstep = getmin(Solverptr->cfl*Parptr->dx/(sqrt(Solverptr->g*Hmax)),Solverptr->InitTstep); //CCS_deletion
			//Solverptr->SGCtmpTstep = getmin(Solverptr->cfl*Parptr->min_dx_dy*Parptr->SGC_m / (sqrt(Solverptr->g*Hmax)), Solverptr->InitTstep);

			//BoundaryFlux(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr);
			//// NOTES: Update Q's handeled within SGC flux equations (SGC_FloodplainQ etc.) time-step calculation intergrated with UpdateH
		}
		else // conventional floodplain models
		{
			// Time step is reset to initial time step at the start of each FloodplainQ calculation
			if (Solverptr->t>C(0.0)) Solverptr->Tstep = Solverptr->InitTstep;

			// Select the aligned step before computing fluxes so both updates use it.
			if (Statesptr->acceleration == ON)
			{
				CalcT(Parptr, Solverptr, Arrptr);
				double remain = (double)Solverptr->Sim_Time - t_exact;
				if (swmm_coupled)
				{
					if (couple_clock.every_step)
					{
						const double r = (TotalDuration - NewRoutingTime) / 1000.0;
						if (r < remain) remain = r;
					}
					else
					{
						const double r1 = couple_clock.interval - couple_clock.elapsed;
						const double r2 = TotalDuration / 1000.0 - t_exact;
						if (r1 < remain) remain = r1;
						if (r2 < remain) remain = r2;
					}
				}
				// Treat nearby rounded rain-frame and exchange times as one target.
				const double rain_remain = dynamic_rain.time_to_next_frame(t_exact);
				if (rain_remain < remain - RAIN_FRAME_TIME_TOL) remain = rain_remain;

                if(!countnum && Parptr->hsfdt>t_exact+1.e-7) remain=(std::min)(remain,(double)Parptr->hsfdt-t_exact);
                Solverptr->Tstep = (NUMERIC_TYPE)coupling_equal_substep((double)Solverptr->Tstep, remain);
			}
			const NUMERIC_TYPE chosen_step = Solverptr->Tstep;

			FloodplainQ(Statesptr, Parptr, Solverptr, Arrptr, SGCptr);
			// CalcT in FloodplainQ only lowers the step, except on a completely dry
			// domain where it resets it to initial_tstep (no fluxes then): keep the target
			if (Statesptr->acceleration == ON && Solverptr->Tstep > chosen_step) Solverptr->Tstep = chosen_step;
			if (Statesptr->acceleration != ON)
			{
			Solverptr->Tstep = dynamic_rain.clip_step(Solverptr->t,
				getmin(Solverptr->Tstep, Solverptr->Sim_Time - Solverptr->t));
			if (swmm_coupled && couple_clock.every_step)
			{
				Solverptr->Tstep = getmin(Solverptr->Tstep,
					(NUMERIC_TYPE)((TotalDuration - NewRoutingTime) / 1000.0));
			}
			else if (swmm_coupled)
			{
				Solverptr->Tstep = (NUMERIC_TYPE)swmm_couple_clip_dt(&couple_clock,
					Solverptr->t, Solverptr->Tstep, Solverptr->Sim_Time,
					TotalDuration / 1000.0);
			}
			}
			if (Solverptr->Tstep <= C(0.0)) break;

			BCs(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);
            drain_nodata_water(Parptr, Solverptr, BCptr, Arrptr);
			if (Statesptr->drychecking == ON) DryCheck(Parptr, Solverptr, Arrptr);

			// Infiltration, evaporation and rainfall routines after time step update (TJF)
			if (Statesptr->calc_infiltration == ON) FPInfiltration(Parptr, Solverptr, Arrptr);
			if (Statesptr->calc_evap == ON) Evaporation(Parptr, Solverptr, Arrptr);
			if (Statesptr->rainfall == ON && Statesptr->routing == OFF) Rainfall(Parptr, Solverptr, Arrptr); // CCS rainfall with routing scheme disabled
			if (Statesptr->routing == ON) Routing(Statesptr, Parptr, Solverptr, Arrptr);	// CCS Routing scheme (controls rainfall if enabled; can be used without rainfall)
            dynamic_rain.update_H(Parptr, Solverptr, Arrptr);

			// Net rain/infiltration first, then inlet exchange.
			// This prevents both processes from removing the same pre-step water.
			if (swmm_coupled)
			{
				BCptr->Qpoint_pos = C(0.0);
				BCptr->Qpoint_neg = C(0.0);
				memset(verlink->flow_neg_count, 0,
					Parptr->xsz * Parptr->ysz * sizeof(int));
				memset(verlink->flow_2D_dem, 0,
					Parptr->xsz * Parptr->ysz * sizeof(double));
				memset(verlink->flow_1D_node, 0,
					verlink->count_index * sizeof(double));
				memset(verlink->flow_2D_node, 0,
					verlink->count_index * sizeof(double));
				memset(Node_Inflow, 0, Nobjects[NODE] * sizeof(double));

				flood::verticalcon::cal_exchange_flow(verlink, Arrptr, Parptr,
					Solverptr, couple_clock.volume);
				if (couple_clock.every_step)
					flood::verticalcon::nodeinflow_cal(verlink);
				else
					swmm_couple_accumulate(&couple_clock, verlink,
						(double)Solverptr->Tstep);

				for (int source = 0; source < verlink->count_dem_index; ++source)
					BCptr->PS_Val[surface_source_start + source] =
						verlink->flow_2D_dem[verlink->dem_Vindex_uni[source]];
				for (int source = 0; source < verlink->outfall_dem_count; ++source)
					BCptr->PS_Val[outfall_source_start + source] = C(0.0);
				// Apply queued returns from the preceding interval.
				if (return_spread.count > 0)
					swmm_return_apply(&return_spread, BCptr->PS_Val + surface_source_start,
						(double)Parptr->dy, (double)Solverptr->Tstep);
				lis::cpu_coupling::apply_coupled_sources(BCptr, Parptr,
					Solverptr, Arrptr, surface_source_start,
					return_spread.count > 0 ? coupled_source_end : outfall_source_start);

				const int exchange_now = swmm_couple_due(&couple_clock,
					t_exact + (double)Solverptr->Tstep, Solverptr->Sim_Time,
					TotalDuration / 1000.0);
				if (exchange_now)
				{
					if (couple_clock.every_step) linkDt = Solverptr->Tstep;
					else
					{
						swmm_couple_write_inflow(&couple_clock, verlink);
						linkDt = couple_clock.elapsed;
					}
					swmm_beginCoupledExchange();
					lis::cpu_coupling::update_outfall_stages(verlink, Parptr,
						Arrptr);
					const auto sw0 = std::chrono::steady_clock::now();
					ErrorCode = swmm_advanceCoupled((double)linkDt, &elapsedTime);
					swmm_wall += std::chrono::duration<double>(std::chrono::steady_clock::now() - sw0).count();
					if (!ErrorCode)
					{
						if (return_spread.count > 0)
						{
							// Queue positive returns; apply reverse outfall flow immediately.
							lis::cpu_coupling::collect_return_sources(verlink, BCptr,
								Parptr, Arrptr, surface_source_start,
								outfall_source_start, (double)linkDt, (double)linkDt);
							const double t_after = t_exact + (double)Solverptr->Tstep;
							double next_interval = couple_clock.interval;
							if ((double)Solverptr->Sim_Time - t_after < next_interval) next_interval = (double)Solverptr->Sim_Time - t_after;
							if (TotalDuration / 1000.0 - t_after < next_interval) next_interval = TotalDuration / 1000.0 - t_after;
							if(next_interval<1.e-6) next_interval=0.0;
							swmm_return_schedule(&return_spread, BCptr->PS_Val + surface_source_start,
								(double)linkDt, next_interval, (double)Parptr->dy, (double)Solverptr->Tstep);
						}
						else
							lis::cpu_coupling::collect_return_sources(verlink, BCptr,
								Parptr, Arrptr, surface_source_start,
								outfall_source_start, (double)linkDt,
								(double)Solverptr->Tstep);
						lis::cpu_coupling::apply_coupled_sources(BCptr, Parptr,
							Solverptr, Arrptr, surface_source_start,
							coupled_source_end);
					}
					if (!couple_clock.every_step)
						swmm_couple_advance(&couple_clock);
				}
			}

			// If Roe == ON used UpdateHRoe else use the normal UpdateH (JN/IV)
			// if(Statesptr->Roe==ON) UpdateHRoe(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
			// else UpdateH(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
			// change to CalcFPQxRoe and CalcFPQyRoe... now return Q's for UpdateH (JCN)
			UpdateH(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr,
				swmm_coupled ? temp_numPS : BCptr->numPS, !swmm_coupled);

			if (Statesptr->hazard == ON) UpdateV(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);

			BoundaryFlux(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr);
			if (Statesptr->acceleration == ON) UpdateQs(Parptr, Arrptr); // Old value of Q updated with current timestep values
			if (Statesptr->Roe == ON) UpdateQsRoe(Parptr, Solverptr, Arrptr); // If Roe==on UpdateQsRoe(JN/IV)
		}

		// Update t with final Tstep calculated
		if (Solverptr->t>C(0.0)) Solverptr->MinTstep = getmin(Solverptr->MinTstep, Solverptr->Tstep);
		t_exact += (double)Solverptr->Tstep;
		Solverptr->t = (NUMERIC_TYPE)t_exact;
		Solverptr->itCount += 1;



		// Update time of initial flood inundation (in hours) and total inundation time (in seconds)
		if ((Statesptr->reset_timeinit == ON) & (Solverptr->t>Parptr->reset_timeinit_time))
		{ //reset the time of initial inundation if called for in parameter file
			Statesptr->reset_timeinit = OFF;

			//#pragma omp parallel for private(j,ptr) // irrelevent benefit JCN
			for (j = 0; j<Parptr->ysz; j++) for (i = 0; i<Parptr->xsz; i++)
			{
				ptr = i + j*Parptr->xsz;
				Arrptr->initHtm[ptr] = (NULLVAL);
			}
			if (verbose == ON) printf("\n Time of initial inundation reset \n");
		}
		// Update maxH, maxHtm, totalHtm and initHtm at the mass interval if mint_hk is specifed in the .par file OR at every time step if not
		if (Solverptr->t >= Parptr->MassTotal || Statesptr->mint_hk == OFF)
		{
#pragma omp parallel for private (i, ptr)
			for (j = 0; j<Parptr->ysz; j++) for (i = 0; i<Parptr->xsz; i++)
			{
				ptr = i + j*Parptr->xsz;
				if (Arrptr->H[ptr] > Solverptr->DepthThresh)
				{
					if (Arrptr->initHtm[ptr] == NULLVAL /*&& (Arrptr->H[ptr]>C(0.01)*/)
						Arrptr->initHtm[ptr] = Solverptr->t / C(3600.0);
					//if (Arrptr->H[ptr]>C(0.01))
					Arrptr->totalHtm[ptr] += Solverptr->Tstep / C(3600.0);
					// Update maximum water depths, and time of maximum (in hours) TDF updated only track maxH when h > DepthThresh
					if (Arrptr->H[ptr] > Arrptr->maxH[ptr])
					{
						Arrptr->maxH[ptr] = Arrptr->H[ptr];
						Arrptr->maxHtm[ptr] = Solverptr->t / C(3600.0);
					}
				}
			}
		}
		// Calculate mass balance error
		if (Solverptr->t >= Parptr->MassTotal)
		{
			Solverptr->vol2 = DomainVol(Statesptr, Parptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr); // CCS

			// calc losses for this mass interval
			loss = (Parptr->InfilTotalLoss - Parptr->InfilLoss) + (Parptr->EvapTotalLoss - Parptr->EvapLoss) - (Parptr->RainTotalLoss - Parptr->RainLoss);

			//Solverptr->Qerror=BCptr->Qin-BCptr->Qout-(Solverptr->vol2+loss-Solverptr->vol1)/Parptr->MassInt;
			// New version using VolInMT and VolOutMT
			// volume error
			Solverptr->Verror = BCptr->VolInMT - BCptr->VolOutMT - (Solverptr->vol2 + loss - Solverptr->vol1);

			// Q error
			Solverptr->Qerror = Solverptr->Verror / Parptr->MassInt;
			// reset to C(0.0)
			BCptr->VolInMT = C(0.0);
			BCptr->VolOutMT = C(0.0);

			// record cumulative loss for next time.
			Parptr->InfilLoss = Parptr->InfilTotalLoss;
			Parptr->EvapLoss = Parptr->EvapTotalLoss;
			Parptr->RainLoss = Parptr->RainTotalLoss;

			// Calculate flood area
			FloodArea = C(0.0);
			dA = Parptr->dA;
#pragma omp parallel for private(i,ptr) firstprivate(dA) reduction ( + : FloodArea)
			for (j = 0; j<Parptr->ysz; j++) for (i = 0; i<Parptr->xsz; i++)
			{
				ptr = i + j*Parptr->xsz;
				if (Statesptr->latlong == ON) dA = Arrptr->dA[ptr]; // if latlong is on change dA to local cell area
				if (Statesptr->porosity == ON)
				{
					if (Arrptr->H[ptr]>C(0.01)) FloodArea += dA*Arrptr->paerial[ptr]; // If porosity used, scale flooded area by porosity (TJF)
				}
				else if (Statesptr->SGC == ON)
				{
					if (Arrptr->H[ptr] - Arrptr->SGCbfH[ptr] > Solverptr->DepthThresh)
						FloodArea += dA; // If sub-grid used remove channel depth
				}
				else
				{
					if (Arrptr->H[ptr] > C(0.01)) FloodArea += dA; // standard area calculation
				}
			}
			Solverptr->FArea = FloodArea;

			// output mass variables line to file
			//fprintf(Fptr->mass_fp,"%-12.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10ld %12.4e %12.4e  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4e %12.4e %12.4e\n",Solverptr->t,Solverptr->Tstep,Solverptr->MinTstep,Solverptr->itCount,Solverptr->FArea,Solverptr->vol2,BCptr->Qin,Solverptr->Hds,BCptr->Qout,Solverptr->Qerror,Solverptr->Verror,Parptr->RainTotalLoss-(Parptr->InfilTotalLoss+Parptr->EvapTotalLoss));
			fprintf(Fptr->mass_fp, "%-12.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10ld %12.4" NUM_FMT" %12.4" NUM_FMT"  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT"\n", Solverptr->t, Solverptr->Tstep, Solverptr->MinTstep, Solverptr->itCount, Solverptr->FArea, Solverptr->vol2, BCptr->Qin, Solverptr->Hds, BCptr->Qout, Solverptr->Qerror, Solverptr->Verror, Parptr->RainTotalLoss - (Parptr->InfilTotalLoss + Parptr->EvapTotalLoss));
			fflush(Fptr->mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.

			Solverptr->vol1 = Solverptr->vol2;
			Parptr->MassTotal += Parptr->MassInt;

			//stage output
			if (Statesptr->save_stages == ON)
			{
				fprintf(Fptr->stage_fp, "%12.3" NUM_FMT"", Solverptr->t);
				for (i = 0; i<Locptr->Nstages; i++)
				{
					int index = Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i] * Parptr->xsz;
					if (Locptr->stage_check[i] == 1) fprintf(Fptr->stage_fp, "%10.4" NUM_FMT"", Arrptr->H[index]);
					else fprintf(Fptr->stage_fp, "-\t");
				}
				fprintf(Fptr->stage_fp, "\n");
				fflush(Fptr->stage_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
				// added to export scalar velocity
				if (Statesptr->voutput_stage == ON)
				{
					fprintf(Fptr->vel_fp, "%12.3" NUM_FMT"", Solverptr->t);
					for (i = 0; i<Locptr->Nstages; i++)
					{
						int index = Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i] * (Parptr->xsz + 1);
						if (Locptr->stage_check[i] == 1) fprintf(Fptr->vel_fp, "%10.4" NUM_FMT"", sqrt(pow(getmax(fabs(Arrptr->Vx[index]), fabs(Arrptr->Vx[index + 1])), 2) + pow(getmax(fabs(Arrptr->Vy[index]), fabs(Arrptr->Vy[index + (Parptr->xsz + 1)])), 2)));
						else fprintf(Fptr->vel_fp, "-\t");
					}
					fprintf(Fptr->vel_fp, "\n");
					fflush(Fptr->vel_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
				}
			}

			//virtual gauge output
			if (Statesptr->gsection == ON)
			{
				fprintf(Fptr->gau_fp, "%12.2" NUM_FMT"", Solverptr->t); // print tiem to file
				for (i = 0; i<Locptr->Ngauges; i++) // loop through each virtual gauge
				{
					// call discharge calculation function
					discharge = CalcVirtualGauge(i, Parptr->xsz + 1, Arrptr->Qx, Arrptr->Qy, Locptr);

					fprintf(Fptr->gau_fp, " %10.3" NUM_FMT"", discharge); // Print discharge to file
				}
				fprintf(Fptr->gau_fp, "\n"); // print end of line
				fflush(Fptr->gau_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
			}

			// Checkpointing
			if (Statesptr->checkpoint == ON)
			{
				//iteration time
				time(&Solverptr->time_check);
				Solverptr->itrn_time_now = Solverptr->itrn_time + (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start);
				if (Solverptr->itrn_time_now >= Parptr->nextcheck)
				{
					WriteCheckpoint(Fnameptr, Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, verbose);
					Parptr->nextcheck = Solverptr->itrn_time_now + (Parptr->checkfreq * 3600);
				}
			}
		}

		// Regular output

		if (Solverptr->t >= Parptr->SaveTotal)
		{

			time(&Solverptr->time_check);
			Comp_time = (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start) / 60;
			if (Comp_time != 0 && Statesptr->comp_out == ON) // only of t is not zero (can't divide by zero)
			{
				Model_Comp_Ratio = ((Solverptr->t / 60) / Comp_time);
				Model_time_left = (Solverptr->Sim_Time - Solverptr->t) / 60;
				Est_Time_Fin = (Model_time_left / Model_Comp_Ratio);
				Est_Time_Tot = Comp_time + Est_Time_Fin;
				printf("T(mins): M: %.1" NUM_FMT", C: %.1" NUM_FMT", M/C: %.2" NUM_FMT", ETot: %.1" NUM_FMT", EFin: %.1" NUM_FMT"\n", (Solverptr->t / C(60.0)), Comp_time, Model_Comp_Ratio, Est_Time_Tot, Est_Time_Fin);
			}

			write_regular_output(Fnameptr, Solverptr, Statesptr, Parptr, Arrptr, SGCptr);

			//regular profile output, if requested in param file
			if (Statesptr->profileoutput == ON)
			{
				write_profile(Fnameptr->resrootname, Parptr->SaveNo, ".profile", Statesptr, ChannelSegments, Arrptr, Parptr, RiversIndexVecPtr, RiversIndexPtr); // CCS
			}

			// update interval counter
			Parptr->SaveTotal += Parptr->SaveInt;
			Parptr->SaveNo += 1;

		}

		// Single overpass
		if (Statesptr->single_op == ON && Solverptr->t >= Parptr->op)
		{
			// write rasters
			if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".opb", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
			else write_ascfile(Fnameptr->resrootname, -1, ".op", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
			// write profiles
			if (Statesptr->profileoutput == ON) write_profile(Fnameptr->resrootname, -1, ".profile", Statesptr, ChannelSegments, Arrptr, Parptr, RiversIndexVecPtr, RiversIndexPtr);
			if (verbose == ON) printf("Writing overpass at %" NUM_FMT" seconds\n", Solverptr->t);
			Statesptr->single_op = OFF;
			// raster elevation output
			if (Statesptr->save_elev == ON)
			{
				if (Statesptr->SGC == ON) // SGC output
				{
					if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".opelevb", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
					else write_ascfile(Fnameptr->resrootname, -1, ".opelev", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
				}
				else // standard model output
				{
					if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".opelevb", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
					else write_ascfile(Fnameptr->resrootname, -1, ".opelev", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
				}
			}
		}

		// Multiple overpasses
		if (Statesptr->multi_op == ON)
		{
			for (i = 0; i<Parptr->op_multinum; i++)
			{
				if (Solverptr->t >= Parptr->op_multisteps[i] && Parptr->op_multiswitch[i] == 0)
				{
					// raster depth output
					if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, i, "-T.opb", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
					else write_ascfile(Fnameptr->resrootname, i, "-T.op", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
					// profiles output
					if (Statesptr->profileoutput == ON) write_profile(Fnameptr->resrootname, i, "-T.profile", Statesptr, ChannelSegments, Arrptr, Parptr, RiversIndexVecPtr, RiversIndexPtr);

					Parptr->op_multiswitch[i] = 1;
					if (verbose == ON) printf("Writing overpass %d at %" NUM_FMT" seconds\n", i, Parptr->op_multisteps[i]);
					// raster elevation output
					if (Statesptr->save_elev == ON)
					{
						if (Statesptr->SGC == ON) // SGC output
						{
							if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, i, "-T.opelevb", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
							else write_ascfile(Fnameptr->resrootname, i, "-T.opelev", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
						}
						else // standard model output
						{
							if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, i, "-T.opelevb", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
							else write_ascfile(Fnameptr->resrootname, i, "-T.opelev", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
						}
					}
				}
			}
		}

		// If requested on command line, check whether we should kill this simulation...
		if (Statesptr->killsim == ON)
		{
			//iteration time
			time(&Solverptr->time_check);
			Solverptr->itrn_time_now = Solverptr->itrn_time + (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start);
			// check if we have reached the kill time
			if (Solverptr->itrn_time_now >= Parptr->killsim_time) {
				if (verbose == ON) printf("Simulation kill time reached... ");
				break;
			}
		}

		// If an auto steady-state has been requested:
		if (Statesptr->steadycheck == ON && Solverptr->t >= Parptr->steadyTotal)
		{
			Parptr->steadyQdiff = BCptr->Qin - BCptr->Qout;
			if (fabs(Parptr->steadyQdiff)<Parptr->steadyQtol) {
				steadyCount += 1; // keep going a bit to make sure...
				if (steadyCount == 10) break;
			}
			else {
				steadyCount = 0;
			}
			Parptr->steadyTotal += Parptr->steadyInt;
		}
        if(!countnum && Parptr->hsfdt>0 && t_exact+1.e-7>=(double)Parptr->hsfdt) {
            if(Fnameptr->hsffloodtemp[0]) save_restart(Fnameptr->hsffloodtemp);
            if(swmm_coupled && Fnameptr->hsfswmmtemp[0]) {
                FILE* legacy=fopen(Fnameptr->hsfswmmtemp,"wb");
                if(!legacy) throw std::runtime_error("Cannot write SWMM warm start");
                saveTempHsf(legacy);
            }
            countnum=1;
        }

	}
	//END main ITERATIONS

	if (timefile != nullptr) fclose(timefile);


    if(Fnameptr->hsffloodend[0]) {
        save_restart(Fnameptr->hsffloodend);
        Fnameptr->hsffloodend[0]='\0'; // Prevent the main driver overwriting the bundled restart.
    }
    time_t loop_end;
	time(&loop_end);

	double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-precise_loop_start).count();
	printf("loop time %lf\n", seconds);
	if (swmm_coupled)
		printf("TIMING loop %.3f s, SWMM %.3f s (%.0f%%), 2D steps %ld, SWMM returns not yet applied %.6g m3\n",
			(double)seconds, swmm_wall, 100.0 * swmm_wall / ((double)seconds > 0 ? (double)seconds : 1.0),
			(long)Solverptr->itCount, swmm_return_pending(&return_spread));
	swmm_return_free(&return_spread);

	if (swmm_coupled)
	{
		swmm_end();
		if (Fout.mode == SCRATCH_FILE) swmm_report();

		// --- close the system
		swmm_close();
		free(Node_Inflow);
		Node_Inflow = nullptr;
	}
	swmm_couple_clock_free(&couple_clock);
	//output regular files
	fileoutput(Fnameptr, Statesptr, Parptr, Arrptr);

	// Must be last because maxH changed to Maximum elevation !!!
	// Write maximum elevation
	if (Statesptr->SGC == ON) // SGC output
	{
		for (i = 0; i<Parptr->xsz; i++) for (j = 0; j<Parptr->ysz; j++)
		{
			ptr = i + j*Parptr->xsz;
			if (Arrptr->maxH[ptr]>1e-3) Arrptr->maxH[ptr] += Arrptr->SGCz[ptr]; else Arrptr->maxH[ptr] = NULLVAL;
		}
		if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".mxeb", Arrptr->maxH, Arrptr->SGCz, 0, Statesptr, Parptr);
		else write_ascfile(Fnameptr->resrootname, -1, ".mxe", Arrptr->maxH, Arrptr->SGCz, 0, Statesptr, Parptr);
	}
	else
	{
		for (i = 0; i<Parptr->xsz; i++) for (j = 0; j<Parptr->ysz; j++)
		{
			ptr = i + j*Parptr->xsz;
			if (Arrptr->maxH[ptr]>1e-3) Arrptr->maxH[ptr] += Arrptr->DEM[ptr]; else Arrptr->maxH[ptr] = NULLVAL;
		}
		if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".mxeb", Arrptr->maxH, Arrptr->DEM, 0, Statesptr, Parptr);
		else write_ascfile(Fnameptr->resrootname, -1, ".mxe", Arrptr->maxH, Arrptr->DEM, 0, Statesptr, Parptr);
	}

	if (swmm_coupled)
	{
		lis::cpu_coupling::free_outfalls(verlink);
		flood::verticalcon::free_vcon(verlink);
	}

	if (verbose == ON) printf("Finished.\n\n");

	return;
}

//----------------------------------------------------------------------------
// SUM Qs INTO A CELL AND UPDATE DEPTH ACCORDINGLY
void UpdateH(States *Statesptr, Pars *Parptr, Solver *Solverptr,
	BoundCs *BCptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr,
	int point_source_end, bool reset_point_stats)
{
	int i, j;
	NUMERIC_TYPE *qxptr0, *qyptr0, *qyptr1, *hptr;
	int *mptr;
	NUMERIC_TYPE dV, himp, qtmp;
	NUMERIC_TYPE dAPorTemp;
	NUMERIC_TYPE Qpnt;

	// Insert point sources ((MT) moved before dV check, in case flow out of cell results in negative H reset and the inflow would prevent this)
	// H point sources moved back to after UpdateH (JCN)
	if (reset_point_stats)
	{
		BCptr->Qpoint_pos = C(0.0);
		BCptr->Qpoint_neg = C(0.0);
	}

	for (int ps_index = 0; ps_index<point_source_end; ps_index++)
	{
		if (BCptr->PS_Ident[ps_index] == QFIX4) // QFIX
		{
			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] += BCptr->PS_Val[ps_index] * Parptr->dx*Solverptr->Tstep / Parptr->dA;
			Qpnt = BCptr->PS_Val[ps_index] * Parptr->dx;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;
		}
		if (BCptr->PS_Ident[ps_index] == QVAR5) // QVAR
		{
			qtmp = InterpolateTimeSeries(BCptr->PS_TimeSeries[ps_index], Solverptr->t);
			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] += qtmp*Parptr->dx*Solverptr->Tstep / Parptr->dA;
			Qpnt = qtmp*Parptr->dx;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;
		}
	}

	// Calculate dV (-ve => outflow) and update NewH
#pragma omp parallel for private( i,qxptr0,qyptr0,qyptr1,hptr,mptr,dV,dAPorTemp)
	for (j = 0; j<Parptr->ysz; j++)
	{
		qxptr0 = Arrptr->Qx + j*(Parptr->xsz + 1);
		qyptr0 = Arrptr->Qy + j*(Parptr->xsz + 1);
		qyptr1 = Arrptr->Qy + (j + 1)*(Parptr->xsz + 1);
		hptr = Arrptr->H + j*Parptr->xsz;
		mptr = Arrptr->ChanMask + j*Parptr->xsz;
		for (i = 0; i<Parptr->xsz; i++)
		{
			if (*mptr == -1)
			{
				if (Statesptr->porosity == ON)
				{
					dV = Solverptr->Tstep*(*qxptr0 - *(qxptr0 + 1) + *qyptr0 - *qyptr1);
					dAPorTemp = PorArea(i, j, Parptr, Arrptr);
					if (dAPorTemp == C(0.0)) (*hptr) += C(0.0);
					else (*hptr) += dV / dAPorTemp;

					if (*hptr<C(0.0)) *hptr = C(0.0);
				}
				else
				{
					dV = Solverptr->Tstep*(*qxptr0 - *(qxptr0 + 1) + *qyptr0 - *qyptr1);
					(*hptr) += dV / Parptr->dA;
					if (*hptr<C(0.0)) *hptr = C(0.0);


					//if (j == 0 && i == 79) {

					//	bool test;
					//	test = 1;
					//}
				}
			}
			qxptr0++;
			qyptr0++;
			qyptr1++;
			hptr++;
			mptr++;
		}
	}

	// Point source HVAR and HFIX
	for (int ps_index = 0; ps_index<point_source_end; ps_index++)
	{
		if (BCptr->PS_Ident[ps_index] == HFIX2) // HFIX
		{
			himp = BCptr->PS_Val[ps_index] - Arrptr->DEM[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz];
			if (himp<C(0.0)) himp = C(0.0);

			Qpnt = (himp - Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz])*Parptr->dA / Solverptr->Tstep;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;

			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] = himp;
		}
		if (BCptr->PS_Ident[ps_index] == HVAR3) // HVAR
		{
			himp = InterpolateTimeSeries(BCptr->PS_TimeSeries[ps_index], Solverptr->t) - Arrptr->DEM[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz];
			if (himp<C(0.0)) himp = C(0.0);

			Qpnt = (himp - Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz])*Parptr->dA / Solverptr->Tstep;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;

			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] = himp;
		}
	}

	return;
}
