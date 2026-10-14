#include "fv1.h"
#include "fv1/modifiedvars.h"
#include "boundary.h"
#include "hll.h"
#include "output.h"
#include "stats.h"
#include <algorithm>
#include <stdio.h>
#include <float.h>
#include "../geometry.h"
extern "C" {
#include "../headers.h"
};
#include "../swmm5.h"

namespace floodfv1 {
	struct verticalcon
	{
		static void load_floodinfo(Fnames* Fnameptr, Vconnect* verlink, Pars* Parptr, Arrays* Arrptr, int verbose);
		static void vcon_init(Vconnect* verlink, Pars* Parptr);
		static void read_vconinfo(Vconnect* verlink, FILE* flood_file, int verbose);
		static void cal_dem_index(Vconnect* verlink, FILE* flood_error, Pars* Parptr, Arrays* Arrptr);
		static void demcountuni_cal(Vconnect* verlink);
		static void check_nodedepth(Vconnect* verlink, FILE* flood_error, Arrays* Arrptr);
		static void cal_exchange_flow(Vconnect* verlink, Arrays* Arrptr, Pars* Parptr, Solver* Solverptr);
		static NUMERIC_TYPE orificeflow_cal(NUMERIC_TYPE c0, NUMERIC_TYPE area, NUMERIC_TYPE H1D, NUMERIC_TYPE H2D);
		static NUMERIC_TYPE weirflow_cal(NUMERIC_TYPE cw, NUMERIC_TYPE w, NUMERIC_TYPE h2d);
		static void free_vcon(Vconnect* verlink);
		static void nodeinflow_cal(Vconnect* verlink);
	};
}

void floodfv1::verticalcon::load_floodinfo
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

void floodfv1::verticalcon::vcon_init
(
	Vconnect* verlink,
	Pars* Parptr
)
{
	//固定值----------------------------------------------------------------
	verlink->node_Vindex = new int[verlink->count_index]();                       //节点配置文件中的节点对应swmm中node数组的index
	verlink->dem_Vindex = new int[verlink->count_index]();                        //节点配置文件中的节点对应dem网格的索引
	verlink->dem_count = new int[Parptr->xsz * Parptr->ysz]();                    //dem每个网格耦合的节点数量
	verlink->node_surfelev = new double[verlink->count_index]();                  //节点配置文件中的节点对应swmm中的自然地面标高
	verlink->node_c0 = new double[verlink->count_index]();                        //节点配置文件中的节点孔口系数
	verlink->node_area = new double[verlink->count_index]();					  //节点配置文件中的节点过流截面积
	verlink->node_qoutmax = new double[verlink->count_index]();                   //节点配置文件中的节点溢流最大流量
	verlink->node_cw = new double[verlink->count_index]();                        //节点配置文件中的节点堰流系数
	verlink->node_w = new double[verlink->count_index]();                         //节点配置文件中的节点周长或边长
	verlink->node_qinmax = new double[verlink->count_index]();                    //节点配置文件中的节点回水最大流量
	verlink->node_c1 = new double[verlink->count_index]();                        //节点配置文件中的节点孔口系数特指入流
	verlink->node_area1 = new double[verlink->count_index]();                     //节点配置文件中的节点过流截面积特指入流

	verlink->count_dem_index = verlink->count_index;                              //节点配置文件中的节点对应dem网格的索引个数(去重之后)
	verlink->dem_Vindex_uni = new int[verlink->count_index]();                    //节点配置文件中的节点对应dem网格的索引列表(去重之后)
	//verlink->flow_type = new int[verlink->count_index]();                       //节点配置文件中的节点对应交换流量计算方法

	verlink->dem_xi = new int[verlink->count_index]();                            //节点配置文件中的节点对应DEM的x索引
	verlink->dem_yi = new int[verlink->count_index]();							  //节点配置文件中的节点对应DEM的y索引
	//----------------------------------------------------------------------
	//每轮计算改变值---------------------------------------------------------
	verlink->flow_1D_node = new double[verlink->count_index]();                   //节点配置文件中的节点 溢流或回流流量数组(对应一维)
	verlink->flow_2D_node = new double[verlink->count_index]();                   //节点配置文件中的节点 溢流或回流流量数组(对应二维)
	verlink->flow_2D_dem = new double[Parptr->xsz * Parptr->ysz]();               //二维流量各个网格的溢流或回流
	verlink->flow_neg_count = new int[Parptr->xsz * Parptr->ysz]();               //二维回流计数
	//----------------------------------------------------------------------
	return;
}

void floodfv1::verticalcon::read_vconinfo(Vconnect* verlink, FILE* flood_file, int verbose)
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

void floodfv1::verticalcon::cal_dem_index
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
		else if (Node[verlink->node_Vindex[j]].type == OUTFALL || Node[verlink->node_Vindex[j]].type == STORAGE || Node[verlink->node_Vindex[j]].type == DIVIDER)
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


void floodfv1::verticalcon::demcountuni_cal
(
	Vconnect* verlink
)
{
	int count_init = verlink->count_index;
	int k_temp = 0; //记录多少个无index的node
	memcpy(verlink->dem_Vindex_uni, verlink->dem_Vindex, verlink->count_index * sizeof(int));
	//去除index==-9999的索引
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

	// 记录每个网格有多少个节点相连
	for (int k = 0; k < verlink->count_index; ++k)
	{
		if (verlink->dem_Vindex[k] != -9999)
			verlink->dem_count[verlink->dem_Vindex[k]] += 1;
	}

	// demindex去重
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


void floodfv1::verticalcon::check_nodedepth
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
	//printf("\n%d %d\n", count, count1);
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "FIX NODE MAX DEPTH END:\n");
	fprintf(flood_error, "====================================================================\n");
	return;
}

void floodfv1::verticalcon::cal_exchange_flow
(
	Vconnect* verlink,
	Arrays* Arrptr,
	Pars* Parptr,
	Solver* Solverptr
)
{
	double H1D, H2D, Z, h2D; //一维水面标高 二维水面标高 分界面标高 二维水深 m
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

		//printf("node:%d,H1D:%lf,H2D:%lf,", i, H1D, H2D);

		if ((H1D - Z > COUPLING_HEAD_TOL) && (H1D - H2D > COUPLING_HEAD_TOL))
		{
			flow_exchange = orificeflow_cal(verlink->node_c0[i],
				verlink->node_area[i],
				H1D,
				H2D);
			flow_exchange = getmin(flow_exchange, verlink->node_qoutmax[i]);
			verlink->flow_1D_node[i] = -1 * flow_exchange;
			verlink->flow_2D_node[i] = flow_exchange;
			//printf("C:%lf,area:%lf,par:%lf\n", verlink->node_c0[i], verlink->node_area[i], 2.0 * SI_GRAVITY * (H1D - H2D));
		}
		else if ((H1D <= Z + COUPLING_HEAD_TOL) && (H2D - Z > COUPLING_HEAD_TOL))
		{
			
			flow_exchange = weirflow_cal(verlink->node_cw[i],
				verlink->node_w[i],
				(H2D - Z));
			verlink->flow_1D_node[i] = flow_exchange;
			verlink->flow_2D_node[i] = flow_exchange * -1;
			/*
			flow_exchange = orificeflow_cal(verlink->node_c1[i],
				verlink->node_area1[i],
				H2D,
				H1D);
			verlink->flow_1D_node[i] = flow_exchange;
			verlink->flow_2D_node[i] = flow_exchange * -1;
			*/
			//printf("C:%lf,area:%lf,par:%lf\n", verlink->node_c1[i], verlink->node_area1[i], 2.0 * SI_GRAVITY * (H2D - H1D));
		}
		else if ((H1D >= Z - COUPLING_HEAD_TOL) && (H2D - H1D > COUPLING_HEAD_TOL))
		{
			flow_exchange = orificeflow_cal(verlink->node_c1[i],
				verlink->node_area1[i],
				H2D,
				H1D);
			verlink->flow_1D_node[i] = flow_exchange;
			verlink->flow_2D_node[i] = flow_exchange * -1;
			//printf("C:%lf,area:%lf,par:%lf\n", verlink->node_c1[i], verlink->node_area1[i], 2.0 * SI_GRAVITY * (H2D - H1D));
		}
		else
		{
			flow_exchange = 0.0;
			verlink->flow_1D_node[i] = 0.0;
			verlink->flow_2D_node[i] = 0.0;
			//printf("\n");
		}
	}
	//printf("--------------------------------------------------------------------------\n");
	//统计各个网格回流数
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
				printf("lf", verlink->flow_2D_node[j]);
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
		//printf("flow1d:%lf,flow2d:%lf,flow2ddem:%f\n", verlink->flow_1D_node[j], verlink->flow_2D_node[j], verlink->flow_2D_dem[verlink->dem_Vindex[j]]);
	}
	//printf("******************************************************************************\n");
	return;
}

NUMERIC_TYPE floodfv1::verticalcon::orificeflow_cal
(
	NUMERIC_TYPE c0,
	NUMERIC_TYPE area,
	NUMERIC_TYPE H1D,
	NUMERIC_TYPE H2D
)
{
	return c0 * area * (sqrt(2.0 * SI_GRAVITY * (H1D - H2D)));
}

NUMERIC_TYPE floodfv1::verticalcon::weirflow_cal
(
	NUMERIC_TYPE cw,
	NUMERIC_TYPE w,
	NUMERIC_TYPE h2d
)
{
	return cw * w * h2d * (sqrt(2.0 * SI_GRAVITY * h2d));
}

void floodfv1::verticalcon::free_vcon
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

void floodfv1::verticalcon::nodeinflow_cal
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

void fv1::solve
(
	Fnames *Fnameptr,
	Files *Fptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Stage *Stageptr,
	Arrays *Arrptr,
	const int verbose
)
{
	//耦合结构体设置
	//--------------------------------------------------------------------------------------------

	double elapsedTime = 0.0;
	ErrorCode = 0;
	const bool swmm_coupled = swmm_coupling_enabled(Fnameptr);

	Vconnect vconnection = {};
	Vconnect* verlink = &vconnection;

	if (swmm_coupled)
	{
		if (swmm_open(Fnameptr->swmminputfilename,
				Fnameptr->swmmreportfilename,
				Fnameptr->swmmoutputfilename) != 0 || swmm_start(1) != 0)
		{
			fprintf(stderr, "ERROR: unable to initialise SWMM coupling.\n");
			swmm_close();
			return;
		}

		floodfv1::verticalcon::load_floodinfo(Fnameptr, verlink, Parptr,
				Arrptr, verbose);
	}
	//--------------------------------------------------------------------------------------------

	Solverptr->Tstep = Solverptr->InitTstep;
	Solverptr->MinTstep = Solverptr->InitTstep;
	Solverptr->vol1 = DomainVol(Statesptr, Parptr, nullptr, Arrptr, nullptr);
	initialise_Zstar(Parptr, Arrptr);
	if (Statesptr->adaptive_ts == ON)
	{
		Solverptr->Tstep = Tstep_from_cfl(Parptr, Solverptr, Arrptr);
	}
    DynamicRain<> rain(Fnameptr->dynamicrainfilename, verbose);

	//--------------------------------------------------------------------
	const int temp_numPS = BCptr->numPS;
	if (swmm_coupled)
	{
	if (temp_numPS > 0)
	{
		int* temp_xpi = new int[temp_numPS + verlink->count_dem_index]();
		int* temp_ypi = new int[temp_numPS + verlink->count_dem_index]();
		ESourceType* temp_PS_Ident = new ESourceType[temp_numPS + verlink->count_dem_index];
		NUMERIC_TYPE* temp_PS_Val = new NUMERIC_TYPE[temp_numPS + verlink->count_dem_index]();

		memcpy(temp_xpi, BCptr->xpi, sizeof(int) * temp_numPS);
		memcpy(temp_ypi, BCptr->ypi, sizeof(int) * temp_numPS);
		memcpy(temp_PS_Ident, BCptr->PS_Ident, sizeof(ESourceType) * temp_numPS);
		memcpy(temp_PS_Val, BCptr->PS_Val, sizeof(NUMERIC_TYPE) * temp_numPS);

		for (size_t tempindex = 0; tempindex < verlink->count_dem_index; tempindex++)
		{
			temp_xpi[temp_numPS + tempindex] = verlink->dem_xi[tempindex];
			temp_ypi[temp_numPS + tempindex] = verlink->dem_yi[tempindex];
			temp_PS_Ident[temp_numPS + tempindex] = QFIX4;
			temp_PS_Val[temp_numPS + tempindex] = verlink->flow_2D_dem[verlink->dem_Vindex_uni[tempindex]];
		}

		delete[] BCptr->xpi;
		delete[] BCptr->ypi;
		delete[] BCptr->PS_Ident;
		memory_free_legacy(&BCptr->PS_Val);

		BCptr->xpi = temp_xpi;
		BCptr->ypi = temp_ypi;
		BCptr->PS_Ident = temp_PS_Ident;
		BCptr->PS_Val = temp_PS_Val;
		BCptr->numPS = temp_numPS + verlink->count_dem_index;
	}
	else
	{
		int* temp_xpi = new int[verlink->count_dem_index](); //used in legacy and read in
		int* temp_ypi = new int[verlink->count_dem_index](); //used in legacy and read in
		ESourceType* temp_PS_Ident = new ESourceType[verlink->count_dem_index];
		NUMERIC_TYPE* temp_PS_Val = new NUMERIC_TYPE[verlink->count_dem_index]();

		for (size_t tempindex = 0; tempindex < verlink->count_dem_index; tempindex++)
		{
			temp_xpi[tempindex] = verlink->dem_xi[tempindex];
			temp_ypi[tempindex] = verlink->dem_yi[tempindex];
			temp_PS_Ident[tempindex] = QFIX4;
			temp_PS_Val[tempindex] = verlink->flow_2D_dem[verlink->dem_Vindex_uni[tempindex]];
		}

		delete[] BCptr->xpi;
		delete[] BCptr->ypi;
		delete[] BCptr->PS_Ident;
		memory_free_legacy(&BCptr->PS_Val);
		delete[] BCptr->PS_Name;

		BCptr->PS_Ident = temp_PS_Ident;
		BCptr->xpi = temp_xpi;
		BCptr->ypi = temp_ypi;
		BCptr->PS_Val = temp_PS_Val;
		BCptr->numPS = verlink->count_dem_index;
	}
	}

	if (swmm_coupled)
	{
		Node_Inflow = (double*)calloc(Nobjects[NODE], sizeof(double));
	}
	//--------------------------------------------------------------------

	time_t loop_start;
	time(&loop_start);

	while (Solverptr->t < Solverptr->Sim_Time
			&& (!swmm_coupled
					|| (NewRoutingTime < TotalDuration && !ErrorCode)))
	{	
		// Freeze one common coupling step for exchange limiting, FV1 and SWMM.
		const NUMERIC_TYPE remaining_2d = Solverptr->Sim_Time - Solverptr->t;
		Solverptr->Tstep = getmin(Solverptr->Tstep, remaining_2d);
		if (swmm_coupled)
		{
			const NUMERIC_TYPE remaining_swmm = static_cast<NUMERIC_TYPE>(
					(TotalDuration - NewRoutingTime) / 1000.0);
			Solverptr->Tstep = getmin(Solverptr->Tstep, remaining_swmm);
		}
		if (Solverptr->Tstep <= C(0.0)) break;

		if (swmm_coupled)
		{
		memset(verlink->flow_neg_count, 0, Parptr->xsz * Parptr->ysz * sizeof(int));
		memset(verlink->flow_2D_dem, 0, Parptr->xsz * Parptr->ysz * sizeof(double));
		memset(verlink->flow_1D_node, 0, verlink->count_index * sizeof(double));
		memset(verlink->flow_2D_node, 0, verlink->count_index * sizeof(double));
		memset(Node_Inflow, 0, Nobjects[NODE] * sizeof(double));

		floodfv1::verticalcon::cal_exchange_flow(verlink, Arrptr, Parptr, Solverptr);
		floodfv1::verticalcon::nodeinflow_cal(verlink);

		if (temp_numPS > 0)
		{
			for (size_t flowindex = 0; flowindex < verlink->count_dem_index; flowindex++)
			{
				BCptr->PS_Val[temp_numPS + flowindex] = verlink->flow_2D_dem[verlink->dem_Vindex_uni[flowindex]];
			}
		}
		else
		{
			for (size_t flowindex = 0; flowindex < verlink->count_dem_index; flowindex++)
			{
				BCptr->PS_Val[flowindex] = verlink->flow_2D_dem[verlink->dem_Vindex_uni[flowindex]];
			}
		}
		}

		zero_flux_stats(BCptr);

        rain.update_H(Parptr, Solverptr, Arrptr);
		update_point_sources(Parptr, Solverptr, BCptr, Arrptr, temp_numPS);
		apply_friction(Parptr, Solverptr, Arrptr);
		update_Hstar(Parptr, Arrptr);
		update_fluxes(Parptr, Solverptr, Arrptr);
		update_fluxes_on_boundaries(Parptr, Solverptr, BCptr, Arrptr);
		update_flow_variables(Parptr, Solverptr, Arrptr);
        fv1::drain_nodata_water(Parptr, Solverptr, BCptr, Arrptr);

		accumulate_boundary_flux_stats(Parptr, Solverptr, BCptr, Arrptr);

		Solverptr->MinTstep = getmin(Solverptr->MinTstep, Solverptr->Tstep);
		int swmm_error = 0;
		if (swmm_coupled)
		{
			linkDt = Solverptr->Tstep;
			swmm_error = swmm_step(&elapsedTime);
		}

		Solverptr->t += Solverptr->Tstep;
		Solverptr->itCount += 1;
		if (swmm_error != 0) break;

		if (Statesptr->adaptive_ts == ON
				&& Solverptr->t < Solverptr->Sim_Time
				&& (!swmm_coupled || NewRoutingTime < TotalDuration))
		{
			Solverptr->Tstep = Tstep_from_cfl(Parptr, Solverptr, Arrptr);
		}

        update_max_field(Parptr, Arrptr);

		if (verbose == ON) printf("t=%f\tdt=%f\n", Solverptr->t, Solverptr->Tstep);

		if (Solverptr->t >= Parptr->MassTotal) {
			update_mass_stats(Statesptr, Parptr, Solverptr, BCptr, Arrptr);
			write_mass_stats(Fptr, Parptr, Solverptr, BCptr);

            if (Statesptr->save_stages == ON)
            {
                write_depth_samples(Fptr, Parptr, Solverptr, Stageptr, Arrptr);
                if (Statesptr->voutput_stage == ON)
                {
                    write_speed_samples(Fptr, Parptr, Solverptr, Stageptr, Arrptr);
                }
            }
		}

		if (Solverptr->t >= Parptr->SaveTotal) {
			if (Statesptr->voutput == ON)
			{
				update_velocity(Parptr, Solverptr, Arrptr);
			}
			write_solution(Fnameptr, Statesptr, Parptr, Solverptr, Arrptr);
		}
	}

	time_t loop_end;
	time(&loop_end);

	double seconds = difftime(loop_end, loop_start);
	printf("loop time %lf\n", seconds);

	if (swmm_coupled)
	{
		swmm_end();
		if (Fout.mode == SCRATCH_FILE) swmm_report();

		// --- close the system
		swmm_close();
	}

    write_max_field(Fnameptr, Statesptr, Parptr, Arrptr);

	if (swmm_coupled) floodfv1::verticalcon::free_vcon(verlink);
	if (verbose == ON) printf("Finished.\n\n");
	return;
}

void fv1::apply_friction
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
	if (Arrptr->Manningsn == nullptr && Parptr->FPn <= C(0.0)) return;

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz + i];
			NUMERIC_TYPE& HU = Arrptr->HU[j*Parptr->xsz + i];
			NUMERIC_TYPE& HV = Arrptr->HV[j*Parptr->xsz + i];

			if (H <= Solverptr->DepthThresh) {
				HU = C(0.0);
				HV = C(0.0);
				continue;
			}

			NUMERIC_TYPE U = HU/H;
			NUMERIC_TYPE V = HV/H;
			if (FABS(U) <= Solverptr->SpeedThresh
					&& FABS(V) <= Solverptr->SpeedThresh)
			{
				HU = C(0.0);
				HV = C(0.0);
				continue;
			}

            NUMERIC_TYPE n = (Arrptr->Manningsn == nullptr)
                ? Parptr->FPn : Arrptr->Manningsn[j*Parptr->xsz + i];

			NUMERIC_TYPE Cf = Solverptr->g*n*n / pow(H, C(1.0)/C(3.0));
			NUMERIC_TYPE speed = SQRT(U*U+V*V);

			NUMERIC_TYPE Sf_x = -Cf*U*speed;
			NUMERIC_TYPE Sf_y = -Cf*V*speed;
			NUMERIC_TYPE D_x = C(1.0) + Solverptr->Tstep*Cf / H
				* (C(2.0)*U*U+V*V)/speed;
			NUMERIC_TYPE D_y = C(1.0) + Solverptr->Tstep*Cf / H
				* (U*U+C(2.0)*V*V)/speed;

			HU += Solverptr->Tstep*Sf_x/D_x;
			HV += Solverptr->Tstep*Sf_y/D_y;
		}
	}
}

void fv1::drain_nodata_water
(
    Pars* Parptr,
    Solver *Solverptr,
    BoundCs *BCptr,
    Arrays *Arrptr
)
{
    if (Parptr->drain_nodata == OFF) return;

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];

            if (FABS(Z - Parptr->nodata_elevation) < 1e-6)
            {
			    NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];
                NUMERIC_TYPE Hold = H;
                H = C(0.0);

				accumulate_point_flux_stats(BCptr,
                        -Hold*Parptr->dA/Solverptr->Tstep);
            }
        }
    }
}

void fv1::update_point_sources
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr,
	const int coupling_start_index
)
{
	for (int ps_i = 0; ps_i < BCptr->numPS; ps_i++)
	{
		const int i = BCptr->xpi[ps_i];
		const int j = BCptr->ypi[ps_i];
		
		NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];

		switch (BCptr->PS_Ident[ps_i])
		{
		case HFIX2:
			{
				NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
				NUMERIC_TYPE H_new = FMIN(C(0.0), BCptr->PS_Val[ps_i] - Z);
				NUMERIC_TYPE Q = (H_new - H)*Parptr->dA / Solverptr->Tstep;
				H = H_new;	

				accumulate_point_flux_stats(BCptr, Q);
			}
			break;
		case HVAR3:
			{
				NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
				NUMERIC_TYPE H_new = FMIN(C(0.0), linear_interpolate(
						BCptr->PS_TimeSeries[ps_i], Solverptr->t) - Z);
				NUMERIC_TYPE Q = (H_new - H)*Parptr->dA / Solverptr->Tstep;
				H = H_new;	
				
				accumulate_point_flux_stats(BCptr, Q);
			}
			break;
		case QFIX4:
			{
				NUMERIC_TYPE Q = BCptr->PS_Val[ps_i];
				NUMERIC_TYPE H_old = H;
				NUMERIC_TYPE H_new = H_old
					+ Q * Parptr->dx * Solverptr->Tstep / Parptr->dA;

				// Coupled surface-to-sewer flow removes the horizontal momentum
				// carried by the drained water. Sewer-to-surface overflow is
				// assumed to enter vertically with zero horizontal momentum.
				if (ps_i >= coupling_start_index && Q < C(0.0))
				{
					const int cell_i = j*Parptr->xsz + i;
					H_new = FMAX(C(0.0), H_new);

					if (H_old > C(0.0))
					{
						NUMERIC_TYPE remaining_depth_fraction = H_new / H_old;
						Arrptr->HU[cell_i] *= remaining_depth_fraction;
						Arrptr->HV[cell_i] *= remaining_depth_fraction;
					}
					else
					{
						Arrptr->HU[cell_i] = C(0.0);
						Arrptr->HV[cell_i] = C(0.0);
					}
				}

				H = H_new;

				accumulate_point_flux_stats(BCptr, Q*Parptr->dx);
			}
			break;
		case QVAR5:
			{
				NUMERIC_TYPE Q = linear_interpolate(
						BCptr->PS_TimeSeries[ps_i], Solverptr->t);
				H += Q * Parptr->dx * Solverptr->Tstep / Parptr->dA;

				accumulate_point_flux_stats(BCptr, Q*Parptr->dx);
			}
			break;
		case NONE0:
		case FREE1:
		default:
			break;
		}
	}
}

void fv1::update_fluxes
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=1; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H_neg = Arrptr->Hstar_neg_x[j*Parptr->xsz + i-1];
			NUMERIC_TYPE HU_neg = HUstar_neg_x(
					Parptr, Solverptr, Arrptr, i-1, j);
			NUMERIC_TYPE HV_neg = HVstar_neg_x(
					Parptr, Solverptr, Arrptr, i-1, j);

			NUMERIC_TYPE H_pos = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
			NUMERIC_TYPE HU_pos = HUstar_pos_x(Parptr, Solverptr, Arrptr, i, j);
			NUMERIC_TYPE HV_pos = HVstar_pos_x(Parptr, Solverptr, Arrptr, i, j);

			NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];

			HLL_x(Solverptr, H_neg, HU_neg, HV_neg, H_pos, HU_pos, HV_pos,
				FHx, FHUx, FHVx);
		}
	}

#pragma omp parallel for
	for (int j=1; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H_neg = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
			NUMERIC_TYPE HU_neg = HUstar_neg_y(Parptr, Solverptr, Arrptr, i, j);
			NUMERIC_TYPE HV_neg = HVstar_neg_y(Parptr, Solverptr, Arrptr, i, j);

			NUMERIC_TYPE H_pos = Arrptr->Hstar_pos_y[(j-1)*Parptr->xsz + i];
			NUMERIC_TYPE HU_pos = HUstar_pos_y(
					Parptr, Solverptr, Arrptr, i, j-1);
			NUMERIC_TYPE HV_pos = HVstar_pos_y(
					Parptr, Solverptr, Arrptr, i, j-1);

			NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

			HLL_y(Solverptr, H_neg, HU_neg, HV_neg, H_pos, HU_pos, HV_pos,
					FHy, FHUy, FHVy);
		}
	}
}

void fv1::update_fluxes_on_boundaries
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
)
{
	// west
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = 0;
		NUMERIC_TYPE H_inside = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = HUstar_pos_x(Parptr, Solverptr, Arrptr, i, j);
		NUMERIC_TYPE HV_inside = HVstar_pos_x(Parptr, Solverptr, Arrptr, i, j);

		NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];
		
		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i, j,
				boundary_index_w(Parptr, i, j), 1,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);

		HLL_x(Solverptr, H_outside, HU_outside, HV_outside,
				H_inside, HU_inside, HV_inside, FHx, FHUx, FHVx);
	}

	// east
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz;
		NUMERIC_TYPE H_inside = Arrptr->Hstar_neg_x[j*Parptr->xsz + i-1];
		NUMERIC_TYPE HU_inside = HUstar_neg_x(
				Parptr, Solverptr, Arrptr, i-1, j);
		NUMERIC_TYPE HV_inside = HVstar_neg_x(
				Parptr, Solverptr, Arrptr, i-1, j);

		NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i-1, j,
				boundary_index_e(Parptr, i, j), -1,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);

		HLL_x(Solverptr, H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside, FHx, FHUx, FHVx);
	}
	
	// north
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		NUMERIC_TYPE H_inside = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = HUstar_neg_y(Parptr, Solverptr, Arrptr, i, j);
		NUMERIC_TYPE HV_inside = HVstar_neg_y(Parptr, Solverptr, Arrptr, i, j);

		NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i, j,
				boundary_index_n(Parptr, i, j), -1,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);

		HLL_y(Solverptr, H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside, FHy, FHUy, FHVy);
	}

	// south
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz;
		NUMERIC_TYPE H_inside = Arrptr->Hstar_pos_y[(j-1)*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = HUstar_pos_y(
				Parptr, Solverptr, Arrptr, i, j-1);
		NUMERIC_TYPE HV_inside = HVstar_pos_y(
				Parptr, Solverptr, Arrptr, i, j-1);

		NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i, j-1,
				boundary_index_s(Parptr, i, j), 1,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);

		HLL_y(Solverptr, H_outside, HU_outside, HV_outside,
				H_inside, HU_inside, HV_inside, FHy, FHUy, FHVy);
	}
}

void fv1::update_flow_variables
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];
			NUMERIC_TYPE H_w = Arrptr->FHx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE H_e = Arrptr->FHx[j*(Parptr->xsz+1) + i+1];
			NUMERIC_TYPE H_n = Arrptr->FHy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE H_s = Arrptr->FHy[(j+1)*(Parptr->xsz+1) + i];

			H = H - Solverptr->Tstep *
				((H_e - H_w)/Parptr->dx + (H_n - H_s)/Parptr->dy);

			NUMERIC_TYPE& HU = Arrptr->HU[j*Parptr->xsz + i];
			NUMERIC_TYPE HU_w = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HU_e = Arrptr->FHUx[j*(Parptr->xsz+1) + i+1];
			NUMERIC_TYPE HU_n = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HU_s = Arrptr->FHUy[(j+1)*(Parptr->xsz+1) + i];

			HU = HU - Solverptr->Tstep *
				((HU_e - HU_w)/Parptr->dx + (HU_n - HU_s)/Parptr->dy
				 - bed_source_x(Parptr, Solverptr, Arrptr, i, j));

			NUMERIC_TYPE& HV = Arrptr->HV[j*Parptr->xsz + i];
			NUMERIC_TYPE HV_w = Arrptr->FHVx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HV_e = Arrptr->FHVx[j*(Parptr->xsz+1) + i+1];
			NUMERIC_TYPE HV_n = Arrptr->FHVy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HV_s = Arrptr->FHVy[(j+1)*(Parptr->xsz+1) + i];

			HV = HV - Solverptr->Tstep *
				((HV_e - HV_w)/Parptr->dx + (HV_n - HV_s)/Parptr->dy
				 - bed_source_y(Parptr, Solverptr, Arrptr, i, j));
		}
	}
}

NUMERIC_TYPE fv1::bed_source_x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg = Arrptr->Hstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
	NUMERIC_TYPE Zdagger_neg = Zdagger_neg_x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger_pos = Zdagger_pos_x(Parptr, Arrptr, i, j);

	return -Solverptr->g * C(0.5)*(Hstar_neg + Hstar_pos)
		* (Zdagger_neg - Zdagger_pos)/Parptr->dx;
}

NUMERIC_TYPE fv1::bed_source_y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos = Arrptr->Hstar_pos_y[j*Parptr->xsz + i];
	NUMERIC_TYPE Zdagger_neg = Zdagger_neg_y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger_pos = Zdagger_pos_y(Parptr, Arrptr, i, j);

	return -Solverptr->g * 0.5*(Hstar_neg + Hstar_pos)
		* (Zdagger_neg - Zdagger_pos)/Parptr->dy;
}

void fv1::set_boundary_values
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr,
	const int z_i,
	const int z_j,
	const int bc_i,
	const int HU_sign,
	NUMERIC_TYPE& H_inside,
	NUMERIC_TYPE& HU_inside,
	NUMERIC_TYPE& HV_inside,
	NUMERIC_TYPE& H_outside,
	NUMERIC_TYPE& HU_outside,
	NUMERIC_TYPE& HV_outside
)
{
	switch (BCptr->BC_Ident[bc_i])
	{
	case FREE1:
		H_outside = H_inside;
		HU_outside = HU_inside;
		HV_outside = HV_outside;
		break;
	case HFIX2:
		{
			NUMERIC_TYPE Z = Arrptr->DEM[z_j*Parptr->xsz + z_i];
			H_outside = FMAX(C(0.0), BCptr->BC_Val[bc_i] - Z);
			//H_inside = H_outside;
			HU_outside = HU_inside;
			HV_outside = HV_inside;
		}
		break;
	case HVAR3:
		{
			NUMERIC_TYPE Z = Arrptr->DEM[z_j*Parptr->xsz + z_i];
			H_outside = FMAX(C(0.0), linear_interpolate(
					BCptr->BC_TimeSeries[bc_i], Solverptr->t) - Z);
			//H_inside = H_outside;
			HU_outside = HU_inside;
			HV_outside = HV_inside;
		}
		break;
	case QFIX4:
		HU_outside = HU_sign * BCptr->BC_Val[bc_i];
		if (FABS(HU_outside) > C(1e-12))
		{
			H_inside = H_outside = FMAX(H_inside, C(1.1)*Solverptr->DepthThresh);
		}
        else
        {
		    H_outside = H_inside;
        }
		HU_inside = HU_outside;
		HV_outside = HV_inside;
		break;
	case QVAR5:
		HU_outside = HU_sign * linear_interpolate(
				BCptr->BC_TimeSeries[bc_i], Solverptr->t);
		if (FABS(HU_outside) > C(1e-12))
		{
			H_inside = H_outside = FMAX(H_inside, C(1.1)*Solverptr->DepthThresh);
		}
        else
        {
		    H_outside = H_inside;
        }
		HU_inside = HU_outside;
		HV_outside = HV_inside;
		break;
	case NONE0:
	default:
		H_outside = H_inside;
		HU_outside = -HU_inside;
		HV_outside = HV_outside;
		break;
	}
}

NUMERIC_TYPE fv1::Tstep_from_cfl
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
	NUMERIC_TYPE dt = Solverptr->InitTstep;
#ifndef _MSC_VER
#pragma omp parallel for reduction(min:dt)
#endif
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz + i];
			if (H > Solverptr->DepthThresh)
			{
				NUMERIC_TYPE HU = Arrptr->HU[j*Parptr->xsz + i];
				NUMERIC_TYPE HV = Arrptr->HV[j*Parptr->xsz + i];
				NUMERIC_TYPE U = HU/H;
				NUMERIC_TYPE V = HV/H;
				
				NUMERIC_TYPE dt_x = Solverptr->cfl*Parptr->dx
					/ (FABS(U)+SQRT(Solverptr->g*H));
				NUMERIC_TYPE dt_y = Solverptr->cfl*Parptr->dy
					/ (FABS(V)+SQRT(Solverptr->g*H));

				dt = (std::min)({dt, dt_x, dt_y});
			}
		}
	}
	return dt;
}
