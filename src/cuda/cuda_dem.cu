#include <unordered_map>
#include <vector>
#include "cuda_dem.cuh"
#include "cuda_coupling_gather.cuh"
#include "cuda_solver.cuh"
#include "io.h"
#include "../swmm5.h"
extern "C" {
#include "../headers.h"
};
#include <float.h>

namespace lis
{
namespace cuda
{

__global__ void gather_coupling_depths_kernel
(
	const NUMERIC_TYPE* H,
	const int* coupling_indices,
	NUMERIC_TYPE* coupling_depths,
	int coupling_count
)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < coupling_count)
	{
		coupling_depths[i] = H[coupling_indices[i]];
	}
}

namespace fv1
{

__global__ void initialise_Zstar_x
(
	NUMERIC_TYPE* __restrict__ Zstar_x,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j+1; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z_neg = DEM[j*cuda::pitch + i];
			NUMERIC_TYPE Z_pos = DEM[j*cuda::pitch + i+1];

			Zstar_x[j*cuda::pitch + i] = FMAX(Z_neg, Z_pos);
		}
	}
}

__global__ void initialise_Zstar_y
(
	NUMERIC_TYPE* __restrict__ Zstar_y,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i+1; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z_neg = DEM[(j+1)*cuda::pitch + i];
			NUMERIC_TYPE Z_pos = DEM[j*cuda::pitch + i];

			Zstar_y[j*cuda::pitch + i] = FMAX(Z_neg, Z_pos);
		}
	}
}

}
}
}

void lis::cuda::gather_coupling_depths
(
	const NUMERIC_TYPE* H,
	const int* coupling_indices,
	NUMERIC_TYPE* coupling_depths,
	int coupling_count
)
{
	if (coupling_count <= 0) return;
	const int blocks = (coupling_count + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;
	lis::cuda::gather_coupling_depths_kernel<<<blocks, CUDA_BLOCK_SIZE>>>(
			H, coupling_indices, coupling_depths, coupling_count);
}

NUMERIC_TYPE* lis::cuda::Topography::load
(
	const char* filename,
	Geometry& geometry,
	int& pitch,
	int& offset,
	NUMERIC_TYPE nodata_elevation,
	int acceleration,
	int verbose
)
{
	FILE* dem_file = fopen_or_die(filename, "rb", "Loading DEM\n", verbose);
	NUMERIC_TYPE original_no_data_value = C(0.0);
	AsciiRaster::read_header(dem_file, geometry, original_no_data_value);
	NUMERIC_TYPE* DEM; 

    if (acceleration == 0){ 
		DEM = lis::GhostRaster::allocate(geometry); 
	    pitch = lis::GhostRaster::pitch(geometry);
	    offset = lis::GhostRaster::offset(geometry);
	} else {                                            
		DEM = lis::GhostRaster::allocate_H(geometry); 
	    pitch = lis::GhostRaster::pitch_ACC(geometry);   
	    offset = lis::GhostRaster::offset_ACC(geometry);	
	}
	AsciiRaster::read(dem_file, DEM, geometry, pitch, offset); 
	AsciiRaster::replace_no_data(DEM, geometry, pitch, offset,
			original_no_data_value, nodata_elevation);

	fclose(dem_file);


	return DEM;
}

void lis::cuda::Topography::initialise_Zstar_x
(
	NUMERIC_TYPE* __restrict__ Zstar_x,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	lis::cuda::fv1::initialise_Zstar_x<<<1, cuda::block_size>>>(Zstar_x, DEM);
}

void lis::cuda::Topography::initialise_Zstar_y
(
	NUMERIC_TYPE* __restrict__ Zstar_y,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	lis::cuda::fv1::initialise_Zstar_y<<<1, cuda::block_size>>>(Zstar_y, DEM);
}

void lis::cuda::Topography::clamp_boundary_values
(
	NUMERIC_TYPE* DEM,
	Geometry& geometry,
	int pitch
)
{
	for (int j=1; j<geometry.ysz+1; j++)
	{
		// west
		{
			const int i=0;
			DEM[j*pitch + i] = DEM[j*pitch + i+1];
		}

		// east
		{
			const int i=geometry.xsz+1;
			//if (j == 161 || j == 162 || j == 163 || j == 164) {
			//	DEM[j * pitch + i] = -C(20.5);
			//}
			//else
			//{
				DEM[j * pitch + i] = DEM[j * pitch + i - 1];
//			}
		}
	}

	for (int i=1; i<geometry.xsz+1; i++)
	{
		// north
		{
			const int j=0;
			DEM[j*pitch + i] = DEM[(j+1)*pitch + i];
		}

		// south
		{
			const int j=geometry.ysz+1;
			DEM[j*pitch + i] = DEM[(j-1)*pitch + i];
		}
	}
}

void lis::cuda::Topography::load_floodconfig
(
	Fnames& filenames,
	Vconnect* verlink,
	Geometry& geometry,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE nodata_elevation,
	int verbose
)
{
	ADDLINK = 1;

	FILE* flood_file = fopen_or_die(filenames.floodfilename, "rb", "Loading FloodConfig\n", verbose);
	FILE* flood_error = fopen_or_die(filenames.flooderrorfile, "w", "Loading FloodError\n", verbose);
	fscanf(flood_file, "%*s %d", &(verlink->count_index));
	flood_init(verlink, geometry);
	fscanf(flood_file, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
	read_floodconfig(flood_file, verlink);
	fclose(flood_file);

	cal_demindex(flood_error, verlink, geometry, DEM, nodata_elevation);
	fix_maxdepth(flood_error, verlink, DEM);
	cal_demcount(verlink);
	map_internal_outfalls(flood_error, verlink, geometry, DEM,
			nodata_elevation, 0);
	fclose(flood_error);

	return;
}


void lis::cuda::Topography::flood_init
(
	Vconnect* verlink,
	Geometry& geometry
)
{

	verlink->node_Vindex = new int[verlink->count_index]();                      
	verlink->dem_Vindex = new int[verlink->count_index]();                       
	verlink->dem_count = new int[verlink->count_index]();
	verlink->node_dem_slot = new int[verlink->count_index]();      
	verlink->node_surfelev = new double[verlink->count_index]();                 
	verlink->node_c0 = new double[verlink->count_index]();                       
	verlink->node_area = new double[verlink->count_index]();					 
	verlink->node_qoutmax = new double[verlink->count_index]();                  
	verlink->node_cw = new double[verlink->count_index]();                       
	verlink->node_w = new double[verlink->count_index]();                        
	verlink->node_qinmax = new double[verlink->count_index]();                   
	//verlink->flow_type = new int[verlink->count_index]();                      
	verlink->node_c1 = new double[verlink->count_index]();                       
	verlink->node_area1 = new double[verlink->count_index]();                    
																				 
	verlink->count_dem_index = verlink->count_index;                             
	verlink->dem_Vindex_uni = new int[verlink->count_index]();                   
																				 
	verlink->dem_xi = new int[verlink->count_index]();                           
	verlink->dem_yi = new int[verlink->count_index]();                           
	//----------------------------------------------------------------------     

	verlink->flow_1D_node = new double[verlink->count_index]();                  
	verlink->flow_2D_node = new double[verlink->count_index]();                  
	verlink->flow_2D_dem = new double[verlink->count_index](); 
	verlink->flow_neg_count = new int[verlink->count_index](); 
	//----------------------------------------------------------------------
	return;
}

void lis::cuda::Topography::read_floodconfig
(
	FILE* flood_file,
	Vconnect* verlink
)
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
	return;
}

void lis::cuda::Topography::cal_demindex
(
	FILE* flood_error,
	Vconnect* verlink,
	Geometry& geometry,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE nodata_elevation
)
{
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "THE CORRESPONDING INDEX CALCULATION START:\n");
	fprintf(flood_error, "====================================================================\n");

	int x_index, y_index, dem_index;
	double x_min = geometry.blx;
	double x_max = geometry.blx + geometry.xsz * geometry.dx;
	double y_min = geometry.bly;
	double y_max = geometry.bly + geometry.ysz * geometry.dy;

	for (int j = 0; j < verlink->count_index; ++j)
	{
		verlink->node_surfelev[j] = (Node[verlink->node_Vindex[j]].invertElev + Node[verlink->node_Vindex[j]].fullDepth) * UCF(LENGTH);

		if ((x_min - NodeXY[verlink->node_Vindex[j]][0]) > 1e-6 || (NodeXY[verlink->node_Vindex[j]][0] - x_max) > 1e-6 ||
			(y_min - NodeXY[verlink->node_Vindex[j]][1]) > 1e-6 || (NodeXY[verlink->node_Vindex[j]][1] - y_max) > 1e-6)
		{
			verlink->dem_Vindex[j] = -9999;
			verlink->dem_xi[j] = -9999;
			verlink->dem_yi[j] = -9999;
			fprintf(flood_error, "node:%-15s index %d is out of DEM range\n", Node[verlink->node_Vindex[j]].ID, verlink->node_Vindex[j]);
		}
		else if (Node[verlink->node_Vindex[j]].type == OUTFALL || Node[verlink->node_Vindex[j]].type == STORAGE || Node[verlink->node_Vindex[j]].type == DIVIDER)
		{
			verlink->dem_Vindex[j] = -9999;
			verlink->dem_xi[j] = -9999;
			verlink->dem_yi[j] = -9999;
			fprintf(flood_error, "node:%-15s index %d is out of DEM range\n", Node[verlink->node_Vindex[j]].ID, verlink->node_Vindex[j]);
		}
		else
		{
			x_index = (int)((NodeXY[verlink->node_Vindex[j]][0] - x_min) / geometry.dx);
			y_index = (int)((y_max - NodeXY[verlink->node_Vindex[j]][1]) / geometry.dy);

			if (x_index >= geometry.xsz) x_index -= 1;
			if (y_index >= geometry.ysz) y_index -= 1;

			dem_index = x_index + y_index * geometry.xsz;

			if (FABS(DEM[dem_index] - nodata_elevation) < 1e-6)
			{
				verlink->dem_Vindex[j] = -9999;
				verlink->dem_xi[j] = -9999;
				verlink->dem_yi[j] = -9999;
				fprintf(flood_error, "node:%-15s index %d is out of DEM range\n", Node[verlink->node_Vindex[j]].ID, verlink->node_Vindex[j]);
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

void lis::cuda::Topography::cal_demcount
(
	Vconnect* verlink
)
{
	std::unordered_map<int, int> slots;
	verlink->count_dem_index = 0;
	for (int i = 0; i < verlink->count_index; ++i)
	{
		const int grid = verlink->dem_Vindex[i];
		if (grid == -9999) { verlink->node_dem_slot[i] = -1; continue; }
		auto found = slots.find(grid);
		int slot;
		if (found == slots.end())
		{
			slot = verlink->count_dem_index++;
			slots.emplace(grid, slot);
			verlink->dem_Vindex_uni[slot] = grid;
			verlink->dem_xi[slot] = verlink->dem_xi[i];
			verlink->dem_yi[slot] = verlink->dem_yi[i];
		}
		else slot = found->second;
		verlink->node_dem_slot[i] = slot;
		++verlink->dem_count[slot];
	}
}

void lis::cuda::Topography::fix_maxdepth
(
	FILE* flood_error,
	Vconnect* verlink,
	NUMERIC_TYPE* DEM
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
			if (FABS(DEM[verlink->dem_Vindex[i]] - verlink->node_surfelev[i]) > SURFLEVELDIFF)
			{
				fprintf(flood_error, "The difference between node: %-15s index %d invertElev:%lf + fullDepth:%lf = %lf and DEM:%lf is over %lf\n",
					Node[verlink->node_Vindex[i]].ID,
					verlink->node_Vindex[i],
					Node[verlink->node_Vindex[i]].invertElev * UCF(LENGTH),
					Node[verlink->node_Vindex[i]].fullDepth * UCF(LENGTH),
					(Node[verlink->node_Vindex[i]].invertElev + Node[verlink->node_Vindex[i]].fullDepth) * UCF(LENGTH),
					DEM[verlink->dem_Vindex[i]],
					SURFLEVELDIFF);
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

void lis::cuda::Topography::cal_flow_exchange
(
	Vconnect* verlink,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H, 
	NUMERIC_TYPE dt,
	Geometry& geometry,
	const double* pending_volumes
)
{
	double H1D, H2D, Z; 
	double flow_exchange;


//#pragma omp parallel for private(i, flow_exchange, H1D, H2D, Z)
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
		H2D = DEM[verlink->dem_Vindex[i]] + H[verlink->node_dem_slot[i]];
		Z = getmax(verlink->node_surfelev[i], DEM[verlink->dem_Vindex[i]]);

		if ((H1D - Z > COUPLING_HEAD_TOL) && (H1D - H2D > COUPLING_HEAD_TOL))
		{
			flow_exchange = cal_orificeflow(verlink->node_c0[i],
												verlink->node_area[i],
												H1D,
												H2D);
			flow_exchange = getmin(flow_exchange, verlink->node_qoutmax[i]);
			const double stored = swmm_getNodeAvailableCoupledVolume(verlink->node_Vindex[i]);
			if (stored >= 0.0)
			{
				const double pending = pending_volumes ? pending_volumes[i] : 0.0;
				flow_exchange = getmin(flow_exchange, getmax(0.0, stored+pending)/dt);
			}
			verlink->flow_1D_node[i] = -1 * flow_exchange;
			verlink->flow_2D_node[i] = flow_exchange;
		}
		else if ((H1D <= Z + COUPLING_HEAD_TOL) && (H2D - Z > COUPLING_HEAD_TOL))
		{
				flow_exchange = cal_weirflow(verlink->node_cw[i],
											 verlink->node_w[i],
											(H2D - Z));
				verlink->flow_1D_node[i] = flow_exchange;
				verlink->flow_2D_node[i] = flow_exchange * -1;
		}
		else if ((H1D >= Z - COUPLING_HEAD_TOL) && (H2D - H1D > COUPLING_HEAD_TOL))
		{
				flow_exchange = cal_orificeflow(verlink->node_c1[i],
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

//#pragma omp parallel for private(k)
	for (size_t k = 0; k < verlink->count_index; k++)
	{
		if (verlink->dem_Vindex[k] == -9999)
		{
			continue;
		}
		if (verlink->flow_2D_node[k] < -DBL_EPSILON)
		{
			verlink->flow_neg_count[verlink->node_dem_slot[k]] += 1;
		}
	}
//#pragma omp parallel for private(j, flow_exchange)
	for (size_t j = 0; j < verlink->count_index; j++)
	{
		if (verlink->dem_Vindex[j] == -9999) continue;
		if (verlink->flow_2D_node[j] > DBL_EPSILON)
		{
			verlink->flow_2D_dem[verlink->node_dem_slot[j]] += verlink->flow_2D_node[j] / geometry.dy;
		}
		else if (verlink->flow_2D_node[j] < -DBL_EPSILON)
		{
			if (verlink->flow_neg_count[verlink->node_dem_slot[j]] == 0)
			{
				printf("wrong-------------------------------------------------------------\n");
				printf("lf", verlink->flow_2D_node[j]);
				continue;
			}
			else
			{
				if (H[verlink->node_dem_slot[j]] < -DBL_EPSILON)
				{
					printf("wrong dem+++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
					continue;
				}
				else
				{
					flow_exchange = getmin(((geometry.dx * geometry.dy * H[verlink->node_dem_slot[j]]) / ((verlink->flow_neg_count[verlink->node_dem_slot[j]]) * dt)),
						verlink->flow_1D_node[j]);
					flow_exchange = getmin(flow_exchange, verlink->node_qinmax[j]);
					verlink->flow_1D_node[j] = flow_exchange;
					verlink->flow_2D_node[j] = flow_exchange * -1;

					verlink->flow_2D_dem[verlink->node_dem_slot[j]] += verlink->flow_2D_node[j] / geometry.dy;
				}
			}
		}
		else
		{
			verlink->flow_2D_dem[verlink->node_dem_slot[j]] += 0;
		}
	}
	return;
}


NUMERIC_TYPE lis::cuda::Topography::cal_orificeflow
(
	NUMERIC_TYPE c0,
	NUMERIC_TYPE area,
	NUMERIC_TYPE H1D,
	NUMERIC_TYPE H2D
)
{
	return c0 * area * (sqrt(2.0 * SI_GRAVITY * (H1D - H2D)));
}

NUMERIC_TYPE lis::cuda::Topography::cal_weirflow
(
	NUMERIC_TYPE cw,
	NUMERIC_TYPE w,
	NUMERIC_TYPE h2d
)
{
	return cw * w * h2d * (sqrt(2.0 * SI_GRAVITY * h2d));
}

void lis::cuda::Topography::free_vconnect
(
	Vconnect* verlink
)
{
	delete[] verlink->node_Vindex;
	delete[] verlink->dem_Vindex;
	delete[] verlink->dem_count;
	delete[] verlink->node_dem_slot;
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
	delete[] verlink->outfall_node_Vindex;
	delete[] verlink->outfall_dem_Vindex;
	delete[] verlink->outfall_dem_Vindex_uni;
	delete[] verlink->outfall_dynamic_stage;
	return;
}


void lis::cuda::Topography::cal_nodeinflow
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

void lis::cuda::Topography::load_floodconfig_fv1
(
	Fnames& filenames,
	Vconnect* verlink,
	Geometry& geometry,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE nodata_elevation,
	int verbose
)
{
	ADDLINK = 1;

	FILE* flood_file = fopen_or_die(filenames.floodfilename, "rb", "Loading FloodConfig\n", verbose);
	FILE* flood_error = fopen_or_die(filenames.flooderrorfile, "w", "Loading FloodError\n", verbose);
	fscanf(flood_file, "%*s %d", &(verlink->count_index));
	flood_init_fv1(verlink, geometry);
	fscanf(flood_file, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
	read_floodconfig(flood_file, verlink);
	fclose(flood_file);

	cal_demindex_fv1(flood_error, verlink, geometry, DEM, nodata_elevation);
	fix_maxdepth(flood_error, verlink, DEM);
	cal_demcount(verlink);
	map_internal_outfalls(flood_error, verlink, geometry, DEM,
			nodata_elevation, 1);
	fclose(flood_error);

	return;
}

void lis::cuda::Topography::flood_init_fv1
(
	Vconnect* verlink,
	Geometry& geometry
)
{

	verlink->node_Vindex = new int[verlink->count_index]();
	verlink->dem_Vindex = new int[verlink->count_index]();
	verlink->dem_count = new int[verlink->count_index]();
	verlink->node_dem_slot = new int[verlink->count_index]();
	verlink->node_surfelev = new double[verlink->count_index]();
	verlink->node_c0 = new double[verlink->count_index]();
	verlink->node_area = new double[verlink->count_index]();
	verlink->node_qoutmax = new double[verlink->count_index]();
	verlink->node_cw = new double[verlink->count_index]();
	verlink->node_w = new double[verlink->count_index]();
	verlink->node_qinmax = new double[verlink->count_index]();
	//verlink->flow_type = new int[verlink->count_index]();                      
	verlink->node_c1 = new double[verlink->count_index]();
	verlink->node_area1 = new double[verlink->count_index]();

	verlink->count_dem_index = verlink->count_index;
	verlink->dem_Vindex_uni = new int[verlink->count_index]();

	verlink->dem_xi = new int[verlink->count_index]();
	verlink->dem_yi = new int[verlink->count_index]();
	//----------------------------------------------------------------------     

	verlink->flow_1D_node = new double[verlink->count_index]();
	verlink->flow_2D_node = new double[verlink->count_index]();
	verlink->flow_2D_dem = new double[verlink->count_index]();
	verlink->flow_neg_count = new int[verlink->count_index]();
	//----------------------------------------------------------------------
	return;
}


void lis::cuda::Topography::cal_demindex_fv1
(
	FILE* flood_error,
	Vconnect* verlink,
	Geometry& geometry,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE nodata_elevation
)
{
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "THE CORRESPONDING INDEX CALCULATION START:\n");
	fprintf(flood_error, "====================================================================\n");

	int x_index, y_index, dem_index;
	double x_min = geometry.blx;
	double x_max = geometry.blx + geometry.xsz * geometry.dx;
	double y_min = geometry.bly;
	double y_max = geometry.bly + geometry.ysz * geometry.dy;

	for (int j = 0; j < verlink->count_index; ++j)
	{
		verlink->node_surfelev[j] = (Node[verlink->node_Vindex[j]].invertElev + Node[verlink->node_Vindex[j]].fullDepth) * UCF(LENGTH);

		if ((x_min - NodeXY[verlink->node_Vindex[j]][0]) > 1e-6 || (NodeXY[verlink->node_Vindex[j]][0] - x_max) > 1e-6 ||
			(y_min - NodeXY[verlink->node_Vindex[j]][1]) > 1e-6 || (NodeXY[verlink->node_Vindex[j]][1] - y_max) > 1e-6)
		{
			verlink->dem_Vindex[j] = -9999;
			verlink->dem_xi[j] = -9999;
			verlink->dem_yi[j] = -9999;
			fprintf(flood_error, "node:%-15s index %d is out of DEM range\n", Node[verlink->node_Vindex[j]].ID, verlink->node_Vindex[j]);
		}
		else if (Node[verlink->node_Vindex[j]].type == OUTFALL || Node[verlink->node_Vindex[j]].type == STORAGE || Node[verlink->node_Vindex[j]].type == DIVIDER)
		{
			verlink->dem_Vindex[j] = -9999;
			verlink->dem_xi[j] = -9999;
			verlink->dem_yi[j] = -9999;
			fprintf(flood_error, "node:%-15s index %d is out of DEM range\n", Node[verlink->node_Vindex[j]].ID, verlink->node_Vindex[j]);
		}
		else
		{
			x_index = (int)((NodeXY[verlink->node_Vindex[j]][0] - x_min) / geometry.dx);
			y_index = (int)((y_max - NodeXY[verlink->node_Vindex[j]][1]) / geometry.dy);

			if (x_index >= geometry.xsz) x_index -= 1;
			if (y_index >= geometry.ysz) y_index -= 1;

			dem_index = (x_index + 1) + (y_index + 1) * (geometry.xsz + 2);

			if (FABS(DEM[dem_index] - nodata_elevation) < 1e-6)
			{
				verlink->dem_Vindex[j] = -9999;
				verlink->dem_xi[j] = -9999;
				verlink->dem_yi[j] = -9999;
				fprintf(flood_error, "node:%-15s index %d is out of DEM range\n", Node[verlink->node_Vindex[j]].ID, verlink->node_Vindex[j]);
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

void lis::cuda::Topography::map_internal_outfalls
(
	FILE* flood_error,
	Vconnect* verlink,
	Geometry& geometry,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE nodata_elevation,
	int ghost_cells
)
{
	verlink->outfall_count = 0;
	verlink->outfall_dem_count = 0;
	verlink->outfall_node_Vindex = nullptr;
	verlink->outfall_dem_Vindex = nullptr;
	verlink->outfall_dem_Vindex_uni = nullptr;
	verlink->outfall_dynamic_stage = nullptr;
	verlink->dynamic_outfall_count = 0;

	int total_outfalls = 0;
	for (int node_index = 0; node_index < Nobjects[NODE]; ++node_index)
	{
		if (Node[node_index].type == OUTFALL) ++total_outfalls;
	}

	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "SWMM OUTFALL TO 2D MAPPING START:\n");
	fprintf(flood_error, "====================================================================\n");

	if (total_outfalls == 0)
	{
		fprintf(flood_error, "No SWMM outfalls were defined.\n");
		return;
	}

	verlink->outfall_node_Vindex = new int[total_outfalls]();
	verlink->outfall_dem_Vindex = new int[total_outfalls]();
	verlink->outfall_dem_Vindex_uni = new int[total_outfalls]();
	verlink->outfall_dynamic_stage = new unsigned char[total_outfalls]();

	const double x_min = geometry.blx;
	const double x_max = geometry.blx + geometry.xsz * geometry.dx;
	const double y_min = geometry.bly;
	const double y_max = geometry.bly + geometry.ysz * geometry.dy;

	for (int node_index = 0; node_index < Nobjects[NODE]; ++node_index)
	{
		if (Node[node_index].type != OUTFALL) continue;

		const double node_x = NodeXY[node_index][0];
		const double node_y = NodeXY[node_index][1];
		if ((x_min - node_x) > 1e-6 || (node_x - x_max) > 1e-6 ||
				(y_min - node_y) > 1e-6 || (node_y - y_max) > 1e-6)
		{
			fprintf(flood_error,
					"outfall:%-20s skipped (outside DEM extent)\n",
					Node[node_index].ID);
			continue;
		}

		int x_index = (int)((node_x - x_min) / geometry.dx);
		int y_index = (int)((y_max - node_y) / geometry.dy);
		if (x_index >= geometry.xsz) x_index = geometry.xsz - 1;
		if (y_index >= geometry.ysz) y_index = geometry.ysz - 1;
		if (x_index < 0 || y_index < 0)
		{
			fprintf(flood_error,
					"outfall:%-20s skipped (outside DEM extent)\n",
					Node[node_index].ID);
			continue;
		}

		const int dem_index = ghost_cells
				? (x_index + 1) + (y_index + 1) * (geometry.xsz + 2)
				: x_index + y_index * geometry.xsz;
		if (FABS(DEM[dem_index] - nodata_elevation) < 1e-6)
		{
			fprintf(flood_error,
					"outfall:%-20s skipped (DEM cell is NoData)\n",
					Node[node_index].ID);
			continue;
		}

		const int outfall_index = verlink->outfall_count++;
		verlink->outfall_node_Vindex[outfall_index] = node_index;
		verlink->outfall_dem_Vindex[outfall_index] = dem_index;

		int unique_index = -1;
		for (int i = 0; i < verlink->outfall_dem_count; ++i)
		{
			if (verlink->outfall_dem_Vindex_uni[i] == dem_index)
			{
				unique_index = i;
				break;
			}
		}
		if (unique_index < 0)
		{
			verlink->outfall_dem_Vindex_uni[verlink->outfall_dem_count++] = dem_index;
		}

		fprintf(flood_error,
				"outfall:%-20s mapped to valid DEM cell (%d, %d), index %d\n",
				Node[node_index].ID, x_index, y_index, dem_index);
	}

	fprintf(flood_error,
			"Mapped %d of %d outfalls to %d unique valid 2D cells; "
			"unmapped outfalls remain external sinks.\n",
			verlink->outfall_count, total_outfalls, verlink->outfall_dem_count);
	// The INP outfall type is the selector. All standard SWMM outfall types
	// retain their ordinary boundary conditions. Every COUPLED outfall
	// must map to a valid 2D cell so that its stage can be updated at runtime.
	int selected = 0;
	for (int node_index = 0; node_index < Nobjects[NODE]; ++node_index)
	{
		if (Node[node_index].type != OUTFALL ||
			Outfall[Node[node_index].subIndex].type != COUPLED_OUTFALL) continue;
		if (RouteModel != DW)
		{
			fprintf(stderr, "COUPLED tailwater requires SWMM DYNWAVE routing\n");
			exit(EXIT_FAILURE);
		}
		int matched = -1;
		for (int i = 0; i < verlink->outfall_count; ++i)
		{
			if (verlink->outfall_node_Vindex[i] == node_index) { matched = i; break; }
		}
		if (matched < 0)
		{
			fprintf(stderr, "COUPLED outfall %s has no valid LISFLOOD DEM cell\n",
				Node[node_index].ID);
			exit(EXIT_FAILURE);
		}
		verlink->outfall_dynamic_stage[matched] = 1;
		++selected;
		fprintf(flood_error, "COUPLED tailwater: %s -> DEM index %d; "
			"2D bed %.6f m; SWMM invert %.6f m\n",
			Node[node_index].ID, verlink->outfall_dem_Vindex[matched],
			(double)DEM[verlink->outfall_dem_Vindex[matched]],
			Node[node_index].invertElev * 0.3048);
	}
	verlink->dynamic_outfall_count = selected;
	fprintf(flood_error, "COUPLED tailwater: %d outfalls\n", selected);
	fprintf(flood_error, "====================================================================\n");
	fprintf(flood_error, "SWMM OUTFALL TO 2D MAPPING END:\n");
	fprintf(flood_error, "====================================================================\n");
}

void lis::cuda::Topography::update_outfall_stages
(
	Vconnect* verlink,
	const NUMERIC_TYPE* DEM,
	const NUMERIC_TYPE* outfall_depths,
	double cell_area
)
{
	for (int i = 0; i < verlink->outfall_count; ++i)
	{
		if (!verlink->outfall_dynamic_stage[i]) continue;
		int cell = -1;
		for (int k = 0; k < verlink->outfall_dem_count; ++k)
			if (verlink->outfall_dem_Vindex_uni[k] == verlink->outfall_dem_Vindex[i]) { cell = k; break; }
		if (cell < 0) { fprintf(stderr, "missing mapped tailwater cell\n"); exit(EXIT_FAILURE); }
		const double bed = DEM[verlink->outfall_dem_Vindex[i]];
		const double depth = outfall_depths[cell];
		if (!isfinite(bed) || !isfinite(depth) || depth < -1.0e-5)
		{
			fprintf(stderr, "invalid river stage at outfall %s (bed %.9g, depth %.9g)\n",
				Node[verlink->outfall_node_Vindex[i]].ID, bed, depth);
			exit(EXIT_FAILURE);
		}
		if (swmm_setOutfallReceivingCell(verlink->outfall_node_Vindex[i],
			verlink->outfall_dem_Vindex[i], bed, depth > 0.0 ? depth : 0.0,
			cell_area) != 0)
		{
			fprintf(stderr, "failed to set SWMM tailwater at %s\n",
				Node[verlink->outfall_node_Vindex[i]].ID);
			exit(EXIT_FAILURE);
		}
	}
}

void lis::cuda::Topography::collect_surface_flood_values(Vconnect* verlink,
	NUMERIC_TYPE* values, int source_offset, NUMERIC_TYPE cell_width,
	double exchange_seconds)
{
	std::vector<int> slots(Nobjects[NODE],-1);
	for (int cell=0;cell<verlink->count_dem_index;++cell) values[source_offset+cell]=C(0.0);
	for (int i=0;i<verlink->count_index;++i)
		slots[verlink->node_Vindex[i]]=verlink->node_dem_slot[i];
	// Include nodes omitted from floodconfig: their flooding is still a real
	// withdrawal from SWMM and must never disappear from the coupled system.
	for (int n=0;n<Nobjects[NODE];++n)
	{
		double v=swmm_getNodeFloodExchangeVolume(n);
		if (v <= 0.0) continue;
		if (slots[n] < 0)
		{
			fprintf(stderr,"flooded node %s has no mapped 2D receiving cell; add a valid floodconfig mapping.\n",Node[n].ID);
			exit(EXIT_FAILURE);
		}
		values[source_offset+slots[n]]+=(NUMERIC_TYPE)(v/(exchange_seconds*cell_width));
	}
}

void lis::cuda::Topography::collect_outfall_source_values
(
	Vconnect* verlink,
	NUMERIC_TYPE* point_source_values,
	int source_offset,
	NUMERIC_TYPE cell_width,
	NUMERIC_TYPE cell_area,
	const NUMERIC_TYPE* tailwater_depths,
	double exchange_seconds
)
{
	for (int i = 0; i < verlink->outfall_dem_count; ++i)
	{
		point_source_values[source_offset + i] = C(0.0);
	}

	for (int i = 0; i < verlink->outfall_count; ++i)
	{
		const int node_index = verlink->outfall_node_Vindex[i];
		// SWMM routing flows are always internal ft3/s, regardless of INP units.
		double outfall_flow = swmm_getOutfallExchangeVolume(node_index) / exchange_seconds;
		if (outfall_flow <= 0.0 && !verlink->outfall_dynamic_stage[i]) continue;
		// The integrated interface volume is already SI (m3).

		for (int cell = 0; cell < verlink->outfall_dem_count; ++cell)
		{
			if (verlink->outfall_dem_Vindex_uni[cell] ==
					verlink->outfall_dem_Vindex[i])
			{
				// LISFLOOD QFIX stores discharge per unit cell width (m2/s).
				point_source_values[source_offset + cell] +=
						(NUMERIC_TYPE)(outfall_flow / cell_width);
				break;
			}
		}
	}
	// SWMM has already admitted any reverse flow. Stop rather than removing
	// more water from the 2D river cell than exists and violating mass balance.
	if (tailwater_depths)
	{
		for (int cell = 0; cell < verlink->outfall_dem_count; ++cell)
		{
			const double q = point_source_values[source_offset + cell] * cell_width;
			const double available = (tailwater_depths[cell] > 0 ? tailwater_depths[cell] : 0) * cell_area;
			const double roundoff = (sizeof(NUMERIC_TYPE) == sizeof(float) ?
				8.0 * FLT_EPSILON : 8.0 * DBL_EPSILON) * getmax(1.0, available);
			if (q < 0.0 && -q * exchange_seconds > available + 1.0e-7 + roundoff)
			{
				fprintf(stderr, "reverse outfall flow exceeds 2D cell water volume at DEM index %d; reduce coupleint or check river cell/datum\n",
					verlink->outfall_dem_Vindex_uni[cell]);
				exit(EXIT_FAILURE);
			}
		}
	}
}
