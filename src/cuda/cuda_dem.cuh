#pragma once
#include "../lisflood.h"
#include "../geometry.h"
#include "stats.h"
#include "ghostraster.h"
#define SURFLEVELDIFF 1e-1
#define depthlimit 1e-6

namespace lis
{
namespace cuda
{

struct Topography
{
	static NUMERIC_TYPE* load
	(
		const char* filename,
		Geometry& geometry,
		int& pitch,
		int& offset,
		NUMERIC_TYPE nodata_elevation,
		int acceleration,
		int verbose
	);

	static void initialise_Zstar_x
	(
		NUMERIC_TYPE* __restrict__ Zstar_x,
		const NUMERIC_TYPE* __restrict__ DEM
	);

	static void initialise_Zstar_y
	(
		NUMERIC_TYPE* __restrict__ Zstar_y,
		const NUMERIC_TYPE* __restrict__ DEM
	);

	static void clamp_boundary_values
	(
		NUMERIC_TYPE* DEM,
		Geometry& geometry,
		int pitch
	);
	static void load_floodconfig
	(
		Fnames& filenames,
		Vconnect* verlink,
		Geometry& geometry,
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE nodata_elevation,
		int verbose
	);

	static void flood_init
	(
		Vconnect* verlink,
		Geometry& geometry
	);

	static void read_floodconfig
	(
		FILE* flood_file,
		Vconnect* verlink
	);

	static void cal_demindex
	(
		FILE* flood_error,
		Vconnect* verlink,
		Geometry& geometry,
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE nodata_elevation
	);

	static void cal_demcount
	(
		Vconnect* verlink
	);

	static void map_internal_outfalls
	(
		FILE* flood_error,
		Vconnect* verlink,
		Geometry& geometry,
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE nodata_elevation,
		int ghost_cells
	);

	static void collect_outfall_source_values
	(
		Vconnect* verlink,
		NUMERIC_TYPE* point_source_values,
		int source_offset,
		NUMERIC_TYPE cell_width,
		NUMERIC_TYPE cell_area,
		const NUMERIC_TYPE* tailwater_depths,
		double exchange_seconds
	);

	static void update_outfall_stages
	(
		Vconnect* verlink,
		const NUMERIC_TYPE* DEM,
		const NUMERIC_TYPE* outfall_depths,
		double cell_area
	);
	static void collect_surface_flood_values(Vconnect* verlink,
		NUMERIC_TYPE* values, int source_offset, NUMERIC_TYPE cell_width,
		double exchange_seconds);

	static void fix_maxdepth
	(
		FILE* flood_error,
		Vconnect* verlink,
		NUMERIC_TYPE* DEM
	);

	static void cal_flow_exchange
	(
		Vconnect* verlink,
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE* H, 
		NUMERIC_TYPE dt,
		Geometry& geometry,
		const double* pending_volumes = nullptr
	);

	static NUMERIC_TYPE cal_orificeflow
	(
		NUMERIC_TYPE c0,
		NUMERIC_TYPE area,
		NUMERIC_TYPE H1D,
		NUMERIC_TYPE H2D
	);

	static NUMERIC_TYPE cal_weirflow
	(
		NUMERIC_TYPE cw,
		NUMERIC_TYPE w,
		NUMERIC_TYPE h2d
	);

	static void free_vconnect
	(
		Vconnect* verlink
	);

	static void cal_nodeinflow
	(
		Vconnect* verlink
	);

	static void load_floodconfig_fv1
	(
		Fnames& filenames,
		Vconnect* verlink,
		Geometry& geometry,
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE nodata_elevation,
		int verbose
	);

	static void flood_init_fv1
	(
		Vconnect* verlink,
		Geometry& geometry
	);

	static void cal_demindex_fv1
	(
		FILE* flood_error,
		Vconnect* verlink,
		Geometry& geometry,
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE nodata_elevation
	);
};

}
}
