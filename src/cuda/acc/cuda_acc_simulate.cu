#include "cuda_dem.cuh"
#include "../cuda_coupling_gather.cuh"
#include "cuda_sample.cuh"
#include "cuda_simulate.cuh"
#include "cuda_acc_flow.cuh"
#include "cuda_acc_simulate.cuh"
#include "cuda_acc_snapshot.cuh"
#include "cuda_acc_solver.cuh"
#include "cuda_rain.cuh"
#include "cuda_max_field.cuh"
#include "cuda_stats.cuh"
#include "cuda_unifiedallocator.cuh"
#include "../rain/rain.h"
#include "io.h"
#include <cstdio>
#include <ctime>
#include <unordered_map>
#include <cfloat>
#include <vector>
#include <chrono>
#include "stats.h"
#include "../../rain/Rainfall_uniform.h"
#include "cuda_boundary.cuh"
extern "C" {
#include "../headers.h"
};
#include "../swmm5.h"
#include "../../swmm_couple.h"
#include "../../coupling_restart.h"
#include "../run_audit.cuh"

//#include "../acc_nugrid/CHECK_CUDA_ERROR.cuh" // use for debug

void lis::cuda::acc::Simulation::run
(
	Fnames& filenames,
	States& states,
	Pars& pars,
	::Solver& solver,
	int verbose
)
{
	print_device_info();
	dim3 grid_size(128, 128);

	Geometry geometry;
	int pitch = 0;
	int offset = 0;

	NUMERIC_TYPE* DEM = Topography::load(filenames.demfilename, geometry,
			pitch, offset, pars.nodata_elevation, states.acceleration, verbose);

	//CHECK_CUDA_ERROR(cuda::peek());
	//CHECK_CUDA_ERROR(cuda::sync());

	update_geometry(pars, geometry);

	cuda::copy_to_symbol(cuda::geometry, &geometry, sizeof(Geometry));

	cuda::copy_to_symbol(cuda::pitch, &pitch, sizeof(int));

	//--------------------------------------------------------------------------------------------
	double elapsedTime = 0.0;
	ErrorCode = 0;
	const bool swmm_coupled = swmm_coupling_enabled(&filenames);

	Vconnect vconnection = {};
	Vconnect* verlink = &vconnection;

	if (swmm_coupled)
	{
		if (swmm_open(filenames.swmminputfilename,
				filenames.swmmreportfilename,
				filenames.swmmoutputfilename) != 0 || swmm_start(1) != 0)
		{
			fprintf(stderr, "ERROR: unable to initialise SWMM coupling.\n");
			swmm_close();
			delete[] DEM;
			throw std::runtime_error("Unable to initialise SWMM coupling; see the SWMM report");
		}

		lis::cuda::Topography::load_floodconfig(filenames, verlink, geometry,
				DEM, pars.nodata_elevation, verbose);
	}

	NUMERIC_TYPE* d_DEM = GhostRaster::allocate_device_H(geometry);

	GhostRaster::copy_H(d_DEM, DEM, geometry);

//	cuda::sync();

    ::DynamicRain<cuda::UnifiedAllocator<NUMERIC_TYPE>> rain(
			filenames.dynamicrainfilename, verbose);
	rain.validate_end(solver.Sim_Time);

	cuda::copy_to_symbol(cuda::rain_geometry, &rain.geometry(), sizeof(Geometry));

	acc::Flow U;

	acc::Flow::allocate_pinned(U, geometry);
	initialise_H(U.H, filenames.startfilename, states, DEM,
			geometry, pitch, offset, verbose);

	nullify_max(U.maxHtm, U.initHtm, U.totalHtm, geometry, pitch, offset, verbose);

	NUMERIC_TYPE* manning = nullptr;
	NUMERIC_TYPE* d_manning = nullptr;
	if (strlen(filenames.nfilename) > 0)
	{
		manning = lis::GhostRaster::allocate_H(geometry);
		initialise_manning(manning, filenames.nfilename, geometry,
				pitch, offset, verbose);
		d_manning = lis::cuda::GhostRaster::allocate_device_H(geometry);
		GhostRaster::copy_H(d_manning, manning, geometry);
	}


	BoundCs boundCs;
	load_boundaries(filenames, states, pars, boundCs, verbose);

	Rainfall_uniform   rainfall_uniform(C(0.0));

//	MaxField max_field(filenames.resrootname, geometry, pitch, offset,
//		grid_size, states.acceleration, pars.output_precision);



	StatsCollector stats_collector(geometry, solver.DepthThresh, states.acceleration);


	//Load checkpointed data if this job has been restarted
	if (states.checkpoint == ON)
	{

		read_checkpoint(filenames, states, pars, solver,
			U.H, U.Qx, U.Qy, U.maxH, U.totalHtm, U.maxHtm, U.initHtm, stats_collector, verbose);

		if (verbose == ON) printf(" - checkpoint output file: %s\n", filenames.checkpointfilename);
	}

	if (strlen(filenames.hsffloodread) != 0 && !coupling_restart_file(filenames.hsffloodread))
	{
		read_hsfflood(filenames.hsffloodread, pars, U.H, U.Qx, U.Qy, verbose);
		if (verbose == ON) printf("hsfflood read done!\n");
	}


	PhysicalParams physical_params(pars, solver);

	cuda::copy_to_symbol(cuda::physical_params, &physical_params,
			sizeof(PhysicalParams));

	SolverParams solver_params(pars, solver, states);

	cuda::copy_to_symbol(cuda::solver_params, &solver_params,
			sizeof(SolverParams));

	//--------------------------------------------------------------------
	/*
	boundCsptr->numPS = verlink->count_dem_index;
	int* temp_xpi = new int[boundCsptr->numPS](); //used in legacy and read in
	int* temp_ypi = new int[boundCsptr->numPS](); //used in legacy and read in
	ESourceType* temp_PS_Ident = new ESourceType[boundCsptr->numPS];
	NUMERIC_TYPE* temp_PS_Val = new NUMERIC_TYPE[boundCsptr->numPS]();
	for (size_t tempindex = 0; tempindex < verlink->count_dem_index; tempindex++)
	{
		temp_xpi[tempindex] = verlink->dem_xi[tempindex];
		temp_ypi[tempindex] = verlink->dem_yi[tempindex];
		temp_PS_Ident[tempindex] = QFIX4;
		temp_PS_Val[tempindex] = verlink->flow_2D_dem[tempindex];
	}
	boundCsptr->PS_Ident = temp_PS_Ident;
	boundCsptr->xpi = temp_xpi;
	boundCsptr->ypi = temp_ypi;
	boundCsptr->PS_Val = temp_PS_Val;
	*/


	if (swmm_coupled)
	{
		Node_Inflow = (double*)calloc(Nobjects[NODE], sizeof(double));
	}
	// Compact gathered depths; no domain-sized host staging raster.


	//--------------------------------------------------------------------


	//cuda::Boundary::initialise(cuda::boundaries, boundCs, pitch, offset);

	auto solve = unique_ptr<cuda::Solver<Flow>>(new Solver(
				U, d_DEM, d_manning, geometry,
				physical_params, grid_size));

	DynamicTimestepACC<acc::Flow> dynamic_dt(cuda::dt, geometry, solver.InitTstep,
			states.adaptive_ts, *solve);


	//start simulation
	time(&solver.time_start);

	solver.itrn_time_now = solver.itrn_time;
    solver.MinTstep = solver.InitTstep;

//	StatsCollector stats_collector(geometry, solver.DepthThresh);

	Stats stats(filenames.resrootname, pars.MassInt, pars.MassTotal, states.checkpoint, solver.t);

	stats.write_header(solver.t);

	Sampler sampler(cuda::sample_points, cuda::sample_buf_idx, verbose);

	if (states.save_stages == ON)
	{
		sampler.load_sample_points(filenames.stagefilename, geometry,
				pitch, offset);

		sampler.open_stage_file(filenames.resrootname, states.checkpoint, solver.t);

		sampler.write_stage_header(DEM, filenames.stagefilename, states.checkpoint, solver.t);

		if (states.voutput_stage == ON)
		{
			sampler.open_gauge_file(filenames.resrootname, states.checkpoint, solver.t);

			sampler.write_gauge_header(DEM, filenames.gaugefilename, solver.t);
		}
	}

	Snapshot snapshot(filenames.resrootname, pars.SaveInt,
			pars.SaveTotal, pars.SaveNo, DEM, U, geometry, pitch, offset,
			verbose, pars.output_precision);

	if (states.save_elev == ON)
	{
		snapshot.enable_elevation_writer();
	}

	if (states.voutput == ON)
	{
		snapshot.enable_velocity_writer(solver.DepthThresh);
	}

	if (states.save_Qs == ON)
	{
		snapshot.enable_discharge_writer();
	}

	Flow& d_U = solve->d_U();
	stats_collector.initialise_volume(d_U.H);
	RunAudit audit(filenames.resrootname,"CUDA_ACC",swmm_coupled,stats_collector.previous_volume,solver.Sim_Time);

	//--------------------------------------------------------------------------------------

	BoundaryConditions temp_BC;

	temp_BC.time_series_count = boundCs.allTimeSeries.size();

	temp_BC.all_time_series = static_cast<TimeSeries*>(cuda::malloc_device(
		temp_BC.time_series_count * sizeof(TimeSeries)));

	// Store boundary time series in contiguous device buffers.
	size_t ts_total = 0;
	for (int i = 0; i < temp_BC.time_series_count; i++)
		ts_total += (size_t)boundCs.allTimeSeries[i].count;
	NUMERIC_TYPE* d_ts_time = static_cast<NUMERIC_TYPE*>(cuda::malloc_device(ts_total * sizeof(NUMERIC_TYPE)));
	NUMERIC_TYPE* d_ts_value = static_cast<NUMERIC_TYPE*>(cuda::malloc_device(ts_total * sizeof(NUMERIC_TYPE)));
	std::unordered_map<const NUMERIC_TYPE*, int> ts_index_by_value;
	ts_index_by_value.reserve((size_t)temp_BC.time_series_count * 2 + 1);
	{
		std::vector<NUMERIC_TYPE> h_ts_time(ts_total), h_ts_value(ts_total);
		TimeSeries* temp_all_time_series1 = new TimeSeries[temp_BC.time_series_count]();
		size_t ts_offset = 0;
		for (int i = 0; i < temp_BC.time_series_count; i++)
		{
			const int count = boundCs.allTimeSeries[i].count;
			temp_all_time_series1[i].count = count;
			temp_all_time_series1[i].time = d_ts_time + ts_offset;
			temp_all_time_series1[i].value = d_ts_value + ts_offset;
			temp_all_time_series1[i].current_value = C(0.0);
			memcpy(h_ts_time.data() + ts_offset, boundCs.allTimeSeries[i].time, count * sizeof(NUMERIC_TYPE));
			memcpy(h_ts_value.data() + ts_offset, boundCs.allTimeSeries[i].value, count * sizeof(NUMERIC_TYPE));
			ts_offset += (size_t)count;
			ts_index_by_value[boundCs.allTimeSeries[i].value] = i;
		}
		cuda::copy(d_ts_time, h_ts_time.data(), ts_total * sizeof(NUMERIC_TYPE));
		cuda::copy(d_ts_value, h_ts_value.data(), ts_total * sizeof(NUMERIC_TYPE));
		cuda::copy(temp_BC.all_time_series, temp_all_time_series1,
			temp_BC.time_series_count * sizeof(TimeSeries));
		delete[] temp_all_time_series1;
	}
	//--------------------------------------------------------------------------------------
	//lis::cuda::Boundary::initialise_all_time_series(temp_BC, &boundCs);

	temp_BC.BC_count = boundCs.numBCs;
	temp_BC.BC_type = static_cast<ESourceType*>(
		cuda::malloc_device(boundCs.numBCs * sizeof(ESourceType)));
	cuda::copy(temp_BC.BC_type, boundCs.BC_Ident, boundCs.numBCs * sizeof(ESourceType));

	temp_BC.BC_value = static_cast<NUMERIC_TYPE*>(
		cuda::malloc_device(boundCs.numBCs * sizeof(NUMERIC_TYPE)));
	cuda::copy(temp_BC.BC_value, boundCs.BC_Val, boundCs.numBCs * sizeof(NUMERIC_TYPE));

	temp_BC.BC_time_series = static_cast<TimeSeries**>(
		cuda::malloc_device(boundCs.numBCs * sizeof(TimeSeries*)));
	TimeSeries** time_series2 = new TimeSeries * [boundCs.numBCs]();

	for (int bc_i = 0; bc_i < boundCs.numBCs; bc_i++)
	{
		if (boundCs.BC_Ident[bc_i] == HVAR3 || boundCs.BC_Ident[bc_i] == QVAR5)
		{
			// Resolve the series through its pointer lookup.
			::TimeSeries* target = boundCs.BC_TimeSeries[bc_i];
			std::unordered_map<const NUMERIC_TYPE*, int>::const_iterator f =
				target != nullptr ? ts_index_by_value.find(target->value) : ts_index_by_value.end();
			if (f == ts_index_by_value.end())
			{
				fprintf(stderr, "ERROR: edge boundary %d has no time series on the GPU.\n", bc_i);
				exit(EXIT_FAILURE);
			}
			time_series2[bc_i] = temp_BC.all_time_series + f->second;
		}
	}
	cuda::copy(temp_BC.BC_time_series,
		time_series2, boundCs.numBCs * sizeof(TimeSeries*));
	delete[] time_series2;

	//lis::cuda::Boundary::initialise_BC(temp_BC, &boundCs);
	//--------------------------------------------------------------------------------------
	const int temp_numPS_acc = MAX(0, boundCs.numPS); // -1 means no BCI point sources
	const int surface_coupling_count = swmm_coupled
			? static_cast<int>(verlink->count_dem_index) : 0;
	const int outfall_source_count = swmm_coupled
			? static_cast<int>(verlink->outfall_dem_count) : 0;
	const int coupling_count = surface_coupling_count + outfall_source_count;
	const int outfall_source_start = temp_numPS_acc + surface_coupling_count;
	NUMERIC_TYPE* coupling_depths = surface_coupling_count > 0
			? new NUMERIC_TYPE[surface_coupling_count]() : nullptr;
	NUMERIC_TYPE* d_coupling_depths = surface_coupling_count > 0
			? static_cast<NUMERIC_TYPE*>(cuda::malloc_device(
				surface_coupling_count * sizeof(NUMERIC_TYPE))) : nullptr;
	const int tailwater_count = swmm_coupled && verlink->dynamic_outfall_count > 0
			? outfall_source_count : 0;
	NUMERIC_TYPE* tailwater_depths = tailwater_count > 0
			? new NUMERIC_TYPE[tailwater_count]() : nullptr;
	NUMERIC_TYPE* d_tailwater_depths = tailwater_count > 0
			? static_cast<NUMERIC_TYPE*>(cuda::malloc_device(
				tailwater_count * sizeof(NUMERIC_TYPE))) : nullptr;
	NUMERIC_TYPE* temp_psvalue;
		if (temp_numPS_acc > 0)
	{
		temp_BC.PS_count = temp_numPS_acc + coupling_count;

		temp_BC.PS_type = static_cast<ESourceType*>(cuda::malloc_device((coupling_count + temp_numPS_acc) * sizeof(ESourceType)));
		ESourceType* temp_sourcetype = new ESourceType[coupling_count + temp_numPS_acc]();
		memcpy(temp_sourcetype, boundCs.PS_Ident, sizeof(ESourceType) * temp_numPS_acc);

		temp_BC.PS_idx = static_cast<int*>(cuda::malloc_device((coupling_count + temp_numPS_acc) * sizeof(int)));
		int* temp_idx = new int[coupling_count + temp_numPS_acc]();
		for (size_t i = 0; i < temp_numPS_acc; i++)
		{
			temp_idx[i] = boundCs.ypi[i] * pitch + boundCs.xpi[i] + offset;
		}

		temp_BC.PS_value = static_cast<NUMERIC_TYPE*>(cuda::malloc_device((coupling_count + temp_numPS_acc) * sizeof(NUMERIC_TYPE)));
		temp_psvalue = new NUMERIC_TYPE[coupling_count + temp_numPS_acc]();
		memcpy(temp_psvalue, boundCs.PS_Val, sizeof(NUMERIC_TYPE)* temp_numPS_acc);

		for (int tempindex = 0; tempindex < surface_coupling_count; tempindex++)
		{
			temp_idx[temp_numPS_acc + tempindex] = verlink->dem_Vindex_uni[tempindex];
			temp_sourcetype[temp_numPS_acc + tempindex] = QFIX4;
			temp_psvalue[temp_numPS_acc + tempindex] = verlink->flow_2D_dem[tempindex];
		}
		for (int tempindex = 0; tempindex < outfall_source_count; tempindex++)
		{
			temp_idx[outfall_source_start + tempindex] =
					verlink->outfall_dem_Vindex_uni[tempindex];
			temp_sourcetype[outfall_source_start + tempindex] = QFIX4;
			temp_psvalue[outfall_source_start + tempindex] = C(0.0);
		}

		cuda::copy(temp_BC.PS_type, temp_sourcetype, (coupling_count + temp_numPS_acc) * sizeof(ESourceType));
		cuda::copy(temp_BC.PS_idx, temp_idx, (coupling_count + temp_numPS_acc) * sizeof(int));
		cuda::copy(temp_BC.PS_value, temp_psvalue, (coupling_count + temp_numPS_acc) * sizeof(NUMERIC_TYPE));

		delete[] temp_sourcetype;
		delete[] temp_idx;

		temp_BC.PS_time_series = static_cast<TimeSeries**>(cuda::malloc_device((coupling_count + temp_numPS_acc) * sizeof(TimeSeries*)));
		TimeSeries** time_series3 = new TimeSeries * [coupling_count + temp_numPS_acc]();

		for (int ps_i = 0; ps_i < temp_numPS_acc; ps_i++)
		{
			if (boundCs.PS_Ident[ps_i] == HVAR3 || boundCs.PS_Ident[ps_i] == QVAR5)
			{
				// Resolve the series through its pointer lookup.
				::TimeSeries* target = boundCs.PS_TimeSeries[ps_i];
				std::unordered_map<const NUMERIC_TYPE*, int>::const_iterator f =
					target != nullptr ? ts_index_by_value.find(target->value) : ts_index_by_value.end();
				if (f == ts_index_by_value.end())
				{
					fprintf(stderr, "ERROR: point boundary %d has no time series on the GPU.\n", ps_i + 1);
					exit(EXIT_FAILURE);
				}
				time_series3[ps_i] = temp_BC.all_time_series + f->second;
			}
		}
		cuda::copy(temp_BC.PS_time_series, time_series3, (coupling_count + temp_numPS_acc) * sizeof(TimeSeries*));
		delete[] time_series3;
	}
	else
	{
		temp_BC.PS_count = coupling_count;
		temp_BC.PS_type = static_cast<ESourceType*>(cuda::malloc_device(coupling_count * sizeof(ESourceType)));
		ESourceType* temp_sourcetype = new ESourceType[coupling_count]();
		for (int i = 0; i < coupling_count; i++)	temp_sourcetype[i] = QFIX4;
		cuda::copy(temp_BC.PS_type, temp_sourcetype, coupling_count * sizeof(ESourceType));
		delete[] temp_sourcetype;

		temp_BC.PS_idx = static_cast<int*>(cuda::malloc_device(coupling_count * sizeof(int)));
		if (surface_coupling_count > 0)
		{
			cuda::copy(temp_BC.PS_idx, verlink->dem_Vindex_uni,
					surface_coupling_count * sizeof(int));
		}
		if (outfall_source_count > 0)
		{
			cuda::copy(temp_BC.PS_idx + surface_coupling_count,
					verlink->outfall_dem_Vindex_uni,
					outfall_source_count * sizeof(int));
		}

		temp_BC.PS_value = static_cast<NUMERIC_TYPE*>(cuda::malloc_device(coupling_count * sizeof(NUMERIC_TYPE)));
		temp_psvalue = new NUMERIC_TYPE[coupling_count]();
		cuda::copy(temp_BC.PS_value, temp_psvalue, coupling_count * sizeof(NUMERIC_TYPE));

		temp_BC.PS_time_series = static_cast<TimeSeries**>(cuda::malloc_device(coupling_count * sizeof(TimeSeries*)));
		TimeSeries** time_series3 = new TimeSeries * [coupling_count]();

		cuda::copy(temp_BC.PS_time_series, time_series3, coupling_count * sizeof(TimeSeries*));
		delete[] time_series3;
	}

	temp_BC.coupling_start = temp_numPS_acc;
	temp_BC.outfall_start = outfall_source_start;
	cuda::copy_to_symbol(cuda::boundaries, &temp_BC, sizeof(BoundaryConditions));
	// Size boundary launches to the actual entry counts.
	cuda::Boundary::set_launch_counts(temp_BC.time_series_count, temp_BC.BC_count, temp_BC.PS_count);
	if (verbose == ON) printf("GPU boundaries: %d bdy series (%zu values), %d edge cells, %d point sources (%d from bci, %d coupling)\n",
		temp_BC.time_series_count, ts_total, temp_BC.BC_count, temp_BC.PS_count, temp_numPS_acc, coupling_count);

	//NUMERIC_TYPE dt_cpu = 0.0;

	time_t loop_start;
	time(&loop_start);
    const auto precise_loop_start=std::chrono::steady_clock::now();
	int countnum = 0;
	SwmmCoupleClock couple_clock;
	swmm_couple_clock_init(&couple_clock,
			!swmm_coupled || pars.swmm_couple_every_step,
			(double)pars.swmm_coupleint,
			swmm_coupled ? verlink->count_index : 0);
	// Queue positive returns for the next fixed exchange interval.
	SwmmReturnSpread return_spread;
	swmm_return_init(&return_spread,
			swmm_coupled && !couple_clock.every_step ? coupling_count : 0);
	// Accumulate time in double; solver.t is a rounded copy.
	double t_exact = (double)solver.t;
	double swmm_wall = 0.0, couple_wall = 0.0;
	if (swmm_coupled)
		printf("SWMM exchange %s, SWMM threads %d\n",
			couple_clock.every_step ? "every 2D step" : "fixed interval (equal sub-steps, returns spread over the next interval)",
			swmm_getNumThreads());
	audit.settings(swmm_coupled ? swmm_getNumThreads() : 0,
			swmm_coupled ? (couple_clock.every_step ? -1.0 : couple_clock.interval) : 0.0);

    if(swmm_coupled && (filenames.hsffloodtemp[0] || filenames.hsffloodend[0] || coupling_restart_file(filenames.hsffloodread)) &&
        !swmm_hydraulicRestartSupported())
        throw std::runtime_error("Hydraulic restart is unavailable with SWMM runoff, RDII, storage exfiltration, routing-interface input, events or active quality routing");
    const std::size_t elements=lis::GhostRaster::elements_H(geometry);
    const std::size_t faces=lis::GhostRaster::elements_Q(geometry);
    const std::vector<CouplingRestartField> restart_fields={{"H",U.H,elements},{"Qx",U.Qx,faces},{"Qy",U.Qy,faces},
        {"Qxold",U.Qxold,faces},{"Qyold",U.Qyold,faces},{"maxH",U.maxH,elements},
        {"maxHtm",U.maxHtm,elements},{"initHtm",U.initHtm,elements},{"totalHtm",U.totalHtm,elements}};
    if(coupling_restart_file(filenames.hsffloodread)) {
        coupling_restart_io(filenames.hsffloodread,true,3,pars,solver,verlink,swmm_coupled,
            couple_clock,return_spread,t_exact,restart_fields);
        Flow::copy(d_U,U,geometry);

        stats.next_save=pars.MassTotal; snapshot.next_save=pars.SaveTotal; snapshot.counter=pars.SaveNo;
        stats_collector.initialise_volume(d_U.H);
        audit.restart_baseline(stats_collector.previous_volume,
            coupling_pending_in(couple_clock)+swmm_return_pending(&return_spread));
    }
    auto save_restart=[&](const char* path) {
        snapshot.wait(); cuda::sync(); Flow::copy(U,d_U,geometry);

        pars.SaveTotal=snapshot.next_save; pars.SaveNo=snapshot.counter; pars.MassTotal=stats.next_save;
        coupling_restart_io(path,false,3,pars,solver,verlink,swmm_coupled,
            couple_clock,return_spread,t_exact,restart_fields);
    };


	while (t_exact < (double)solver.Sim_Time-1.e-9
			&& (!swmm_coupled
					|| (NewRoutingTime < TotalDuration && !ErrorCode)))
	{
		//---------------------------------------------------------------------------------------------------------------------
		dynamic_dt.update_dt_ACC();

        cuda::dt=(NUMERIC_TYPE)coupling_step(couple_clock,swmm_coupled,t_exact,
            (double)cuda::dt,(double)solver.Sim_Time,NewRoutingTime/1000.0,
            TotalDuration/1000.0,rain.time_to_next_frame(t_exact),
            countnum?0.0:(double)pars.hsfdt);
        if(cuda::dt<=C(0.0)) break;

		// Actual net rain/infiltration is applied before the inlet budget.
		stats_collector.zero_instantaneous_mass();
		if (rain.enabled()) {
			rain.update_time(solver.t);
			lis::cuda::DynamicRain::updateACC<<<grid_size, CUDA_BLOCK_SIZE>>>(
				d_DEM, d_U.H, rain.data(), stats_collector.instantaneous_mass());
		}
		if (states.rainfall == ON) {
			rainfall_uniform.update_rainfall(filenames.rainfilename, solver.t);
			solve->update_uniform_rain(rainfall_uniform.value);
		}

		if (swmm_coupled)
		{
			// Reset only cells that can exchange with SWMM.  The historical
			// implementation cleared two complete domain-sized arrays here.
			for (int coupling_i = 0; coupling_i < surface_coupling_count; ++coupling_i)
			{
				const int dem_i = coupling_i;
				verlink->flow_neg_count[dem_i] = 0;
				verlink->flow_2D_dem[dem_i] = 0.0;
			}
			memset(verlink->flow_1D_node, 0, verlink->count_index * sizeof(double));
			memset(verlink->flow_2D_node, 0, verlink->count_index * sizeof(double));
			memset(Node_Inflow, 0, Nobjects[NODE] * sizeof(double));

			const auto cw0 = std::chrono::steady_clock::now();
			// Gather and transfer only the unique coupled-cell depths (693 or
			// fewer for the San Francisco case), rather than the 158M-cell H grid.
			lis::cuda::gather_coupling_depths(
				d_U.H, temp_BC.PS_idx + temp_numPS_acc,
				d_coupling_depths, surface_coupling_count);
			cuda::copy(coupling_depths, d_coupling_depths,
				surface_coupling_count * sizeof(NUMERIC_TYPE));
			lis::cuda::Topography::cal_flow_exchange(verlink, DEM, coupling_depths, cuda::dt, geometry, couple_clock.volume);
			for(int k=0;k<surface_coupling_count;++k)
				audit.exchange(verlink->flow_2D_dem[k]*geometry.dy*(double)cuda::dt);
		if (couple_clock.every_step)
		{
			lis::cuda::Topography::cal_nodeinflow(verlink);
		}
		else
		{
			swmm_couple_accumulate(&couple_clock, verlink, (double)cuda::dt);
		}

		if (temp_numPS_acc > 0)
		{
			for (size_t flowindex = 0; flowindex < verlink->count_dem_index; flowindex++)
			{
				temp_psvalue[temp_numPS_acc + flowindex] = verlink->flow_2D_dem[flowindex];
			}
			for (int flowindex = 0; flowindex < outfall_source_count; ++flowindex)
				temp_psvalue[outfall_source_start + flowindex] = C(0.0);
			// Apply queued returns from the preceding interval.
			if (return_spread.count > 0)
				audit.exchange(swmm_return_apply(&return_spread, temp_psvalue + temp_numPS_acc,
						(double)geometry.dy, (double)cuda::dt));
			// Copy coupling sources only; user boundaries are maintained on device.
			cuda::copy(temp_BC.PS_value + temp_numPS_acc, temp_psvalue + temp_numPS_acc,
					coupling_count * sizeof(NUMERIC_TYPE));
		}
		else
		{
			for (size_t flowindex = 0; flowindex < verlink->count_dem_index; flowindex++)
			{
				temp_psvalue[flowindex] = verlink->flow_2D_dem[flowindex];
			}
			for (int flowindex = 0; flowindex < outfall_source_count; ++flowindex)
				temp_psvalue[surface_coupling_count + flowindex] = C(0.0);
			// Apply queued returns from the preceding interval.
			if (return_spread.count > 0)
				audit.exchange(swmm_return_apply(&return_spread, temp_psvalue,
						(double)geometry.dy, (double)cuda::dt));
			cuda::copy(temp_BC.PS_value, temp_psvalue,
					coupling_count * sizeof(NUMERIC_TYPE));
		}
		couple_wall += std::chrono::duration<double>(std::chrono::steady_clock::now() - cw0).count();
		}


		// The timestep remains frozen after net rain and inlet exchange.

		solve->FloodplainQ();

		cuda::Boundary::update_time_series(solver.t);

        solve->update_ghost_cells();

		if (pars.drain_nodata)
		{

			cuda::Boundary::drain_nodata_waterACC(d_U.H, d_DEM,
					stats_collector.instantaneous_mass(), grid_size);
		}

		// Keep drain_nodata and point-source/coupling fluxes in the same
		// instantaneous Qin/Qout, matching CPU LISFLOOD.
		cuda::Boundary::update_H(d_U.H, d_DEM, d_U.Qx, d_U.Qy,
				stats_collector.instantaneous_mass(), grid_size, states.drychecking);

		int swmm_error = 0;
		if (swmm_coupled)
		{
			// update_H is asynchronous.  Synchronize before the host reads the
			// managed time-step value and advances SWMM (required on WDDM too).
			cuda::sync();
			const double swmm_end = TotalDuration / 1000.0;
			const int exchange_now = swmm_couple_due(&couple_clock,
					t_exact + (double)cuda::dt, solver.Sim_Time, swmm_end);
			if (exchange_now)
			{
				if (couple_clock.every_step) linkDt = cuda::dt;
				else
				{
					swmm_couple_write_inflow(&couple_clock, verlink);
					linkDt = couple_clock.elapsed;
				}
				swmm_beginCoupledExchange();
				if (tailwater_count > 0)
				{
					lis::cuda::gather_coupling_depths(d_U.H,
						temp_BC.PS_idx + outfall_source_start,
						d_tailwater_depths, tailwater_count);
					cuda::copy(tailwater_depths, d_tailwater_depths,
						tailwater_count * sizeof(NUMERIC_TYPE));
					lis::cuda::Topography::update_outfall_stages(
						verlink, DEM, tailwater_depths, geometry.dx * geometry.dy);
				}
				const auto sw0 = std::chrono::steady_clock::now();
				swmm_error = swmm_advanceCoupled((double)linkDt, &elapsedTime);
				swmm_wall += std::chrono::duration<double>(std::chrono::steady_clock::now() - sw0).count();
				cuda::sync();
				if (swmm_error == 0 && coupling_count > 0)
				{
					lis::cuda::Topography::collect_outfall_source_values(
							verlink, temp_psvalue, outfall_source_start, geometry.dy,
							geometry.dx * geometry.dy, tailwater_depths,
							(double)linkDt);
					lis::cuda::Topography::collect_surface_flood_values(
						verlink, temp_psvalue, temp_numPS_acc, geometry.dy, (double)linkDt);
                    double network_return=0.0, network_inlet=0.0;
                    for(int k=0;k<coupling_count;++k)
                        network_return+=(double)temp_psvalue[temp_numPS_acc+k]*(double)geometry.dy*(double)linkDt;
                    for(int n=0;n<Nobjects[NODE];++n) network_inlet+=Node_Inflow[n]*(double)linkDt;
                    audit.network_exchange(network_return-network_inlet);
					if (!couple_clock.every_step && cuda::dt > C(0.0))
					{
						// Queue positive returns and apply reverse flow now; flush at the end.
						const double t_after = t_exact + (double)cuda::dt;
						double next_interval = couple_clock.interval;
						if ((double)solver.Sim_Time - t_after < next_interval) next_interval = (double)solver.Sim_Time - t_after;
						if (swmm_end - t_after < next_interval) next_interval = swmm_end - t_after;
						if(next_interval<1.e-6) next_interval=0.0;
						swmm_return_schedule(&return_spread, temp_psvalue + temp_numPS_acc,
								couple_clock.elapsed, next_interval, (double)geometry.dy, (double)cuda::dt);
					}
					for(int k=0;k<coupling_count;++k)
						audit.exchange(temp_psvalue[temp_numPS_acc+k]*geometry.dy*(double)cuda::dt);
					cudaError_t copy_status = cudaMemcpy(temp_BC.PS_value + temp_numPS_acc,
							temp_psvalue + temp_numPS_acc,
							coupling_count * sizeof(NUMERIC_TYPE),cudaMemcpyHostToDevice);
					if (copy_status != cudaSuccess) { fprintf(stderr,"source copy failed: %s\n",cudaGetErrorString(copy_status)); exit(EXIT_FAILURE); }
					cuda::Boundary::update_point_sources_range(
							d_U.H, nullptr, nullptr, d_DEM,
							stats_collector.instantaneous_mass(),
							temp_numPS_acc);
				}
				if (!couple_clock.every_step) swmm_couple_advance(&couple_clock);
			}
		}

		stats_collector.accumulate_mass();

		d_U = solve->update_flow_variables(
				stats_collector.instantaneous_mass());

		cuda::sync();

		solver.MinTstep = FMIN(solver.MinTstep, cuda::dt);

		if (verbose) printf("t=%f\tdt=%f\n", solver.t, cuda::dt);
		t_exact += (double)cuda::dt;
		solver.t = (NUMERIC_TYPE)t_exact;
		solver.itCount++;
		if (swmm_error != 0) break;




		if (stats.need_to_write(solver.t) || states.mint_hk == OFF)
		{
			solve->updateMaxFieldACC(solver.t);
		}

		snapshot.write_if_needed(d_U, solver.t, states.call_gzip);
		pars.SaveTotal = snapshot.next_save;
		pars.SaveNo = snapshot.counter;

		if (stats.need_to_write(solver.t))
		{
//			solve->updateMaxFieldACC(solver.t);

			sampler.sample_ACC(d_U.H, d_U.Vx, d_U.Vy, solver.t);

//			lis::StatsEntry Stats_Entry;

			lis::StatsEntry Stats_Entry = stats_collector.create_entry(solver, d_U.H);

			stats.write(Stats_Entry);
			double pending=0;
			if(swmm_coupled && !couple_clock.every_step)
				for(int k=0;k<verlink->count_index;++k) pending+=couple_clock.volume[k];
			audit.sample(Stats_Entry,stats_collector,pending,swmm_return_pending(&return_spread),swmm_wall,couple_wall);

			pars.MassTotal = stats.next_save;

			sampler.write_if_buffer_full();

			// Checkpointing
			if (states.checkpoint == ON)
			{
				//iteration time
				time(&solver.time_check);
				solver.itrn_time_now = solver.itrn_time + (NUMERIC_TYPE)difftime(solver.time_check, solver.time_start);

				snapshot.wait();
				Flow::copy(U, d_U, geometry);

				if (solver.itrn_time_now >= pars.nextcheck)
				{
//					NUMERIC_TYPE SaveTotal, MassTotal;
//					int SaveNo;


					write_checkpoint(filenames, states, pars, solver,
						U.H, U.Qx, U.Qy, U.maxH, U.totalHtm, U.maxHtm, U.initHtm, Stats_Entry, cuda::dt, verbose);

					pars.nextcheck = solver.itrn_time_now + (pars.checkfreq * 3600);
				}
			}
		}

		// If requested on command line, check whether we should kill this simulation...
		if (states.killsim == ON)
		{
			//iteration time
			time(&solver.time_check);
			solver.itrn_time_now = solver.itrn_time + (NUMERIC_TYPE)difftime(solver.time_check, solver.time_start);
			// check if we have reached the kill time
			if (solver.itrn_time_now >= pars.killsim_time) {
				if (verbose == ON) printf("Simulation kill time reached... ");
				break;
			}
		}
		//printf("done222\n");
        if(!countnum && pars.hsfdt>0 && t_exact+1.e-7>=(double)pars.hsfdt) {
            if(filenames.hsffloodtemp[0]) save_restart(filenames.hsffloodtemp);
            if(swmm_coupled && filenames.hsfswmmtemp[0]) {
                FILE* legacy=fopen(filenames.hsfswmmtemp,"wb");
                if(!legacy) throw std::runtime_error("Cannot write SWMM warm start");
                saveTempHsf(legacy);
            }
            countnum=1;
        }

	}
	//---------------------------------------------------------------------------
	//fclose(dUfiletest2);
	//---------------------------------------------------------------------------
	if(filenames.hsffloodend[0]) save_restart(filenames.hsffloodend);
	double final_pending=0;
	if(swmm_coupled && !couple_clock.every_step)
		for(int k=0;k<verlink->count_index;++k) final_pending+=couple_clock.volume[k];
	audit.sample(stats_collector.create_entry(solver,d_U.H),stats_collector,final_pending,swmm_return_pending(&return_spread),swmm_wall,couple_wall);
	swmm_couple_clock_free(&couple_clock);
	swmm_return_free(&return_spread);
	if (swmm_coupled)
	{
		swmm_end();
		if (Fout.mode == SCRATCH_FILE) swmm_report();

		// --- close the system
		swmm_close();
		free(Node_Inflow); Node_Inflow = nullptr;
	}

	time_t loop_end;
	time(&loop_end);
	double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-precise_loop_start).count();
	printf("loop time %lf\n", seconds);
	if (swmm_coupled) audit.timing((double)seconds, swmm_wall, couple_wall, (long)solver.itCount);

	//output max files
	snapshot.write_max_files(d_U);


	if (swmm_coupled) lis::cuda::Topography::free_vconnect(verlink);

	snapshot.wait();
	audit.finish(solver.t,solver.Sim_Time,ErrorCode);
	delete[] DEM;
	delete[] manning;

	delete[] coupling_depths;
	delete[] tailwater_depths;
	delete[] temp_psvalue;
	cuda::free_device(d_coupling_depths);
	cuda::free_device(d_tailwater_depths);
	cuda::free_device(d_DEM);
	cuda::free_device(d_manning);
	Flow::free_pinned(U);
	// Release the allocated boundary buffers.
	cuda::free_device(d_ts_time);
	cuda::free_device(d_ts_value);
	cuda::free_device(temp_BC.all_time_series);
	cuda::free_device(temp_BC.BC_type);
	cuda::free_device(temp_BC.BC_value);
	cuda::free_device(temp_BC.BC_time_series);
	cuda::free_device(temp_BC.PS_type);
	cuda::free_device(temp_BC.PS_value);
	cuda::free_device(temp_BC.PS_time_series);
	cuda::free_device(temp_BC.PS_idx);
}
