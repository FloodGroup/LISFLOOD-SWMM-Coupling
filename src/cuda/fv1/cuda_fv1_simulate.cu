#include "cuda_dem.cuh"
#include "../cuda_coupling_gather.cuh"
#include "cuda_sample.cuh"
#include "cuda_simulate.cuh"
#include "cuda_fv1_flow.cuh"
#include "cuda_fv1_simulate.cuh"
#include "cuda_fv1_snapshot.cuh"
#include "cuda_fv1_solver.cuh"
#include "cuda_rain.cuh"
#include "cuda_max_field.cuh"
#include "cuda_stats.cuh"
#include "cuda_unifiedallocator.cuh"
#include "../rain/rain.h"
#include "io.h"
#include <cstdio>
#include <ctime>
#include <chrono>
#include <vector>
#include "../../rain/Rainfall_uniform.h"
#include "cuda_boundary.cuh"
extern "C" {
#include "../headers.h"
};
#include "../swmm5.h"
#include "../../swmm_couple.h"
#include "../../coupling_restart.h"
#include "../run_audit.cuh"
#include "stats.h"

void lis::cuda::fv1::Simulation::run
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

	NUMERIC_TYPE* DEM = Topography::load(filenames.demfilename, geometry, pitch, offset, pars.nodata_elevation, states.acceleration, verbose);
	Topography::clamp_boundary_values(DEM, geometry, pitch);

	update_geometry(pars, geometry);
	cuda::copy_to_symbol(cuda::geometry, &geometry, sizeof(Geometry));
	cuda::copy_to_symbol(cuda::pitch, &pitch, sizeof(int));

	//??????????
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

		lis::cuda::Topography::load_floodconfig_fv1(filenames, verlink,
				geometry, DEM, pars.nodata_elevation, verbose);
	}

	//--------------------------------------------------------------------------------------------

	NUMERIC_TYPE* d_DEM = GhostRaster::allocate_device(geometry);
	GhostRaster::copy(d_DEM, DEM, geometry);

	NUMERIC_TYPE* Zstar_x = GhostRaster::allocate_device(geometry);
	Topography::initialise_Zstar_x(Zstar_x, d_DEM);

	NUMERIC_TYPE* Zstar_y = GhostRaster::allocate_device(geometry);
	Topography::initialise_Zstar_y(Zstar_y, d_DEM);

	cuda::sync();

    ::DynamicRain<cuda::UnifiedAllocator<NUMERIC_TYPE>> rain(filenames.dynamicrainfilename, verbose);
	rain.validate_end(solver.Sim_Time);
	cuda::copy_to_symbol(cuda::rain_geometry, &rain.geometry(), sizeof(Geometry));

	fv1::Flow U; 
	fv1::Flow::allocate_pinned(U, geometry);
	initialise_H(U.H, filenames.startfilename, states, DEM, geometry, pitch, offset, verbose);
	initialise_discharge(U.HU, U.HV, filenames.startfilename, states, geometry, pitch, offset, verbose);

	NUMERIC_TYPE* manning = nullptr;
	NUMERIC_TYPE* d_manning = nullptr;
	if (strlen(filenames.nfilename) > 0)
	{
	    manning = lis::GhostRaster::allocate(geometry);
		initialise_manning(manning, filenames.nfilename, geometry, pitch, offset, verbose);
		d_manning = lis::cuda::GhostRaster::allocate_device(geometry);
		GhostRaster::copy(d_manning, manning, geometry);
	}

	PhysicalParams physical_params(pars, solver);
	cuda::copy_to_symbol(cuda::physical_params, &physical_params, sizeof(PhysicalParams));

	SolverParams solver_params(pars, solver, states);
	cuda::copy_to_symbol(cuda::solver_params, &solver_params, sizeof(SolverParams));

	Rainfall_uniform   rainfall_uniform(C(0.0));

	BoundCs boundCs;
	load_boundaries(filenames, states, pars, boundCs, verbose);
	//cuda::Boundary::initialise(cuda::boundaries, boundCs, pitch, offset);

	auto solve = unique_ptr<cuda::Solver<Flow>>(new Solver(U, d_DEM, Zstar_x, Zstar_y, d_manning, geometry, physical_params, grid_size));

	DynamicTimestep<fv1::Flow> dynamic_dt(cuda::dt, geometry, solver.InitTstep, states.adaptive_ts, *solve);
	solver.MinTstep = solver.InitTstep;

	StatsCollector stats_collector(geometry, solver.DepthThresh, states.acceleration);
	MaxField max_field(filenames.resrootname, geometry, pitch, offset, grid_size, states.acceleration, pars.output_precision);

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

	Snapshot snapshot(filenames.resrootname, pars.SaveInt, pars.SaveTotal, pars.SaveNo, DEM, U, geometry, pitch, offset, verbose, pars.output_precision);
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
	RunAudit audit(filenames.resrootname,"CUDA_FV1",swmm_coupled,stats_collector.previous_volume,solver.Sim_Time);

	//------------------------------------------------------------------------------------------------------------
	if (swmm_coupled)
	{
		Node_Inflow = (double*)calloc(Nobjects[NODE], sizeof(double));
	}


	BoundaryConditions temp_BC;

	temp_BC.time_series_count = boundCs.allTimeSeries.size();

	temp_BC.all_time_series = static_cast<TimeSeries*>(cuda::malloc_device(temp_BC.time_series_count * sizeof(TimeSeries)));

	TimeSeries* temp_all_time_series1 = new TimeSeries[temp_BC.time_series_count];
	for (int i = 0; i < temp_BC.time_series_count; i++)
	{
		int count = boundCs.allTimeSeries[i].count;

		temp_all_time_series1[i].count = count;

		temp_all_time_series1[i].time = static_cast<NUMERIC_TYPE*>(cuda::malloc_device(count * sizeof(NUMERIC_TYPE)));
		temp_all_time_series1[i].value = static_cast<NUMERIC_TYPE*>(cuda::malloc_device(count * sizeof(NUMERIC_TYPE)));

		cuda::copy(temp_all_time_series1[i].time, boundCs.allTimeSeries[i].time, count * sizeof(NUMERIC_TYPE));
		cuda::copy(temp_all_time_series1[i].value, boundCs.allTimeSeries[i].value, count * sizeof(NUMERIC_TYPE));
	}
	cuda::copy(temp_BC.all_time_series, temp_all_time_series1, temp_BC.time_series_count * sizeof(TimeSeries));
	delete[] temp_all_time_series1;
	//--------------------------------------------------------------------------------------
	//lis::cuda::Boundary::initialise_all_time_series(temp_BC, &boundCs);

	temp_BC.BC_count = boundCs.numBCs;
	temp_BC.BC_type = static_cast<ESourceType*>(cuda::malloc_device(boundCs.numBCs * sizeof(ESourceType)));
	cuda::copy(temp_BC.BC_type, boundCs.BC_Ident, boundCs.numBCs * sizeof(ESourceType));

	temp_BC.BC_value = static_cast<NUMERIC_TYPE*>(cuda::malloc_device(boundCs.numBCs * sizeof(NUMERIC_TYPE)));
	cuda::copy(temp_BC.BC_value, boundCs.BC_Val, boundCs.numBCs * sizeof(NUMERIC_TYPE));

	temp_BC.BC_time_series = static_cast<TimeSeries**>(cuda::malloc_device(boundCs.numBCs * sizeof(TimeSeries*)));
	TimeSeries** time_series2 = new TimeSeries * [boundCs.numBCs];

	for (int bc_i = 0; bc_i < boundCs.numBCs; bc_i++)
	{
		if (boundCs.BC_Ident[bc_i] == HVAR3 || boundCs.BC_Ident[bc_i] == QVAR5)
		{

			::TimeSeries* target = boundCs.BC_TimeSeries[bc_i];

			for (std::vector<::TimeSeries>::size_type ts_i = 0,
				e = boundCs.allTimeSeries.size(); ts_i != e; ts_i++)
			{
				::TimeSeries* candidate = &(boundCs.allTimeSeries[ts_i]);
				if (candidate->value == target->value)
				{
					time_series2[bc_i] = &(temp_BC.all_time_series[ts_i]);
					//+ ts_i*sizeof(TimeSeries*);
				}
			}
		}
	}
	cuda::copy(temp_BC.BC_time_series, time_series2, boundCs.numBCs * sizeof(TimeSeries*));
	delete[] time_series2;

	//lis::cuda::Boundary::initialise_BC(temp_BC, &boundCs);
	//--------------------------------------------------------------------------------------
	const int temp_numPS_fv1 = MAX(0, boundCs.numPS); // -1 means no BCI point sources
	const int surface_coupling_count = swmm_coupled
			? static_cast<int>(verlink->count_dem_index) : 0;
	const int outfall_source_count = swmm_coupled
			? static_cast<int>(verlink->outfall_dem_count) : 0;
	const int coupling_count = surface_coupling_count + outfall_source_count;
	const int outfall_source_start = temp_numPS_fv1 + surface_coupling_count;
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
	if (temp_numPS_fv1 > 0)
	{
		temp_BC.PS_count = temp_numPS_fv1 + coupling_count;

		temp_BC.PS_type = static_cast<ESourceType*>(cuda::malloc_device((coupling_count + temp_numPS_fv1) * sizeof(ESourceType)));
		ESourceType* temp_sourcetype = new ESourceType[coupling_count + temp_numPS_fv1]();
		memcpy(temp_sourcetype, boundCs.PS_Ident, sizeof(ESourceType) * temp_numPS_fv1);

		temp_BC.PS_idx = static_cast<int*>(cuda::malloc_device((coupling_count + temp_numPS_fv1) * sizeof(int)));
		int* temp_idx = new int[coupling_count + temp_numPS_fv1]();
		for (size_t i = 0; i < temp_numPS_fv1; i++)
		{
			temp_idx[i] = boundCs.ypi[i] * pitch + boundCs.xpi[i] + offset;
		}

		temp_BC.PS_value = static_cast<NUMERIC_TYPE*>(cuda::malloc_device((coupling_count + temp_numPS_fv1) * sizeof(NUMERIC_TYPE)));
		temp_psvalue = new NUMERIC_TYPE[coupling_count + temp_numPS_fv1]();
		memcpy(temp_psvalue, boundCs.PS_Val, sizeof(NUMERIC_TYPE)* temp_numPS_fv1);

		for (int tempindex = 0; tempindex < surface_coupling_count; tempindex++)
		{
			temp_idx[temp_numPS_fv1 + tempindex] = verlink->dem_Vindex_uni[tempindex];
			temp_sourcetype[temp_numPS_fv1 + tempindex] = QFIX4;
			temp_psvalue[temp_numPS_fv1 + tempindex] = verlink->flow_2D_dem[tempindex];
		}
		for (int tempindex = 0; tempindex < outfall_source_count; tempindex++)
		{
			temp_idx[outfall_source_start + tempindex] =
					verlink->outfall_dem_Vindex_uni[tempindex];
			temp_sourcetype[outfall_source_start + tempindex] = QFIX4;
			temp_psvalue[outfall_source_start + tempindex] = C(0.0);
		}

		cuda::copy(temp_BC.PS_type, temp_sourcetype, (coupling_count + temp_numPS_fv1) * sizeof(ESourceType));
		cuda::copy(temp_BC.PS_idx, temp_idx, (coupling_count + temp_numPS_fv1) * sizeof(int));
		cuda::copy(temp_BC.PS_value, temp_psvalue, (coupling_count + temp_numPS_fv1) * sizeof(NUMERIC_TYPE));

		delete[] temp_sourcetype;
		delete[] temp_idx;

		temp_BC.PS_time_series = static_cast<TimeSeries**>(cuda::malloc_device((coupling_count + temp_numPS_fv1) * sizeof(TimeSeries*)));
		TimeSeries** time_series3 = new TimeSeries * [coupling_count + temp_numPS_fv1]();

		for (int ps_i = 0; ps_i < temp_numPS_fv1; ps_i++)
		{
			if (boundCs.PS_Ident[ps_i] == HVAR3 || boundCs.PS_Ident[ps_i] == QVAR5)
			{
				::TimeSeries* target = boundCs.PS_TimeSeries[ps_i];

				for (std::vector<::TimeSeries>::size_type ts_i = 0,
					e = boundCs.allTimeSeries.size(); ts_i != e; ts_i++)
				{
					::TimeSeries* candidate = &(boundCs.allTimeSeries[ts_i]);
					if (candidate->value == target->value)
					{
						time_series3[ps_i] = temp_BC.all_time_series + ts_i;
					}
				}
			}
		}
		cuda::copy(temp_BC.PS_time_series, time_series3, (coupling_count + temp_numPS_fv1) * sizeof(TimeSeries*));
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
	temp_BC.coupling_start = temp_numPS_fv1;
	temp_BC.outfall_start = outfall_source_start;
	cuda::copy_to_symbol(cuda::boundaries, &temp_BC, sizeof(BoundaryConditions));

	//--------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------
	time_t loop_start;                                                          
	time(&loop_start);
    const auto precise_loop_start=std::chrono::steady_clock::now();
	SwmmCoupleClock couple_clock;
	swmm_couple_clock_init(&couple_clock,
			!swmm_coupled || pars.swmm_couple_every_step,
			(double)pars.swmm_coupleint,
			swmm_coupled ? verlink->count_index : 0);


    SwmmReturnSpread return_spread;
    swmm_return_init(&return_spread,swmm_coupled && !couple_clock.every_step?coupling_count:0);
    int countnum=0;
    double t_exact=(double)solver.t,swmm_wall=0.0,couple_wall=0.0;
    if(swmm_coupled && (filenames.hsffloodtemp[0] || filenames.hsffloodend[0] || coupling_restart_file(filenames.hsffloodread)) &&
        !swmm_hydraulicRestartSupported())
        throw std::runtime_error("Hydraulic restart is unavailable with SWMM runoff, RDII, storage exfiltration, routing-interface input, events or active quality routing");
    const std::size_t elements=lis::GhostRaster::elements(geometry);
    std::vector<NUMERIC_TYPE> restart_max(elements,0);
    const std::vector<CouplingRestartField> restart_fields={
        {"H",U.H,elements},{"HU",U.HU,elements},{"HV",U.HV,elements},{"maxH",restart_max.data(),elements}};
    if(coupling_restart_file(filenames.hsffloodread)) {
        coupling_restart_io(filenames.hsffloodread,true,4,pars,solver,verlink,swmm_coupled,
            couple_clock,return_spread,t_exact,restart_fields);
        Flow::copy(d_U,U,geometry);
        cuda::copy(max_field.restart_depth(),restart_max.data(),elements*sizeof(NUMERIC_TYPE));
        stats.next_save=pars.MassTotal; snapshot.next_save=pars.SaveTotal; snapshot.counter=pars.SaveNo;
        solve->zero_ghost_cells(); cuda::sync();
        stats_collector.initialise_volume(d_U.H);
        audit.restart_baseline(stats_collector.previous_volume,
            coupling_pending_in(couple_clock)+swmm_return_pending(&return_spread));
    } else if(filenames.hsffloodread[0]) {
        throw std::runtime_error("FV1 requires a bundled hydraulic restart, not an ACC warm start");
    }
    auto save_restart=[&](const char* path) {
        snapshot.wait(); cuda::sync(); Flow::copy(U,d_U,geometry);
        cuda::copy(restart_max.data(),max_field.restart_depth(),elements*sizeof(NUMERIC_TYPE));
        pars.SaveTotal=snapshot.next_save; pars.SaveNo=snapshot.counter; pars.MassTotal=stats.next_save;
        coupling_restart_io(path,false,4,pars,solver,verlink,swmm_coupled,
            couple_clock,return_spread,t_exact,restart_fields);
    };
    audit.settings(swmm_coupled?swmm_getNumThreads():0,
        swmm_coupled?(couple_clock.every_step?-1.0:couple_clock.interval):0.0);
    while (t_exact < (double)solver.Sim_Time-1.e-9
			&& (!swmm_coupled
					|| (NewRoutingTime < TotalDuration && !ErrorCode)))
	{
		//---------------------------------------------------------------------------------------------------------------------
		// Determine and freeze one common step before calculating exchange.
		// The second ghost update below includes rainfall/coupling changes for
		// the FV1 flux calculation without changing this step.
		solve->update_ghost_cells();
		dynamic_dt.update_dt();

        cuda::dt=(NUMERIC_TYPE)coupling_step(couple_clock,swmm_coupled,t_exact,
            (double)FMIN(cuda::dt,solver.InitTstep),(double)solver.Sim_Time,
            NewRoutingTime/1000.0,TotalDuration/1000.0,rain.time_to_next_frame(t_exact),
            countnum?0.0:(double)pars.hsfdt);
		if (cuda::dt <= C(0.0)) break;
		solver.MinTstep = FMIN(solver.MinTstep, cuda::dt);
		// Use remaining water after net rain for the inlet budget.
		stats_collector.zero_instantaneous_mass();
		if (rain.enabled()) {
			rain.update_time(solver.t);
			lis::cuda::DynamicRain::update<<<grid_size, CUDA_BLOCK_SIZE>>>(
				d_DEM, d_U.H, rain.data(), stats_collector.instantaneous_mass());
		}
		if (states.rainfall == ON) {
			rainfall_uniform.update_rainfall(filenames.rainfilename, solver.t);
			solve->update_uniform_rain(rainfall_uniform.value);
		}

		if (swmm_coupled)
		{
		for (int coupling_i = 0; coupling_i < surface_coupling_count; ++coupling_i)
		{
			const int dem_i = coupling_i;
			verlink->flow_neg_count[dem_i] = 0;
			verlink->flow_2D_dem[dem_i] = 0.0;
		}
		memset(verlink->flow_1D_node, 0, verlink->count_index * sizeof(double));
		memset(verlink->flow_2D_node, 0, verlink->count_index * sizeof(double));
		memset(Node_Inflow, 0, Nobjects[NODE] * sizeof(double));

		const auto cw0=std::chrono::steady_clock::now();
		lis::cuda::gather_coupling_depths(
			d_U.H, temp_BC.PS_idx + temp_numPS_fv1,
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

		if (temp_numPS_fv1 > 0)
		{
			for (size_t flowindex = 0; flowindex < verlink->count_dem_index; flowindex++)
			{
				temp_psvalue[temp_numPS_fv1 + flowindex] = verlink->flow_2D_dem[flowindex];
			}
			for (int flowindex = 0; flowindex < outfall_source_count; ++flowindex)
				temp_psvalue[outfall_source_start + flowindex] = C(0.0);

            if(return_spread.count>0) audit.exchange(swmm_return_apply(&return_spread,
                temp_psvalue+temp_numPS_fv1,(double)geometry.dy,(double)cuda::dt));
            cuda::copy(temp_BC.PS_value+temp_numPS_fv1,temp_psvalue+temp_numPS_fv1,coupling_count*sizeof(NUMERIC_TYPE));
		}
		else
		{
			for (size_t flowindex = 0; flowindex < verlink->count_dem_index; flowindex++)
			{
				temp_psvalue[flowindex] = verlink->flow_2D_dem[flowindex];
			}
			for (int flowindex = 0; flowindex < outfall_source_count; ++flowindex)
				temp_psvalue[surface_coupling_count + flowindex] = C(0.0);

            if(return_spread.count>0) audit.exchange(swmm_return_apply(&return_spread,
                temp_psvalue,(double)geometry.dy,(double)cuda::dt));
            cuda::copy(temp_BC.PS_value,temp_psvalue,coupling_count*sizeof(NUMERIC_TYPE));
		}
        couple_wall+=std::chrono::duration<double>(std::chrono::steady_clock::now()-cw0).count();
		}

		cuda::Boundary::update_time_series(solver.t);
		cuda::Boundary::update_point_sources_fv1(
				d_U.H, d_U.HU, d_U.HV, d_DEM,
				stats_collector.instantaneous_mass(), temp_numPS_fv1);
		// Do not clear instantaneous mass here. CPU LISFLOOD writes Qin/Qout
		// as this step's total domain exchange, including point sources and
		// SWMM coupling. Accumulate once after all fluxes to avoid double
		// counting in Qerror.
		solve->update_ghost_cells();
		// Windows WDDM does not allow host access to __managed__ dt
		// while any kernel is still running.
		cuda::sync();

		//--------------------------------------------------------------------------------
		int swmm_error = 0;
		if (swmm_coupled)
		{
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

                const auto sw0=std::chrono::steady_clock::now();
                swmm_error = swmm_advanceCoupled((double)linkDt, &elapsedTime);
                swmm_wall+=std::chrono::duration<double>(std::chrono::steady_clock::now()-sw0).count();
				if (swmm_error == 0 && coupling_count > 0)
				{
					lis::cuda::Topography::collect_outfall_source_values(
							verlink, temp_psvalue, outfall_source_start, geometry.dy,
							geometry.dx * geometry.dy, tailwater_depths,
							(double)linkDt);
					lis::cuda::Topography::collect_surface_flood_values(
						verlink, temp_psvalue, temp_numPS_fv1, geometry.dy, (double)linkDt);
                    double network_return=0.0, network_inlet=0.0;
                    for(int k=0;k<coupling_count;++k)
                        network_return+=(double)temp_psvalue[temp_numPS_fv1+k]*(double)geometry.dy*(double)linkDt;
                    for(int n=0;n<Nobjects[NODE];++n) network_inlet+=Node_Inflow[n]*(double)linkDt;
                    audit.network_exchange(network_return-network_inlet);
					if (!couple_clock.every_step && cuda::dt > C(0.0))
					{

                        double next=(std::min)(couple_clock.interval,
                            (std::min)((double)solver.Sim_Time,swmm_end)-t_exact-(double)cuda::dt);
                        if(next<1.e-6) next=0.0;
                        swmm_return_schedule(&return_spread,temp_psvalue+temp_numPS_fv1,
                            couple_clock.elapsed,next,(double)geometry.dy,(double)cuda::dt);
					}
					for(int k=0;k<coupling_count;++k)
						audit.exchange(temp_psvalue[temp_numPS_fv1+k]*geometry.dy*(double)cuda::dt);
					cuda::copy(temp_BC.PS_value + temp_numPS_fv1,
							temp_psvalue + temp_numPS_fv1,
							coupling_count * sizeof(NUMERIC_TYPE));
					cuda::Boundary::update_point_sources_range(
							d_U.H, d_U.HU, d_U.HV, d_DEM,
							stats_collector.instantaneous_mass(),
							temp_numPS_fv1);
					// cuda::dt is managed memory. Under Windows WDDM, synchronize
					// before the following host-side timestep access/printing.
					cuda::sync();
				}
				if (!couple_clock.every_step) swmm_couple_advance(&couple_clock);
			}
		}
		//--------------------------------------------------------------------------------

		if (verbose) printf("t=%f\tdt=%f\n", solver.t, cuda::dt);
		t_exact+=(double)cuda::dt;
        solver.t=(NUMERIC_TYPE)t_exact;
		solver.itCount++;

		d_U = solve->update_flow_variables(
				stats_collector.instantaneous_mass());
		if (pars.drain_nodata)
		{
			cuda::Boundary::drain_nodata_water(d_U.H, d_DEM, 
					stats_collector.instantaneous_mass(), grid_size);
		}
		stats_collector.accumulate_mass();
		if (swmm_error != 0) break;

		max_field.update(d_U.H);

		if (stats.need_to_write(solver.t))
		{
			sampler.sample(d_U.H, d_U.HU, d_U.HV, solver.t);
			solve->zero_ghost_cells();
			cuda::sync();

			lis::StatsEntry entry=stats_collector.create_entry(solver,d_U.H);
			stats.write(entry);
			double pending=0;
			if(swmm_coupled && !couple_clock.every_step)
				for(int k=0;k<verlink->count_index;++k) pending+=couple_clock.volume[k];
			audit.sample(entry,stats_collector,pending,swmm_return_pending(&return_spread),swmm_wall,couple_wall);
			sampler.write_if_buffer_full();
		}
		snapshot.write_if_needed(d_U, solver.t, states.call_gzip);
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
	if(filenames.hsffloodend[0]) save_restart(filenames.hsffloodend);
	double final_pending=0;
	if(swmm_coupled && !couple_clock.every_step)
		for(int k=0;k<verlink->count_index;++k) final_pending+=couple_clock.volume[k];
	solve->zero_ghost_cells(); cuda::sync();
	audit.sample(stats_collector.create_entry(solver,d_U.H),stats_collector,final_pending,swmm_return_pending(&return_spread),swmm_wall,couple_wall);
	swmm_return_free(&return_spread);
	swmm_couple_clock_free(&couple_clock);
	if (swmm_coupled)
	{
		swmm_end();
		if (Fout.mode == SCRATCH_FILE) swmm_report();
		// --- close the system
		swmm_close();
		free(Node_Inflow); Node_Inflow = nullptr;
	}
	//---------------------------------------------------------------------------

	time_t loop_end;
	time(&loop_end);
	double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-precise_loop_start).count();
	printf("loop time %lf\n", seconds);
    if(swmm_coupled) audit.timing((double)seconds,swmm_wall,couple_wall,(long)solver.itCount);

	sampler.write();
	snapshot.wait();
	max_field.write();
	audit.finish(solver.t,solver.Sim_Time,ErrorCode);

	if (swmm_coupled) lis::cuda::Topography::free_vconnect(verlink);

	delete[] DEM;
	delete[] manning;

	delete[] coupling_depths;
	delete[] tailwater_depths;
	delete[] temp_psvalue;
	cuda::free_device(d_coupling_depths);
	cuda::free_device(d_tailwater_depths);
	cuda::free_device(d_DEM);
	cuda::free_device(d_manning);
	cuda::free_device(Zstar_x);
	cuda::free_device(Zstar_y);
	Flow::free_pinned(U);
	cuda::Boundary::free_device(boundaries);
}
