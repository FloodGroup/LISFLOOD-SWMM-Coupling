#include "cuda_stats.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include <cub/cub.cuh>
#include <limits>
#include <stdexcept>

namespace lis
{
namespace cuda
{

__global__ void accumulate_mass
(
	MassStats* cumulative,
	MassStats* instantaneous,
	double* audit
)
{
	cumulative->in += instantaneous->in * cuda::dt;
	cumulative->out += instantaneous->out * cuda::dt;
	// Separate long-run diagnostic totals, in double precision. Hydraulics unchanged.
	const NUMERIC_TYPE rates[9] = {instantaneous->in,instantaneous->out,
		instantaneous->rain_in,instantaneous->rain_out,instantaneous->nodata_out,
		instantaneous->inlet_in,instantaneous->inlet_out,instantaneous->outfall_in,instantaneous->outfall_out};
	for (int k=0;k<9;++k) audit[k] += (double)rates[k]*(double)cuda::dt;
}

}
}

lis::cuda::StatsCollector::StatsCollector
(
	Geometry& geometry,
	NUMERIC_TYPE DepthThresh,
	int acceleration
)
:
geometry(geometry),
elements((acceleration == 0) ? lis::GhostRaster::elements(geometry) : lis::GhostRaster::elements_H(geometry))
{
	area_temp = nullptr;
	volume_temp = nullptr;
	NUMERIC_TYPE* dummy_in = nullptr;
	NUMERIC_TYPE* dummy_reduce_out = nullptr;
	int* dummy_histogram_out = nullptr;

	histogram_levels = static_cast<NUMERIC_TYPE*>(
			malloc_unified(3*sizeof(NUMERIC_TYPE)));

	histogram_levels[0] = C(0.0);
	histogram_levels[1] = DepthThresh;
	histogram_levels[2] = std::numeric_limits<NUMERIC_TYPE>::max();

	histogram_counts = static_cast<int*>(malloc_unified(2*sizeof(int)));

	d_histogram_levels = static_cast<NUMERIC_TYPE*>(
			malloc_device(3*sizeof(NUMERIC_TYPE)));
	cuda::copy(d_histogram_levels, histogram_levels, 3*sizeof(NUMERIC_TYPE));
	d_histogram_counts = static_cast<int*>(malloc_device(2*sizeof(int)));

	checkCudaErrors(cub::DeviceHistogram::HistogramRange(area_temp, area_bytes,
				dummy_in, dummy_histogram_out, 3, histogram_levels, elements));
	area_temp = malloc_device(area_bytes);

	checkCudaErrors(cub::DeviceReduce::Sum(volume_temp, volume_bytes,
				dummy_in, dummy_reduce_out, elements));
	volume_temp = malloc_device(volume_bytes);

	d_volume = static_cast<NUMERIC_TYPE*>(malloc_device(sizeof(NUMERIC_TYPE)));

	d_instantaneous_mass = static_cast<MassStats*>(
			malloc_device(sizeof(MassStats)));
	d_cumulative_mass = static_cast<MassStats*>(
			malloc_device(sizeof(MassStats)));
	zero_cumulative_mass();
	d_audit_totals = static_cast<double*>(malloc_device(9*sizeof(double)));
	double audit_zero[9] = {};
	cuda::copy(d_audit_totals,audit_zero,sizeof(audit_zero));
}

void lis::cuda::StatsCollector::zero_instantaneous_mass()
{
	MassStats temp = {};
	cuda::copy(d_instantaneous_mass, &temp, sizeof(MassStats));
}

lis::MassStats* lis::cuda::StatsCollector::instantaneous_mass()
{
	return d_instantaneous_mass;
}

void lis::cuda::StatsCollector::accumulate_mass()
{
	lis::cuda::accumulate_mass<<<1, 1>>>(d_cumulative_mass,
			d_instantaneous_mass, d_audit_totals);
	cuda::sync();
	cumulative_time += cuda::dt;
}

void lis::cuda::StatsCollector::audit_totals(double* totals)
{
	cuda::copy(totals,d_audit_totals,9*sizeof(double));
}

void lis::cuda::StatsCollector::initialise_volume(NUMERIC_TYPE* H)
{
	previous_volume = volume(H);
}

lis::StatsEntry lis::cuda::StatsCollector::create_entry
(
	::Solver& solver,
	NUMERIC_TYPE* H
)
{
	NUMERIC_TYPE current_volume = volume(H);
	if (!std::isfinite(current_volume)) {
		fprintf(stderr,"V11 ERROR nonfinite 2D storage at time_s=%.12g\n",(double)solver.t);
		throw std::runtime_error("V11 nonfinite 2D storage");
	}
	NUMERIC_TYPE volume_change = current_volume - previous_volume;
	previous_volume = current_volume;

	MassStats h_instantaneous_mass;
	cuda::copy(&h_instantaneous_mass, d_instantaneous_mass, sizeof(MassStats));

	MassStats h_cumulative_mass;
	cuda::copy(&h_cumulative_mass, d_cumulative_mass, sizeof(MassStats));

	NUMERIC_TYPE time = cumulative_time;

	NUMERIC_TYPE net_source_flux =
		h_cumulative_mass.in - h_cumulative_mass.out;

	zero_cumulative_mass();

	return { 
		solver.t,
		cuda::dt,
		solver.itCount,
		solver.MinTstep,
		area(H),
		current_volume,
		h_instantaneous_mass,
		net_source_flux - volume_change,
		time > C(0.0) ? (net_source_flux - volume_change)/time : C(0.0)
	};
}

lis::cuda::StatsCollector::~StatsCollector()
{
	free_unified(histogram_levels);
	free_unified(histogram_counts);
	free_device(d_histogram_levels);
	free_device(d_histogram_counts);
	free_device(area_temp);
	free_device(volume_temp);
	free_device(d_volume);
	free_device(d_instantaneous_mass);
	free_device(d_cumulative_mass);
	free_device(d_audit_totals);
}

NUMERIC_TYPE lis::cuda::StatsCollector::area
(
	NUMERIC_TYPE* H
)
{
	checkCudaErrors(cub::DeviceHistogram::HistogramRange(area_temp, area_bytes,
				H, d_histogram_counts, 3, d_histogram_levels, elements));
	cuda::copy(histogram_counts, d_histogram_counts, 2*sizeof(int));
	return histogram_counts[1] * geometry.dx * geometry.dy;
}

NUMERIC_TYPE lis::cuda::StatsCollector::volume
(
	NUMERIC_TYPE* H
)
{

	checkCudaErrors(cub::DeviceReduce::Sum(volume_temp, volume_bytes,
				H, d_volume, elements));

	NUMERIC_TYPE h_volume;

	cuda::copy(&h_volume, d_volume, sizeof(NUMERIC_TYPE));

	return h_volume * geometry.dx * geometry.dy;
}

void lis::cuda::StatsCollector::zero_cumulative_mass()
{
	MassStats temp_stats = {};
	cuda::copy(d_cumulative_mass, &temp_stats, sizeof(MassStats));

	cumulative_time = C(0.0);
}
