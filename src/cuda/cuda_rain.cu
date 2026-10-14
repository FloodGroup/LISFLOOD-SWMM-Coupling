#include "cuda_rain.cuh"
#include "cuda_atomic.cuh"
#include "../rain/rain.h"
#include "cuda_unifiedallocator.cuh"

#if _NETCDF == 1
	#include "../rain/rain.tpp"
#else
	#include "../rain/rain_stub.tpp"
#endif

template class DynamicRain<lis::cuda::UnifiedAllocator<NUMERIC_TYPE>>;

__device__ 
NUMERIC_TYPE lis::cuda::rate_at_cell
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j
)
{
    int tile_i = (i-1) * cuda::geometry.dx / cuda::rain_geometry.dx;

    NUMERIC_TYPE top_gap = cuda::geometry.tly - cuda::rain_geometry.tly;
    int tile_j = (j-1 - top_gap/cuda::geometry.dy) *
		cuda::geometry.dy / cuda::rain_geometry.dy;

    if (tile_i >= 0 && tile_i < cuda::rain_geometry.xsz &&
            tile_j >= 0 && tile_j < cuda::rain_geometry.ysz)
    {
        return rain_data[tile_j*cuda::rain_geometry.xsz + tile_i];
    }
    else
    {
        return C(0.0);
    }
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
lis::cuda::DynamicRain::update
	(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* rain_data,
	MassStats* mass_stats
)
{
#if _NETCDF == 1
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j+1; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i+1; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z = DEM[j*cuda::pitch + i];
            if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6)) continue;

			NUMERIC_TYPE& Hval = H[j*cuda::pitch + i];
			NUMERIC_TYPE rain_rate = rate_at_cell(rain_data, i, j);
			if (!isfinite(rain_rate) || !isfinite(Hval)) {
				if (atomicAdd(&mass_stats->invalid_count,1U)==0)
					printf("V11 ERROR rain/water nonfinite at valid FV1 cell x=%d y=%d H=%g rate=%g; check rain frame log\n",i,j,(double)Hval,(double)rain_rate);
				asm("trap;");
			}
			NUMERIC_TYPE H_old = FMAX(Hval, C(0.0));
			Hval = H_old;
			Hval = FMAX(C(0.0), H_old + rain_rate * cuda::dt);
			if (mass_stats != nullptr)
			{
				// Report the water actually added/removed.  A negative
				// rainfall value remains an infiltration demand, but it
				// cannot remove more water than the cell currently stores.
				NUMERIC_TYPE discharge = (Hval - H_old)
					/ cuda::dt * cuda::geometry.dx * cuda::geometry.dy;
				if (discharge >= C(0.0)) {
					atomicAdd(&(mass_stats->in), discharge);
					atomicAdd(&(mass_stats->rain_in), discharge);
				} else {
					atomicAdd(&(mass_stats->out), -discharge);
					atomicAdd(&(mass_stats->rain_out), -discharge);
				}
			}
		}
	}
#endif
}


__device__ 
NUMERIC_TYPE lis::cuda::rate_at_cell_ACC
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j
)
{
    int tile_i = (i) * cuda::geometry.dx / cuda::rain_geometry.dx;

    NUMERIC_TYPE top_gap = cuda::geometry.tly - cuda::rain_geometry.tly;
    int tile_j = (j - top_gap/cuda::geometry.dy) *
		cuda::geometry.dy / cuda::rain_geometry.dy;

    if (tile_i >= 0 && tile_i < cuda::rain_geometry.xsz &&
            tile_j >= 0 && tile_j < cuda::rain_geometry.ysz)
    {
        return rain_data[tile_j*cuda::rain_geometry.xsz + tile_i];
    }
    else
    {
        return C(0.0);
    }
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
lis::cuda::DynamicRain::updateACC
(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* rain_data,
	MassStats* mass_stats
)
{
#if _NETCDF == 1
	__shared__ NUMERIC_TYPE block_in[CUDA_BLOCK_SIZE];
	__shared__ NUMERIC_TYPE block_out[CUDA_BLOCK_SIZE];

	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;
	int local_idx = threadIdx.y * blockDim.x + threadIdx.x;
	NUMERIC_TYPE local_in = C(0.0);
	NUMERIC_TYPE local_out = C(0.0);

	for (int j=global_j; j<cuda::geometry.ysz; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::geometry.xsz; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z = DEM[j*cuda::geometry.xsz + i];
            if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6)) continue;

			NUMERIC_TYPE& Hval = H[j*cuda::geometry.xsz + i];
			NUMERIC_TYPE rain_rate = rate_at_cell_ACC(rain_data, i, j);
			if (!isfinite(rain_rate) || !isfinite(Hval)) {
				if (atomicAdd(&mass_stats->invalid_count,1U)==0)
					printf("V11 ERROR rain/water nonfinite at valid ACC cell x=%d y=%d H=%g rate=%g; check rain frame log\n",i,j,(double)Hval,(double)rain_rate);
				asm("trap;");
			}
			NUMERIC_TYPE H_old = FMAX(Hval, C(0.0));
			Hval = H_old;
			Hval = FMAX(C(0.0), H_old
				+ rain_rate * cuda::dt);

			// Accumulate the water actually added or removed.  Negative rain
			// remains an infiltration demand, limited only by available water.
			NUMERIC_TYPE discharge = (Hval - H_old)
				/ cuda::dt * cuda::geometry.dx * cuda::geometry.dy;
			if (discharge >= C(0.0)) local_in += discharge;
			else local_out -= discharge;
		}
	}

	block_in[local_idx] = local_in;
	block_out[local_idx] = local_out;
	__syncthreads();

	for (int stride = CUDA_BLOCK_SIZE / 2; stride > 0; stride >>= 1)
	{
		if (local_idx < stride)
		{
			block_in[local_idx] += block_in[local_idx + stride];
			block_out[local_idx] += block_out[local_idx + stride];
		}
		__syncthreads();
	}

	if (local_idx == 0 && mass_stats != nullptr)
	{
		atomicAdd(&(mass_stats->in), block_in[0]);
		atomicAdd(&(mass_stats->out), block_out[0]);
		atomicAdd(&(mass_stats->rain_in), block_in[0]);
		atomicAdd(&(mass_stats->rain_out), block_out[0]);
	}
#endif
}
