#include "get_level.cuh"

__device__ int get_level(index_1D idx)
{
	return log( C(3.0) * idx + 1) / log( C(4.0) );
}