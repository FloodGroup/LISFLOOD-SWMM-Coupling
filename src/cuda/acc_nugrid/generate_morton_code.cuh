#pragma once

#include "cuda_runtime.h"
#include "dilate.cuh"

__host__ __device__ __forceinline__ MortonCode generate_morton_code
(
	int x,
	int y
)
{
	// please grep 'BITWISE OPERATORS AND HEXADECIMAL' for explanation
	return dilate(x) | (dilate(y) << 1);
}