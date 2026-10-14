#pragma once

#include "../lisflood.h"

namespace lis
{
namespace cuda
{

// Gather only the 2D depths used by SWMM coupling. The index and output
// arrays reside on the GPU, so only this compact buffer is copied to the CPU.
void gather_coupling_depths
(
	const NUMERIC_TYPE* H,
	const int* coupling_indices,
	NUMERIC_TYPE* coupling_depths,
	int coupling_count
);

}
}
