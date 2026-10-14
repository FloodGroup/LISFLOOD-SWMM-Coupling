#pragma once

#include "cuda_runtime.h"
#include "cub/cub.cuh"
#include "cuda_utils.cuh"
#include "AssembledSolution.h"
#include "MortonCode.h"

void reverse_z_order_act_idcs
(
	MortonCode*       d_reverse_z_order,
	MortonCode*       d_indices,
	AssembledSolution d_buf_sol,
	AssembledSolution d_sol,
	int               array_length
);