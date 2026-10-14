#pragma once

#include <stdio.h>
#include "AssembledSolution.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__
void update_uniform_rain
(
	AssembledSolution    d_assem_sol,
	Pars pars,
	NUMERIC_TYPE                 dt,
	NUMERIC_TYPE  rain_rate
);





