#pragma once

#include "generate_morton_code.cuh"
#include "Directions.h"
#include "Boundary.h"

void gen_bdy_morton_codes
(
	const Boundary& boundary,
	const Pars& pars,
	const int& direction
);