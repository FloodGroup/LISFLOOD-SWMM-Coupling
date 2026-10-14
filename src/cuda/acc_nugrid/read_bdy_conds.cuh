#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Boundary.h"
#include "gen_bdy_morton_codes.cuh"
#include "../../lisflood.h"

void read_bdy_conds
(
	const char*                 bcifilename,
	const int                   direction,
	Boundary&                   boundary,
	const Pars&                 pars
);