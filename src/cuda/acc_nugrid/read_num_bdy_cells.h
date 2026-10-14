#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Directions.h"
#include "InletTypes.h"
#include "../../lisflood.h"

int read_num_bdy_cells
(
	const char* bcifilename,
	const Pars& pars,
	const int                   direction
);