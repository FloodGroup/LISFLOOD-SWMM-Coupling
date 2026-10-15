#pragma once

#include "geometry.h"

namespace lis { namespace cpu_coupling {

void map_internal_outfalls(Fnames* filenames, Vconnect* coupling,
    Pars* pars, Arrays* arrays, int verbose);

void free_outfalls(Vconnect* coupling);

void update_outfall_stages(Vconnect* coupling, Pars* pars, Arrays* arrays);

void collect_return_sources(Vconnect* coupling, BoundCs* boundaries,
    Pars* pars, Arrays* arrays, int surface_source_start,
    int outfall_source_start, double exchange_seconds,
    double application_seconds);

void apply_coupled_sources(BoundCs* boundaries, Pars* pars, Solver* solver,
    Arrays* arrays, int source_start, int source_end);

} }
