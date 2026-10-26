# Changelog

## 1.3.0 - 2026-10-26

- Synchronise fixed-interval coupling and hydraulic restart across CPU ACC, CPU FV1, CUDA ACC and CUDA FV1, retaining each solver's own CFL condition.
- Choose the CPU ACC time step before flux computation. Recompute CUDA ACC's CFL limit from current depth instead of retaining the smallest historical limit.
- Align exchange, rain-frame, restart and end times with equal substeps bounded by the current CFL limit, and accumulate simulation time in double precision.
- Spread positive outfall and node-flooding returns over the following fixed exchange interval. Keep reverse outfall withdrawals immediate and release remaining returns at the end; every-step coupling remains immediate. This changes local water-depth histories and can change peak depths.
- Save and restore surface fields, maxima, SWMM hydraulic history, exchange progress, pending inlet and return volumes, and return rates in one bundled restart file. Retain legacy ACC warm-start reading and the separate surface-only checkpoint mechanism.
- Validate restart format, solver/backend, precision, grid, mapping, exchange settings and record completeness. Reject unrecognised coupled formats instead of falling back to a fresh run. Bundled restart requires the same inputs and build, with hydraulic-only SWMM routing; reject unsupported runoff, RDII, storage exfiltration, routing-interface input, events and active quality configurations before running.
- Compile SWMM C routing sources with OpenMP, balance conduit and node workloads dynamically, and protect the shared flooded-volume statistic.
- Replace repeated receiving-cell scans with hashed cell lookup and precomputed source slots on CPU and CUDA.
- Include both pending-water reservoirs in the CUDA combined water balance, remove internal transfers from external budgets, and log an independent interface residual and restart baseline.
- Report SWMM thread count and measured solver/coupling wall time for both surface solvers. Return a nonzero exit code for SWMM initialisation and restart errors.
- Add reproducible fixed-interval, every-step and restart regressions with positive and negative net rain, alongside the existing wet-coupling and coincident-source checks.
- Update the user guide and source map. Redraw the exchange and timing schematics with signed-flow paths, pending-water reservoirs, delayed-return timelines, CFL alignment and defined units; keep implementation names and runtime messages independent of internal development labels.
- Add named software authors, their institutional affiliation and matching README, CFF and BibTeX citations.

## 1.2.0

- Preserve SWMM's conduit midpoint flow-area history when accounting for coupled storage, avoiding spurious inertial forcing.
- Safeguard the volume-to-head inversion with bisection when Newton progress stalls, and return the iterate with the smallest residual.
- Make the non-OpenMP thread-count fallback file-local so it cannot override the surface solver's OpenMP runtime.
- Validate point-boundary coordinates and time-series name lengths; report and disable out-of-domain points and missing time-series references. For multiple prescribed stages in one cell, keep the last entry.
- Match boundary time series through hash tables instead of repeated full scans.
- Store CUDA ACC boundary time-series data in contiguous device buffers, parallelise interpolation, size boundary launches to the point count, and transfer only coupling source values each step.
- Use atomic accumulation for coincident CUDA ACC discharge sources, correct boundary-kernel strides, and release the allocated ACC boundary buffers.
- Add a reproducible smoke-test script for wet CPU/CUDA coupling and coincident ACC boundary sources.
- Redraw the six documentation figures with fuller labels, typeset equations and a consistent colour scheme, matching the current coupling implementation.

## 1.1.0

- Raise the `.bci` point-boundary capacity from 50,000 to 200,000 entries, so large urban cases can load without hitting the previous limit. Multiple entries may refer to the same cell.

## 1.0.0 - first public release

Two-way coupling of LISFLOOD-FP and EPA SWMM in a single executable, on CPU and GPU.

- Based on LISFLOOD-FP 8.1.0 (git commit `637106b`) and EPA SWMM 5.1.015.
- Coupled 2-D solvers: ACC and FV1, each on the CPU (OpenMP) or on the GPU (CUDA); the GPU is optional.
- Manhole exchange in both directions: orifice surcharge, weir free inflow, orifice submerged inflow, with supply limits on both sides.
- `COUPLED` outfall type into finite-volume 2-D receiving cells; SWMM node flooding returned to the surface.
- Exchange every 2-D step or at a fixed interval (`coupleint`), with exact time landing.
- GPL-2.0-only. Bundled Intel / MSVC runtimes in binary archives keep their own licences (see `THIRD_PARTY_NOTICES.txt` in the archive).
- Documentation, six coupling figures, and two runnable examples under `examples/`.
