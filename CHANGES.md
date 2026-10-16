# Changelog

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
