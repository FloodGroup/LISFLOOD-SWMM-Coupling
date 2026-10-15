## 1.1.0

- Raise the `.bci` point-boundary capacity from 50,000 to 200,000 unique cells (`maxpi` in `src/input.cpp`), so large urban cases such as Tokyo T01¨CT08 can load without hitting the previous limit.
# Changelog

## 1.0.0 - first public release

Two-way coupling of LISFLOOD-FP and EPA SWMM in a single executable, on CPU and GPU.

- Based on LISFLOOD-FP 8.1.0 (git commit `637106b`) and EPA SWMM 5.1.015.
- Coupled 2-D solvers: ACC and FV1, each on the CPU (OpenMP) or on the GPU (CUDA); the GPU is optional.
- Manhole exchange in both directions: orifice surcharge, weir free inflow, orifice submerged inflow, with supply limits on both sides.
- `COUPLED` outfall type into finite-volume 2-D receiving cells; SWMM node flooding returned to the surface.
- Exchange every 2-D step or at a fixed interval (`coupleint`), with exact time landing.
- GPL-2.0-only. Bundled Intel / MSVC runtimes in binary archives keep their own licences (see `THIRD_PARTY_NOTICES.txt` in the archive).
- Documentation, six coupling figures, and two runnable examples under `examples/`.
