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
