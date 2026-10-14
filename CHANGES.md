# Changelog

## 1.0.0 — first public release

Two-way coupling of LISFLOOD-FP and EPA SWMM in a single Windows executable.

- Based on LISFLOOD-FP 8.1.0 (git commit `637106b`) and EPA SWMM 5.1.015.
- Manhole exchange (orifice / weir) in both directions; COUPLED outfall type into 2-D receiving cells; SWMM node flooding returned to the surface.
- CUDA ACC and CUDA FV1 coupled runs; GPL-2.0-only. Bundled Intel/MSVC runtimes retain their own licences (see THIRD_PARTY_NOTICES if present).
- Documentation and six coupling figures; synthetic examples under `examples/`.
