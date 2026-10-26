# Short regression validation

Validated on 2026-10-26 on Windows x64 with MSVC 14.39, CUDA 12.8 and a Release build (surface float, SWMM double). C and C++ sources were built with OpenMP. The CUDA binary includes native sm_120 and sm_75/compute_75 compatibility code.

The executable was rebuilt from this public source tree. Its SHA-256 is:

```text
9773cef45542e26e0abd18326833cbf79dad7e5897b460992313d86e189c825b
```

All 34 checks passed: 18 successful simulations and 16 intentional restart rejections. Tests used the staged executable and its nine runtime DLLs, with the child-process PATH limited to that package and Windows system directories. CUDA runs confirmed actual GPU selection. Public runtime messages contain no internal development labels.

## Coupling and boundary smoke checks

`tools/test_coupling_smoke.py` runs four wet coupled simulations and two boundary checks. Every successful run reached its requested end time and produced finite, nonnegative depth rasters.

| Test | Simulated time | Result |
|---|---:|---|
| Wet CPU ACC | 60 s | Passed; final surface volume 994.9329965 m3 |
| Wet CPU FV1 | 60 s | Passed; final surface volume 994.9328549 m3 |
| Wet CUDA ACC | 60 s | Passed; maximum absolute combined residual 0.0000868568 m3 |
| Wet CUDA FV1 | 60 s | Passed; maximum absolute combined residual 0.000279445 m3 |
| CPU ACC coincident boundaries | 2 s | Passed; injected volume 0.6000000 m3 |
| CUDA ACC coincident boundaries | 2 s | Passed; injected volume 0.6000001 m3 |

The wet cases start with 0.1 m surface water. CUDA checks require nonzero manhole and outfall exchange. Boundary checks combine two discharge sources in one cell with an out-of-domain point and an unresolved time-series name; the invalid entries must be reported and disabled. Expected injected volume is 0.6 m3.

## Fixed-interval coupling and restart

`tools/test_coupled_restart.py` runs each solver/backend continuously, resumes it from a saved state, and separately checks every-step exchange: 12 successful simulations. The fixture has a 10 x 10 grid, 10 m cells, 16 independent conduits with shared coupling cells, and four actual SWMM routing threads. Positive and negative NetCDF net-rain intervals are included.

The 12-second event uses a 2-second exchange interval. A restart at 5.25 seconds retains 1.25 seconds of exchange progress, with both pending-water reservoirs nonzero.

| Solver/backend | Saved pending to 1-D (m3) | Saved pending to 2-D (m3) | Resumed final-depth difference | Resumed maximum-depth difference |
|---|---:|---:|---:|---:|
| CPU ACC | 10.558336608 | 15.094086528 | 0 m | 0 m |
| CPU FV1 | 7.747125526 | 14.750303328 | 0 m | 0 m |
| CUDA ACC | 10.558422183 | 15.095468760 | 0 m | 0 m |
| CUDA FV1 | 7.747125526 | 14.750300646 | 0 m | 0 m |

Saved and restored queue diagnostics match exactly at the printed precision. Every-step runs also reach 12 seconds. Across these CUDA fixed-interval, resumed and every-step runs, the largest absolute combined residual is 0.000245005 m3 and the largest independent interface residual is 0.0000245722 m3. Final pending return volume is zero.

Each of the four solver/backend paths also rejects a wrong exchange interval, wrong solver, truncated file and unsupported coupled-file signature. All 16 rejection cases return exit code 1, including CUDA paths, instead of completing from a fresh state.

## Reproduction and limits

From the repository root, after building:

```bash
python tools/test_coupling_smoke.py --exe build/src/lisflood --output build/smoke-check --cuda
python tools/test_coupled_restart.py --exe build/src/lisflood --output build/restart-check --cuda --rain
```

On Windows, use `build/release/lisflood.exe` or the packaged executable. Omit `--cuda` for CPU-only builds; `--rain` requires the Python `netCDF4` package. Output directories must not already exist. Each script retains inputs, console logs, outputs and `validation.json`; the binary archive also includes the final validation records and a source/executable manifest.

These tests establish short-run coupling and same-build restart consistency. They do not establish full-event stability, field accuracy, Linux multi-thread correctness, cross-platform restart portability or production speed. The two regression scripts were run concurrently, so their wall times are not performance benchmarks. No full production event was run. See the [hotstart scope](user_guide.md#8-hotstart) and [exchange-interval limitations](coupling_methodology.md#8-time-synchronisation).
