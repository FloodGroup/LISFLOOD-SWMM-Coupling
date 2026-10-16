# Short regression validation

Validated on 2026-10-16 on Windows with MSVC 14.39, CUDA 12.8 and a Release build (surface float, SWMM double). CPU tests used four OpenMP threads. CUDA tests confirmed GPU selection in the console log.

The executable was built from this public source tree. Its SHA-256 was:

```text
b5afa750a93575a451b291462696b67001ad5deb142814e7012de5f07d396ae4
```

## Results

All six tests in `tools/test_coupling_smoke.py` completed with exit code zero, reached their requested simulation time and produced finite, nonnegative depth outputs.

| Test | Simulated time | Result |
|---|---:|---|
| Wet CPU ACC coupling | 60 s | Passed; final surface volume 994.9329977 m3 |
| Wet CPU FV1 coupling | 60 s | Passed; final surface volume 994.9328580 m3 |
| Wet CUDA ACC coupling | 60 s | Passed; maximum absolute combined balance residual 0.000105808 m3 |
| Wet CUDA FV1 coupling | 60 s | Passed; maximum absolute combined balance residual 0.000280992 m3 |
| CPU ACC boundary handling | 2 s | Passed; injected volume 0.5999999 m3 versus 0.6 m3 expected |
| CUDA ACC boundary handling | 2 s | Passed; injected volume 0.6 m3 versus 0.6 m3 expected |

The wet cases start with 0.1 m of surface water. GPU audit checks require nonzero surface-to-manhole and outfall exchange. The boundary cases combine two time-varying discharge points in one cell, an out-of-domain point and an unresolved time-series name; the latter two are reported and disabled.

All 34 depth snapshots also matched the supplied reference executable's output numerically, with a maximum difference of 0 m at the saved precision. Public startup logs contained no internal engine revision labels.

## Reproduction and limits

See the [user guide](user_guide.md#short-regression-checks) for commands. The script creates a new output directory containing inputs, console logs, model outputs and `validation.json`.

These short tests check the source synchronisation and basic runtime behaviour. They do not establish full-event stability, large-network accuracy, Linux multi-thread correctness or performance on other GPUs. No full production event was run for this validation.
