# Examples

Two small, self-contained cases for checking a build.

| Case | Solver | What it checks |
|---|---|---|
| [`01_surface_only_drycheck`](01_surface_only_drycheck/) | CPU ACC, no SWMM | the executable starts and writes output |
| [`02_coupled_synthetic`](02_coupled_synthetic/) | CPU ACC + SWMM; switch to FV1 or add `cuda` for the GPU | coupled input wiring: one junction inlet, one `COUPLED` outfall |

Both cases run without a GPU.

## Running a case

1. Build `lisflood` ([user guide §2](../docs/user_guide.md#2-build)).
2. Put the executable on your `PATH`, or copy it (and on Windows its DLLs) into the case folder.
3. Run from the case folder so that relative paths resolve:

```bash
cd examples/02_coupled_synthetic
lisflood case.par
```

Results are written to `results/` inside the case folder. Each case README lists what to expect.
