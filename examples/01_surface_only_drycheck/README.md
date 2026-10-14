# Example 01 鈥?surface-only dry check

Tiny DEM/start pair adapted from the in-tree CUDA ACC dry-check inputs.
No SWMM coupling. Useful as a smoke test that the executable starts and writes mass/output files.

## Run

```
lisflood case.par
```

## Expected

- Completes `sim_time = 1 s`
- Creates `results/` with mass and depth outputs prefixed `drycheck`
