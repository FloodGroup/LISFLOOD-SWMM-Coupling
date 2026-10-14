# Example 01: surface-only dry check

A tiny DEM and initial-depth pair adapted from the upstream LISFLOOD-FP dry-check inputs. There is no SWMM network, so this case runs on any build, with or without a GPU.

| File | Role |
|---|---|
| `case.par` | LISFLOOD-FP parameters, CPU ACC solver |
| `dem.dem` | DEM |
| `start.start` | initial water depth |

## Run

```bash
lisflood case.par
```

## Expected result

- The run completes `sim_time = 1 s` without errors.
- `results/` contains the mass-balance and depth files prefixed `drycheck`.
