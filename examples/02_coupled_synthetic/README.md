# Example 02: synthetic coupled case

A minimal coupled case for checking that a build reads and maps all coupling inputs. It is synthetic, not a calibrated study site.

```text
     100 m × 100 m flat DEM, 10 m cells, z = 10 m

             J1 (45, 45)        OF1 (55, 45)
   surface ────■──────────────────■────────
               │ invert 9.0 m     │ COUPLED outfall, invert 8.5 m
               └────── C1 ────────┘
```

| File | Role |
|---|---|
| `case.par` | LISFLOOD-FP parameters (`cuda`, `acceleration`) and coupling keywords |
| `dem.asc` | 10 × 10 flat DEM, elevation 10 m |
| `start.asc` | initial depth, 0 everywhere |
| `network.inp` | SWMM network: junction `J1`, `COUPLED` outfall `OF1`, conduit `C1`, `DYNWAVE` routing |
| `floodconfig.txt` | manhole table with one row for `J1` (node index 0) |
| `expected.md` | expected behaviour |

## Run

Requires a CUDA build and an NVIDIA GPU, because coupled runs use the CUDA solvers.

```bash
lisflood case.par
```

To try the FV1 solver, replace `acceleration` with `fv1` in `case.par`. To exchange at a fixed interval, set `coupleint` to a number of seconds (for example `coupleint 5`).

## Inspect

- `flood_error.txt` has no `out of DEM range` line for `J1`, and has one `V8 COUPLED tailwater: OF1 -> DEM index …` line.
- `results/syn.water_balance.csv` lists every exchange term over time; `combined_balance_residual_m3` should stay near zero.
- `network.rpt` is the standard SWMM report.
