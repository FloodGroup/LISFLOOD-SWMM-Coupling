# Example 02: synthetic coupled case

A minimal coupled case for checking that a build reads and maps all coupling inputs, on the CPU or on the GPU. It is synthetic, not a calibrated study site.

```text
     100 m × 100 m flat DEM, 10 m cells, z = 10 m

             J1 (45, 45)        OF1 (55, 45)
   surface ────■──────────────────■────────
               │ invert 9.0 m     │ COUPLED outfall, invert 8.5 m
               └────── C1 ────────┘
```

| File | Role |
|---|---|
| `case.par` | LISFLOOD-FP parameters (`acceleration`) and coupling keywords |
| `dem.asc` | 10 × 10 flat DEM, elevation 10 m |
| `start.asc` | initial depth, 0 everywhere |
| `network.inp` | SWMM network: junction `J1`, `COUPLED` outfall `OF1`, conduit `C1`, `DYNWAVE` routing |
| `floodconfig.txt` | manhole table with one row for `J1` (node index 0) |
| `expected.md` | expected behaviour |

## Run

```bash
lisflood case.par
```

As shipped, the case runs on the CPU and needs no GPU. The same case covers all four coupled drivers:

| Driver | Change to `case.par` |
|---|---|
| CPU ACC | none |
| CPU FV1 | replace `acceleration` with `fv1` |
| CUDA ACC | add a line `cuda` (needs a CUDA build and an NVIDIA GPU) |
| CUDA FV1 | add `cuda` and replace `acceleration` with `fv1` |

To exchange at a fixed interval instead of every step, set `coupleint` to a number of seconds (for example `coupleint 5`).

## Inspect

- `flood_error.txt` has no `out of DEM range` or `wrong node type` line for `J1`, and maps `OF1` to DEM index 55 (`outfall:OF1 cell=55 coupled_stage=1` on the CPU, `COUPLED tailwater: OF1 -> DEM index 55` on the GPU).
- `network.rpt` is the standard SWMM report; its continuity errors should be near zero.
- On the GPU, `results/syn.water_balance.csv` lists every exchange term over time, and `combined_balance_residual_m3` should stay near zero.
