# Example 02 — synthetic coupled case (**synthetic**)

This case is **constructed for the open-source layout**. It is not a calibrated study site.

## Files

| File | Role |
|---|---|
| `case.par` | LISFLOOD-FP controls + `swmminp` / `floodconfig` / `coupleint` |
| `dem.asc` | 10×10 flat DEM (10 m cells), elevation 10 m |
| `start.asc` | Initial depth 0 |
| `network.inp` | SWMM network: junction J1, COUPLED outfall OF1, one conduit |
| `floodconfig.txt` | Vertical connection table for J1 |
| `expected.md` | Qualitative expected behaviour |

## Run

```text
lisflood case.par
```

Requires a build with SWMM coupling enabled (default in this repository). CUDA is optional for this tiny grid; CPU acceleration path may be used depending on flags in `case.par`.

## Note

If your binary was built with CUDA-only coupling intervals, keep `coupleint step` as shipped.
