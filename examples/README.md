# Examples

Each subdirectory is a small, self-contained case (< 5 MB).

| Case | Type | Notes |
|---|---|---|
| `01_surface_only_drycheck` | Surface only | Tiny DEM from the upstream-style dry-check suite; no SWMM |
| `02_coupled_synthetic` | **Synthetic** coupled | Minimal DEM + SWMM `.inp` with one `COUPLED` outfall + one junction inlet; for layout demonstration |

## How to run

1. Build `lisflood` (see `docs/user_guide.md`).
2. Copy the executable (and NetCDF DLLs on Windows) next to the case, or put them on `PATH`.
3. From the case directory:

```text
lisflood case.par
```

Expected outputs are described in each case `README.md`.
