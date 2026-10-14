# CHANGES — modifications relative to upstream

Generated for the public release of the LISFLOOD-FP × SWMM coupling (internal label V11).

## Contents

- [Upstream baselines](#upstream-baselines)
- [Summary](#summary)
- [Added coupling sources](#added-coupling-sources)
- [Modified areas (high level)](#modified-areas-high-level)
- [Licence obligations](#licence-obligations)

---

## Upstream baselines

| Component | Baseline | Licence |
|---|---|---|
| LISFLOOD-FP | Development commit **`637106b`** (2022-04-20), reported as **8.1.0** in `VersionHistory.h`. Not the Zenodo 8.1/8.2 zip contents. | GPL-2.0 ([Zenodo 8.1](https://zenodo.org/records/6912932)) |
| EPA SWMM | **5.1.015** (`VERSION 51015`) | Public domain |
| NVIDIA CUB | Bundled under `src/cuda/cub/`, unmodified | BSD-3-Clause |

---

## Summary

| Metric | Count |
|---|---|
| **Modified** | 64 files (47 LISFLOOD-FP, 17 SWMM) |
| **Added** | Core coupling sources listed below (development changelogs were removed from this release tree) |

---

## Added coupling sources

| File | Role |
|---|---|
| `src/coupled_boundary.c` / `.h` | Finite-volume COUPLED receiving cells, flux limits, exchange volumes |
| `src/swmm_couple.h` | Exchange clock (`coupleint`) |
| `src/cuda/cuda_coupling_gather.cuh` | Sparse GPU depth gather for coupling cells |

---

## Modified areas (high level)

| Area | Files (examples) | Change |
|---|---|---|
| Build | `CMakeLists.txt`, `config.default.cmake` | Compile SWMM + coupling into `lisflood` |
| Parameters | `pars.cpp`, `lisflood.h` | `swmminp` / `swmmrpt` / `swmmout` / `coupleint` / `floodconfig` |
| CPU loops | `iterateq.cpp`, `swe/fv1.cpp`, `swe/dg2.cpp` | Vertical connection orifice exchange; SWMM advance |
| CUDA loops | `cuda/acc/cuda_acc_simulate.cu`, `cuda/fv1/cuda_fv1_simulate.cu`, `cuda/cuda_dem.cu`, `cuda/cuda_boundary.cu` | Bidirectional exchange, mapping, point sources, audits |
| SWMM core | `swmm5.c` / `.h`, `dynwave.c`, `node.c`, `routing.c`, `enums.h`, `keywords.c`, `text.h`, … | `COUPLED` outfall type; coupling API; volume hooks |
| Hotstart | `chkpnt.cpp`, `hotstart.c` | Coupled surface / SWMM checkpoints |

> A full per-file `+/-` line table from the automated upstream diff is available from the maintainers on request; this release document intentionally stays concise.

---

## Licence obligations

- Entire repository: **GPL-2.0** (because LISFLOOD-FP is GPL-2.0 and was modified).
- Keep `LICENSE`, `LICENCE-3RD-PARTY.txt`, and upstream copyright notices.
- Do not commit Microsoft / Intel redistributable DLLs into git; ship them only in a separate binary package with vendor notices.
