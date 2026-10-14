# User guide

## 1. Obtain the code

Clone or download this repository. All compilable sources are under `src/`. Documentation is under `docs/`. Small runnable cases are under `examples/`.

Optional Windows binaries (not in git) are packaged separately for GitHub Releases; see the repository README.

## 2. Build from source

### Prerequisites

| Dependency | Notes |
|---|---|
| CMake ≥ 3.13 | |
| MSVC 2019/2022 or Linux GCC/Clang | C++14 |
| CUDA Toolkit | Optional; needed for GPU ACC / FV1 coupling |
| NetCDF + HDF5 | Windows: headers/libs under `src/windep/` (bundled). Linux: install `libnetcdf-dev` |
| OpenMP | Usually with the compiler |

### Windows (script)

Edit paths, then:

```powershell
.\tools\build.ps1 `
  -VcVars "C:\BuildTools\VC\Auxiliary\Build\vcvars64.bat" `
  -CMake "C:\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" `
  -CudaArchitecture 86 `
  -RedistDirectory "C:\BuildTools\VC\Redist\MSVC\<version>\x64"
```

The script configures NMake Makefiles in a build directory, compiles target `lisflood`, and copies `lisflood.exe` plus NetCDF/HDF5/zlib DLLs to the output directory. MSVC CRT / OpenMP redistributables are copied from `-RedistDirectory` when provided.

### Windows (Visual Studio CMake)

Open the **repository root** (or `src/`) as a CMake folder, select a Release x64 preset, build target `lisflood`.

### Linux

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=75
cmake --build build --target lisflood -j
```

## 3. Run a simulation

```text
lisflood path/to/case.par
```

Working directory is typically the case folder so relative paths in `.par` resolve.

### Essential `.par` keywords

| Keyword | Example | Meaning |
|---|---|---|
| `DEMfile` | `dem.asc` | ASCII DEM |
| `resroot` | `out` | Output name prefix |
| `dirroot` | `results` | Output directory |
| `sim_time` | `3600` | Simulation length (s) |
| `acceleration` / `cuda` / `fv1` | flags | Solver selection |
| `swmminp` | `network.inp` | SWMM input |
| `swmmrpt` | `network.rpt` | SWMM report |
| `swmmout` | `network.out` | SWMM binary out |
| `coupleint` | `step` or `1.0` | Exchange every 2-D step, or every *N* seconds |

Without `swmminp`, the model runs as standalone LISFLOOD-FP (no sewer coupling).

Upstream LISFLOOD-FP parameters (friction, boundary files, rainfall, …) follow the official manuals:

- LISFLOOD-FP 8.1 on Zenodo: <https://zenodo.org/records/6912932>
- EPA SWMM: <https://www.epa.gov/water-research/storm-water-management-model-swmm>

## 4. SWMM input for coupling

### `COUPLED` outfalls

In the `[OUTFALLS]` section, use type **`COUPLED`** instead of `FREE` / `NORMAL` / `FIXED` / … for nodes that exchange volume with the 2-D surface. The coupler maps each such outfall to a DEM cell and updates an external stage from the cell water surface (`coupled_set_cell` in `coupled_boundary.c`).

Coupled runs require **continuous dynamic-wave** routing (`ROUTING DYNWAVE`), no routing `EVENTS`, and hydraulic-only routing when pollutants would otherwise be active (`coupled_open`).

### Manhole / inlet table (`floodconfig`)

Vertical connections (surface ↔ junction) are listed in a text table named by the `.par` keyword `floodconfig`. The CPU path reads it with `flood::verticalcon::read_vconinfo` (`iterateq.cpp`) and the CUDA path with `lis::cuda::Topography::read_floodconfig` (`cuda_dem.cu`); both use the same format:

- line 1: `node_num <N>`;
- line 2: a header of exactly **10** words (skipped by the reader, so the names are free but the word count must be 10);
- then *N* rows of 10 fields: an ID string followed by 9 numbers.

```text
node_num      <N>
ID      Index      Orifice_QCoef      Pond_Area      Max_Qout      Weir_QCoef      Length_W      Max_Qin      Orifice_QCoef_In      Orifice_Area_In
```

| # | Column | Code variable | Role |
|---|---|---|---|
| 1 | `ID` | (skipped) | Row label; not used by the code |
| 2 | `Index` | `node_Vindex` | 0-based index of the SWMM node in the `.inp` node list |
| 3 | `Orifice_QCoef` | `node_c0` | Orifice coefficient $C_0$ for sewer → surface flow |
| 4 | `Pond_Area` | `node_area` | Orifice area $A_0$ (m²) for sewer → surface flow. Despite the name, this is **not** SWMM's `pondedArea` |
| 5 | `Max_Qout` | `node_qoutmax` | Cap on sewer → surface discharge (m³/s) |
| 6 | `Weir_QCoef` | `node_cw` | Weir coefficient $C_w$; read but currently unused (weir branch commented out) |
| 7 | `Length_W` | `node_w` | Weir crest length (m); read but currently unused |
| 8 | `Max_Qin` | `node_qinmax` | Cap on surface → sewer discharge (m³/s) |
| 9 | `Orifice_QCoef_In` | `node_c1` | Orifice coefficient $C_1$ for surface → sewer flow |
| 10 | `Orifice_Area_In` | `node_area1` | Orifice area $A_1$ (m²) for surface → sewer flow |

Exchange discharge is computed in `cal_exchange_flow` (`iterateq.cpp`) from the SWMM node head, the surface water level and the crest elevation. Every active branch uses the orifice equation: sewer → surface with $C_0, A_0$ (capped by `Max_Qout`), and surface → sewer with $C_1, A_1$ in both the free and the submerged case (capped by `Max_Qin` and by the water in the cell). The weir equation is in the source but commented out. The three branches are written out in [`coupling_methodology.md`](coupling_methodology.md). The signed result is stored in `Vconnect::flow_1D_node` and passed to SWMM as `Node_Inflow`.

## 5. Outputs

Standard LISFLOOD-FP rasters / mass files appear under `dirroot`. SWMM writes `swmmrpt` / `swmmout`. Coupling builds may also emit mass-audit / exchange diagnostics depending on solver flags (see `cuda/cuda_stats.*`).

## 6. Hotstart

Surface and SWMM hotstart paths were extended for coupled states (`chkpnt.cpp`, `hotstart.c`). Prefer rebuilding hotstarts with the same coupling version; mixing with pre-coupling hotstarts can disagree on pipe storage definitions.

## 7. Licence reminders when redistributing binaries

- Source and GPL-covered binaries: GPL-2.0.
- MSVC CRT / OpenMP and Intel OpenMP / math DLLs (if present in a release zip) are under their vendors’ redistributable EULAs — **not** GPL. Ship them only with `THIRD_PARTY_NOTICES.txt`, never as a substitute for providing corresponding source.
