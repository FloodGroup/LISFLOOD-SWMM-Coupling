# Source map

`src/` contains LISFLOOD-FP, SWMM and the coupling layer. All of it is compiled into one executable, `lisflood`, by `CMakeLists.txt`. The CUDA sources are compiled only when CMake finds a CUDA compiler; the CPU coupling is always built. For build instructions, see the [user guide](../docs/user_guide.md#2-build).

## Upstream baselines

| Component | Baseline |
|---|---|
| LISFLOOD-FP | development commit `637106b` (2022-04-20), reported as 8.1.0 in `VersionHistory.h` |
| EPA SWMM | 5.1.015 (`VERSION 51015` in `consts.h`) |

The tables below locate the coupling implementation. Public release changes are recorded only in [CHANGES.md](../CHANGES.md).

## Coupling layer

| File | Back end | Role |
|---|---|---|
| `swmm_couple.h` | shared | exchange clock, aligned substeps, inlet accumulation and pending return reservoirs |
| `coupling_restart.h`, `hydraulic_restart_io.h` | shared | bundled surface, SWMM and pending-water restart I/O and compatibility checks |
| `coupled_boundary.c` / `.h` | shared (SWMM side) | `COUPLED` outfalls: receiving-cell volumes, donor limits, backward-Euler node update, exchange volumes |
| `cpu_coupling.cpp` / `.h` | CPU | outfall mapping, `COUPLED` outfall stages, returned outfall and flooding volumes, point-source application |
| `cuda/cuda_coupling_gather.cuh` | CUDA | gather of coupled-cell depths from GPU to host |
| `cuda/run_audit.cuh` | CUDA | `<resroot>.audit.log` and `<resroot>.water_balance.csv` |

## Solver and engine areas

| Area | Files | Role |
|---|---|---|
| Build | `CMakeLists.txt`, `config.default.cmake` | SWMM and the coupling layer are compiled into `lisflood`; CUDA is optional |
| Parameters | `pars.cpp`, `lisflood.h`, `lisflood.cpp` | coupling keywords; checks that coupled runs use ACC or FV1 (CPU or CUDA) |
| CPU drivers | `iterateq.cpp` (ACC), `swe/fv1.cpp`, `swe/fv1.h` (FV1) | `floodconfig` reader, node-to-cell mapping, manhole discharge (`cal_exchange_flow`), exchange cycle inside the 2-D time loop |
| CUDA drivers | `cuda/acc/cuda_acc_simulate.cu`, `cuda/fv1/cuda_fv1_simulate.cu` | exchange cycle inside the 2-D time loop |
| CUDA mapping and exchange | `cuda/cuda_dem.cu`, `cuda/cuda_dem.cuh` | `floodconfig` reader, node-to-cell mapping, manhole discharge (`cal_flow_exchange`), outfall and flooding sources |
| CUDA point sources | `cuda/cuda_boundary.cu` | exchange discharge injected as 2-D sources |
| SWMM core | `swmm5.c` / `.h`, `routing.c`, `dynwave.c`, `node.c`, `enums.h`, `keywords.c`, `text.h`, … | `COUPLED` outfall type, coupling API, volume hooks |
| Restart | `coupling_restart.h`, `hydraulic_restart_io.h`, `swmm5.c`, `routing.c`, `dynwave.c` | complete hydraulic restart for CPU/CUDA ACC and FV1 |
| Legacy warm starts | `chkpnt.cpp`, `hotstart.c` | surface-only checkpoints and conventional SWMM warm starts |
| CPU DG2 | `swe/dg2.cpp` | contains exchange routines that are not used: DG2 runs standalone only |

## Upstream layout

| Path | Content |
|---|---|
| `*.cpp`, `lisflood.h` | LISFLOOD-FP CPU core: input, output, boundaries, ACC solver |
| `*.c`, `*.h` (lower-case SWMM names) | EPA SWMM 5.1.015 engine |
| `cuda/` | LISFLOOD-FP GPU solvers (`acc/`, `fv1/`, `fv2/`, `dg2/`, `acc_nugrid/`) and shared kernels |
| `swe/` | CPU FV1 / DG2 shallow-water solvers |
| `lisflood2/` | sub-grid channel model |
| `rain/` | rainfall input, including NetCDF |
| `preprocess/`, `postprocess/` | auxiliary tools (`generateDG2DEM`, `generateDG2start`, …) |
| `test/` | unit tests (Catch) |
| `windep/` | NetCDF headers and `netcdf.lib` for Windows builds |
| `cmake/` | CMake helper modules |
