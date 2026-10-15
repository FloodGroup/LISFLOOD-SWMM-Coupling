# User guide

How to build the coupled model, prepare a case, run it, and read the results. For the physics and the exchange algorithm, see [coupling_methodology.md](coupling_methodology.md).

**Contents**

1. [Requirements](#1-requirements)
2. [Build](#2-build)
3. [Prepare a coupled case](#3-prepare-a-coupled-case)
4. [The `floodconfig` table](#4-the-floodconfig-table)
5. [Parameter file reference](#5-parameter-file-reference)
6. [Run](#6-run)
7. [Outputs](#7-outputs)
8. [Hotstart](#8-hotstart)
9. [Troubleshooting](#9-troubleshooting)
10. [Redistributing binaries](#10-redistributing-binaries)

---

## 1. Requirements

| Item | Notes |
|---|---|
| CMake ≥ 3.13 | |
| C++14 compiler | MSVC 2019 / 2022 on Windows; GCC or Clang on Linux |
| OpenMP | Usually ships with the compiler |
| NetCDF + HDF5 + zlib | For NetCDF output and dynamic rainfall; can be disabled (see below) |
| `libnuma` | Linux only (`libnuma-dev` on Debian / Ubuntu) |
| NVIDIA GPU + CUDA Toolkit ≥ 11 | Optional, for GPU acceleration. CUB is taken from the toolkit. |

Coupled runs work on the CPU without any GPU. When CMake finds a CUDA compiler, the same executable also contains the CUDA solvers, and the `cuda` keyword selects them at run time.

## 2. Build

The build is defined in `src/CMakeLists.txt`; the root `CMakeLists.txt` forwards to it. The only target you need is `lisflood`.

CMake checks for a CUDA compiler (`check_language(CUDA)`):

| CUDA compiler found | Result |
|---|---|
| no | CPU build: ACC and FV1 coupling on the CPU, parallelised with OpenMP |
| yes | CPU + GPU build: the CPU solvers as above, plus CUDA ACC and CUDA FV1, selected with the `cuda` keyword |

To build without the CUDA solvers on a machine that has CUDA installed, configure in a fresh build folder from a shell in which `nvcc` is not on the `PATH` and `CUDACXX` is not set.

### Linux

```bash
sudo apt install cmake libnuma-dev libnetcdf-dev     # Debian / Ubuntu
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86
cmake --build build --target lisflood -j
```

The executable is `build/src/lisflood`. `CMAKE_CUDA_ARCHITECTURES` is ignored in a CPU build.

### Windows, scripted

`tools/build.ps1` loads the MSVC environment, configures with NMake and builds `lisflood`:

```powershell
.\tools\build.ps1 `
  -VcVars   "C:\BuildTools\VC\Auxiliary\Build\vcvars64.bat" `
  -CMake    "C:\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" `
  -CudaArchitecture 86 `
  -RedistDirectory "C:\BuildTools\VC\Redist\MSVC\<version>\x64"
```

| Parameter | Default | Meaning |
|---|---|---|
| `-VcVars` | required | path to `vcvars64.bat` |
| `-CMake` | required | path to `cmake.exe` |
| `-CudaArchitecture` | `75` | GPU compute capability without the dot; ignored in a CPU build |
| `-MsvcToolset` | `14.39` | MSVC toolset passed to `vcvars64.bat` |
| `-RedistDirectory` | empty | if set, the MSVC runtime and OpenMP DLLs are copied next to the executable |
| `-OutputDirectory` | `build\release` | where `lisflood.exe` and the DLLs are collected |

You can also open the repository folder in Visual Studio as a CMake project and build the `lisflood` target in a Release x64 configuration.

### NetCDF on Windows

CMake links Windows builds against `src/windep/netCDF4-64/lib/netcdf.lib`, which is included in the repository. The runtime DLLs (`netcdf.dll`, `hdf5.dll`, `hdf5_hl.dll`, `zlib.dll`) are not in git; place them next to the executable. To build without NetCDF, copy `src/config.default.cmake`, change `set(_NETCDF 1)` to `set(_NETCDF 0)`, and configure with `-D_CONFIG=<your copy>`.

### Build options

`src/config.default.cmake` sets numeric precision and feature switches. To change them, copy the file and pass it with `-D_CONFIG=<file>`. In a CUDA build the file also adds `compute_75` code generation, so the executable runs on any GPU of compute capability 7.5 or newer, while `CMAKE_CUDA_ARCHITECTURES` adds native code for your GPU.

## 3. Prepare a coupled case

A coupled case is a normal LISFLOOD-FP case plus a SWMM network and a manhole table:

```text
my_case/
├── case.par          LISFLOOD-FP parameters + coupling keywords
├── dem.asc           DEM (ASCII grid)
├── start.asc         initial depth (optional)
├── network.inp       SWMM network
└── floodconfig.txt   manhole <-> cell table
```

Checklist:

- [ ] The DEM and the SWMM `[COORDINATES]` use the **same projected coordinate system**, in metres.
- [ ] `FLOW_ROUTING DYNWAVE` in the SWMM `[OPTIONS]`.
- [ ] No `[EVENTS]` section; routing must be continuous.
- [ ] No pollutants, or `IGNORE_QUALITY YES`.
- [ ] Outfalls that should exchange with the surface have type `COUPLED` and lie on a valid DEM cell.
- [ ] Every junction that can exchange with the surface has a row in `floodconfig`.
- [ ] Junction rim elevations (invert + max depth) are close to the DEM; differences above 0.1 m are reported in the `flooderror` log.
- [ ] The `.par` file contains `acceleration` or `fv1`, plus `cuda` if the run should use the GPU.

### `COUPLED` outfalls

In the SWMM `[OUTFALLS]` section, use type `COUPLED` in place of `FREE`, `NORMAL`, `FIXED`, `TIDAL` or `TIMESERIES`. It takes no stage data:

```text
[OUTFALLS]
;;Name   Elevation  Type      Gated
OF1      8.5        COUPLED   NO
```

The outfall's stage is then set every exchange from the water level of the 2-D cell that contains it ([methodology §6](coupling_methodology.md#6-coupled-outfalls)). Ordinary outfalls inside the DEM keep their SWMM boundary condition, and their discharge is added to the cell they fall in.

## 4. The `floodconfig` table

The file named by the `floodconfig` keyword lists the SWMM nodes that exchange with the surface through a manhole or inlet. It is read by `verticalcon::load_floodinfo` on the CPU (`iterateq.cpp`, `swe/fv1.cpp`) and by `Topography::read_floodconfig` on the GPU (`cuda_dem.cu`), with the same format:

- line 1: `node_num <N>`;
- line 2: a header of exactly **10** words, which is skipped (the words are free, but there must be 10);
- then *N* rows of 10 whitespace-separated fields: a label followed by 9 numbers.

```text
node_num      1
ID   Index   Orifice_QCoef   Pond_Area   Max_Qout   Weir_QCoef   Length_W   Max_Qin   Orifice_QCoef_In   Orifice_Area_In
J1   0       0.6             1.0         1.0        0.5          1.0        1.0       0.6                1.0
```

| # | Column | Variable | Meaning | Used in |
|:-:|---|---|---|---|
| 1 | `ID` | — | row label, ignored by the code | — |
| 2 | `Index` | `node_Vindex` | 0-based index of the node in SWMM's node list | mapping |
| 3 | `Orifice_QCoef` | `node_c0` | orifice coefficient \(C_0\), sewer → surface | surcharge |
| 4 | `Pond_Area` | `node_area` | orifice area \(A_0\) (m²), sewer → surface. Not SWMM's ponded area. | surcharge |
| 5 | `Max_Qout` | `node_qoutmax` | maximum discharge sewer → surface (m³/s) | surcharge |
| 6 | `Weir_QCoef` | `node_cw` | weir coefficient \(C_w\) | free inflow |
| 7 | `Length_W` | `node_w` | weir crest length \(w\) (m) | free inflow |
| 8 | `Max_Qin` | `node_qinmax` | maximum discharge surface → sewer (m³/s) | free and submerged inflow |
| 9 | `Orifice_QCoef_In` | `node_c1` | orifice coefficient \(C_1\), surface → sewer | submerged inflow |
| 10 | `Orifice_Area_In` | `node_area1` | orifice area \(A_1\) (m²), surface → sewer | submerged inflow |

**Node index.** SWMM numbers nodes in the order their lines appear in the input file, across the `[JUNCTIONS]`, `[OUTFALLS]`, `[DIVIDERS]` and `[STORAGE]` sections (`input.c`). The first node line is index 0. In a file written by the SWMM GUI, `[JUNCTIONS]` comes first, so the junctions are numbered 0, 1, 2, … in the order listed.

Rows that point to an outfall, storage unit or divider, lie outside the DEM, or fall on a NoData cell are kept but do not exchange. Each one is reported in the `flooderror` log.

## 5. Parameter file reference

### Coupling keywords

| Keyword | Value | Required | Meaning |
|---|---|:-:|---|
| `swmminp` | file | yes | SWMM input file |
| `floodconfig` | file | yes | manhole table ([§4](#4-the-floodconfig-table)); must be given together with `swmminp` |
| `swmmrpt` | file | | SWMM report file |
| `swmmout` | file | | SWMM binary results file |
| `flooderror` | file | | log of node and outfall mapping decisions |
| `coupleint` | `step` or seconds | | exchange after every 2-D step (`step`), or every *T* seconds |
| `hsfdt` | seconds | | time at which to write hotstart files (ACC only, [§8](#8-hotstart)) |
| `hsfswmmtemp` | file | | SWMM hotstart written at `hsfdt` |
| `hsffloodtemp` | file | | 2-D hotstart written at `hsfdt` |
| `hsffloodend` | file | | 2-D hotstart written at the end of the run |
| `hsffloodread` | file | | 2-D hotstart to start from |

Leaving out both `swmminp` and `floodconfig` gives a standalone LISFLOOD-FP run.

### Solver keywords for coupled runs

| Keyword | Meaning |
|---|---|
| `acceleration` | local-inertia (ACC) solver |
| `fv1` | first-order finite-volume (Godunov) solver |
| `cuda` | optional; run the chosen solver on the GPU. Needs an executable built with CUDA. |

Without `cuda`, the solver runs on the CPU with OpenMP threads. The four combinations (CPU ACC, CPU FV1, CUDA ACC, CUDA FV1) all support two-way coupling and read the same SWMM and `floodconfig` files, so a case can be moved between CPU and GPU by adding or removing one line.

All other keywords (`DEMfile`, `startfile`, `sim_time`, `initial_tstep`, `massint`, `saveint`, `fpfric`, boundary and rainfall files, …) follow the LISFLOOD-FP 8.1 manual on [Zenodo](https://zenodo.org/records/6912932).

### Example

```text
DEMfile         dem.asc
startfile       start.asc
resroot         syn
dirroot         results
sim_time        60.0
initial_tstep   1.0
massint         10.0
saveint         60.0
acceleration
# cuda          uncomment to run on the GPU

swmminp         network.inp
swmmrpt         network.rpt
swmmout         network.out
floodconfig     floodconfig.txt
flooderror      flood_error.txt
coupleint       step
```

## 6. Run

Run from the case folder so that relative paths in the `.par` file resolve:

```bash
lisflood case.par
```

On Windows the binary package includes `run.ps1`, which runs the `lisflood.exe` next to it:

```powershell
.\run.ps1 path\to\case.par
```

At start-up the console prints `LISFLOOD-SWMM | numeric=… | built=…` and the coupling interval (`SWMM coupling interval: step` or `… : 30.000 s`); a GPU run also prints `Using CUDA device: …`. A GPU run ends with a line beginning `END state=`, which gives the final state and the simulated time; the same line is written to `<resroot>.audit.log`.

On the CPU, set the number of OpenMP threads with the `OMP_NUM_THREADS` environment variable.

## 7. Outputs

| File | Written by | Content |
|---|---|---|
| `dirroot/resroot-*.wd`, `*.elev`, … | LISFLOOD-FP | depth and elevation rasters at every `saveint` |
| `dirroot/resroot.max`, `*.mxe`, … | LISFLOOD-FP | maximum depth and related maxima |
| `dirroot/resroot.mass` | LISFLOOD-FP | surface mass balance at every `massint` |
| `dirroot/resroot.audit.log` | coupler, GPU runs | run settings, progress, coupled balance residual |
| `dirroot/resroot.water_balance.csv` | coupler, GPU runs | 2-D and 1-D storage and every exchange term over time (m³) |
| `swmmrpt` | SWMM | standard SWMM report, including its continuity errors |
| `swmmout` | SWMM | binary time series, readable by SWMM GUI and `swmm-toolkit` |
| `flooderror` | coupler | which nodes and outfalls were mapped to which cells, and why others were skipped |

The columns of `water_balance.csv` are described in [methodology §10](coupling_methodology.md#10-water-balance-audit). For CPU runs, check the balance with `resroot.mass` and the continuity errors in the SWMM report.

## 8. Hotstart

Hotstart is available for the ACC solver, on the CPU and the GPU. FV1 runs always start from `startfile`.

1. In the first run, set `hsfdt`, `hsfswmmtemp` and `hsffloodtemp` (or `hsffloodend`).
2. In the next run, start the surface from the saved state with `hsffloodread`. Start SWMM from its saved state with a `USE HOTSTART` line in the SWMM `[FILES]` section.

Hotstart files written by an uncoupled LISFLOOD-FP or SWMM may attribute pipe storage differently. Regenerate them with this model.

## 9. Troubleshooting

| Message | Cause | Fix |
|---|---|---|
| `SWMM coupling requires both 'swmminp' and 'floodconfig'` | only one of the two keywords is set | give both, or neither for a standalone run |
| `SWMM coupling requires acceleration or fv1 (CPU or CUDA)` | coupled run with a solver other than ACC or FV1 | use `acceleration` or `fv1` |
| `COUPLED requires continuous DYNWAVE and hydraulic-only routing`, `COUPLED outfall requires DYNWAVE` | SWMM routing is not `DYNWAVE`, has `[EVENTS]`, or routes pollutants | fix the SWMM `[OPTIONS]` |
| `COUPLED outfall … has no valid … DEM cell` | outfall outside the DEM or on NoData | move the outfall, extend the DEM, or use an ordinary outfall type |
| `flooded node … has no mapped 2D receiving cell` | a node floods but has no valid `floodconfig` row | add the node to `floodconfig`, or set its surcharge depth so it cannot flood |
| `reverse outfall flow exceeds receiving-cell water at cell …` (CPU) | a `COUPLED` outfall took more water from its cell than the cell held | reduce SWMM `ROUTING_STEP` or `coupleint` |
| `invalid receiving-cell stage at …` (CPU) | the receiving cell's depth or bed is not finite | check the DEM and the 2-D run for instability at that cell |
| `finite receiving cell … did not converge` | receiving-cell volume equation missed its tolerance | reduce SWMM `ROUTING_STEP` or `coupleint` |
| `node … volume equation failed` | node volume equation missed its tolerance | reduce SWMM `ROUTING_STEP` |
| `SWMM nonfinite state node=…` | SWMM produced NaN or Inf at a node | check that node's geometry and offsets; reduce `ROUTING_STEP` |
| `remove swmmriveroutfalls` | parameter from an older version | delete the line; `COUPLED` outfalls are read from the SWMM input |

If a run behaves unexpectedly but does not stop, check the `flooderror` log first. Most problems come from nodes that were not mapped to the cells you expected.

## 10. Redistributing binaries

- Source code and binaries built from it are GPL-2.0. A binary release must be accompanied by the corresponding source, or a written offer of it.
- A CPU build runs on any machine, with or without a GPU.
- Microsoft (MSVC CRT, OpenMP) and Intel (OpenMP, math library) runtime DLLs in a binary archive are covered by their vendors' redistribution licences, not by the GPL. Include `THIRD_PARTY_NOTICES.txt` and do not commit these DLLs to git.
