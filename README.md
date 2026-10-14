<div align="center">

# LISFLOOD-FP × SWMM

**Two-way coupling of a GPU 2-D surface flood model and a 1-D urban drainage model**

[![Licence: GPL-2.0](https://img.shields.io/badge/licence-GPL--2.0-2C5577.svg)](LICENSE)
[![Release v1.0.0](https://img.shields.io/badge/release-v1.0.0-2E7D5B.svg)](https://github.com/FloodGroup/LISFLOOD-SWMM-Coupling/releases/tag/v1.0.0)
[![LISFLOOD-FP 8.1](https://img.shields.io/badge/LISFLOOD--FP-8.1.0-3F6F9A.svg)](https://zenodo.org/records/6912932)
[![EPA SWMM 5.1.015](https://img.shields.io/badge/EPA%20SWMM-5.1.015-3F7F6E.svg)](https://github.com/USEPA/Stormwater-Management-Model)
[![CUDA](https://img.shields.io/badge/CUDA-required%20for%20coupling-B4643A.svg)](https://developer.nvidia.com/cuda-toolkit)
[![CMake](https://img.shields.io/badge/build-CMake%20%E2%89%A5%203.13-7A7A7A.svg)](https://cmake.org)

[Overview](#overview) ·
[How it works](#how-the-coupling-works) ·
[Quick start](#quick-start) ·
[Documentation](#documentation) ·
[Citation](#citation)

</div>

<p align="center">
  <img src="docs/figures/surface.svg" width="100%" alt="Cross-section of a street with a 2-D surface grid above a buried SWMM sewer. Water drains into one manhole and surcharges from another; the sewer discharges through a COUPLED outfall into a river cell.">
</p>

## Overview

This repository couples the 2-D hydrodynamic model [LISFLOOD-FP](https://www.bristol.ac.uk/geography/research/hydrology/models/lisflood/) with the 1-D sewer model [EPA SWMM](https://www.epa.gov/water-research/storm-water-management-model-swmm) inside a single executable. Surface water and pipe water exchange in both directions while both models run:

- **Manholes** drain surface water into the sewer, or spill sewer water onto the street when the network surcharges. Discharge follows orifice and weir laws driven by the head difference between the two models.
- **`COUPLED` outfalls**, a new SWMM outfall type, let the sewer discharge into a river or coastal cell of the 2-D grid and take water back when that cell is higher. The receiving cell is treated as a finite volume, so the outfall cannot create or lose water.
- **SWMM node flooding** is returned to the 2-D cell above the node instead of disappearing from the model.

Both physics engines are unchanged: LISFLOOD-FP still solves the shallow-water equations, and SWMM still solves the dynamic-wave equations in its conduits. The coupling layer only decides what water crosses the interface, when, and how much.

| | |
|---|---|
| **Surface solvers** | CUDA ACC (local inertia) and CUDA FV1 (first-order Godunov) |
| **Exchange interval** | every 2-D step, or every *T* seconds with exact time landing |
| **Data transfer** | only the depths of coupled cells move between GPU and host |
| **Mass balance** | volume-based continuity at coupled outfalls, plus a per-run water-balance audit |

## How the coupling works

<p align="center">
  <img src="docs/figures/layers.svg" width="100%" alt="Three-column architecture: LISFLOOD-FP modules on the left, the coupling layer in the middle, SWMM modules on the right, with arrows for depth, point sources, Node_Inflow, receiving-cell geometry, exchange volume and stored volume.">
</p>

Each 2-D time step, the coupler reads the water depth in every cell that contains a manhole, computes the manhole discharge, and adds it to the surface as a point source. The same discharge is accumulated as inflow for SWMM. When an exchange is due, SWMM advances by exactly the elapsed interval using its own sub-steps. It receives the accumulated inflow as a constant rate and the current water level of every receiving cell. SWMM then returns the volume that actually crossed each `COUPLED` outfall and each flooded node, and LISFLOOD-FP applies those volumes to the surface.

The full derivation, every equation with the function that implements it, and the time-synchronisation rules are in [**docs/coupling_methodology.md**](docs/coupling_methodology.md).

## Quick start

### Requirements

| | |
|---|---|
| GPU | NVIDIA GPU with CUDA Toolkit ≥ 11 (coupled runs are GPU-only) |
| Compiler | MSVC 2019 / 2022 on Windows, or GCC / Clang on Linux (C++14) |
| Build | CMake ≥ 3.13, OpenMP |
| Libraries | NetCDF + HDF5 + zlib (can be disabled), `libnuma` on Linux |

### Build

```bash
# Linux
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86
cmake --build build --target lisflood -j
```

```powershell
# Windows (MSVC + NMake)
.\tools\build.ps1 -VcVars "C:\...\vcvars64.bat" -CMake "C:\...\cmake.exe" -CudaArchitecture 86
```

Set `CMAKE_CUDA_ARCHITECTURES` to your GPU's compute capability (86 for RTX 30xx, 89 for RTX 40xx). Windows NetCDF libraries and other build options are covered in the [user guide](docs/user_guide.md#2-build).

### Run the coupled example

```bash
cd examples/02_coupled_synthetic
lisflood case.par
```

A coupled `.par` file adds five keywords to a normal LISFLOOD-FP run:

```text
cuda
acceleration                    # or fv1
swmminp      network.inp        # SWMM network with COUPLED outfalls
floodconfig  floodconfig.txt    # manhole <-> cell exchange table
coupleint    step               # or a positive number of seconds
```

## Documentation

| Document | What it covers |
|---|---|
| [User guide](docs/user_guide.md) | Building, preparing a coupled case, every coupling keyword, the `floodconfig` format, outputs, hotstart, troubleshooting |
| [Coupling methodology](docs/coupling_methodology.md) | Architecture, the step-by-step exchange cycle, manhole and outfall equations, time synchronisation, mass-balance audit |
| [Examples](examples/README.md) | Two small runnable cases: surface-only smoke test and a synthetic coupled network |
| [Source map](src/README.md) | Where the coupling code lives inside `src/`, and what was changed relative to upstream |

## Repository layout

```text
.
├── docs/
│   ├── coupling_methodology.md   how the two models exchange water
│   ├── user_guide.md             build, inputs, outputs, troubleshooting
│   └── figures/                  SVG schematics used by the docs
├── examples/
│   ├── 01_surface_only_drycheck/ standalone LISFLOOD-FP smoke test
│   └── 02_coupled_synthetic/     minimal coupled case (one junction, one COUPLED outfall)
├── src/                          LISFLOOD-FP 8.1 + SWMM 5.1.015 + coupling layer
├── tools/
│   ├── build.ps1                 Windows build script
│   └── run.ps1                   Windows launcher used in the binary package
├── CITATION.cff
└── LICENSE                       GPL-2.0
```

## Upstream versions and licence

| Component | Version | Licence |
|---|---|---|
| LISFLOOD-FP | **8.1.0** development snapshot, git commit `637106b` (2022-04-20). This commit is in the Zenodo 8.2 source history but is not the Zenodo 8.1 / 8.2 release archive. | GPL-2.0 ([Zenodo 8.1](https://zenodo.org/records/6912932)) |
| EPA SWMM | **5.1.015** (`VERSION 51015` in `src/consts.h`; cf. [USEPA v5.1.15](https://github.com/USEPA/Stormwater-Management-Model)) | Public domain |

**The whole repository is released under GPL-2.0** ([`LICENSE`](LICENSE)). SWMM is public domain and is distributed here under the same licence; its use does not imply endorsement by the U.S. EPA.

Pre-built Windows binaries are published on the **Releases** page, not in the git tree. They include Microsoft and Intel runtime libraries under those vendors' redistribution terms; see `THIRD_PARTY_NOTICES.txt` inside the archive.

## Citation

If you use this software, please cite it using [`CITATION.cff`](CITATION.cff) (GitHub shows a *Cite this repository* button), and also cite LISFLOOD-FP and SWMM:

- Shaw, J., Kesserwani, G., Neal, J., Bates, P. and Sharifian, M. K. (2021). LISFLOOD-FP 8.0: the new discontinuous Galerkin shallow-water solver for multi-core CPUs and GPUs. *Geoscientific Model Development*, 14, 3577–3602. <https://doi.org/10.5194/gmd-14-3577-2021>
- Rossman, L. A. (2015). *Storm Water Management Model User's Manual, Version 5.1*. EPA/600/R-14/413b, U.S. Environmental Protection Agency.

## Acknowledgements

- LISFLOOD-FP: University of Bristol Hydrology Research Group and contributors
- EPA SWMM: U.S. Environmental Protection Agency
- Coupling layer: FloodGroup
