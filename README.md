# LISFLOOD-FP × SWMM Bidirectional Coupling

Two-way coupling between the 2-D surface flood model [LISFLOOD-FP](https://www.bristol.ac.uk/geography/research/hydrology/models/lisflood/) (University of Bristol) and the 1-D urban drainage model [EPA SWMM](https://www.epa.gov/water-research/storm-water-management-model-swmm).

Surface–sewer exchange occurs at manholes (orifice inlet exchange) and at `COUPLED` outfalls (volume-conserving bidirectional exchange with the overlying 2-D cell). GPU paths (CUDA ACC and CUDA FV1) and CPU paths are supported.

| | |
|---|---|
| 📐 Coupling | Methodology, equations, and diagrams → [`docs/`](docs/) |
| 📋 Changes | File-level differences from upstream → [`CHANGES.md`](CHANGES.md) |
| 📖 Guide | Build, run, inputs → [`docs/user_guide.md`](docs/user_guide.md) |

---

## Contents

- [Upstream versions and licence](#upstream-versions-and-licence)
- [Repository layout](#repository-layout)
- [Dependencies](#dependencies)
- [Build](#build)
- [Run](#run)
- [Documentation](#documentation)
- [Examples](#examples)
- [Binary release package](#binary-release-package)
- [Citation](#citation)
- [Acknowledgements](#acknowledgements)

---

## Upstream versions and licence

| Component | Version | Licence |
|---|---|---|
| LISFLOOD-FP | **8.1.0** development snapshot, git commit `637106b` (2022-04-20). This commit appears in the Zenodo 8.2 source history but is **not** the Zenodo 8.1 / 8.2 release zip. | **GPL-2.0** ([Zenodo 8.1 record](https://zenodo.org/records/6912932)) |
| EPA SWMM | **5.1.015** (`VERSION 51015` in `consts.h`; cf. [USEPA tag v5.1.15](https://github.com/USEPA/Stormwater-Management-Model)) | Public domain |
| NVIDIA CUB (`src/cuda/cub/`) | Bundled with upstream, unmodified | BSD-3-Clause (`LICENCE-3RD-PARTY.txt`) |

> **Licence.** Because this repository contains and modifies GPL-2.0 code, **the whole repository is released under GPL-2.0**. See root [`LICENSE`](LICENSE). Do not imply EPA endorsement. Microsoft / Intel redistributable runtimes are **not** included in the source tree; they ship only in the optional binary release package (see below).

---

## Repository layout

```text
.
├── README.md                 # this file
├── LICENSE                   # GPL-2.0 full text
├── LICENCE-3RD-PARTY.txt     # NVIDIA CUB and other third-party notices
├── CHANGES.md                # modifications relative to upstream
├── CITATION.cff
├── CMakeLists.txt            # thin wrapper → src/
├── src/                      # full buildable source tree
│   ├── CMakeLists.txt
│   ├── coupled_boundary.*    # coupling volume / continuity core
│   ├── swmm_couple.h         # exchange clock
│   ├── cuda/                 # GPU solvers + coupling gather
│   ├── swe/                  # CPU SWE solvers
│   └── …                     # LISFLOOD-FP + SWMM sources (mixed upstream layout)
├── docs/                     # English documentation + figures
├── examples/                 # small runnable cases
└── tools/                    # user-facing build script
```

The flat mix of LISFLOOD-FP and SWMM files under `src/` matches the upstream build layout (shared include paths). Coupling-specific additions are named explicitly above and documented in [`docs/coupling_methodology.md`](docs/coupling_methodology.md).

---

## Dependencies

- CMake ≥ 3.13
- C/C++ compiler (MSVC 2019/2022 recommended on Windows)
- CUDA Toolkit (optional; required for GPU ACC / FV1 coupling)
- NetCDF + HDF5 (+ zlib) for NetCDF output / dynamic rainfall
- OpenMP

---

## Build

### Windows

```powershell
.\tools\build.ps1 -VcVars "C:\Path\To\vcvars64.bat" -CMake "C:\Path\To\cmake.exe" -CudaArchitecture 86
```

See comments in `tools/build.ps1` for all parameters. Or open `src/` as a CMake project in Visual Studio.

### Linux

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target lisflood
```

> 📦 **Tip.** Full prerequisite tables and Visual Studio CMake steps are in [`docs/user_guide.md`](docs/user_guide.md).

---

## Run

```text
lisflood.exe path\to\case.par
```

Essential `.par` keywords for coupling:

| Keyword | Meaning |
|---|---|
| `swmminp` | Path to SWMM `.inp` |
| `swmmrpt` | SWMM report file |
| `swmmout` | SWMM binary output |
| `coupleint` | `step` (exchange every 2-D step) or positive seconds |

In the SWMM `.inp`, mark outfalls that exchange with the 2-D surface as type **`COUPLED`**. Manhole inlet parameters are supplied via a vertical-connection table named by the `.par` keyword `floodconfig` (see [`docs/user_guide.md`](docs/user_guide.md) and [`examples/`](examples/)).

---

## Documentation

| Document | Contents |
|---|---|
| [`docs/user_guide.md`](docs/user_guide.md) | Build, run, inputs, coupling parameters |
| [`docs/coupling_methodology.md`](docs/coupling_methodology.md) | Exchange mode, equations (from source), time synchronisation, GPU gather |
| [`docs/figures/`](docs/figures/) | Flowcharts and schematics (Mermaid sources + PNG/SVG) |

---

## Examples

See [`examples/README.md`](examples/README.md). Each case is self-contained and under 5 MB.

---

## Binary release package

Pre-built Windows binaries (with the redistributable Intel and Microsoft runtime DLLs) are published as a zip on the repository's **GitHub Releases** page, not in this git tree.

> ⚠️ **Do not** commit `.exe` / proprietary DLLs to git. See `THIRD_PARTY_NOTICES.txt` inside the zip.

---

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite the upstream LISFLOOD-FP and SWMM publications.

---

## Acknowledgements

- LISFLOOD-FP: University of Bristol Hydrology Research Group
- SWMM: U.S. Environmental Protection Agency (public domain; use does not imply endorsement)
- Coupling modifications: FloodGroup / Ruiyi Yang

This README is technical documentation only and is **not legal advice**.
