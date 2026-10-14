# Coupling methodology

This document describes the bidirectional surface–sewer coupling implemented in this repository. Every formula and behavioural claim cites the implementing file and function. Internal SWMM storage uses US customary units (feet, ft³, s); the public exchange API converts to SI (m, m³) with the factor `0.028316846592` ft³ → m³ (`coupled_boundary.c`).

## Overview

```mermaid
flowchart LR
  subgraph Surface["2-D surface (LISFLOOD-FP)"]
    DEM[DEM / water depth H]
    SOLVER["ACC / FV1 / other SWE"]
  end
  subgraph Exchange["Coupling layer"]
    CLOCK["SwmmCoupleClock\n(swmm_couple.h)"]
    GATHER["gather_coupling_depths\n(cuda_coupling_gather.cuh)"]
    CB["coupled_boundary.c\nvolume CV + limits"]
    VIN["Manhole inlet Q\n(weir / orifice)"]
  end
  subgraph Sewer["1-D sewer (SWMM 5.1.015)"]
    DW["Dynamic wave\n(dynwave.c)"]
    OUT["COUPLED outfalls"]
    NODES["Junctions / links"]
  end
  SOLVER --> CLOCK
  DEM --> GATHER
  GATHER --> CB
  CLOCK --> CB
  VIN --> NODES
  CB --> OUT
  OUT --> DW
  DW --> CB
  CB --> SOLVER
```

![Coupling overview](figures/coupling_overview.png)

## What was changed relative to upstream

1. **New files**: `coupled_boundary.c` / `.h`, `swmm_couple.h`, `cuda/cuda_coupling_gather.cuh` (plus packaging / build helpers).
2. **LISFLOOD-FP (~47 files)**: main loops (CPU and CUDA ACC/FV1) call SWMM advance + exchange; DEM–node mapping; point-source injection of exchange volumes; `.par` keywords `swmminp`, `swmmrpt`, `swmmout`, `coupleint`.
3. **SWMM (~17 files)**: new outfall type `COUPLED_OUTFALL`; coupling APIs in `swmm5.c`; dynamic-wave hooks; volume accounting.

Full file list: [`CHANGES.md`](../CHANGES.md).

## Coupling modes

### A. `COUPLED` outfalls ↔ receiving 2-D cells (volume control volumes)

Outfalls declared `COUPLED` in the SWMM `.inp` (`text.h`: `w_COUPLED`; `enums.h`: `COUPLED_OUTFALL`) are mapped to one or more surface cells. The outfall is treated as an **exchange face**, not an extra storage node (`coupled_boundary.c` file header and `end_volume()`).

Receiving-cell geometry is supplied each exchange interval by `coupled_set_cell(node, cell, bed, depth, area)`:

- `bed`, `depth`, `area` are SI (m, m, m²).
- Internally: `bed_ft = bed/0.3048`, `area_ft2 = area/(0.3048²)`, `water_ft3 = depth*area/0.028316846592`.
- External stage applied to the outfall:  
  \( H_{\mathrm{ext}} = z_{\mathrm{bed}} + V_{\mathrm{cell}} / A_{\mathrm{cell}} \)  
  (`coupled_set_cell`, lines setting `Outfall[].externalStage`).

**Cell water volume** (`cell_volume`):

\[
V_{\mathrm{cell}}(H) = A_{\mathrm{cell}} \max(0,\, H - z_{\mathrm{bed}})
\]

**Nodal storage volume** for marked (coupled-subgraph) nodes (`node_volume` = `own_volume` + `end_volume`):

- Junction / storage body: SWMM `node_getVolume`, plus ponding or minimum surface area above full depth (`own_volume`).
- Conduit end storage attributed to endpoints with weights: interior–interior \(1/2\)–\(1/2\); conduit ending at an outfall attributes the **full** pipe storage to the interior node and **zero** to the outfall face (`end_volume`).

**Backward-Euler volume update** for a non-outfall coupled node (`coupled_solve_node`):

\[
V^{n+1}_{\mathrm{target}} = V^{n}_{\mathrm{old}} + \Delta t\,(Q_{\mathrm{in}} - Q_{\mathrm{out}})
\]

Depth is recovered by safeguarded Newton inversion of `node_volume` (`inverse_volume`). Overflow (if no ponding) is the residual volume rate above the capped head.

**Receiving-cell head** after summing net nodal flow into the cell (`coupled_solve_cells`):

\[
H = z_{\mathrm{bed}} + \max\!\left(0,\, \frac{V_{\mathrm{old}} + \Delta t\,\sum(Q_{\mathrm{in}}-Q_{\mathrm{out}})}{A_{\mathrm{cell}}}\right)
\]

**Finite-supply flux limiting** (`coupled_limit_flows` / `scale_donor`): withdrawals from a cell or junction cannot exceed available volume over \(\Delta t\) (including water reserved for 2-D manhole outflow over the remaining exchange interval). Donor link flows, evaporation and seepage are scaled by a common factor.

**Exchange volume credited at a COUPLED outfall** over a SWMM sub-step (`coupled_commit`):

\[
\Delta V_{\mathrm{ex}} = Q_{\mathrm{net}}\Delta t - \big(V_{\mathrm{end}}^{n+1} - V_{\mathrm{end}}^{n}\big)
\]

returned to LISFLOOD in m³ via `coupled_exchange()` (= `exchange * 0.028316846592`).

Requirements (`coupled_open`): continuous **DYNWAVE** routing, no routing `EVENTS`, hydraulic-only (no quality) when pollutants are present.

### B. Manhole inlets (surface → sewer)

Vertical connections are described by a table (`vconnect.config` format in production models; see examples). Columns include orifice and weir coefficients, ponded area, max inflow/outflow and weir crest length. Runtime structures are `Vconnect` in `geometry.h` (`node_Vindex`, `flow_1D_node`, …).

Inlet discharges computed for the 2-D → 1-D direction are accumulated into `Node_Inflow[]` and applied through `coupled_sync_inlets()` when `ADDLINK == 1`. The interval-average rate is held constant on both ends of SWMM’s trapezoidal lateral-flow integration.

### Manhole exchange formulas (implemented)

Defined in `iterateq.cpp`:

```text
orificeflow_cal(c, A, H_up, H_down) = c * A * sqrt(2 * g * (H_up - H_down))
weirflow_cal(cw, w, h)             = cw * w * h * sqrt(2 * g * h)
```

with \(g =\) `SI_GRAVITY`. Heads are in metres (`UCF(LENGTH)` converts SWMM node head to SI). Crest elevation \(Z = \max(z_{\mathrm{node\,surface}},\, z_{\mathrm{DEM}})\).

Logic in `cal_exchange_flow` (tolerance `COUPLING_HEAD_TOL`):

1. **Sewer → surface** when \(H_{1D}-Z > arepsilon\) and \(H_{1D}-H_{2D} > arepsilon\):  
   \(Q = \minig(\mathrm{orifice}(C_0,A,H_{1D},H_{2D}),\, Q_{\mathrm{out,max}}ig)\);  
   `flow_1D_node = -Q`, `flow_2D_node = +Q` (`node_c0`, `node_area`, `node_qoutmax`).
2. **Surface → sewer** when the surface is above the crest and above (or equal within tolerance to) the sewer head:  
   \(Q = \mathrm{orifice}(C_1, A_1, H_{2D}, H_{1D})\) using `node_c1` / `node_area1`.  
   A weir branch calling `weirflow_cal(node_cw, node_w, H2D-Z)` exists in source but is **commented out**; the active path uses the orifice form above.
3. Surface→sewer rates are then limited by available cell volume  
   \(\min(Q,\, \Delta x\Delta y\, h / (n_{\mathrm{neg}}\Delta t),\, Q_{\mathrm{in,max}})\).

The `.par` keyword **`floodconfig`** points to the vertical-connection table (`pars.cpp`).

## Time synchronisation

![Exchange clock](figures/time_sync.png)


`SwmmCoupleClock` (`swmm_couple.h`):

- `coupleint step` → `every_step = 1`: one SWMM exchange per 2-D step.
- `coupleint <seconds>`: accumulate 2-D volumes until elapsed time reaches the interval (or sim / SWMM end), then exchange once.
- `swmm_couple_clip_dt` shortens the 2-D `dt` so exchange times land exactly on the interval.
- Inside each exchange window SWMM advances with its native dynamic-wave sub-steps; volumes passed to 2-D are **integrated**, not instantaneous \(Q\times\Delta t_{\mathrm{couple}}\) approximations (`coupled_commit` / SWMM `swmm_advanceCoupled` API in `swmm5.c`).

## GPU path

For CUDA ACC / FV1:

1. Only depths at coupling cell indices are gathered: `lis::cuda::gather_coupling_depths` (`cuda/cuda_coupling_gather.cuh`) — compact GPU buffer, then copy to host.
2. Host runs SWMM + `coupled_boundary` for the interval.
3. Exchange volumes / stages are written back as point sources through `cuda/cuda_boundary.cu` (ordered updates when multiple sources share a cell).

`coupleint` with a positive interval is implemented for **CUDA acceleration** and **CUDA fv1** only (`pars.cpp` validation).

## Governing equations (unchanged physics engines)

| Domain | Equations | Where |
|---|---|---|
| 2-D surface | Local-inertia / shallow-water solvers already in LISFLOOD-FP (ACC, FV1, DG2, …) | `fp_acc.cpp`, `swe/*`, `cuda/acc/*`, `cuda/fv1/*` |
| 1-D sewer | SWMM dynamic-wave momentum + continuity | `dynwave.c`, `dwflow.c`, … |
| Coupled subgraph continuity | Backward Euler finite-volume update above | `coupled_boundary.c` |

Momentum in conduits remains SWMM’s dynamic wave; only the **continuity / storage attribution** for components that contain `COUPLED` outfalls is replaced.

## Schematic: outfall receiving cell

```mermaid
flowchart TB
  subgraph Cell["2-D receiving cell"]
    Z["bed elevation z_bed"]
    H2["surface head H = z_bed + h"]
    A["cell area A"]
  end
  subgraph Pipe["Terminal conduit"]
    N["interior node storage"]
    F["COUPLED outfall face\n(no half-pipe storage)"]
  end
  H2 <-->|"ΔV_ex integrated"| F
  N --- F
```

![Outfall schematic](figures/outfall_schematic.png)

## References in tree

| Topic | Primary source |
|---|---|
| Volume CVs, limits, exchange | `src/coupled_boundary.c` |
| Exchange clock | `src/swmm_couple.h` |
| GPU depth gather | `src/cuda/cuda_coupling_gather.cuh` |
| `.par` keywords | `src/pars.cpp` |
| `COUPLED` keyword | `src/text.h`, `src/enums.h`, `src/keywords.c` |
| SWMM coupling API | `src/swmm5.c`, `src/swmm5.h` |
