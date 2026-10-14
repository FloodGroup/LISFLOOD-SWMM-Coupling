# Coupling methodology

This document explains how LISFLOOD-FP and SWMM exchange water in this repository. Each equation is followed by the function that implements it, so the text can be checked against the source.

**Contents**

1. [Design principles](#1-design-principles)
2. [Architecture](#2-architecture)
3. [One coupling cycle](#3-one-coupling-cycle)
4. [Linking SWMM nodes to 2-D cells](#4-linking-swmm-nodes-to-2-d-cells)
5. [Manhole exchange](#5-manhole-exchange)
6. [COUPLED outfalls](#6-coupled-outfalls)
7. [Ordinary outfalls and node flooding](#7-ordinary-outfalls-and-node-flooding)
8. [Time synchronisation](#8-time-synchronisation)
9. [Units](#9-units)
10. [Water-balance audit](#10-water-balance-audit)
11. [Requirements and limitations](#11-requirements-and-limitations)
12. [Source map](#12-source-map)

---

## 1. Design principles

- **Keep both physics engines.** LISFLOOD-FP still solves the 2-D shallow-water equations and SWMM still solves the 1-D dynamic-wave equations. The coupling layer only adds exchange terms and, in network components that contain a `COUPLED` outfall, replaces SWMM's continuity update with a volume-based one.
- **Exchange volumes, not instantaneous rates.** SWMM is advanced across the whole exchange interval with its own sub-steps, and the volume that actually crossed each interface is returned to the surface (`swmm_advanceCoupled` in `swmm5.c`).
- **Never withdraw water that is not there.** Every withdrawal is capped by the water available in the donor: the 2-D cell for inlet drainage, the sewer storage for surcharge, and the receiving cell or junction for outfall exchange.
- **Move as little data as possible.** Only the depths of coupled cells are copied from the GPU to the host each step (`gather_coupling_depths`).

## 2. Architecture

<p align="center">
  <img src="figures/layers.svg" width="100%" alt="Architecture of the coupled model">
</p>

The coupled executable has three layers. LISFLOOD-FP owns the surface grid and the main time loop. The coupling layer (new in this repository) owns the exchange clock, the depth gather, the manhole discharge law and the receiving-cell bookkeeping. SWMM owns the sewer network and exposes a small coupling API.

| SWMM API (`swmm5.h`) | Purpose |
|---|---|
| `swmm_beginCoupledExchange()` | reset per-interval exchange totals |
| `swmm_setOutfallReceivingCell(node, cell, bed, depth, area)` | give a `COUPLED` outfall the current state of its 2-D cell (SI units) |
| `swmm_advanceCoupled(T, &elapsed)` | route the network forward by exactly `T` seconds |
| `swmm_getOutfallExchangeVolume(node)` | volume that crossed an outfall during the interval (m³) |
| `swmm_getNodeFloodExchangeVolume(node)` | node flooding volume during the interval (m³) |
| `swmm_getNodeAvailableCoupledVolume(node)` | water stored at a node, used to cap surcharge (m³) |

## 3. One coupling cycle

<p align="center">
  <img src="figures/exchange_cycle.svg" width="100%" alt="Nine-step coupling cycle across three swim lanes">
</p>

The cycle runs inside the main loop of `cuda/acc/cuda_acc_simulate.cu` (the FV1 driver `cuda/fv1/cuda_fv1_simulate.cu` follows the same structure).

| Step | What happens | Code |
|:---:|---|---|
| 1 | The 2-D step `dt` is set by the CFL condition and rainfall clipping. In interval mode it is also shortened so the exchange time lands exactly on the interval boundary. | `update_dt_ACC`, `swmm_couple_clip_dt` |
| 2 | Water depths in the coupled cells are gathered into a compact GPU buffer and copied to the host. | `gather_coupling_depths` |
| 3 | The discharge through every manhole is computed from the SWMM node head and the surface water level ([§5](#5-manhole-exchange)). | `Topography::cal_flow_exchange` |
| 4 | The discharges are written as point sources, and the 2-D fluxes and depths are updated. | `FloodplainQ`, `Boundary::update_H` |
| 5 | Inlet volumes are accumulated: \(V_i \mathrel{+}= Q_i\,dt\). | `swmm_couple_accumulate` |
| — | If the interval is not complete, the loop continues with the next 2-D step. | `swmm_couple_due` |
| 6 | The accumulated volume is converted to a constant inflow rate `Node_Inflow = V / T`, and each `COUPLED` outfall receives the bed, depth and area of its cell. | `swmm_couple_write_inflow`, `update_outfall_stages` |
| 7 | SWMM routes the network forward by `T`. Inside every dynamic-wave sub-step the coupled components use the volume-based update of [§6](#6-coupled-outfalls). | `swmm_advanceCoupled` |
| 8 | Exchanged volumes at outfalls and flooded nodes are read back and converted to point-source fluxes. | `collect_outfall_source_values`, `collect_surface_flood_values` |
| 9 | The sources are added to the surface depths and the clocks are reset. | `update_point_sources_range`, `swmm_couple_advance` |

With `coupleint step`, steps 5 and 6 collapse: the inflow of the current step is passed directly and `T` equals `dt`.

## 4. Linking SWMM nodes to 2-D cells

Manholes are listed explicitly in the `floodconfig` table ([user guide §4](user_guide.md#4-the-floodconfig-table)). Outfalls are found automatically from the SWMM input.

**Manholes.** For each row of `floodconfig`, the node's SWMM coordinates are converted to a DEM cell (`cal_demindex` in `cuda_dem.cu`):

$$
i_x = \left\lfloor \frac{x - x_{\min}}{\Delta x} \right\rfloor, \qquad
i_y = \left\lfloor \frac{y_{\max} - y}{\Delta y} \right\rfloor
$$

A node is excluded from exchange (index `-9999`) when it lies outside the DEM, falls on a NoData cell, or is an `OUTFALL`, `STORAGE` or `DIVIDER` node. The rim elevation \(z_{\text{rim}} = z_{\text{invert}} + d_{\text{full}}\) is compared with the DEM, and differences above 0.1 m are reported (`fix_maxdepth`). When several manholes share one cell, they are grouped so the cell's water is shared fairly (`cal_demcount`).

**Outfalls.** Every SWMM outfall whose coordinates fall inside the DEM is mapped to a cell (`map_internal_outfalls`). Outfalls outside the DEM remain ordinary external sinks. A `COUPLED` outfall that cannot be mapped stops the run.

All mapping decisions are written to the file named by the `flooderror` keyword.

## 5. Manhole exchange

<p align="center">
  <img src="figures/manhole_exchange.svg" width="100%" alt="Three manhole regimes: surcharge, free inflow, submerged inflow">
</p>

Three heads are compared at every coupled manhole (`Topography::cal_flow_exchange`, `cuda_dem.cu`):

| Symbol | Definition |
|---|---|
| \(H_{1D}\) | SWMM node head, invert + depth, converted to metres |
| \(H_{2D}\) | water surface in the linked cell, DEM + depth |
| \(Z\) | crest elevation, \(\max(z_{\text{rim}},\, z_{\text{DEM}})\) |

The branches are tested in order with tolerance \(\varepsilon = 10^{-3}\) m (`COUPLING_HEAD_TOL`, `consts.h`):

**(a) Surcharge, sewer → surface**, when \(H_{1D} - Z > \varepsilon\) and \(H_{1D} - H_{2D} > \varepsilon\):

$$
Q = \min\!\Big( C_0 A_0 \sqrt{2g\,(H_{1D} - H_{2D})},\;\; Q_{\text{out,max}},\;\; \frac{V_{\text{stored}} + V_{\text{pending}}}{dt} \Big)
$$

The last term prevents the manhole from releasing more water than SWMM currently holds at that node, including inflow already accumulated for the next exchange.

**(b) Free inflow, surface → sewer**, when \(H_{1D} \le Z + \varepsilon\) and \(H_{2D} - Z > \varepsilon\). The sewer head is below the crest, so the inlet behaves as a weir:

$$
Q = C_w\, w\, (H_{2D} - Z)\sqrt{2g\,(H_{2D} - Z)}
$$

**(c) Submerged inflow, surface → sewer**, when \(H_{1D} \ge Z - \varepsilon\) and \(H_{2D} - H_{1D} > \varepsilon\):

$$
Q = C_1 A_1 \sqrt{2g\,(H_{2D} - H_{1D})}
$$

Otherwise \(Q = 0\).

**Supply limit for inflow.** After branch selection, inflow in (b) and (c) is limited by the water in the cell, shared among the \(n\) manholes draining that cell in this step, and by the inlet capacity:

$$
Q_{\text{in}} = \min\!\Big( Q,\;\; \frac{h\,\Delta x\,\Delta y}{n\,dt},\;\; Q_{\text{in,max}} \Big)
$$

**Sign convention.** `flow_1D_node` is positive into the sewer; `flow_2D_node` has the opposite sign. The surface receives `flow_2D_node / Δy` as a per-width point-source flux (m²/s), which is how LISFLOOD-FP stores point sources.

> [!NOTE]
> The CPU routines in `iterateq.cpp`, `swe/fv1.cpp` and `swe/dg2.cpp` contain an older version of this exchange. They are not reachable, because `lisflood.cpp` rejects coupled runs that do not use CUDA ACC or CUDA FV1. In the CPU ACC routine, branch (b) uses the orifice law instead of the weir law.

## 6. COUPLED outfalls

<p align="center">
  <img src="figures/outfall_cv.svg" width="100%" alt="Longitudinal section of a COUPLED outfall and its receiving 2-D cell, with the volume equations">
</p>

A standard SWMM outfall imposes a stage, such as free outfall, normal depth or a fixed level, and its receiving water is infinite. A `COUPLED` outfall instead connects the last conduit to a finite 2-D cell, so water leaving the sewer raises that cell and a high cell pushes water back into the pipe. The implementation is in `coupled_boundary.c`.

### 6.1 Which nodes use the volume update

At start-up `coupled_open` finds the connected components of the network with a union–find pass over all links. Every node in a component that contains at least one `COUPLED` outfall is *marked*. Marked nodes use the volume-based continuity below; all other nodes keep SWMM's native depth update. Conduit momentum is solved by SWMM's dynamic wave everywhere.

### 6.2 Storage attribution

Conduit storage is split between the two end nodes (`end_volume`):

| Conduit between | Weight at each end |
|---|---|
| two interior nodes | ½ and ½ |
| an interior node and an outfall | 1 at the interior node, 0 at the outfall |

The outfall is an **exchange face**, not a storage node. Giving it zero weight means a dry pipe cannot gain water just because the receiving cell has a bed elevation above the outfall invert. The volume of a marked node is

$$
V_{\text{node}}(H) = V_{\text{own}}(H) + V_{\text{end}}(H)
$$

where \(V_{\text{own}}\) is the SWMM junction or storage volume, extended above full depth by the ponded area (if ponding is allowed) or by the minimum surface area (`own_volume`).

### 6.3 Receiving cell

Each exchange, the 2-D model passes the cell's bed \(z_{\text{bed}}\), depth \(h\) and area \(A\) (`coupled_set_cell`). The cell holds \(V_{\text{cell}} = A\,h\), and the outfall sees the stage

$$
H_{\text{ext}} = z_{\text{bed}} + \frac{V_{\text{cell}}}{A}
$$

Several outfalls may share one cell; their net flows are summed.

### 6.4 Volume update inside each SWMM sub-step

For each dynamic-wave sub-step \(\delta t\) (hooks in `dynwave.c` and `routing.c`):

1. **Limit donors** (`coupled_limit_flows`). The sum of withdrawals from a receiving cell or a marked junction may not exceed what is available. Water already promised to the 2-D surface for the rest of the interval is reserved:

   $$
   \sum Q_{\text{out}} \;\le\; \frac{V^{n} - V_{\text{reserved}}}{\delta t} + \sum Q_{\text{in}}
   $$

   If the limit is exceeded, all donor link flows, evaporation and seepage at that node are scaled by one common factor (`scale_donor`). Up to eight passes resolve chains of donors.

2. **Solve node heads** (`coupled_solve_node`). Backward-Euler continuity gives a target volume, and the head is found by a safeguarded Newton inversion of \(V_{\text{node}}(H)\) (`inverse_volume`):

   $$
   V^{n+1} = V^{n} + \delta t\,(Q_{\text{in}} - Q_{\text{out}})
   $$

   If ponding is off and the head exceeds the rim, the excess becomes the node overflow rate.

3. **Solve receiving cells** (`coupled_solve_cells`):

   $$
   H = z_{\text{bed}} + \frac{\max\!\big(0,\; V_{\text{cell}}^{n} + \delta t \sum (Q_{\text{in}} - Q_{\text{out}})\big)}{A}
   $$

4. **Commit** (`coupled_commit`). The volume equations are checked against their targets (relative tolerance about \(10^{-8}\)). The run stops with an error rather than continuing with a mass error. The volume credited at each outfall face is

   $$
   \Delta V_{\text{ex}} \mathrel{+}= Q_{\text{net}}\,\delta t - \big(V_{\text{end}}^{n+1} - V_{\text{end}}^{n}\big) = Q_{\text{net}}\,\delta t
   $$

   because the end-volume weight at the face is zero. The general form is kept in the code for clarity.

5. **Account** (`coupled_account`). SWMM's routing totals are corrected so its own continuity report stays consistent with the backward-Euler update.

After the interval, \(\sum \Delta V_{\text{ex}}\) is returned to LISFLOOD-FP in m³ (`coupled_exchange`) and applied to the receiving cell as a point source.

### 6.5 Requirements checked at start-up

`coupled_open` and `map_internal_outfalls` stop the run unless all of these hold:

- `FLOW_ROUTING DYNWAVE`, with continuous routing (no `[EVENTS]` section);
- no water-quality routing, or `IGNORE_QUALITY YES`;
- every `COUPLED` outfall lies on a valid (non-NoData) DEM cell.

## 7. Ordinary outfalls and node flooding

| Situation | Treatment |
|---|---|
| Ordinary outfall (`FREE`, `NORMAL`, `FIXED`, …) **inside** the DEM | keeps its SWMM boundary condition; its positive discharge is added to the cell as a one-way source |
| Ordinary outfall **outside** the DEM | external sink, as in standalone SWMM |
| SWMM node flooding (overflow above the rim) | the flooded volume is returned to the mapped 2-D cell (`collect_surface_flood_values`); a flooded node with no mapping stops the run |

## 8. Time synchronisation

<p align="center">
  <img src="figures/time_sync.svg" width="100%" alt="Timeline of 2-D steps, accumulated volume and SWMM sub-steps for interval and step coupling">
</p>

The exchange clock is `SwmmCoupleClock` in `swmm_couple.h`:

| Mode | Behaviour |
|---|---|
| `coupleint step` | one SWMM exchange after every 2-D step; `T = dt` |
| `coupleint T` (seconds) | 2-D steps accumulate inlet volume until their durations sum to `T`, then SWMM advances once by `T` |

In interval mode:

- `swmm_couple_clip_dt` shortens the last 2-D step so the interval ends exactly on `T`, or at the end of the 2-D run or the SWMM run, whichever comes first.
- `swmm_couple_write_inflow` passes the interval-mean inflow `V / T` to SWMM. `coupled_sync_inlets` holds this rate constant at both ends of SWMM's trapezoidal lateral-inflow integration, so the volume SWMM receives equals the volume removed from the surface.
- `swmm_advanceCoupled` caps each SWMM routing step so the network lands exactly on `t + T`, and fails if it does not.
- Outfall and flooding volumes returned for the interval are applied to the surface in the step that closes the interval.

A shorter interval couples the two models more tightly. A longer interval reduces the number of SWMM calls and host–device transfers.

## 9. Units

SWMM works internally in feet, cubic feet and seconds, whatever `FLOW_UNITS` the input uses. The coupling API converts at the boundary:

| Quantity | Conversion |
|---|---|
| volume | 1 ft³ = 0.028316846592 m³ |
| length | 1 ft = 0.3048 m |
| node head to the 2-D model | `UCF(LENGTH)` |

DEM coordinates and SWMM `[COORDINATES]` must use the same projected coordinate system, in metres.

## 10. Water-balance audit

The CUDA drivers write two files next to the other outputs (`cuda/run_audit.cuh`):

| File | Content |
|---|---|
| `<resroot>.audit.log` | run settings, periodic progress lines with the combined residual, and the final state |
| `<resroot>.water_balance.csv` | time series of 2-D and 1-D storage, rain, inlet and outfall exchange in both directions, pending exchange, and the combined balance residual (m³) |

Coupling entries are signed from the 2-D point of view. Inflow accumulated but not yet passed to SWMM is reported as `pending_to_1d_m3` and included in the combined balance.

## 11. Requirements and limitations

- Coupled runs require `cuda` together with `acceleration` or `fv1`. Other solvers, including CUDA FV2 and DG2, run standalone only.
- `coupleint` with a positive interval is validated for CUDA ACC and CUDA FV1 (`pars.cpp`).
- A manhole exchanges with exactly one 2-D cell, the cell that contains its coordinates.
- The `floodconfig` row index is the node's position in the SWMM input file, not its name ([user guide §4](user_guide.md#4-the-floodconfig-table)). Reordering the input file changes the mapping.
- Hotstart files written by an uncoupled version may define pipe storage differently. Rebuild hotstarts with this version.

## 12. Source map

| Topic | Files |
|---|---|
| Exchange clock | `src/swmm_couple.h` |
| Receiving-cell volumes, donor limits, exchange volumes | `src/coupled_boundary.c`, `src/coupled_boundary.h` |
| SWMM coupling API | `src/swmm5.c`, `src/swmm5.h` |
| Hooks in SWMM routing | `src/dynwave.c`, `src/routing.c`, `src/node.c` |
| `COUPLED` keyword | `src/text.h`, `src/enums.h`, `src/keywords.c` |
| Node and outfall mapping, manhole discharge | `src/cuda/cuda_dem.cu` |
| GPU depth gather | `src/cuda/cuda_coupling_gather.cuh` |
| Point sources | `src/cuda/cuda_boundary.cu` |
| Main loops | `src/cuda/acc/cuda_acc_simulate.cu`, `src/cuda/fv1/cuda_fv1_simulate.cu` |
| Water-balance audit | `src/cuda/run_audit.cuh` |
| `.par` keywords and validation | `src/pars.cpp`, `src/lisflood.cpp` |
