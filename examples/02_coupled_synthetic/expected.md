# Expected behaviour

- The model starts, reads the DEM and the SWMM input, maps `J1` through `floodconfig.txt`, and maps `OF1` as a `COUPLED` outfall.
- With no initial surface water, rainfall or inflow, all exchange volumes stay at about zero.
- The run completes `sim_time = 60 s` and prints a final line beginning `V11 END state=`.
- `combined_residual_m3` in the `PROGRESS` lines of `results/syn.audit.log` stays at round-off level.

The case checks input wiring, not hydrological realism.
