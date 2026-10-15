# Expected behaviour

- The model starts, prints `SWMM coupling interval: step`, reads the DEM and the SWMM input, maps `J1` through `floodconfig.txt`, and maps `OF1` as a `COUPLED` outfall in DEM cell 55.
- With no initial surface water, rainfall or inflow, all exchange volumes stay at about zero.
- The run completes `sim_time = 60 s` and exits with code 0. The output folder `results/` holds `syn-0000.wd`, `syn-0001.wd`, `syn.max` and `syn.mass`.
- GPU runs (`cuda`) also print a final line containing `END state=COMPLETED time_s=60`, and write `results/syn.audit.log` and `results/syn.water_balance.csv`; `combined_residual_m3` in the `PROGRESS` lines of the audit log stays at round-off level.

The case checks input wiring, not hydrological realism.
