# Expected behaviour (qualitative)

- Model should start, read DEM + SWMM INP, and map `OF1` as a `COUPLED` outfall and `J1` via `floodconfig.txt`.
- With zero initial surface water and no rainfall/inflow, exchange volumes remain ~0 and the run should complete `sim_time`.
- This case validates **file wiring**, not hydrologic realism.
