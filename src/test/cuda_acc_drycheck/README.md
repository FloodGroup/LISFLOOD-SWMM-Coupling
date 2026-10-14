# CUDA-ACC dry-cell conservation regression

This two-by-two closed-domain case starts with one metre of water in each
left-hand cell and a 50 m bed drop toward the dry right-hand cells.  The first
ACC step therefore requests much more discharge than the donor cells store.

Expected V4 behaviour after one second:

- no negative depths;
- total volume remains 200 m3 (within floating-point tolerance);
- the limited shared face flux is used by both donor and receiver cells.

The pre-V4 CUDA dry check scales only thread-local flux copies.  It removes
200 m3 from the donors but adds the original, much larger flux to the receivers,
so this case is also an A/B regression against the historical implementation.

Measured Release results on 2026-09-14:

- historical V3: 200.0 m3 -> 1000.3 m3, volume error -800.27 m3;
- V4: 200.0 m3 -> 200.0 m3, volume error 0.0 m3.

The raw outputs are retained in `results_v3_historical` and `results`.

## CPU-net-flux direction test

The three-cell chain cases exercise a wet middle cell with simultaneous inflow
and outflow.  The original serial CPU `DryCheck` conserves 400 m3 when flow is
aligned with its row-major scan, but the reversed chain ends at approximately
1000.27 m3.  A two-phase CUDA translation of that net-flux formula also ends at
approximately 1000.27 m3 because it credits inflow that an upstream donor later
limits.

The retained strict V4 limiter conserves 400 m3 in both chain directions with
zero reported volume error.  Raw results are in `results_chain_cpu_serial`,
`results_chain_reverse_cpu_serial`, `results_chain_cpu_formula`,
`results_chain_strict`, and `results_chain_reverse_strict`.
