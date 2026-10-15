# CUDA-ACC dry-cell conservation regression

This two-by-two closed-domain case starts with one metre of water in each
left-hand cell and a 50 m bed drop toward the dry right-hand cells.  The first
ACC step therefore requests much more discharge than the donor cells store.

Expected behaviour after one second:

- no negative depths;
- total volume remains 200 m3 (within floating-point tolerance);
- the limited shared face flux is used by both donor and receiver cells.

Scaling only thread-local flux copies removes 200 m3 from the donors but adds
the original, much larger flux to the receivers. This case checks that the
shared-face limit is applied to both sides.

## CPU-net-flux direction test

The three-cell chain cases exercise a wet middle cell with simultaneous inflow
and outflow.  The original serial CPU `DryCheck` conserves 400 m3 when flow is
aligned with its row-major scan, but the reversed chain ends at approximately
1000.27 m3.  A two-phase CUDA translation of that net-flux formula also ends at
approximately 1000.27 m3 because it credits inflow that an upstream donor later
limits.

The strict limiter conserves 400 m3 in both chain directions with zero reported
volume error.  The parameter files are `chain_cpu_serial.par`,
`chain_reverse_cpu_serial.par`, `chain_cpu_formula.par`, `chain_strict.par`,
and `chain_reverse_strict.par`.
