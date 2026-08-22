# B20 random-scan common-reference preflight (issue #725)

Exact algebraic feasibility preflight mandated as the issue's first
step: before any fresh simulation, decide on retained source data
whether random-scan mixtures of overlapping conditional-resampling
projectors admit a nonconstant protected observable.

- `preflight_random_scan.py`: producer.  Recounts the pinned run
  `b12_prereg_16k_20260806` exactly (integers, then
  `fractions.Fraction`), forms the faithful visit-count reference, and
  for every subset of at least two committed packet fields builds the
  heat-bath projectors and their random-scan mixture under two declared
  schedulers (uniform and block-count-proportional). For the primary uniform
  scheduler it records row-stochasticity, exact stationarity, a
  non-idempotence witness, join component count, and exact fixed-space
  dimension. The second-scheduler packet records its rule and weights, exact
  stationarity, and fixed-space dimension. A second arena
  runs the same battery on the realized 256-state record/companion
  conditional-resampling structure.  Every retained step field beyond
  the packet fields is inventoried and must be constant, so the
  certified grammar covers all retained step fields.
- `runtime/b20_preflight_certificate.json`: the committed certificate.
  Vendors the state labels, visit counts, and arena-2 class table, so
  the full battery is recomputable offline; SHA-256 digests pin the
  five run inputs.
- `validate_random_scan.py`: independent validator with two layers:
  an offline algebra layer that recomputes every battery entry from the
  vendored objects through its own code path (independent connectivity
  routine, opposite-order elimination pivoting), and a custody layer
  that re-derives the vendored objects from the pinned run when the run
  directory is present (`--require-custody` makes its absence fatal).
- `test_b20_preflight.py`: semantic mutation suite (tampered pass
  flags, component counts, dimensions, witnesses, visit counts, state
  labels, arena-2 masses, verdict, reference, extra-field inventory all
  rejected) plus synthetic-instance oracle cross-checks of the
  connectivity and elimination routines.

Result: **negative** across both arenas. Every computed uniform-scheduler
mixture is row-stochastic, exactly stationary for the shared reference, and
non-idempotent. Under both declared schedulers every computed fixed space is
one-dimensional. The analytically covered constant-field mixtures inherit
only this fixed-space conclusion; no non-idempotence claim is made for them.
Thus the only protected observables in the certified grammar are constants.
The certificate states the no-go
grammar precisely and records the next viable route (an enriched export
with a disconnected join, or a dilated construction) rather than
declaring anything beyond the grammar impossible.  Per the issue
contract, a fresh prospectively frozen run is not triggered.
