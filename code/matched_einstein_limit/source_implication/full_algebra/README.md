# An entropy maximum on the complete two-head state

On the registered 80-cell radar experiment, a faithful preparation is the unique global maximum of stored-head entropy plus calibrated tetrahedral area at the measured comparison volume. The twelve directed ranges come from independent active-law probe effects. The gain comes from a separate ordinary exchange standard. The maximum covers all 31 parity-block coordinates, and pinching extends it to arbitrary coherences on the same eight-dimensional support. A positive mixture of native product preparations reaches the state. [Proof](PROOF.md) · [Certified bounds](certificate.json) · [State, source weights and held measurements](results.json)

Every other state at the same survey volume loses at least 0.14 times its squared trace-norm distance in the objective, measured in nats. A second faithful native preparation at that exact volume has a strictly positive first derivative between 7.47 and 7.49 nats per unit odd fraction. The source lotteries select which state the observers survey. A feedback calculation starts from the tracial state and uses the ordinary calibration and volume target to set those lotteries. For any number of disjoint copies, the product preparation uniquely maximizes the joint objective under both separate and pooled survey readings. The range cycle ratios at that equilibrium differ from one even after adjusting each station's clock scale. The same state and gain also fail stationarity when the clock timeout changes from six to four units.

On the same connected board, unconditioned, uniformly launched short probes relax the head coins to sector traces. One of these stationary states has the calibrated volume and a nonzero fixed-volume entropy-plus-area slope. Reusing the variational preparation for one full capped radar cycle decreases its left-coin reading on every charge and route. [Connected equilibration](PROOF.md#equilibration-on-one-connected-board) · [Complete radar readback](reuse/README.md)

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.station
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.certificate
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra
```

From the radar timestamps alone, 31 pre-cutoff return probabilities identify the full two-head state, and calibrated clock weights predict the independently measured signed force on both carriers. The twelve means and both force summaries together miss an entropy-increasing change. [Complete clock measurement and proof](PROOF.md#recovering-the-complete-source-from-the-original-clocks)
