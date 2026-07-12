# OPH Code

This directory is the canonical code surface for OPH derivation work.

Primary derivation surfaces:

- `P_derivation/`: compressed `P -> alpha` closure and fixed-point artifacts.
- `particles/`: particle-spectrum builders, status surfaces, gap campaigns, and
  the electroweak hierarchy certificate bundle under `particles/hierarchy/`.
- `consensus/`: packet-net, consensus-protocol, and fixed-cutoff Z2/S3 reference-architecture benchmark artifacts.
- `geometry/`: machine receipts for the Einstein branch-entry packets
  (GitHub #503, #523-#528). The suite covers the quotient-intrinsic
  incidence/cap readout with its topology countermodels, the null-net
  standardness and Markov modular-locality witnesses (including the Gibbs
  nonlocality counterexample), the synthetic event reconstruction with its
  dimension and Hausdorff countermodels, the realized cyclic cap-net repair
  run, the free-fermion modular-clock instrumentation (120-digit arc
  entanglement Hamiltonians, resummed velocity profiles), the null-net
  receipt instrumentation (half-sided compression, modular Lie closure),
  the realized-event receipts, and the bulk-depth channel with its
  Lorentzian (1,3) cone verdict and (2,2) strong-coupling countermodel.
  Machine state: sphere/mesh/naturality receipts pass on a genuine repair
  run under explicit branch selection; boundary-collar modular receipts
  pass; null-net receipts pass at one-particle level; realized records
  produce (1,2) and (1,3) causal cones; cone-margin convergence, limit
  clauses, cap-interior data, and the physical-identification receipts are
  work in progress, and the realized geometric branch is not certified
  nonempty. Reports live in `code/geometry/runs/`. Run with
  `python3 -m pytest code/geometry/` (needs numpy, scipy, mpmath).
- `dark_matter/`: pre-likelihood dark-sector simulation and likelihood scaffolds
  for the dark-matter supplement.
- `ibm_quantum_cloud/`: completed IBM hardware benchmark scripts and redacted
  raw-result JSON for OPH-inspired reduced-sector states. This is a benchmark
  bundle, not a theorem input to the particle or cosmology papers.

For the reference-architecture benchmark suite:

```bash
python3 code/consensus/reference_architecture_benchmark_suite.py
python3 -m pytest code/consensus/test_reference_architecture_benchmark_suite.py
```

The emitted artifact is `code/consensus/runs/reference_architecture_benchmark_suite_current.json`.

For the particle program, start with:

```bash
cd particles
python3 compute_current_output_table.py --no-print-table --show-paths
python3 scripts/build_derivation_gap_ledger.py
```

The particle gap campaign lives at `particles/campaigns/gap_bundle/`.

## License And Patent Policy

This code surface is part of the OPH public repository. See the main
[LICENSE](../LICENSE) and [OPH Open Use And Anti-Patent Covenant](../PATENTS.md).
OPH-derived software, simulations, algorithms, hardware methods, applications,
and implementation patterns may be studied, tested, implemented, modified,
deployed, manufactured, and shared, but may not be used to create patent claims
that restrict others from practicing OPH-derived work.
