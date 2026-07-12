# OPH Code

This directory is the canonical code surface for OPH derivation work.

Primary derivation surfaces:

- `P_derivation/`: compressed `P -> alpha` closure and fixed-point artifacts.
- `particles/`: particle-spectrum builders, status surfaces, gap campaigns, and
  the electroweak hierarchy certificate bundle under `particles/hierarchy/`.
- `consensus/`: packet-net, consensus-protocol, and fixed-cutoff Z2/S3 reference-architecture benchmark artifacts.
- `geometry/`: machine receipts for the Einstein branch-entry packets (GitHub
  #523/#524/#525): the quotient-intrinsic incidence/cap readout with its
  topology countermodels, the null-net standardness and Markov modular-locality
  witnesses (including the Gibbs nonlocality counterexample), and the synthetic
  end-to-end event-manifold reconstruction with its dimension/Hausdorff
  countermodels, plus the Einstein-closure receipts (#526-#528, #503): null
  tomography of the stress tensor, the exact first-law split with edge term,
  the MaxEnt multiplier identity, baseline countermodels, and the
  realized-branch receipt evaluation
  (`python3 code/geometry/realized_branch_receipts.py` writes
  `code/geometry/runs/realized_branch_receipt_report.json`), and the realized
  cyclic cap-net repair run (`python3 code/geometry/cyclic_cap_net_run.py`):
  a genuine transactional repair tower whose repaired output passes the
  spherical-incidence/mesh/naturality receipts under explicit branch
  selection, and the free-fermion modular-clock instrumentation
  (`python3 code/geometry/modular_clock_instrumentation.py`, needs mpmath):
  120-digit arc entanglement Hamiltonians on the boundary collars whose
  resummed velocity profiles witness the 2-pi-KMS and modular cross-ratio
  receipts, and the null-net instrumentation
  (`python3 code/geometry/null_net_receipts.py`): NTI, weak additivity,
  separating moduli, mixed-GNS Cauchy, one-particle HSM compression, and
  percent-level modular Lie closure (rate open), and the realized-event
  instrumentation (`python3 code/geometry/realized_event_receipts.py`):
  E1/E2/E4 plus an intrinsic Lorentzian (1,2) ancestry cone on records
  from the actual repair dynamics, with the realized (1+2) screen-sheet
  dimension measured honestly, and the bulk-depth channel
  (`python3 code/geometry/bulk_depth_receipts.py`): multi-scale records
  with the produced depth dictionary witness E3 rank four, bulk PCA four,
  and a seed-stable Lorentzian (1,3) ancestry cone on the
  corpus-consistent inverse-system dynamics (narrow-margin caveats
  explicit; the strong-coupling variant measures (2,2) and is kept as a
  countermodel). Current verdict: all structural receipt families
  realized at finite/one-particle level; margin convergence, Cyc limit,
  Lie-closure rate, second-quantized MI, cap-interior, and physical
  families pending; realized geometric branch not certified nonempty.
  Run with `python3 -m pytest code/geometry/`.

Supplemental and benchmark surfaces:

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
