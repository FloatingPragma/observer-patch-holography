# OPH Executable Evidence

This directory contains the certificates, simulations, exact finite
calculations, public receipts, and experimental adapters behind the OPH paper
stack. The papers state the claims. These artifacts let a reader recompute the
finite and numerical parts, inspect their inputs, and see which physical bridge
each result requires.

The strongest artifacts are small and exact: finite repair checks, the
twelve-port coefficient algebra, exterior Standard Model branching, anomaly
arithmetic, and interval uniqueness certificates. Larger simulations test
whether the same architecture survives scale, noise, refinement, and physical
readout. A simulation is evidence for its declared branch. It does not replace
the theorem or physical carrier named by that branch.

## Main Routes

| Surface | What to inspect |
| --- | --- |
| [`consensus/`](consensus/) | Finite patch-net repair, schedule independence, and reference-architecture benchmarks |
| [`geometry/`](geometry/) | Screen incidence, modular clocks, event reconstruction, causal cones, and Einstein-branch receipts |
| [`a5_closure/`](a5_closure/) | Exact twelve-port $A_5$ coefficient algebra and gauge-closure checks |
| [`particles/`](particles/) | Particle carriers, hierarchy, flavor, neutrino, hadron, and current status surfaces |
| [`P_derivation/`](P_derivation/) | Pixel fixed-point maps, interval contraction, uniqueness, and provenance |
| [`capacity_readback/`](capacity_readback/) | Finite record-capacity readback and conditional global closure |
| [`regulator_gluing/`](regulator_gluing/) | Regulator compatibility and refinement evidence |
| [`ibm_quantum_cloud/`](ibm_quantum_cloud/) | Redacted hardware runs and reproducible benchmark analysis |

## Supporting Modules

| Surface | Role |
| --- | --- |
| [`collar_alignment/`](collar_alignment/) | Boundary-collar alignment checks |
| [`maxent/`](maxent/) | Finite maximum-entropy reconstruction |
| [`edge_sectors/`](edge_sectors/) | Edge-sector and transport calculations |
| [`e8_triality_claim_statement/`](e8_triality_claim_statement/) | Exact exceptional-symmetry certificate specification |
| [`precision_ledgers/`](precision_ledgers/) | Machine-readable numerical provenance ledgers |

## Quick Reproduction

The repository-level [reproduction guide](../REPRODUCE.md) gives the clean-clone
environment and mandatory checks. The commands below are short entry points
into the main evidence families.

Run the finite consensus benchmark from the repository root:

```bash
python3 code/consensus/reference_architecture_benchmark_suite.py
python3 -m pytest code/consensus/test_reference_architecture_benchmark_suite.py
```

Inspect the current particle output surface:

```bash
cd code/particles
python3 compute_current_output_table.py --no-print-table --show-paths
python3 scripts/build_derivation_gap_ledger.py
```

The geometry suite is available through:

```bash
python3 -m pytest code/geometry/
```

Generated JSON receipts under `runs/` and `runtime/` are evidence artifacts when they are tracked beside the producing script and test. Local caches, virtual environments, logs, and archives are not repository material.

## License

The code in this tree is licensed under [Apache-2.0](LICENSE). The main [LICENSE](../LICENSE) gives the repository-wide licensing map.
