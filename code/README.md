# OPH Executable Evidence

This directory contains the certificates, simulations, exact finite
calculations, public receipts, and experimental adapters behind the OPH paper
stack and runnable model. The papers and model contracts state their respective
claims. These artifacts let a reader recompute the
finite and numerical parts, inspect their inputs, and see which physical bridge
each result requires.

The candidate executable model is [A-contrast-memory-spherical-v1](../model/candidate_laws/contrast_memory_spherical_v1/README.md),
whose contrast-memory law is a candidate native repair law. This directory
contains supporting component calculations; their source-specific hypotheses
are stated individually.

The strongest artifacts are small and exact: finite repair checks, the
twelve-port coefficient algebra, exterior Standard Model branching, anomaly
arithmetic, and interval uniqueness certificates. Larger simulations test
whether the same architecture survives scale, noise, refinement, and physical
readout. A simulation is evidence for its declared branch. It does not replace
the theorem or physical carrier named by that branch.

## Main Routes

| Surface | What to inspect |
| --- | --- |
| [`../model/candidate_laws/contrast_memory_spherical_v1/`](../model/candidate_laws/contrast_memory_spherical_v1/README.md) | The executable model: the contrast-memory law, the binding theorem, particle stability on agreed records, bound objects, port response and the unrecovered collar |
| [`consensus/`](consensus/) | Finite patch-net repair, schedule independence, and reference-architecture benchmarks |
| [`geometry/`](geometry/) | Screen incidence, modular clocks, event reconstruction, causal cones, and Einstein-branch receipts |
| [`m1_operational_clocks/`](m1_operational_clocks/) | Massive charged flights, all-band speed bounds, resolved local reads and an observable moving clock |
| [`a5_closure/`](a5_closure/) | Exact twelve-port $A_5$ coefficient algebra and gauge-closure checks |
| [`P_derivation/`](P_derivation/) | Pixel fixed-point maps, interval contraction, uniqueness, and provenance |
| [`capacity_readback/`](capacity_readback/) | Finite record-capacity readback and conditional global closure |
| [`cosmology/`](cosmology/) | Radial-lift, clock-certificate, posterior-chain and cosmology postdiction calculations |

## Supporting Modules

| Surface | Role |
| --- | --- |
| [`collar_alignment/`](collar_alignment/) | Boundary-collar alignment checks |
| [`maxent/`](maxent/) | Finite maximum-entropy reconstruction |
| [`edge_sectors/`](edge_sectors/) | Edge-sector and transport calculations |

## Quick Reproduction

The repository-level [reproduction guide](../REPRODUCE.md) gives the clean-clone
environment and mandatory checks. The commands below are short entry points
into the main evidence families.

Run the finite consensus benchmark from the repository root:

```bash
python3 code/consensus/reference_architecture_benchmark_suite.py
python3 -m pytest code/consensus/test_reference_architecture_benchmark_suite.py
```

Run the executable model and its tests:

```bash
python -m model
python -m pytest -q -W error model/candidate_laws/contrast_memory_spherical_v1
```

The geometry suite is available through:

```bash
python3 -m pytest code/geometry/
```

Generated JSON receipts under `runs/` and `runtime/` are evidence artifacts when they are tracked beside the producing script and test. Local caches, virtual environments, logs, and archives are not repository material.

## License

The code in this tree is licensed under [Apache-2.0](LICENSE). The main [LICENSE](../LICENSE) gives the repository-wide licensing map.
