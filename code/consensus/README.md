# Consensus Code

This directory contains OPH packet-net and consensus-protocol artifacts.

Canonical runners:

- `export_verified_tree_packet_net.py`: exports the verified tree packet-net repair domain.
- `reference_architecture_benchmark_suite.py`: runs the fixed-cutoff Z2/S3 reference-architecture benchmark suite for issue #237.

Run the benchmark suite from the repo root:

```bash
python3 code/consensus/reference_architecture_benchmark_suite.py
python3 -m pytest code/consensus/test_reference_architecture_benchmark_suite.py
```

The current emitted suite artifact is `code/consensus/runs/reference_architecture_benchmark_suite_current.json`.
It is a fixed-cutoff analytic benchmark surface only; it does not claim continuum/gravity closure or uniqueness of the microscopic UV completion.

## Finite Repair-Projection Receipt

For the finite conditional-expectation claim, use the identifications

```text
protected observation B = repaired datum rho_C
full-support weight mu = stationary weight pi_r
fiber-resampling projector P_B = E_C
```

The tested transition matrix must be extracted independently from the declared
local transition table or from frozen transition counts. It must not be built
from the target conditional-expectation formula. A valid receipt records the
ordered quotient-state list, `rho_C` value for every state, exact or
outward-rounded stationary weights, the raw transition source, the extracted
row-stochastic matrix, tolerances, and these recomputed checks:

- `R1_fiber_support`: a positive transition never changes `rho_C`;
- `R2_equal_rows_in_fiber`: starting states with the same `rho_C` have the
  same transition row;
- `R3_weighted_detailed_balance`: `pi_r[x] * K[x,y]` equals
  `pi_r[y] * K[y,x]` within the declared exact or interval arithmetic.

On finite full support, R1--R3 recognize exactly the stationary-weighted
fiber-resampling kernel. The resulting operator is idempotent, self-adjoint in
`L2(pi_r)`, contractive, and has the `rho_C`-measurable functions as its fixed
space. The receipt does not establish a spectral gap, identify a different
repair dynamics with conditional expectation, or supply a continuum/GNS
transfer theorem.

A simulator artifact carrying this claim should expose at least:

```json
{
  "claim": "finite_repair_is_conditional_expectation",
  "state_order": [],
  "rho_C_by_state": {},
  "stationary_weights": {},
  "transition_source": {
    "kind": "declared_local_transitions_or_frozen_counts",
    "sha256": ""
  },
  "extracted_transition_matrix": [],
  "arithmetic": {"mode": "exact_or_outward_rounded", "tolerance": 0.0},
  "checks": {
    "row_stochastic": false,
    "full_support": false,
    "R1_fiber_support": false,
    "R2_equal_rows_in_fiber": false,
    "R3_weighted_detailed_balance": false
  },
  "target_formula_used_to_construct_matrix": false
}
```

The claim fails closed unless every check is recomputed from the attached raw
objects and the final anti-circularity field is `false`.

## License

This code surface is part of the OPH public repository. See the main
[LICENSE](../../LICENSE).
