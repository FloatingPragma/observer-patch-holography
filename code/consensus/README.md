# Consensus Code

This directory contains OPH packet-net and consensus-protocol artifacts.

Canonical runners:

- `export_verified_tree_packet_net.py`: exports the verified tree packet-net repair domain.
- `reference_architecture_benchmark_suite.py`: runs the fixed-cutoff Z2/S3 reference-architecture benchmark suite for issue #237.
- `verify_issue_517_proof_obligations.py`: recomputes the finite
  conflict-component, prepared-lock BFT, refinement-modulus, \(\ell^p\), and
  selector/separation receipts used by proof-audit issue #517.

Run the benchmark suite from the repo root:

```bash
python3 code/consensus/reference_architecture_benchmark_suite.py
python3 -m pytest code/consensus/test_reference_architecture_benchmark_suite.py
```

The current emitted suite artifact is `code/consensus/runs/reference_architecture_benchmark_suite_current.json`.
It is a fixed-cutoff analytic benchmark surface only; it does not claim continuum/gravity closure or uniqueness of the microscopic UV completion.

Run and verify the issue #517 receipt from the repo root:

```bash
python3 code/consensus/verify_issue_517_proof_obligations.py emit \
  --output code/consensus/runs/issue_517_proof_obligations.json
python3 code/consensus/verify_issue_517_proof_obligations.py verify \
  --receipt code/consensus/runs/issue_517_proof_obligations.json
python3 -m pytest code/consensus/test_issue_517_proof_obligations.py
```

`TXN-DIAMOND-1` exhausts a finite reference engine; it is evidence that the
declared read/write, support-reclosed component-merge, prepared-batch,
protected-record, and descent contract is realizable, not a substitute for checking that contract in another
engine. `BFT-LOCK-1` exhausts quorum intersections for the displayed
\((n,f,q)\) instances, parses finite certificate/view-change traces, and
executes reference honest-progress phases; the arbitrary-view and wall-clock
arguments remain the paper theorem. The refinement controls are finite
descriptions of parametric infinite counterfamilies, not fixed-depth
surrogates for a limit; the arbitrary-depth positive telescope remains the
mathematical/Lean theorem. The receipt
also recomputes the independent twelve-port selector certificate and keeps
coefficient reconstruction, physical currents, global-form descent, and
matter realization in their independent receipt classes.

Formalization status is deliberately split:

- Transactional local confluence has a paper proof and an exhaustive Python
  reference-engine receipt. Its general support-closure and prepared-batch
  theorem is not formalized in Lean, TLA+, or a protocol model checker.
- Prepared-lock BFT has a paper proof plus finite Python certificate,
  next-view, negative-control, and honest-progress traces at three parameter
  points. Its arbitrary-view protocol theorem is not Lean/TLA+/model-checker
  verified.
- `Lean/ObservableNormalForms/Refinement.lean` proves the metric telescope.
  The family-uniform inverse/residual moduli and cofinal-limit quantifiers are
  paper proofs with executable witness schemas, not separate Lean theorems.

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
