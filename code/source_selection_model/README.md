# Scalar seams and M1 source selection

This package constructs a mathematical source with primitive central ports,
a complete observable twelve-dimensional response, same-response proper
chart transports, an all-level spherical refinement and an explicit A3
state problem. Distinct carriers share scalar seams and cannot transmit
their private records. The [derivation](DERIVATION.md) gives the analytic
model-membership argument and checks each clause of the current typed
axioms. The [contract](CONTRACT.md) makes that argument necessary for the
non-entailment result; passing finite tests is insufficient.

The resulting no-communication theorem covers every local program, retained
local memory and stopping rule, not only the finite sample histories. It is
a non-entailment construction for the stated A1--A3 data-access clauses,
not a no-go for M1 in every source, the canonical scalar-mean branch or the
registered amplitude simulator. The exact central-channel decoder criterion names
the missing transmission capability and distinguishes it from the separate
population, demand, preservation and execution choices.

The companion [RG sufficiency theorem](RECORD_GLUING.md) supplies the positive
exit: one explicit proposed source law gives the operational causal/count
limit, with M1 as an alternative regulator. RG is a substantial new transport
law with several explicit requirements, not a consequence of the original
axioms or a property already verified in the simulator. The
[reference execution](../source_record_gluing/README.md) checks its generated
process theory and keeps that boundary explicit.

From the repository root, put `code` on `PYTHONPATH` and run:

```sh
python -m source_selection_model.verify
python -m pytest -q code/source_selection_model
```

To regenerate the two compact exact source packets and their receipt:

```sh
python -m source_selection_model.response
python -m source_selection_model.records
python -m source_selection_model.verify --write-receipt
```

The response producer executes ordered jets and Cayley factors. Its
independent verifier reconstructs the generators by a different formula
and checks observable rank, all 66 brackets, 60 closed paths, 720 covariance
identities and 3,600 group compositions. It imports no source constructor.
The record verifier independently replays all 6,912 local events in 144
paired interventions. Deduplication retains twelve identical local tapes
once each without discarding any intervention. These are exact software
executions, not laboratory observations.

The support checker covers four levels with 12, 42, 162 and 642 carriers;
it verifies the five-valent seed and reconstructs every oriented midpoint
subdivision from incidence. Equal-count non-icosahedral spheres and altered
fine triangulations are rejected.
The all-level chain and mesh proofs are analytic. All 431 independent state
coordinates have positive-state omission witnesses. Scalar, noiseless,
redundant and invertible-but-noisy channels test exact record decoding.
Hostile-input tests include actual CLI receipt and source-byte tampering,
missing cases, forged paths, response coefficients, local events and
noncanonical numbers. The package also re-executes both producers and
compares their complete output with the retained packets.

From `Lean/`, run:

```sh
lake build Geometry.SourceSelectionLocalityAxiomAudit
```

Fifteen audited declarations cover all-program locality, transcript/stopping
obstructions, exact channel support and composition, and population
arithmetic. Full axiom-model membership, complete quantum instruments,
Umegaki minimization and geometric mesh convergence are analytic arguments;
they are not purported Lean results. CI runs the finite evidence and hostile
tests on Linux, with Lean checked by its separate workflow.
