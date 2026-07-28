# OPH-native W/Z source packet

This package owns the producer-side frontier of GitHub issue #594. It is
separate from the external Standard-Model calculation in issue #593.

The committed inventory binds the finite source parents that can be consumed
without a physical W/Z claim:

| Parent | Bound result | Boundary retained |
| --- | --- | --- |
| #565 | Twelve-port carrier, incidence, inverse pairing, A5 frame, and finite refinement maps | No universal carrier selection or electroweak load identification |
| #566 | Finite source-model current algebra and response representation | No continuum or laboratory current, coupling, or kinetic normalization |
| #314 | Rank-fifteen chiral matter pair, finite Spin typing, anomaly checks, and conditional Yukawa channel lines | No scalar source, numerical Yukawa matrices, families, or pole masses |
| #567 | Finite diagonal Z6 form, common kernel, lattices, and screen sector transport | No continuum instanton, theta, monopole, or laboratory attachment |

The rank-three family-band result from #569 is carried as
`conditional_context_only`. Its finite screen-response resolvent is not a
physical matter pole. Matter-pole identification, family chirality and
locality, refinement, and laboratory-current attachment remain open.

The strongest machine status is:

```text
FINITE_SOURCE_PARENTS_BOUND__NATIVE_ACTION_AND_PHYSICAL_ATTACHMENT_OPEN
```

`promotion_allowed` is fixed to `false`. The inventory emits no pole
coordinate and no physical-unit value.

## Producer work queue

| ID | Work item | Status | Main owners |
| --- | --- | --- | --- |
| P594-01 | Hash-bind and independently resolve the finite source-parent inventory | Complete | #565, #566, #314, #567 |
| P594-02 | Define the target firewall, closed input allowlist, separate comparison boundary, and sealed-input replay | Partial: path/content closure and sealed-input replay pass; runtime transcript and human-selection ancestry remain open | #594 |
| P594-03 | Produce a typed event-base and source action with a scoped complete field/operator census | Blocked | #503, #547, #594; the #609 sterile countermodel must remain explicit |
| P594-04 | Complete the physical family, seam-action, and matter-pole attachment | Blocked | #569 |
| P594-05 | Prove the source-Higgs coordinate to FJ-coordinate map | Blocked | #34, #547 |
| P594-06 | Emit full target-clean Yukawa matrices and coefficient law | Blocked | #34, #569, #594 |
| P594-07 | Emit target-clean running, thresholds, finite maps, Jacobians, masks, and remainders | Blocked | #32 |
| P594-08 | Emit a unique point law or target-independent joint law and pole covariance | Blocked | #545, #594 |
| P594-09 | Substitute the native packet into the frozen QFT calculator without changing its algorithms and recompute the common subject digest | Blocked | #593 |
| P594-10 | Attach physical units through a source-derived operational clock | Blocked | #334 |
| P594-11 | Run clean-room replay, then mount known W/Z values in a separate read-only comparison process | Blocked by P594-03 through P594-10 | #594 |

The GitHub issue is the live status surface. Acceptance boxes stay unchecked
until their full clauses are discharged. Progress comments should cite the
inventory digest and distinguish complete, partial, and blocked rows.

## Files

- `data/source_parent_policy_v1.json`: trusted role, path, scope, and firewall
  policy.
- `schemas/source_parent_inventory_v1.schema.json`: closed output schema.
- `build_source_parent_inventory.py`: deterministic producer.
- `check_source_parent_inventory.py`: independent resolver. It imports no
  producer code and can replay the five native parent verifiers.
- `outputs/source_parent_inventory.json`: committed non-promoting frontier
  receipt.
- `tests/test_source_parent_inventory.py`: deterministic, tamper, target-leak,
  coordinate, unit, family-boundary, and DAG mutations.

## Run

Fast local replay:

```bash
python3 build_source_parent_inventory.py --check-byte-exact
python3 check_source_parent_inventory.py --skip-native-verifiers
python3 -m pytest -q tests/test_source_parent_inventory.py
```

Full parent replay:

```bash
python3 run_all.py
```

The full checker executes the production verifiers for #565, #566, #314,
#567, and the conditional #569 family-band certificate. It does not invoke or
modify the #593 calculator.

## Coordinate and unit boundary

`v_chart` and `v_F` are different typed coordinates. The equality receipt is
absent and relabeling is rejected. Before the source clock closes, a native
packet may expose only dimensionless coordinates normalized by `E_star`.

The current external-validation schemas are pinned only as provisional
consumer shapes. They are not a frozen native interface, and their GeV fields
cannot be populated by this producer.

