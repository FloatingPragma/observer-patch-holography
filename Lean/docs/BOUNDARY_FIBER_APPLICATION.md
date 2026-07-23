# Boundary-fiber identifiability: the OPH application bridge (issue #304)

This note connects the substrate-neutral identifiability theorem in
`ObservableNormalForms/` to the concrete OPH objects, hypothesis by
hypothesis, and records exactly what is proved, on which domain, and what
stays a named premise. The machine-checked side lives in
`ObserverPatchHolography/BoundaryFiber.lean`.

## The gate statement

Issue #304 leaves one application obligation after the generic layer:

> Define the declared physical boundary/sector map `B_OPH` on the intended
> consistent-state domain and prove: `u, v ∈ C_OPH` and
> `B_OPH(u) = B_OPH(v)` imply `u ~gauge v`.

## The three OPH objects

| Gate object | Concrete definition | Source |
|---|---|---|
| Intended domain | `TreePacketNet` / `TreePacketNet.carrier`: the verified rooted-tree packet-net domain, per-patch state `A × K_i`, interfaces reading the packet component only | *Reality as a Consensus Protocol*, Definition `def:tree-packet-domain`; `BoundaryFiber.lean` |
| `B_OPH` | `TreePacketNet.BOPH`: the root-packet readback `x ↦ (x root).1`, the "root-packet map" entry of the protected-data taxonomy | Definition `def:finite-quotient-repair-presentation`; `BoundaryFiber.lean` |
| `~gauge` | `gaugeEquiv` on `Records`, the kernel of `obsMap`; the hidden labels `K_i` are the declared redundancy data and move inside one class | `Primitives.lean`; non-triviality witness `fourVertexNet_gauge_nontrivial` |

## The generic theorem and its hypotheses, checked one by one

The neutral theorem
`ObservableNormalForms.boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo`
states: under (G1) and (G2), boundary identification modulo `E` on the
consistent set is equivalent to cross-source normal-endpoint uniqueness
modulo `E`. Instantiation: `Q := Records T.carrier`,
`C := {x | Consistent T.carrier x}`, `r := T.treeStep` (the paper's tree
repair `T_i` as an accepted-step relation), `B := T.BOPH`,
`E := gaugeEquiv T.carrier`.

| Hypothesis | Content | Status on the declared domain |
|---|---|---|
| G1 `ObservationPreserving r B` | every accepted step preserves the boundary readback | **Proved**: `TreePacketNet.treeStep_observationPreserving` (tree repairs never write the root) |
| G2 `CompleteFor r C` | normal forms are exactly the consistent records | **Proved**: `TreePacketNet.treeStep_completeFor` (each edge is the parent edge of one non-root vertex; the copy move at `i` fires iff that edge is broken) |
| G3 `BoundaryIdentifiesModulo C B E` | `B_OPH` injective modulo gauge on the consistent set | **Proved**: `TreePacketNet.BOPH_boundaryIdentifiesModulo`, from the gate theorem `TreePacketNet.BOPH_injective_modulo_gauge` |

Since G1, G2, and G3 all hold, the equivalence yields the payoff form
`TreePacketNet.BOPH_observerEndpointUnique`: any two records with equal
`B_OPH` readback settle, along any maximal accepted tree-repair schedules, to
gauge-equal normal forms. Confluence enters nowhere. Weak normalization is
needed only for endpoint existence; termination of the tree repair is the
paper's Theorem `thm:tree-packet-domain` (Lyapunov weights) and is not
re-proved in the Lean module.

## The gate theorem as proved

`TreePacketNet.BOPH_injective_modulo_gauge`: for every rooted-tree packet net
`T` (every finite rooted tree, every packet alphabet, every hidden-label
family) and all `u, v : Records T.carrier`,

```
Consistent u  ∧  Consistent v  ∧  B_OPH(u) = B_OPH(v)   ⟹   gaugeEquiv u v.
```

Exact premise list: (1) `T` is a `TreePacketNet` (finite rooted tree with
strictly depth-decreasing parent map; per-patch state `A × K_i`; every
interface projection reads the packet component); (2) `u` and `v` are
consistent (`Φ = 0`, equivalently edge-consistent); (3) equal root-packet
readback. Nothing else: no repair law, no confluence, no schedule input.
The proof consumes consistency through
`packet_eq_root_of_edgeConsistent` (strong induction on tree depth
propagates the root packet through the bulk), and the conclusion is
observable equality, not record equality: `fourVertexNet_gauge_nontrivial`
shows two distinct same-boundary consistent records that the theorem
correctly identifies only modulo gauge.

## Domain boundary: what stays a named premise

The unconditional statement "for every finite patch net and every declared
boundary map, `B` is injective modulo gauge on the consistent set" is
**false**, with machine-checked countermodels in this tree:
`demoCarrier_Hfib_fails` (trivial boundary), `rule90_Hfib_bad_fails`
(deficient information set), and `fourVHiddenReadback_not_identifying`
(readback missing the protected packet, on the verified domain itself).

Consequences for claim surfaces:

- On the declared rooted-tree packet-net domain, physical same-boundary
  uniqueness is a theorem; no premise remains there.
- The paper-level generalization to single-parent deterministic extension
  maps is Theorem `thm:functional-selected-fiber`; its Lean transcription is
  not in the module (the packet-copy rule is the realized instance). The
  multi-parent layered class is Theorem `thm:layered-carrier-HB-Hfib`,
  likewise paper-level.
- For any patch net outside the declared class, or any boundary map other
  than the declared protected readback, boundary identifiability
  (`BoundaryIdentifiesModulo`, the `Hfib` binder of
  `boundary_fiber_observer_unique`) is a named per-net premise, and the
  countermodels above show it can genuinely fail. Claim surfaces that need
  same-boundary uniqueness beyond the declared class must carry that premise
  explicitly.

## Verification

```
cd <repository>/Lean
LC_ALL=C LANG=C lake build
```

The build passes with zero admissions; the former `Primitives.lean`
admissions (`localRepair`, `Repair`, `repair_respects_gauge`) are discharged,
and `BoundaryFiber.lean` adds none. The in-file `#print axioms` block reports
only `[propext, Classical.choice, Quot.sound]` for all eight audited
declarations; no `sorryAx`.
