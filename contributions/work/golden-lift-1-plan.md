---
title: GOLDEN-LIFT-1 — abstract PSL2(F5) and canonical SL2(F5) cover
status: work-plan
claim_level: planning
oph_dependency: true
physical_claim: false
---

# GOLDEN-LIFT-1

## Goal

Link the repository's explicit six-axis action to Mathlib's abstract projective special linear group, and expose the canonical special-linear quotient interface. The typed target of this lot is deliberately limited to the abstract group and the six-axis subgroup:

```text
SL(2,F5) → PSL(2,F5) ≃ SixAxisGroup ≤ Equiv.Perm (Fin 6)
```

The first arrow is the canonical quotient by the center. The isomorphism must identify the concrete action, not be inferred from cardinality or group order. The existing port bridge remains a separate, pointwise comparison of the sixty indexed rows; it is not yet a group homomorphism that can be appended to this chain.

## Existing exact ingredients

OPH already has:

- `Lean/Screen/A5SixAxes.lean`: the explicit 60-row action on the six points conventionally labelled `{0,1,2,3,4,infinity}`, with designated permutations `t : z ↦ z+1` and `s : z ↦ -1/z`; it also proves list membership, distinctness, identity membership, and closure under multiplication and inverse.
- `Lean/Screen/A5PortSixAxesBridge.lean`: an exact pointwise conjugacy between the antipodal quotient of the 12-port action and that six-axis model, plus injectivity of the indexed quotient rows. It provides `quotientAxis`, `rowEquiv`, `bridged_axis_eq_six_axis`, and `quotient_action_faithful`; it does not currently package the port rows as a group or provide a group isomorphism.
- `Lean/Screen/GoldenSectorCharacters.lean`: the golden `3+3'` character sector evaluated on the same committed port-action rows.

On the repository's pinned Mathlib revision, the usable ingredients are:

- `Matrix.ProjectiveSpecialLinearGroup`, defined as `SL / center(SL)`;
- the generic projectivization action in `Mathlib.LinearAlgebra.Projectivization.Action`;
- `Matrix.SpecialLinearGroup.toLin'`, `MulAction.compHom`, and `MulAction.toPermHom` for constructing the `SL` permutation action;
- `Projectivization.map` and `Projectivization.map_comp` for reasoning about that action;
- `OnePoint.equivProjectivization` and its affine/infinity evaluation lemmas, together with `LinearEquiv.finTwoArrow`, for the coordinate relabelling;
- `QuotientGroup.mk'`, `QuotientGroup.mk'_surjective`, `QuotientGroup.ker_mk'`, `QuotientGroup.lift`, and `QuotientGroup.quotientKerEquivRange` for the quotient and faithful descended action;
- `Equiv.permCongrHom` for conjugating a permutation action;
- `Matrix.SpecialLinearGroup.mem_center_iff` for the center calculation.

There is no pre-packaged, specialized faithful action of `PSL2F5` with the required kernel theorem on the pinned revision. GOLDEN-LIFT-1 must construct the action, prove its kernel, and descend it locally.

## Preferred Lean surface

Suggested file:

```text
Lean/Screen/PSL2F5SixAxesBridge.lean
```

Suggested namespace and required setup:

```lean
open scoped LinearAlgebra.Projectivization

namespace OPH.PSL2F5SixAxesBridge

local instance : Fact (Nat.Prime 5) := ⟨by decide⟩

abbrev F5 := ZMod 5
abbrev SL2F5 := Matrix.SpecialLinearGroup (Fin 2) F5
abbrev PSL2F5 := Matrix.ProjectiveSpecialLinearGroup (Fin 2) F5
abbrev P1F5 := ℙ F5 (Fin 2 → F5)
```

The prime fact is needed before the field/division-ring instances for `ZMod 5` can be synthesized. Use Lean's Unicode homomorphism arrow `→*`, not ASCII text that merely resembles it.

## Theorem targets

### 1. Canonical quotient interface

Expose the canonical homomorphism

```lean
def slToPsl : SL2F5 →* PSL2F5 :=
  QuotientGroup.mk' (Subgroup.center SL2F5)
```

and prove or wrap:

```lean
slToPsl_surjective
slToPsl_ker_center
```

This is a typed specialization of `QuotientGroup.mk'`, `QuotientGroup.mk'_surjective`, and `QuotientGroup.ker_mk'`.

### 2. Explicit projective-line relabelling

Construct an explicit equivalence

```lean
p1EquivSix : P1F5 ≃ Fin 6
```

using the OPH convention

```text
[z:1] ↔ z,  z ∈ F5
[1:0] ↔ infinity = 5.
```

Build this from `OnePoint.equivProjectivization`, `LinearEquiv.finTwoArrow`, the finite equivalence for `ZMod 5`, and the equivalence which places the extra point last in `Fin 6`. Prove explicit evaluation lemmas for every affine `z` and for infinity. Do not use an arbitrary `Fintype.equivFin`; the coordinate convention must be load-bearing.

### 3. Construct and descend the projective action

First construct the `SL2F5` action on `P1F5` from `Matrix.SpecialLinearGroup.toLin'` and the generic projectivization action, for example through `MulAction.compHom`, then expose its permutation homomorphism with `MulAction.toPermHom`:

```lean
slProjectiveAction : SL2F5 →* Equiv.Perm P1F5
```

Prove locally that its kernel is the center:

```lean
slProjectiveAction_ker : slProjectiveAction.ker = Subgroup.center SL2F5
```

Both inclusions must be proved. For the nontrivial inclusion, constrain a matrix acting trivially on the projective line by its action on the explicit points infinity, zero, and one; these force it to be scalar. Either work directly with `Projectivization.mk_eq_mk_iff'` on representatives `[1,0]`, `[0,1]`, and `[1,1]`, or first prove that `p1EquivSix` intertwines the action induced by `Matrix.SpecialLinearGroup.toLin'` with the corresponding `OnePoint` action; the existing `OnePoint` action lemmas do not provide that bridge automatically.

Descend the action through `QuotientGroup.lift`, using the center-to-kernel proof, and prove the descended homomorphism injective with `QuotientGroup.ker_lift`; alternatively use `QuotientGroup.quotientKerEquivRange` and explicitly transport its domain along `slProjectiveAction_ker`, for example with `QuotientGroup.quotientMulEquivOfEq`. Finally conjugate it through `p1EquivSix` using `Equiv.permCongrHom`:

```lean
pslProjectiveAction : PSL2F5 →* Equiv.Perm P1F5
pslToSix : PSL2F5 →* Equiv.Perm (Fin 6)
pslToSix_injective : Function.Injective pslToSix
```

Do not assume an unverified specialized PSL action or faithfulness theorem.

### 4. Match the OPH generators and action

Verify directly that the classes of the standard determinant-one lifts

```text
T = [[1, 1], [0, 1]]
S = [[0, -1], [1, 0]]
```

act, under the explicit relabelling, as the existing OPH generators

```text
A5SixAxes.t = [1,2,3,4,0,5]
A5SixAxes.s = [5,4,2,3,1,0].
```

Create a typed subgroup

```lean
SixAxisGroup : Subgroup (Equiv.Perm (Fin 6))
```

whose carrier is proved equivalent to membership in `A5SixAxes.L60`. This construction is mandatory: `L60` is currently a `List`, not a group type. Reuse the existing receipts `one_mem`, `inv_closed`, `mul_closed`, `length_L60`, `el_apply`, and the raw multiplication/inverse tables rather than rechecking expensive equality of `Equiv.Perm` values.

Then prove the image of `pslToSix` is exactly `SixAxisGroup`. Either of these proof routes is acceptable:

- match the explicit lifts of `t` and `s`, prove `Subgroup.closure {A5SixAxes.t, A5SixAxes.s} = SixAxisGroup`, **and** prove that the corresponding classes generate `PSL2F5` (or otherwise prove the reverse image inclusion); or
- give an exhaustive, kernel-checked image equality with witnesses in both directions.

The current corpus proves only `t_mem` and `s_mem`; it does not yet prove that these two elements generate all of `L60`. A cardinality-60 comparison by itself is not an identification proof.

Headline statement:

```lean
psl_equiv_six_axis_group : PSL2F5 ≃* SixAxisGroup
```

### 5. Center of SL(2,F5)

Use `Matrix.SpecialLinearGroup.mem_center_iff` specialized to dimension two over `ZMod 5` to prove that the center consists exactly of scalar `+I` and `-I`.

Desired interface:

```lean
center_eq_plus_minus_one
center_card_two
```

Then record the exact short-cover statement

```text
1 → {+I,-I} → SL(2,F5) → PSL(2,F5) → 1
```

at the theorem/interface level supported by Mathlib's quotient.

For abstract cardinalities, prefer `Nat.card` and quotient-cardinality theorems. `PSL2F5` has a `Finite` instance on the pinned revision but not a default computational `Fintype`/`DecidableEq` pair. If exhaustive quotient computation is used, install a local decidable center-membership instance via `Matrix.SpecialLinearGroup.mem_center_iff`; do not silently rely on unavailable quotient instances.

It is acceptable to prove `|SL(2,F5)| = 120` by finite enumeration or formula and derive `|PSL(2,F5)| = 60`; cardinality must remain a consistency check, not the identification proof.

## Scope boundary

This work MAY claim, if proved:

- a concrete isomorphism between Mathlib's `PSL(2,ZMod 5)` and the committed OPH six-axis subgroup;
- that the canonical `SL(2,ZMod 5) → PSL(2,ZMod 5)` quotient has kernel equal to the center;
- that in this specialization the center is `{+I,-I}`.

This lot may also cite the existing pointwise port comparison after proving the typed six-axis isomorphism, but it may not present that comparison as another group-homomorphism arrow. A later lot may package the port rows as a group and prove a typed port/six-axis isomorphism.

This work MUST NOT claim without additional proof:

- `PSL(2,5) ≅ A5` as an abstract classification theorem;
- `SL(2,5) ≅ 2I` or a binary icosahedral identification;
- McKay `E8 ↔ 2I`;
- that the golden OPH `3` and `3'` sectors are the McKay defining representation;
- a derivation or selection of `phi`;
- `27^phi` or any mass law;
- any physical rotation, family, coupling, or observable identification.

## Important representation-theoretic warning

The current Golden modules evaluate matrices on raw committed row indices; they do not yet define a representation of `PSL2F5`, prove compatibility with multiplication, or prove a pullback theorem along `slToPsl`. In particular, `GoldenSectorComplexIrreducibility` records the convention `gR(comp p q) = gR q * gR p`, so the direction of any future representation transport must be handled explicitly.

Once a typed Golden representation of the six-axis or port group and its transport through `PSL2F5` have been proved, its pullback along `SL2F5 → PSL2F5` will have trivial action of the central `-I`. That conclusion is conditional on the missing typed transport and is not a theorem delivered by this lot.

The McKay `E8` correspondence for the binary icosahedral group is driven by a faithful two-dimensional representation on which the center is nontrivial. Therefore a future McKay bridge cannot simply identify the OPH golden triplets with the McKay defining doublet; it must prove how the triplets arise inside the representation ring generated by that doublet.

This distinction is part of the no-cheating boundary for GOLDEN-LIFT-2/3.

## Validation contract

Before any PR:

1. add `PSL2F5SixAxesBridge` to the explicit `OPHScreen` roots in `Lean/lakefile.lean`; import it from `Lean/Screen/OPHScreen.lean` only if it is intended to be part of that umbrella surface;
2. build the new module directly, then build `OPHScreen`;
3. add `#print axioms` for the quotient-kernel theorem, descended-action injectivity, concrete image equality, center theorem, and final group isomorphism;
4. extend `tools/test_third_wave_surfaces.py`, the existing owner of the Golden/port bridge surface, so it checks the module path, headline declarations, axiom prints, Lake root, proof-index entry, and the stated no-cheating boundary;
5. add or retain negative fixtures proving that the gate fails when each required surface is removed or weakened;
6. update `Lean/docs/PROOF_INDEX.md` only after the kernel build passes;
7. update the standalone statement in the owner paper and synchronize the relevant claim-registry `statement`, `oph_specific_delta`, `evidence`, `status`, and `falsifier` fields only to the exact theorem level proved; synchronize the novelty and falsification matrices whenever that claim's scope, delta, or falsifier changes, and synchronize the dependency graph whenever an ID or dependency changes, as required by `claims/README.md`;
8. make `tools/test_third_wave_surfaces.py` check those paper/registry boundary updates rather than accepting an abstract-group promotion on Lean evidence alone;
9. regenerate the inventory with `python3 tools/check_axiom_consistency.py --inventory` and commit the resulting `claims/active_surface_inventory.json` delta;
10. run the standard `python tools/run_mandatory_suite.py`. Do not trigger opt-in heavy certificates or full pools unless this implementation actually changes their owned inputs.

No `sorry`, `admit`, or new axioms.

## Suggested commit sequence

1. `Add abstract PSL2(F5) six-axis bridge`
2. `Document PSL2(F5) cover interface`

The present file is a planning seed only and is not scientific evidence.
