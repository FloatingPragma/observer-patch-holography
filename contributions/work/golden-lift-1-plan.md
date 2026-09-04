---
title: GOLDEN-LIFT-1 — abstract PSL2(F5) and canonical SL2(F5) cover
status: work-plan
claim_level: planning
oph_dependency: true
physical_claim: false
---

# GOLDEN-LIFT-1

## Goal

Close the exact mathematical boundary left by `A5PortSixAxesBridge` by linking the repository's explicit six-axis `PSL(2,F5)` action to Mathlib's abstract projective special linear group, and expose its canonical special-linear double-cover interface.

Target chain:

```text
SL(2,F5) -> PSL(2,F5) -> OPH six-axis action -> antipodal port quotient
```

The first arrow is the canonical quotient by the center. The second arrow must be an explicit action equivalence, not an inference from cardinality or group order.

## Existing exact ingredients

OPH already has:

- `Lean/Screen/A5SixAxes.lean`: the explicit 60-row action on the six points conventionally labelled `{0,1,2,3,4,infinity}`, with generators `t : z -> z+1` and `s : z -> -1/z`.
- `Lean/Screen/A5PortSixAxesBridge.lean`: exact conjugacy between the antipodal quotient of the 12-port action and that six-axis model, including faithful quotient action.
- `Lean/Screen/GoldenSectorCharacters.lean`: the golden 3+3' character sector evaluated on the same committed port-action rows.

Mathlib already provides:

- `Matrix.ProjectiveSpecialLinearGroup`: by definition `SL / center(SL)`.
- the projective action of `SL` and `PSL` on projectivization;
- `Projectivization.SL_mulAction_ker`: the kernel of the projective `SL` action is exactly the center;
- `Projectivization.PSLAction.toPermHom`;
- `Matrix.ProjectiveSpecialLinearGroup.toPermHom_injective`;
- the canonical quotient map `QuotientGroup.mk'` and its surjectivity/kernel theorems.

## Preferred Lean surface

Suggested file:

```text
Lean/Screen/PSL2F5SixAxesBridge.lean
```

Suggested namespace:

```lean
namespace OPH.PSL2F5SixAxesBridge
```

Suggested type aliases:

```lean
abbrev F5 := ZMod 5
abbrev SL2F5 := Matrix.SpecialLinearGroup (Fin 2) F5
abbrev PSL2F5 := Matrix.ProjectiveSpecialLinearGroup (Fin 2) F5
abbrev P1F5 := ℙ F5 (Fin 2 -> F5)
```

## Theorem targets

### 1. Canonical quotient interface

Expose the canonical homomorphism

```lean
slToPsl : SL2F5 ->* PSL2F5
```

and prove/wrap:

```lean
slToPsl_surjective
slToPsl_ker_center
```

No new mathematics is required here: this should be a typed specialization of Mathlib's quotient construction.

### 2. Explicit projective-line relabelling

Construct an explicit equivalence

```lean
p1EquivSix : P1F5 ≃ Fin 6
```

using the OPH convention

```text
[z:1] <-> z,  z in F5
[1:0] <-> infinity = 5.
```

Do not use an arbitrary `Fintype.equivFin`; the coordinate convention must be load-bearing.

### 3. Transport the Mathlib PSL action to six points

Define

```lean
pslToSix : PSL2F5 ->* Equiv.Perm (Fin 6)
```

by conjugating `Projectivization.PSLAction.toPermHom` through `p1EquivSix`.

Prove it is injective by transporting Mathlib's faithful PSL action.

### 4. Match the OPH generators/action

Verify directly that suitable lifts/classes in `PSL2F5` act as the existing OPH generators

```text
A5SixAxes.t : z -> z+1
A5SixAxes.s : z -> -1/z.
```

Then prove the image of `pslToSix` is exactly the subgroup/action represented by `A5SixAxes.L60` / `rowF`.

Preferred headline statement:

```lean
psl_equiv_six_axis_group : PSL2F5 ≃* SixAxisGroup
```

where `SixAxisGroup` is either an existing subgroup type extracted from `L60` or a small new typed subgroup whose carrier is proved exactly the sixty committed rows.

Do not prove equality merely because both sides have cardinality 60. Equality of the concrete actions or generator closure must carry the proof.

### 5. Center of SL(2,F5)

Use `Matrix.SpecialLinearGroup.mem_center_iff` specialized to dimension two over `ZMod 5` to prove that the center consists exactly of scalar `+I` and `-I`.

Desired interface:

```lean
center_eq_plus_minus_one
center_card_two
```

Then record the exact short-cover statement

```text
1 -> {+I,-I} -> SL(2,F5) -> PSL(2,F5) -> 1
```

at the theorem/interface level supported by Mathlib's quotient.

If cardinality is useful, it is acceptable to prove `|SL(2,F5)| = 120` by finite enumeration/formula and derive `|PSL(2,F5)| = 60`; cardinality must remain a consistency check, not the identification proof.

## Scope boundary

This work MAY claim, if proved:

- a concrete isomorphism between Mathlib's `PSL(2,ZMod 5)` and the committed OPH six-axis action;
- that the canonical `SL(2,ZMod 5) -> PSL(2,ZMod 5)` quotient has kernel equal to the center;
- that in this specialization the center is `{+I,-I}`.

This work MUST NOT claim without additional proof:

- `PSL(2,5) ≅ A5` as an abstract classification theorem;
- `SL(2,5) ≅ 2I` / binary icosahedral identification;
- McKay `E8 <-> 2I`;
- that the golden OPH `3` and `3'` sectors are the McKay defining representation;
- a derivation/selection of `phi`;
- `27^phi` or any mass law;
- any physical rotation, family, coupling, or observable identification.

## Important representation-theoretic warning

The OPH golden sectors `3` and `3'` already descend to the projective quotient. Their pullbacks to `SL(2,5)` therefore have trivial action of the central `-I`.

The McKay `E8` correspondence for the binary icosahedral group is driven by a faithful two-dimensional representation on which the center is nontrivial. Therefore a future McKay bridge cannot simply identify the OPH golden triplets with the McKay defining doublet; it must prove how the triplets arise inside the representation ring generated by that doublet.

This distinction is part of the no-cheating boundary for GOLDEN-LIFT-2/3.

## Validation contract

Before any PR:

1. build the new Lean module;
2. build `OPHScreen`;
3. add `#print axioms` for headline theorems;
4. add the structural gate in the appropriate screen/golden test surface;
5. update `Lean/docs/PROOF_INDEX.md` only after the kernel build passes;
6. regenerate `claims/active_surface_inventory.json`;
7. update the golden/icosahedral claim boundary only to the exact level actually proved;
8. run the mandatory suite.

No `sorry`, `admit`, or new axioms.

## Suggested commit sequence

1. `Add abstract PSL2(F5) six-axis bridge`
2. `Document PSL2(F5) cover interface`

The present file is a planning seed only and is not scientific evidence.