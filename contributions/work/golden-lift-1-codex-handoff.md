---
title: GOLDEN-LIFT-1 — Codex handoff
status: work-plan
claim_level: planning
physical_claim: false
---

# GOLDEN-LIFT-1 — Codex handoff

Resume work on branch `arithmon/golden-lift-1`.

Read `contributions/work/golden-lift-1-plan.md` first. That file is the scientific and API contract; do not broaden its scope.

## Work in two revisions

### Revision 1 — proof surface only

Implement the actual Lean bridge before promoting any scientific claim.

Primary target:

```text
SL(2,F5) → PSL(2,F5) ≃ SixAxisGroup ≤ Equiv.Perm (Fin 6)
```

Expected new Lean surface:

```text
Lean/Screen/PSL2F5SixAxesBridge.lean
```

Required ingredients/theorem interfaces, modulo better names if the implementation suggests them:

- `F5 := ZMod 5`, `SL2F5`, `PSL2F5`, `P1F5`;
- canonical quotient `slToPsl`, with surjectivity and kernel=center;
- explicit coordinate-preserving `p1EquivSix : P1F5 ≃ Fin 6`, with `[z:1] ↔ z` and `[1:0] ↔ 5` proved by evaluation lemmas;
- locally constructed projective `SL2F5` action on `P1F5` using only APIs available on the repository's pinned Mathlib revision;
- local theorem that the kernel of that `SL2F5` projective action is exactly the center;
- descended faithful `PSL2F5` action and `pslToSix`;
- a typed `SixAxisGroup : Subgroup (Equiv.Perm (Fin 6))` whose carrier is exactly the committed `A5SixAxes.L60` action;
- direct identification of the standard determinant-one lifts of `T=[[1,1],[0,1]]` and `S=[[0,-1],[1,0]]` with the existing OPH `A5SixAxes.t` and `A5SixAxes.s` under the explicit relabelling;
- image equality, not a cardinality-only argument;
- headline `psl_equiv_six_axis_group : PSL2F5 ≃* SixAxisGroup`;
- center theorem specialized to `SL(2,F5)`: center exactly `{+I,-I}`, plus `center_card_two` if practical;
- `#print axioms` on the headline quotient-kernel, faithful descended action, center, image-equality, and final isomorphism theorems.

Wire the new module into `OPHScreen`/`Lean/lakefile.lean` and extend the existing structural test owner (`tools/test_third_wave_surfaces.py`) only as needed to gate the new proof surface.

Before Revision 1 is considered complete, run at minimum:

```text
cd Lean
lake build PSL2F5SixAxesBridge
lake build OPHScreen
cd ..
python -m pytest -q tools/test_third_wave_surfaces.py
```

If repository conventions require equivalent commands, use them and record what was run.

Do not add paper/claim promotion merely to make the proof look more complete. If the theorem cannot be proved at the planned strength on the pinned Mathlib revision, stop and report the exact obstruction rather than weakening the statement silently.

Suggested Revision-1 commit message:

```text
Add abstract PSL2(F5) six-axis bridge
```

### Revision 2 — documentation / claim synchronization

Only start this revision after Revision 1 builds cleanly.

Then update the repository surfaces that own the result, to exactly the theorem level actually proved:

- `Lean/docs/PROOF_INDEX.md`;
- the standalone statement in the relevant owner paper, if the repository's claim ownership requires it;
- the relevant claim-registry row(s): synchronize `statement`, `oph_specific_delta`, `evidence`, `status`, and `falsifier` only where the new theorem genuinely changes them;
- novelty/falsification matrices if required by `claims/README.md` for a changed claim scope/delta/falsifier;
- dependency graph only if an ID/dependency actually changes;
- `tools/test_third_wave_surfaces.py` should gate the exact promoted wording and retain the no-cheating boundary;
- regenerate and commit `claims/active_surface_inventory.json` using the canonical inventory procedure.

Then run:

```text
python3 tools/check_axiom_consistency.py --inventory
python3 tools/check_axiom_consistency.py --check-inventory
python tools/run_mandatory_suite.py
```

Also rerun the direct Lean/OPHScreen builds if Revision 2 touches any Lean-owned surface.

Suggested Revision-2 commit message:

```text
Document PSL2(F5) cover interface
```

## Hard scope boundary for both revisions

Do not claim or implement in GOLDEN-LIFT-1 without a separate explicit decision:

- abstract `PSL(2,5) ≅ A5` classification;
- `SL(2,5) ≅ 2I` / binary icosahedral identification;
- McKay `E8 ↔ 2I`;
- transport of the OPH golden `3,3'` sector as a typed `PSL2F5` representation unless separately proved;
- identification of those triplets with the McKay defining doublet;
- derivation/selection of `phi`;
- `27^phi` or any charged-lepton mass law;
- any physical rotation, family, coupling, or observable claim.

The existing `A5PortSixAxesBridge` remains a pointwise indexed-row bridge, not yet a group-homomorphism arrow. Do not append it to the typed group chain unless a separate group packaging/isomorphism is proved.

## Stop condition

After both revisions are complete, report:

1. final branch SHA;
2. exact changed files by revision;
3. theorem names actually delivered;
4. direct Lean build results;
5. structural-test result;
6. mandatory-suite result;
7. any deviation from the planned theorem strength;
8. whether the branch is ready for a Draft PR to `FloatingPragma/observer-patch-holography:main`.

Do not open the upstream PR automatically unless explicitly asked.
