# Physical W/Z upstream mathematical completion package

This directory preserves the scientific theorem/specification package between
the structural OPH branch and the already constructed strict W/Z pole map. The
integrated v4.2 text is a corrected **draft sufficiency specification**, not a
production scientific receipt system.

## Main files

- `OPH_SM_EFT_WZ_COMPLETION_THEOREMS.md`: corrected non-entailment results, one explicit sufficient augmented branch, action/BRST/FJ/renormalization/pole theorems, proofs, and scientific boundaries.
- `candidate_axioms/OPH_PLUS_SM_EFT_FJ_V1.md`: explicit sufficient candidate source augmentation for human review; no minimality claim.
- `docs/OPH_SOURCE_PARENT_CORRECTIONS.md`: H3/event-base, mixed-GNS, cone, causality and clock parent repairs.
- `schemas/`: the baseline Draft 2020-12 checklist schemas plus the
  fail-closed diagnostic and certified-contour receipt schemas. The baseline
  templates remain specification-only.
- `templates/`: deliberately incomplete, non-promotable packet templates.
- `proofs/symbolic_completion_proofs.py`: exact algebraic checks.
- `checkers/check_completion_bundle.py`: fail-closed fixed-template linter. It is unconditionally non-promoting because external artifact resolution is not implemented.
- `producers/wz_pole_receipts.py`: a target-free sampled W/Z boundary diagnostic. Its fixed false certification flags prevent the sampled winding from being used as a root, Laurent, current-pole, or OPH-native receipt.
- `checkers/check_wz_pole_diagnostic.py`: independent schema, source-pin, exact-correction, contour, self-digest, and fail-closed flag replay for that diagnostic.
- `producers/certified_wz_contours.py`: directed complex-interval
  principal-sheet zero exclusion on declared upper-half-plane boxes. It
  includes the finite `d = 4 - 2 epsilon` prefactor terms and replays one
  exact rational boundary partition at 128/192/256 bits.
- `checkers/check_certified_wz_contours.py`: fail-closed artifact and
  enclosure-evidence checker for the principal-sheet receipt. It checks
  immutable inputs, exact fixture/correction data, partition identity,
  interval nesting, residual gates, claim scope, and the self-digest. It is
  an arithmetic/shape validator: it does not re-evaluate the loop functions,
  authenticate the producer, or supply an independent clean-room numerical
  third verifier.
- `producers/certified_second_sheet_poles.py`: scalar W/Z pole certificates
  on an explicitly declared channel-by-channel continuation chart. The
  receipt uses a mass-exchange-symmetric one-mass `B0` chart and records
  which cuts are crossed, winding one, interval-Newton enclosures,
  derivative-denominator and scalar-residue balls, and the same precision
  evidence ladder.
- `checkers/check_certified_second_sheet_poles.py`: fail-closed artifact and
  evidence checker for that declared scalar chart. Matrix-rank Laurent data,
  a BRST-invariant current amplitude, and independent numerical
  re-evaluation remain outside this receipt.
- `tests/test_completion.py`: schema/dependency tests plus a regression proving that forged self-attestation can satisfy the candidate predicate but can never promote.
- `data/nonlinear_gauge_grid_v1.json`: frozen 45-point gauge stress grid.
- `data/receipt_dependency_dag_v4.json`: acyclic receipt dependency graph.

## Run

```bash
python3 run_all.py
```

Expected final status:

```text
DRAFT_SUFFICIENCY_STACK_DEFINED__SIMULATION_RECEIPTS_OPEN__NO_OPH_NATIVE_POLE_PROMOTION
```

The templates are not physics inputs. They define a preliminary data shape and
intentionally keep evidence gates false. Even a forged all-true template remains
non-promoting. A future verifier must resolve immutable artifacts, recompute
digests and equations, and derive every evidence result itself.

The sampled W/Z contour output remains non-promoting: its mpmath radii are
roundoff heuristics without directed rounding, and its finite regulator does not
certify an analytic continuation sheet. It is superseded for contour claims by
the two directed-interval receipts above.

The principal receipt proves only zero exclusion on its declared first-sheet
boxes. Its coefficient-polynomial denominator checks are not Laurent
denominator balls. The pole receipt proves a scalar subclaim on one declared
multi-channel continuation, not the full neutral matrix-rank or
physical-current evidence boundary. Both receipts keep the corresponding
full-stack completion flags false, and all OPH-native and unit claims remain
false. The scalar pole receipt fixes its engine convention as
`G(s)=s-m_tree^2-Pi_engine(s)` and does not certify a sign bridge to the
separately written theorem convention. A receipt self-hash detects accidental
corruption but does not authenticate numerical evidence; the bundled
validation therefore performs byte-exact producer regeneration before running
either structural checker.

`UPSTREAM_PACKAGE_MANIFEST.json` is an immutable record of the pristine
35-file source archive; its historical file pins are retained even when a
source-archive file is not part of the current package surface. Current
integrated file identities are recorded separately in
`INTEGRATION_MANIFEST.json`; running validation does not rewrite either
manifest.
