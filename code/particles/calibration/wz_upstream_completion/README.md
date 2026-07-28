# Physical W/Z upstream mathematical completion package

This directory preserves and hardens Pro's upstream theorem/specification
package between the structural OPH branch and the already constructed strict
W/Z pole map. The integrated v4.2 text is a corrected **draft sufficiency
specification**, not a production scientific receipt system.

## Main files

- `OPH_SM_EFT_WZ_COMPLETION_THEOREMS.md` — corrected non-entailment results, one explicit sufficient augmented branch, action/BRST/FJ/renormalization/pole theorems, proofs, and scientific boundaries.
- `INDEPENDENT_AUDIT.md` — archive verification, theorem audit, checker exploit, required production objects, and closure consequences.
- `candidate_axioms/OPH_PLUS_SM_EFT_FJ_V1.md` — explicit sufficient candidate source augmentation for human review; no minimality claim.
- `AGENT_WORK_ORDER_V4.md` — path-by-path coding-agent instructions.
- `docs/OPH_SOURCE_PARENT_CORRECTIONS.md` — H3/event-base, mixed-GNS, cone, causality and clock parent repairs.
- `schemas/` — nine Draft 2020-12 checklist schemas used for ten receipt instances (W and Z share one pole schema). They are not yet proof-bearing production schemas.
- `templates/` — deliberately incomplete, non-promotable packet templates.
- `proofs/symbolic_completion_proofs.py` — exact algebraic checks.
- `checkers/check_completion_bundle.py` — fail-closed fixed-template linter. It is unconditionally non-promoting because external artifact resolution is not implemented.
- `tests/test_completion.py` — schema/dependency tests plus a regression proving that forged self-attestation can satisfy the candidate predicate but can never promote.
- `data/nonlinear_gauge_grid_v1.json` — frozen 45-point gauge stress grid.
- `data/receipt_dependency_dag_v4.json` — acyclic receipt dependency graph.

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

`UPSTREAM_PACKAGE_MANIFEST.json` records Pro's pristine 35-file archive. The
post-audit integrated file identities are recorded separately in
`INTEGRATION_MANIFEST.json`; running validation does not rewrite either source
manifest.
