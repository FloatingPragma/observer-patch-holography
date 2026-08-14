# OPH Claim Registry

The papers are the standalone source for every theorem, assumption, falsifier, and claim boundary.
This directory contains the machine-readable scientific registry that keeps those standalone
statements synchronized across the public stack.
It is not public reading-path material. Do not link the registry
from the top-level public README files; point readers to the papers, falsifiability map, and public
explainers instead.
README numeric summaries should distinguish source-only rows, empirical closures, compare-only
rows, and SI convention/display rows.

The registry is part of the working process:

- `claims/axiom_registry.yaml` records the normative three-axiom identities,
  interfaces, exclusions, and scientific realization boundaries.
- `claims/claim_registry.yaml` records top-level claim IDs, owner papers, claim tiers, imported
  mathematics, OPH-specific deltas, assumptions, evidence, falsifiers, and survival rules.
- `claims/novelty_matrix.csv` maps each claim against prior work.
- `claims/falsification_matrix.csv` records mathematical, physical-identification, and
  phenomenological failure modes.
- `claims/dependency_graph.json` records cross-claim dependencies.
- `claims/assumption_dictionary.md` gives stable names to recurring assumptions.
- `claims/frozen_prediction_register.json` records frozen and pending
  prediction contracts; `tools/build_fz_registry.py` validates it and renders
  `docs/FROZEN_PREDICTION_LADDER.md`.
- `claims/emergent_instrument_register.json` records scientific simulation and
  measurement instruments; `tools/build_instrument_register.py` validates it
  and renders `docs/INSTRUMENT_REGISTER_V3.md`.
- `claims/selection_ledger.json` and
  `claims/physical_identification_registry.json` record scientific selection
  classes and physical-identification boundaries; `tools/build_selection_ledger.py`
  validates them and renders `docs/SELECTION_LEDGER.md`.
- `claims/gravity_premise_ladder.json` records the gravity premise-elimination
  rungs; `tools/build_gravity_ladder.py` validates it and renders
  `docs/GRAVITY_PREMISE_LADDER.md`.
- `claims/public_surface_quantitative_claims.json` controls quantitative
  statements on public summary surfaces and is checked by
  `tools/check_public_surface_claims.py`.
- `claims/active_surface_inventory.json` is a generated reachability
  projection written by `tools/check_axiom_consistency.py --inventory`; it is
  not a hand-edited registry or project-status surface.

The validator is:

```bash
python3 tools/check_claim_registry.py
```

It checks that the registry release ID matches `paper/release_info.tex`, that every claim has an
owner paper and falsifier, that the novelty/falsification matrices and dependency graph contain
every canonical claim ID with no unknown IDs, that the one-row-per-claim novelty and DAG node
projections have no duplicates, and that paper sources do not depend on direct paths to this
registry. The falsification matrix may keep several independently scoped rows for one claim.

The GitHub workflow runs the validator on registry changes and on public claim-surface changes.
When a pull request changes paper TeX or the README claim narrative, it must also touch this
registry/check surface. That rule keeps the registry from becoming a stale snapshot.
