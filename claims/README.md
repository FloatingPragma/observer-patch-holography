# OPH Claim Registry

The papers are the standalone source for every theorem, assumption, falsifier, and claim boundary.
This directory is the development mirror used to keep those standalone statements synchronized
across the public stack.
It is an internal release-audit surface, not public reading-path material. Do not link the registry
from the top-level public README files; point readers to the papers, falsifiability map, and public
explainers instead.
README numeric summaries should distinguish source-only rows, empirical closures, compare-only
rows, and SI convention/display rows.

The registry is part of the working process:

- `claims/claim_registry.yaml` records top-level claim IDs, owner papers, claim tiers, imported
  mathematics, OPH-specific deltas, assumptions, evidence, falsifiers, and survival rules.
- `claims/novelty_matrix.csv` maps each claim against prior work.
- `claims/falsification_matrix.csv` records mathematical, physical-identification, and
  phenomenological failure modes.
- `claims/dependency_graph.json` records cross-claim dependencies.
- `claims/assumption_dictionary.md` gives stable names to recurring assumptions.
- `bibliography/doi_audit.csv`, `predictions/prospective_predictions.json`, and
  `negative_results/README.md` are adjacent audit surfaces.

The validator is:

```bash
python3 tools/check_claim_registry.py
```

It checks that the registry release ID matches `paper/release_info.tex`, that every claim has an
owner paper and falsifier, that the novelty/falsification matrices and dependency graph use known
claim IDs, and that paper sources do not depend on direct paths to this registry.

The GitHub workflow runs the validator on registry changes and on public claim-surface changes.
When a pull request changes paper TeX or the README claim narrative, it must also touch this
registry/check surface. That rule keeps the registry from becoming a one-time audit artifact.
