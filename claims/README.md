# OPH Claim Registry

The candidate executable model is [A-contrast-memory-spherical-v1](../model/contrast_memory_spherical_v1/README.md).
Its contrast-memory law is a candidate native repair law of the theory, and
the [scalar channel](../model/scalar_channel/README.md) on the same carrier is a second.
Claims about the model are owned by the
[executable-model paper](../paper/executable_observer_construction.tex), and
`python -m model` recomputes their evidence. Other claims keep their stated
source, preparation, reference, clock and observation domains; a common name
supplies no transfer.

The papers are the standalone source for their theorems, assumptions,
falsifiers and claim boundaries. Executed constructions may instead name a
self-contained protocol record as their owner. This directory contains the
machine-readable scientific registry that keeps those statements synchronized
across the public stack.
It is not public reading-path material. Do not link the registry
from the top-level public README files; point readers to the papers, falsifiability map, and public
explainers instead.
README numeric summaries should distinguish source-only rows, empirical closures, compare-only
rows, and SI convention/display rows.

`claims/claim_registry.yaml` is written as JSON, which YAML 1.2 also accepts.
The validator parses it as strict JSON and rejects duplicate keys; a YAML loader
reads the same data but lets a duplicate key through silently, so consumers
should parse it as JSON.

## The five ledgers

The repository keeps its scientific bookkeeping in five ledgers. Each ledger
has one machine source; the other files are views rendered from it or
sub-registers validated with it, and the checks fail when a view drifts from
its source.

| Ledger | Machine source | Views and sub-registers |
| --- | --- | --- |
| Claims | `claims/claim_registry.yaml`: the claim rows and a `claim_context` section with each claim's prior work and citations, owner section, falsifiers and dependencies; `claims/prior_work.json`, the verified bibliography of that prior work; `claims/assumption_dictionary.md`, the assumption names | `claims/novelty_matrix.csv`, `claims/falsification_matrix.csv` and `claims/dependency_graph.json` (rendered by `tools/build_claim_views.py`); the credited-prior-work section of each paper (`tools/build_prior_work.py`); the selection and physical-identification registers `claims/selection_ledger.json` and `claims/physical_identification_registry.json` (the selection ledger is rendered to `docs/registers/SELECTION_LEDGER.md` by `tools/build_selection_ledger.py`, which validates the physical-identification registry without rendering it); the public quantitative-display policy `claims/public_surface_quantitative_claims.json` |
| Premises | `tracking/premise_register.json`: the premise rows, with the model-supply crosswalk, the Standard Model coverage and the gravity rungs as its `model_supply`, `sm_coverage` and `gravity_rungs` sections; the theorem-library premises A1-A3 in `claims/axiom_registry.yaml`, stated in `docs/AXIOM_REFERENCE.md` | `docs/registers/PREMISE_REGISTER_V3.md`, whose closing sections show the Standard Model coverage and the gravity rungs |
| Predictions and tests | `claims/frozen_prediction_register.json`, the custody record of frozen predictions, and `tracking/observation_ledger.json`, whose `instruments` and `constants` sections hold the emergent-instrument register and the constants ancestry | `docs/FROZEN_PREDICTION_LADDER.md`, `docs/registers/OBSERVATION_LEDGER_V3.md`, `tracking/null_model_scorecard.md`; `docs/OPH_FALSIFICATION_PROGRAM.md` summarizes the falsification rules for readers |
| Evidence | the per-package manifests and verifiers under `evidence/`, `code/` and `model/`, the receipt-promotion and external-data registries under `code/audit/`, and `Lean/docs/PROOF_INDEX.md` | `evidence/README.md` |
| Papers | `paper/content_inventory.yml` and the paper list in `tools/build_tex_papers.py` | `paper/paper_release_manifest.json`, `docs/PAPER_SET.md` |

The executable model uses the existing claim identities and premises through
its [specification](../model/contrast_memory_spherical_v1/SPEC.md) and the
[premise register](../docs/registers/PREMISE_REGISTER_V3.md). Candidate
hypotheses and issue closure do not change a registered claim's scientific
status. Change this registry when the established claim or its evidence
changes.

Finite theorem validity, realization of its hypotheses by an executable
instance, observational equivalence to an effective field description, and
agreement with physical data are separate claims. An `emitted_artifact` row
can establish the second without establishing the last two. Its evidence
must identify the actual state sector, instruments, scheduler, observer access
and retained outputs. The [executable-model hypothesis](../docs/OPH_EXECUTABLE_HYPOTHESIS.md#how-the-model-is-tested)
states the operational comparison: distributions of admissible observer
records against what observers in our world record; counts of theorem references do not
discharge its physical hypotheses. Retrospective conditional fits belong in
the observation/postdiction ledgers and cannot become prospective prediction
successes through a status or wording change.

`claims/active_surface_inventory.json` is a generated reachability projection
written by `tools/check_axiom_consistency.py --inventory`; it is not a ledger.

The validator is:

```bash
python3 tools/check_claim_registry.py
```

It checks that the registry release ID matches `paper/release_info.tex`, that every claim has an
owner file and falsifier, that the novelty/falsification matrices and dependency graph contain
every canonical claim ID with no unknown IDs, that the one-row-per-claim novelty and DAG node
projections have no duplicates, that those three views match their rendering from the registry,
that every established result a claim uses is cited from a verified entry of
`claims/prior_work.json` or marked as a result of the project itself, that the papers'
credited-prior-work sections are current, and that paper sources do not depend on direct paths
to this registry. The falsification matrix may keep several independently scoped rows for one claim.

Three row-level contracts carry their own machine-checked declaration:

- Evidence medium. The artifact medium of every evidence path is read from its suffix, and an
  unlisted suffix fails closed. A theorem-asserting row (`conditional_implication`,
  `branch_entry`, `empirical_implementation`, `physical_establishment`) whose entire evidence
  list is prose declares `proof_medium: paper_prose`; the validator rejects that field on every
  other row, so a Lean-backed or run-backed row cannot acquire the prose label. An empty
  evidence list fails.
- Historical lanes and live owners. `gates` records the lane issues under which a claim's open
  work was scoped. The V3 and V4 lanes are closed; every gate is either an open issue or a
  historical lane that `tools/closed_lanes.py` maps to the open issue owning its remaining work,
  and generated pages render each historical number as labelled history beside that live owner.
  A gate never names a lane without a live owner. Each row that the historical topical-lane policy in
  `tools/check_claim_registry.py` names declares `required_topical_gate_owners`, a sorted subset
  of that row's own `gates`. The validator compares the declaration against the policy in both
  directions and rejects a policy key that names no registered claim, so a one-sided edit to
  either surface fails. This field is the one gate-named key admitted outside `gates`; every
  other gate-named key is rejected as a side channel.
- Owner medium. An owner file under `paper/`, `extra/`, `cosmology/`, `flagship/` or
  `docs/research/` is a paper source and declares no medium. An owner outside those roots declares
  `owner_medium: protocol_record`, which is admissible only for
  `claim_class: emitted_artifact`.

The GitHub workflow runs the validator on registry changes and on public claim-surface changes.
When a pull request changes paper TeX or the README claim narrative, it must also touch this
registry/check surface. That rule keeps the registry from becoming a stale snapshot.

## Executable model

The executable-model paper owns the claims about
[A-contrast-memory-spherical-v1](../model/contrast_memory_spherical_v1/README.md): the exact reduction of the
private seam memory to an exchange plus an agreement energy, on every state
of the memory; the [binding theorem](../model/contrast_memory_spherical_v1/PROOF.md) on every
federation of complete cells; particle stability on the observers' agreed
records, with the controls that lose it; the bound objects of one cell; the
port response; and the unrecovered collar. Their supplied parts are listed in
the [specification](../model/contrast_memory_spherical_v1/SPEC.md). The model reproduces
particle stability; its other observables are work in progress, and no
registered claim promotes them.

Every test of the model compares distributions of admissible observer records
with what observers in our world record. The state of the simulation is no
physical space, and the theory proves that observers who inhabit the
simulation can describe their experiences as relativistic 3+1D continuum
spacetimes. Finite theorem validity, realization of a theorem's hypotheses by
the model, and agreement with physical data stay separate claims.
