# Invariant-mining generation campaign

This package is the source-only control surface for the first bounded invariant
mining campaign. The campaign vocabulary was frozen before any candidate was
generated or evaluated, and the package contains no public comparison
payload, measured target, or comparison reader.

The committed state is:

`GENERATION_ENABLED__COMPARISON_SEALED`

That state has a narrow meaning. The source-feature registry is final for the
frozen campaign scope; the producer slots, observable grammar, nuisance
taxonomy, ranking rule, skip rule, exposed-data registry, campaign comparison
budget, and direct-\(N\) exclusion are byte-pinned and unchanged since
finalization; the deterministic candidate generator and evaluator are
enabled; and public-data access, target registration, and scoring stay
sealed. Candidates accumulate in `outputs/candidate_registry.json` under the
frozen grammar and ranking, and scoring waits for the complete registry and
the single issue-639 comparison.

Completeness is relative to the frozen scope: the seven registered producer
slots and the frozen grammar. The closure checks recorded in the registry are
machine-enforced by the builder and the independent verifier: every slot's
required features and grammar classes resolve, every feature pins live
artifacts by exact bytes, every non-direct-\(N\) slot has an exposure
surface, and the source-admissible completion nuisance covers every slot.

## Files

- `data/source_feature_registry.json` inventories the registered source
  surfaces with their proof or receipt ancestry and the finalized
  completeness block.
- `data/producer_slots.json` registers the search families without executing
  them.
- `data/observable_grammar.json` fixes the typed expression grammar and
  complexity budget.
- `data/nuisance_registry.json` fixes the nuisance taxonomy and the unresolved
  directions that every later producer must address.
- `data/ranking_policy.json` fixes numerical weights, eligibility, tie
  breaking, and skip behavior.
- `data/exposed_data_registry.json` fixes the exposure classes, the public
  surfaces presumed exposed, and the one-comparison scoring rule.
- `policy/pregeneration_policy.json` binds the documents into one fail-closed
  state contract with the campaign comparison budget.
- `outputs/source_projection.json` pins every registered source and control
  file by exact bytes.
- `outputs/pregeneration_freeze.json` binds the policy and projection into one
  content-addressed freeze.
- `tools/generate_candidates.py` is the deterministic generator; it reads
  only the frozen registries, certifies every relation by finite exhaustion
  or a pinned index, scores with the frozen weights, and rebuilds
  byte-exactly.
- `outputs/candidate_registry.json` holds the generated candidates, the
  per-slot generation states, and the sealed scoring boundary.

The direct-\(N\) lane is registered only to prove non-reentry. It is owned by
the dedicated capacity and cosmology chain, carries the locked
non-identifiability verdict of that chain, and is ineligible for fallback
ranking.

## Verification

From the repository root:

```bash
python code/invariant_mining/tools/build_source_projection.py --check
python code/invariant_mining/tools/build_pregeneration_freeze.py --check
python code/invariant_mining/tools/generate_candidates.py --check
python code/invariant_mining/tools/verify_pregeneration_freeze_independent.py
python -m pytest -q code/invariant_mining/tests
```

The independent verifier imports no builder or project helper. It recomputes
the source-file set, byte hashes, policy bindings, freeze identifier, required
registry rows, exposure coverage, comparison budget, disabled execution
state, and direct-\(N\) exclusion.

## Promotion boundary

Public-data access, target registration, and scoring stay sealed. The
grammar, ranking, nuisance taxonomy, producer list, exposure registry,
comparison budget, and direct-\(N\) rule cannot change while generation is
enabled. Scoring requires the complete registry, a separately reviewed
unsealing change, and issue #639 custody; the campaign terminates after its
first physical comparison.
