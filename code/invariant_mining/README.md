# Invariant-mining finalized registry lock

This package is the source-only control surface for the first bounded invariant
mining campaign. It freezes the campaign vocabulary before any candidate is
generated or evaluated. It contains no public comparison payload, measured
target, candidate table, scoring code, or comparison reader.

The committed state is:

`REGISTRY_FINALIZED__GENERATOR_DISABLED_PENDING_ENABLEMENT_REVIEW`

That state has a narrow meaning. The source-feature registry is final for the
frozen campaign scope, and the producer slots, observable grammar, nuisance
taxonomy, ranking rule, skip rule, exposed-data registry, campaign comparison
budget, and direct-\(N\) exclusion are byte-pinned. Candidate production and
evaluation remain disabled until a separately reviewed enablement change
flips the generator flag without touching any registry content.

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

The direct-\(N\) lane is registered only to prove non-reentry. It is owned by
the dedicated capacity and cosmology chain, carries the locked
non-identifiability verdict of that chain, and is ineligible for fallback
ranking.

## Verification

From the repository root:

```bash
python code/invariant_mining/tools/build_source_projection.py --check
python code/invariant_mining/tools/build_pregeneration_freeze.py --check
python code/invariant_mining/tools/verify_pregeneration_freeze_independent.py
python -m pytest -q code/invariant_mining/tests
```

The independent verifier imports no builder or project helper. It recomputes
the source-file set, byte hashes, policy bindings, freeze identifier, required
registry rows, exposure coverage, comparison budget, disabled execution
state, and direct-\(N\) exclusion.

## Promotion boundary

Generator and evaluator enablement is a separate reviewed change. It must
flip only the enablement flags and add producer code, and it cannot change
the grammar, ranking, nuisance taxonomy, producer list, exposure registry,
comparison budget, or direct-\(N\) rule in the same mutation. The campaign
terminates after its first physical comparison, which issue #639 owns.
