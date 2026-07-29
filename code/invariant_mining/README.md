# Invariant-mining pre-generation lock

This package is the source-only control surface for the first bounded invariant
mining campaign. It freezes the campaign vocabulary before any candidate is
generated or evaluated. It contains no public comparison payload, measured
target, candidate table, scoring code, or comparison reader.

The committed state is:

`PREGENERATION_POLICY_LOCKED__GENERATOR_DISABLED_PENDING_REGISTRY_FINALIZATION`

That state has a narrow meaning. The seed source-feature registry, producer
slots, observable grammar, nuisance taxonomy, ranking rule, skip rule, and
direct-\(N\) exclusion are byte-pinned. Candidate production and evaluation
remain disabled until a separately reviewed registry-finalization change
replaces this state.

## Files

- `data/source_feature_registry.json` inventories the bounded seed source
  surfaces and their proof or receipt ancestry.
- `data/producer_slots.json` registers the first search families without
  executing them.
- `data/observable_grammar.json` fixes the typed expression grammar and
  complexity budget.
- `data/nuisance_registry.json` fixes the nuisance taxonomy and the unresolved
  directions that every later producer must address.
- `data/ranking_policy.json` fixes numerical weights, eligibility, tie
  breaking, and skip behavior.
- `policy/pregeneration_policy.json` binds the documents into one fail-closed
  state contract.
- `outputs/source_projection.json` pins every registered source and control
  file by exact bytes.
- `outputs/pregeneration_freeze.json` binds the policy and projection into one
  content-addressed freeze.

The direct-\(N\) lane is registered only to prove non-reentry. It is owned by
the dedicated capacity and cosmology chain and is ineligible for fallback
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
registry rows, disabled execution state, and direct-\(N\) exclusion.

## Promotion boundary

Registry finalization is a separate reviewed change. It must replace the
status, prove that the declared source inventory is complete for the frozen
campaign scope, and preserve a pre-generation commit boundary. Candidate
generation cannot be enabled in the same unreviewed mutation that changes the
grammar, ranking, nuisance taxonomy, producer list, or direct-\(N\) rule.
