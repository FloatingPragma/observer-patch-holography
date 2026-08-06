# H0 source-only preflight exit evaluation (v1)

Mirror of `source_only_preflight_exit_v1.json`. Issue #708. All five preflight-completion
gates are frozen as content (decision rule in the contract, the other four in
`source_only_preflight_gate_freeze_v1.json`).

## Exit classification

**NOT_EVALUABLE_MINIMAL_PREMISE_LIST.**

Against the declared menu:

- `SOURCE_ONLY_CANDIDATE_FROZEN` is unearned: zero of thirty registered rows attain all
  seven source-only row gates; the strongest nine rows attain four of seven.
- `BOUNDED_NO_GO` is unearned: the audit fails no strong instance structurally on all
  seven gates, the grammar admits successor-version extension by its own closure rule,
  and exhaustion over unexecuted source laws is a declared open proof obligation
  (`future_source_laws_exhausted = false`, `candidate_inventory_exhaustiveness =
  NOT_ESTABLISHED`). The fail-closed rule bars the no-go while the qualification space
  is unexhausted.

## Production space at the frozen grammar version

Finite at depth budget 5: nine listed productions (eight closed-form, one interval
certified), three barred productions, 197 expression instances (4 CAP-K, 4 CAP-P power,
2 CAP-P linear, 180 CAP-L sublattice, 1 coupled, 1 common load, 2 reserve, 1 RC-LOAD,
2 P interval modes), audited as 30 registry rows (22 grammar rows after CAP-L
aggregation, 3 barred or control entries, 5 hierarchy packets). Instance coverage is
complete for this version; the space bounds the frozen version only, so it does not
bound the H0 question.

## Gate failures per family

Three gates are open on every row: `target_independent_candidate_selection`,
`same_quantity_constructed`, `source_return_map_complete`. Their failure mode is
evidentiary absence under the fail-closed rule.

| Family | Rows | Attained | Classification |
|---|---|---|---|
| P_INTERVAL_MAP | 2 | 4/7 | contingent |
| CAP-K linear | 4 | 2/7 | structural family-local (no positive fixed point) |
| CAP-P power | 4 | 4/7 | contingent (seed self-selection barred in Lean) |
| CAP-P linear | 2 | 2/7 | structural family-local (no positive fixed point) |
| CAP-L | 5 (180) | 2/7 | contingent (uniqueness and stability uncertified) |
| CAP-B | 1 | 1/7 | structural as registered (barred target bridge) |
| coupled | 1 | 4/7 | ineligible as registered, contingent under extension |
| direct control | 1 | 1/7 | structural as registered (fixed-cutoff control) |
| common load | 1 | 1/7 | contingent |
| reserve branches | 2 | 1/7 | contingent (branch selection unsourced) |
| RC-LOAD | 1 | 4/7 | contingent |
| archive 190 | 1 | 2/7 | structural as registered (enumeration is metadata) |
| hierarchy packets | 5 | 1-4/7 | mixed (four target exposed, one open witness) |

Earned family-local no-gos: linear families with s != 1 admit no positive fixed point;
seed self-selection is barred (`seededReadback_fixedPt_iff`); codomain sharing is
insufficient as a bridge (`sameCodomain_noncommuting_witness`).

## Minimal not-evaluable premise list

1. **P1** A frozen target-independent selection rule (selection rule status is
   DRAFT_UNFROZEN; every candidate freeze_evidence is null); the selector must be
   external to the candidate fixed-point equation.
2. **P2** A constructed same-quantity bridge (same_quantity_constructed is null for
   every candidate); the bridge must exhibit the commuting square.
3. **P3** A complete source return map (P consumer-policy flags false; no N family
   carries a return leg; archive enumeration is barred as a return map).
4. **P4** N-side certification or branch selection: certified uniqueness and stability
   for a CAP-L family, or a source rule selecting one reserve branch.
5. **P5** The exhaustion premise over unexecuted source laws, required before any
   bounded no-go could be earned at a successor version.

## Replay

`python3 code/closure_preflight/verify_source_only_closure_preflight_independent.py`
prints `SOURCE_ONLY_CLOSURE_PREFLIGHT_INDEPENDENT_VALID` and exits 0. This document is
additive: no audit surface pins it, and issue closure stays unauthorized by it.
