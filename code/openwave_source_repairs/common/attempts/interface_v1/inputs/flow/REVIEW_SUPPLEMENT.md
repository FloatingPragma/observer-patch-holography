# Source-word error accounting

Reviewer: Codex energy-repair worker, same group, separate from the stage2
producer. This review covers the two proof points identified by root. It does
not confer other-group G6 or a physical clock interpretation. The stage2
source and its preregistration are unchanged.

## Retained source-padding outcomes

Let U_j be the ideal rank-one word on the valid240-dimensional source chart,
extended by identity on its16 padding coordinates. Let A_j be the retained
success product of its stored M6 primitives. Set
delta_j >= ||A_j-U_j|| and g_j >= product_{k<j} ||A_k||. Since every U_j
preserves the valid chart, the j-th first-padding-abort map has norm at most
delta_j g_j. Its effect is positive, so the total padding-abort effect obeys

    ||sum_j B_j* B_j|| <= sum_j (delta_j g_j)^2
                         <= (sum_j delta_j g_j)^2.

The scalar telescoping estimate used by the implementation bounds the final
sum by e=(N b+q)/(1-N b): all raw-pulse gain factors are included, and the
target preparation/phase errors q compare unitary products. Thus e² bounds
the sum of all padding-abort effects. This is one bound on the sum, not a
separate allocation of the entire e² to each word. The same first-failure
factorization retains processor complements before each source export.

## Full M6 targets and proper-code dimensions

The retained compiler's _source routine checks exact Q(sqrt5) linear
combinations of all twelve source generators before any pulse is used. Its
X,Y,Z,H and I targets act on the first two coordinates; J is the phase on
the first three coordinates. These statements concern all six coordinates,
including their complements. The numerical native_pulse evaluates the full
exponential of the reconstructed generator, and primitive_target encloses
the corresponding exact full6x6 target entry by entry.

For a one-qubit request d=2. For J, d=4 and the logical target is
diag(exp(i theta),exp(i theta),exp(i theta),1). Every target leaves its proper
code invariant. Its stored complement is full[d:,:d], which apply_word
retains before subsequent pulses. The exact source identity, the chosen
code embedding and the supplied load/export service are distinct premises.
The comparison does not assume that a three-dimensional block is a qubit,
or omit its third coordinate when J uses the four-dimensional code.

The existing prepare/inverse/phase words are finite compositions of these
targets. Their binary64 state propagation is diagnostic; the rational
full-matrix enclosure and complete-map bound carry the arbitrary-reference
claim. Reviewers should not read the small observed padding probability as
a replacement for the operator estimate.

## Replay custody

The original attempt stage2_001 is retained with its partial files and log.
The identical frozen program was restarted as stage2_002 because the former
worker session was interrupted. Each attempt has a separate output path.
No existing receipt or scientific gate is overwritten. Completion and actual
runtime are taken from stage2_002/RECEIPT.json when that run terminates.
