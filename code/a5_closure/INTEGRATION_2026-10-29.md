# Reviewed source custody after integration

The operational attachment probes were reviewed together with the current
native-model changes. Their original work baselines remain in each receipt.
The original four receipts, producers and independent verifiers are preserved
at commit `78c4eb137cf3bf20e7c1d907561f2f2abe22cbf4`; for example
`git show 78c4eb137cf3bf20e7c1d907561f2f2abe22cbf4:code/a5_closure/receipts/proper_code_attachment_probe.json`
recovers the exact earlier receipt. Current source and prerequisite pins are
advanced separately from those historical work baselines.

The complete imported-source comparison found only these changes relative to
the probes' scientific source pins:

- `model/native_clock/configuration.py`: `success_probabilities` now returns
  fresh `Fraction` objects with the same exact numerator and denominator.
  This prevents callers changing the cache through mutable protected slots.
  Endpoint maps, probability values, source controls and preparation laws
  are unchanged. Four binary/full78 and scan/retry cases were checked for
  exact value equality and independence from mutated returned objects.
- `model/common_source/QUANTUM_NATIVE.md`, used only by the coherent-access
  probe: a paragraph clarifies that the old expanded pulse totals are
  historical, identifies the existing saved 9,738 pulses per flip and links
  the already current receipt. No source code, pulse schedule or number was
  changed by that paragraph.

All other pinned scientific inputs were byte identical. In particular the
native branch, compiler, tomography, pointer and coherent routing code were
unchanged. Current scientific inputs are pinned to
`3f650db7a625f0df66495bdaea4332861507f5d4`. Dependent prerequisite pins identify
the subsequent immutable commits containing their reviewed current receipts.
The original Galois drift diagnostic retains its original source label and
diff hash; integration does not erase that older unresolved provenance check.

On 2026-11-01 the source-realization receipt was rebuilt twice: for the
dependency and action upgrades, and for the move of CI to Ubuntu only, which
changed the workflow file it binds. Its topology is unchanged and its
floating diagnostics move by at most 3.3e-16. The source pin and the
port-support base advanced to `83113a9cc33d82a56c84a8634445c8edd1e599a2`, where every other pinned input is byte
identical. The C12 prerequisite pin advanced to `e5a59e228540703224f9d1dcc2c32e6d25de87d4` and the coherent
native-access prerequisite pin to `564a0b6eb13332e02f27c07e4870a660e52efea6`. Each refreshed receipt differs from
its predecessor only in custody fields.

Each refreshed receipt is compared with its original after removing only
source/prerequisite/implementation custody fields. Every mathematical result,
matrix-unit map, resource count, output distribution, obstruction scope and
verdict must remain exactly equal. Producers and independent verifiers retain
their strict byte comparisons against the new pinned commits. An arbitrary
future source modification is still rejected; this is no hash bypass.

The current source inventory also refreshes its complete directory snapshot.
The only added path is this note; the sixteen changed paths are the four
probes, their independent verifiers, receipts and explanatory pages reviewed
here. Every candidate, classification, qualification, scope, source pin and
implementation pin in that inventory remains exactly unchanged. Only its
directory snapshot and the resulting inventory digest change.

The retained conclusions remain bounded: proper-code export alone does not
identify spatial loads with quantum states; coherent native access does not
provide a source-derived spatial decoder; and the clean-register decoder
obstruction applies to the fixed finite template grammar. No spectral-gap
operator map, physical particle or calibrated mass is added by integration.

## A word-length-independent corollary inside the fixed grammar

The clean-register decoder proof uses field closure, not a search through
the 31-block resource envelope. Every exact ideal template in the declared
alphabet has matrix entries in `F = Q(i,sqrt(2))`. Any finite product of
these templates, and its adjoint, still has entries in `F`. If such a unitary
decoded the GV code into the declared clean receiver, its code projector
would be `U* P_receiver U`, with all entries in `F`, since the clean receiver
projector has rational entries. The independently reconstructed GV projector
has entry `1/20 + sqrt(5)/100`, outside `F` because
`Q(sqrt(5)) intersect Q(i,sqrt(2)) = Q`.

Thus the same exact clean-register obstruction holds for every finite word
length in this fixed template alphabet. The receipt's 31-block admitted
resource envelope and its replay are unchanged. This analytic corollary does
not cover arbitrary angles, approximate decoders, nontrivial residual states,
postselection, a changed code or another alphabet. It does not exclude the
broader executable program class. No new decoder search or source execution
is claimed by the corollary.
