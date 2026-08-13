# OL-A1 follow-up: carrier count by absolute support size

Status: design only. This document does not freeze a campaign, draw a seed,
authorize a run, or change the `FAILED` verdict of INS-01. Any execution needs
a new immutable preregistration and a separate freeze event under issues #737
and #738.

## What the first campaign did and did not show

The conformant INS-01 campaign returned `FAILED` overall, but its component
results should not be collapsed into one undifferentiated statement:

- The 65,536-carrier A2 arm reproduced the threshold `(1,3)` form in five of
  five replicates and matched its preregistered robust reference in five of
  five replicates.
- The 16,384-carrier A1 arm reproduced the threshold form in two of five
  replicates and matched its robust reference in one of five. This is direct
  evidence of small-rung or seed sensitivity, not evidence that the small rung
  is known to be the cause.
- The A2/A1 cone-margin ratio missed the declared band in all five paired
  replicates.
- The ancestry-permutation null did not discriminate reliably. The current
  estimator therefore did not attribute the form to ancestry structure.
- The low-support control greatly reduced cross edges, but the frozen firing
  rule did not fire. That control did not isolate a clean support-size
  mechanism.

It is therefore plausible that the smallest rung was too small, but INS-01
cannot establish that explanation. Carrier count changed together with
observer count, and only one low-support control was included. A factorial
follow-up is needed to separate finite-size, support, and ancestry effects.

## Candidate factorial ladder

The next freeze should vary carrier count and absolute observer support size
independently, while holding observer count and every other producer and
analysis choice fixed within each paired seed block. "Support size" here means
the integer number of carriers read by one observer. The derived support
fraction is support size divided by carrier count; because that fraction
changes along the carrier-count axis, this document does not call the second
axis "support density." A reasonable design envelope is:

| Axis | Candidate levels | Purpose |
| --- | --- | --- |
| carrier count | 16,384; 65,536; 131,072; 262,144 | retains the seed-sensitive rung and adds three larger rungs, so the small-rung hypothesis is directly testable |
| absolute support size | 48; 96; 192 | low, reference, and high absolute support at every carrier count; support fraction is reported as a derived quantity |
| observer count | one value fixed at the future preregistration, provisionally 256 | prevents observer-count scaling from being confused with carrier scaling |
| fresh replicate blocks | determined prospectively before freeze | chosen from declared minimum effects, confidence-interval widths, multiplicity, and attrition rather than assuming that eight or any other convenient count is sufficient |

These numbers are design candidates, not frozen values. Before a freeze, a
runtime-only pilot may measure cost but must not expose signature, margin, or
control outcomes. A target-blind prospective calculation must set the number
of replicate blocks using predeclared minimum scientifically relevant main
effects, interactions, ancestry contrasts, equivalence widths, familywise
error or false-discovery control, desired power, and an attrition allowance.
It must publish the calculation and its inputs; it may not use INS-02 outcome
estimates. The final preregistration must either keep the complete factorial
grid or declare a target-blind resource reduction before any outcome is
inspected. There is no authorization to execute this ladder now.

Each fresh master seed forms one block spanning every retained
carrier-count-by-support-size cell and every matched control. The freeze must
pin a deterministic, collision-free derivation from the master seed to each
cell seed, randomize execution order within each block, and analyze the paired
contrasts at block level. Seed redraw, cell substitution, optional stopping,
and analysis that treats paired cells as independent are forbidden.

## Stronger ancestry controls

A mere label permutation may preserve the structures used by the estimator.
The next campaign should preregister all of the following:

1. An ancestry-destroying, degree-preserving rewiring that preserves carrier
   count, degree sequence, support counts, and one-point record marginals while
   breaking parent-child and shared-lineage correlations.
2. A depth-stratified lineage shuffle that preserves the population at every
   declared depth while breaking cross-depth family identity.
3. A sham relabeling that preserves the complete ancestry graph. The analysis
   must remain invariant under this negative control.
4. A sensitivity control with a synthetic ancestry-dependent perturbation of
   known sign and predeclared magnitude. It checks whether the fixed estimator
   can recover an ancestry effect when one is present; it has no promotion
   authority.

Every control must run through the same observable code path and the same
fixed estimator as its matched main cell. A separate structural verifier may
check the invariants of a transform, but it may not replace the estimator in
an inferential gate. Before freeze, the preregistration must define one
continuous ancestry-sensitive estimator output and quantitative thresholds for
all of the following gates:

- **destruction:** each destructive transform reduces the declared lineage
  association metric by at least a fixed amount while its preserved graph and
  marginal quantities remain within fixed tolerances;
- **estimator contrast:** the paired actual-minus-destroyed contrast has the
  declared sign and its multiplicity-adjusted interval clears a minimum effect;
- **sham equivalence:** the paired actual-minus-sham interval lies wholly
  inside a predeclared equivalence band;
- **synthetic sensitivity:** the same estimator recovers the injected effect's
  sign and clears a predeclared minimum fraction of its known magnitude.

The thresholds, interval construction, multiplicity family, estimator code,
and missing-data rule must all be frozen. Failure of destruction, sham, or
synthetic-sensitivity gates prevents an ancestry-specific positive conclusion,
even if a threshold signature appears in the main cells.

## Separate questions and outcomes

The follow-up must report the following separately before applying any overall
verdict:

- threshold-signature stability at each predeclared threshold;
- the continuous normalized eigenvalue spectrum and degeneracy margins;
- carrier-count main effect, absolute-support-size main effect, and their
  interaction, with the derived support fraction reported but not relabeled as
  an independently varied axis;
- ancestry-destroying effects versus sham effects;
- cone-margin scaling and the old ratio-band result as a distinct endpoint;
- cell-level pass, fail, or inconclusive status with multiplicity and missing-
  data rules fixed in advance.

The future frozen decision rule should distinguish at least these scientific
outcomes: size-supported stabilization, support-sensitive stabilization,
ancestry-specific stabilization, no discrimination, and unresolved. A large
rung passing cannot overwrite a failed ancestry control or failed margin
endpoint. Conversely, a failed small rung cannot erase a reproducible large-
rung component result. Only a separately frozen overall rule may return
`REPLICATED`, `FAILED`, or `INCONCLUSIVE` for the emergent OL-A1 claim. Even
`REPLICATED` is not an automatic ledger promotion: it can support or qualify
OL-A1 only under the anchored architecture version used by the run, after an
origin-anchored independent audit qualifies OL-A1, and after the register
explicitly selects the successor as controlling. A controlling `FAILED`
verdict blocks or demotes OL-A1 to `owed`.

## Custody and stopping rule

The final campaign must retain raw feature matrices or sufficient source
captures for independent observable recomputation, not only derived receipt
fields. It must pin the exact simulator commit, manifest inventory, code,
configuration, seeds, environment, captures, receipts, and independent
verifier before comparison. No adaptive seed redraw, optional stopping,
post-hoc threshold, or selective cell omission is allowed.

Until that freeze and execution happen, INS-02 remains `SPECIFIED`, OL-A1
remains `owed`, and INS-01 remains the controlling completed verdict. INS-02
does not supersede INS-01 merely by existing as a design.
