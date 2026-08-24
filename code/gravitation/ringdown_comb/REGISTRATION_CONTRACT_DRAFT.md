# Registration Contract Draft: Integer-k Kerr Ringdown Comb

Status: DRAFT. This document is not a registration. It becomes a
registered contract only when the owner anchors it (hash commit of this
text plus the completed likelihood artifacts, before any decision data
is opened). Until that owner action, nothing in this directory carries
verdict weight, and no numerical statement here is a frozen prediction.
No comparison access has occurred in the preparation of this draft: no
gravitational-wave event posterior, strain segment, catalog entry, or
published measurement was fetched, opened, or numerically evaluated.

## 1. What the erratum demands

`falsification/frozen_targets/fz01_2026-07-17/SCIENTIFIC_STATUS_ERRATUM_2026-07-29.md`
records that no ringdown value was frozen under FZ-01 and that the
former FZ-06 standing-freeze description is void. It states the
registration obligation this draft answers: a future integer-k
transition test requires a new source-derived transition, nuisance
model, event selection, likelihood, and decision rule registered before
comparison access. Each of the following sections supplies one of those
components in draft form.

## 2. Source-derived transition

The transition law is the frozen ratio law of
`falsification/frozen_targets/fz01_2026-07-17/frozen_target_integer_k_comb_2026-07-17.md`
(companion derivation in
`proof/epic_wins/ringdown_comb/INTEGER_K_COMB_STATEMENT.md`), taken as
the source of record and not modified here. Candidate spectral features
above the rotation line satisfy, under the comb hypothesis,

    (f_a - m*Omega_H/(2*pi)) / (f_b - m*Omega_H/(2*pi)) = ln(k_a)/ln(k_b)

for integers k >= 2, equivalently universal-coordinate positions
x_k = ln(k)/(8*pi) with x = (G*M/(c^3*g(chi))) * (omega - m*Omega_H)
and g(chi) = 2*sqrt(1-chi^2)/(1+sqrt(1-chi^2)). Secondary structure:
the (k-1)/k KMS weight hierarchy and the mass-independent fractional
linewidth 64*pi^2*p_0/(a*ln(k)) with declared a in [1, 10]. The
build-stage numeric instrument for this law is
`integer_k_comb_template.py` in this directory, with the independently
verified receipt `runtime/integer_k_comb_template_receipt.json` and the
exact invariance theorems in `Lean/Geometry/IntegerKCombInvariance.lean`
(mass, spin, and redshift invariance of the offset-subtracted ratio;
strict tooth monotonicity; rational bracketing of the reference ladder;
KMS hierarchy).

## 3. Event-selection rule

- Eligible events: ringdown observations published AFTER the owner
  anchors this contract. The anchoring timestamp partitions the
  literature; nothing published on or before it can enter the decision
  dataset.
- Analyst inputs per event: the published remnant mass and spin
  posteriors (M, chi) and any published secondary spectral feature
  candidates (frequency posteriors of post-merger features beyond the
  dominant mode), at the published-posterior level. No reprocessing of
  strain below the published-posterior level is part of this contract.
- Seen data, audit-only: GW150914 and GW250114 are named here as
  provenance references for the historical alpha = 4 audits recorded in
  the frozen target's supersession section. Both are already-seen data
  and are excluded from any future comparison under this contract; they
  can never contribute to a verdict.
- Catalog versions and data cuts: owner-fixed cells, to be enumerated
  at anchoring [OWNER-FIX: catalog list], [OWNER-FIX: data cuts].

## 4. Nuisance model

- Linewidth nuisance: a in [1, 10], declared, from the frozen
  statement; marginalized with a prior fixed at anchoring
  [OWNER-FIX: prior density on a].
- Page coefficient: p_0 = 2e-4, declared, statement-pinned.
- Remnant marginalization, stated abstractly: per event, the comb
  hypothesis is evaluated marginally over the published joint (M, chi)
  posterior, mapping each posterior draw to tooth positions through the
  template; frequency-error-only propagation is inadmissible per the
  frozen target. The concrete marginalization estimator is an
  owner-fixed cell [OWNER-FIX: estimator and convergence criterion].
- Redshift: offset-subtracted ratios and universal-coordinate positions
  are redshift-free; no redshift nuisance enters the ratio observable.
- Selection effects and trials accounting: owner-fixed cells
  [OWNER-FIX: selection model], [OWNER-FIX: trials correction].

## 5. Likelihood specification (structure only)

Per event, the contract compares two hypotheses at the
published-posterior level:

- H_comb: the secondary feature frequencies sit at the template teeth
  f_{k,m} = m*Omega_H/(2*pi) + c^3*g(chi)*ln(k)/(16*pi^2*G*M) for one
  integer assignment per feature, k in {2, ..., 12}, weighted across k
  by the normalized (k-1)/k hierarchy times GR greybody factors
  [OWNER-FIX: greybody normalization artifact], with linewidths from
  the nuisance model.
- H_no_comb: the published no-comb description of the same features
  (owner-fixed reference model [OWNER-FIX: no-comb reference]).

The per-event statistic is the Bayes factor BF = Z(H_comb)/Z(H_no_comb)
with both evidences computed over the same published posteriors and the
nuisance model of section 4. Stacking across events is performed in the
universal coordinate x. This section fixes structure only: no numeric
evaluation of any likelihood, evidence, or Bayes factor has been
performed, and none may be performed before anchoring. The frozen
target additionally requires a derived strain/asymptotic-readout
likelihood artifact before any verdict; that artifact is open and is a
precondition for anchoring, not replaced by this draft.

## 6. Proposed decision rule and kill bands (PROPOSAL, pending freeze)

Every threshold in this section is a proposal. The owner fixes each
cell at anchoring; until then no decision rule exists.

| Cell | Draft value | Status |
|---|---|---|
| Confidence level (kill and detection) | 99% | [OWNER-FIX] |
| Ringdown SNR gate per contributing event | >= 30 | [OWNER-FIX] |
| Resolved features required per event | >= 2 above rotation line | [OWNER-FIX] |
| Frequency-uncertainty gate | < 1/3 of local k = 2 tooth spacing | [OWNER-FIX] |
| Stacked evidence threshold | abs(ln BF) [OWNER-FIX] | [OWNER-FIX] |
| Integer ladder set for exclusion | k_a, k_b in {2, ..., 12} | [OWNER-FIX] |

Proposed verdict structure, mirroring the frozen target:

- KILL band: a resolved pair clearing the gates whose offset-subtracted
  ratio excludes every ratio ln(k_a)/ln(k_b) with k_a, k_b in the
  ladder set at the frozen confidence, under the frozen likelihood and
  covariance. This kills the integer-division continuation template; it
  does not by itself kill the derived area-spectrum statement.
- DETECTION-CANDIDATE band: coherent stacking at the predicted x_k
  positions with the strongest-pair ratio consistent with a single
  integer pair at the frozen confidence with trials correction; a
  single-event pass is reportable as a candidate only.
- INCONCLUSIVE: everything else, including all datasets in which only
  one feature clears the gate.

## 7. Boundary

This contract is a draft, registered only when the owner anchors it.
No comparison access has occurred; no event likelihood was evaluated;
no event was selected. The template instrument in this directory is
target-blind and its receipt contains no event-derived number. The
derived strain likelihood, the normalized pre-data prior over k, the
event list, and the trials accounting are open owner actions listed in
the frozen target. Nothing in this directory may be cited as a
registered, frozen, or scored prediction.
