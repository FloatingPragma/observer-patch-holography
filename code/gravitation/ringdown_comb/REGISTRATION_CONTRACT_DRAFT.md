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

## 2. Imported continuation transition; source derivation open

The transition law is imported from the attested but explicitly unanchored
FZ-01 draft
`falsification/frozen_targets/fz01_2026-07-17/frozen_target_integer_k_comb_2026-07-17.md`
(companion derivation in
`proof/epic_wins/ringdown_comb/INTEGER_K_COMB_STATEMENT.md`). It is a
Bekenstein--Mukhanov-style continuation premise, not a transition derived from
the OPH source. The FZ-01 erratum's source-derivation obligation therefore
remains open.

For an emission, the rule assumes positive integer record dimensions with
`d_before = k*d_after`, so `k >= 2` must divide `d_before`. The signed
black-hole entropy change is

    ln(d_after) - ln(d_before) = -ln(k),

while the positive entropy loss entering the emitted-energy formula is
`ln(k)`. The tooth formula uses the leading small-transition Kerr first law;
a finite-step model must account for the changing background. Conditional on
that imported selection, candidate spectral features
above the rotation line satisfy

    (f_a - m*Omega_H/(2*pi)) / (f_b - m*Omega_H/(2*pi)) = ln(k_a)/ln(k_b)

for eligible integer divisors `k >= 2`, equivalently universal-coordinate positions
x_k = ln(k)/(8*pi) with x = (G*M/(c^3*g(chi))) * (omega - m*Omega_H)
and g(chi) = 2*sqrt(1-chi^2)/(1+sqrt(1-chi^2)). Secondary structure:
the within-line (k-1)/k KMS net-response factor and the mass-independent but
spin-dependent linewidth-to-spacing ratio
64*pi^2*p_0/(a*g(chi)^2*ln(k)) with declared a in
[1, 10]. The earlier draft omitted the g(chi)^(-2) factor; its k=2
1.8--18 percent band is only the Schwarzschild limit. The constant-p_0 model
is not controlled near extremality, so a spin-dependent Page/greybody power
model is required before such remnants are eligible. A fixed record dimension
admits only its actual divisors; the all-integer display is a union of
candidate continuations, not a source-produced transition set. The
build-stage numeric instrument for this law is
`integer_k_comb_template.py` in this directory, with the independently
verified receipt `runtime/integer_k_comb_template_receipt.json` and the
exact invariance theorems in `Lean/Geometry/IntegerKCombInvariance.lean`
(abstract scale-and-offset invariance of the offset-subtracted ratio;
strict tooth monotonicity; rational bracketing of the reference ladder; and
algebraic properties of the declared KMS factor). The KMS factor is not a
transition probability or prior across different `k`.

## 3. Event-selection rule

- Eligible events: ringdown observations published AFTER the owner
  anchors this contract. The anchoring timestamp partitions the
  literature; nothing published on or before it can enter the decision
  dataset.
- Required comparison input per event: detector strain/readout data together
  with the prospectively specified likelihood, or a released sufficient
  likelihood product that supports both hypotheses and preserves the needed
  normalization. Published posterior samples may be nuisance-proposal samples
  only when their sampling prior, likelihood values, evidence normalization,
  support, and reweighting diagnostics are available. Posterior samples alone
  are not a likelihood or model evidence.
- Frame inputs per event: the joint detector-frame remnant mass and spin
  distribution and detector-frame secondary-feature frequencies. If a release
  supplies source-frame mass, the contract must transform it with the same
  redshift inference and covariance used by the likelihood:
  `M_det=(1+z)M_source`.
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
- Remnant marginalization, stated abstractly: per event, the comb likelihood
  is marginalized over the joint detector-frame `(M_det, chi)` nuisance model.
  A pre-existing posterior may be used as a proposal only under the validity
  conditions above; treating posterior density as likelihood is forbidden.
  Frequency-error-only propagation is inadmissible. The concrete estimator,
  proposal correction, and convergence rule are owner-fixed cells
  [OWNER-FIX: estimator, proposal correction, and convergence criterion].
- Redshift and frame: observed frequencies use `M_det=(1+z)M_source` in both
  `Omega_H` and the tooth offset. The offset-subtracted ratio is redshift-free
  only when the offset and teeth use the same frame. Any source-frame
  conversion and its covariance belong to the frozen nuisance model.
- Selection effects and trials accounting: owner-fixed cells
  [OWNER-FIX: selection model], [OWNER-FIX: trials correction].

## 5. Likelihood specification (structure only)

Per event, the contract must compare two hypotheses on the same detector data
`D` through a common strain/readout likelihood interface:

- H_comb: the secondary feature frequencies sit at the template teeth
  f_{k,m} = m*Omega_H/(2*pi) + c^3*g(chi)*ln(k)/(16*pi^2*G*M_det) for one
  integer assignment per feature, k in {2, ..., 12}. A separate normalized
  pre-data transition prior across k is required [OWNER-FIX: transition-rate
  model and normalized prior]. The within-line `(k-1)/k` KMS net-response
  factor may multiply a GR greybody response but does not supply that prior
  [OWNER-FIX: greybody normalization artifact]. Linewidths come from the
  nuisance model.
- H_no_comb: the published no-comb description of the same features
  (owner-fixed reference model [OWNER-FIX: no-comb reference]).

For each hypothesis, the evidence is

    Z(H) = integral L(D | theta, H) pi(theta | H) d theta.

The per-event statistic is the Bayes factor BF = Z(H_comb)/Z(H_no_comb).
It cannot be obtained merely by integrating over a published posterior from a
different model; any posterior-recycling estimator must divide out the known
sampling prior, retain the required likelihood/evidence normalization, cover
the target support, and pass frozen importance-weight diagnostics. Stacking
across events is performed in the universal coordinate computed consistently
with `M_det`. This section fixes structure only: no numeric
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

- BOUNDED-SUBMODEL KILL band: a resolved pair clearing the gates whose offset-subtracted
  ratio excludes every ratio ln(k_a)/ln(k_b) with k_a, k_b in the
  ladder set at the frozen confidence, under the frozen likelihood and
  covariance. This kills only the prospectively normalized finite submodel
  with k in {2,...,12}. It cannot kill the unrestricted integer-k family:
  ratios ln(m)/ln(n) with unbounded integers are dense on the positive line,
  so a finite-resolution interval admits sufficiently large integer pairs.
  A claim about all k therefore requires a normalized prior with a declared
  tail and a likelihood-level exclusion of its prior-predictive distribution.
  If a source theory later supplies a fixed `d_before`, its divisor support
  replaces the unrestricted family and must be frozen before comparison.
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
target-blind and its receipt contains no event-derived number. The derived
source transition, strain/readout likelihood or sufficient likelihood product,
detector-frame nuisance convention, normalized pre-data prior over k
(including a source/declared transition-rate model and the finite-submodel or
tail decision), finite-step background treatment, event list and cuts,
selection model, trials accounting, greybody/Page normalization, no-comb
reference, and estimator convergence receipt are open prerequisites. This
draft is not anchorable until those artifacts exist. Nothing in this directory
may be cited as a registered, frozen, or scored prediction.
