# FZ-13 registration proposal: fixed-capacity dark-energy w-law, pre-DESI-DR3

Status: PROPOSAL. This document is referenced by a
`registered_pending_freeze` index row, but that listing is not a scientific
freeze and claims no frozen or scored prediction. Adoption of the completed
cells as a frozen target, anchoring, the numeric kill-band cells, the
`content_sha256`, the `frozen_utc` stamp, and the custody and attestation
records are owner actions, per the custody discipline of the DK-01 draft
target (see Custody below). Until the owner freezes and anchors a row, this
document binds no one and carries no verdict weight.

## Sources of record

- Theorem file (source of record for the w-law, its premises P1 through P4,
  the drift map, and the sign audit):
  `proof/epic_wins/dk01_wlaw/DK01_FIXED_N_WLAW.md` in the oph-meta
  metarepository.
- FZ-01 custody path (the anchored 2026-07-17 registration set carries a
  dark-energy w-law target file):
  `falsification/frozen_targets/fz01_2026-07-17/frozen_target_dark_energy_wlaw_2026-07-17.md`
  in the oph-meta metarepository. Custody, anchoring state, and any
  supersession decision for that file are owner matters; this proposal does
  not alter FZ-01.
- Machine-checked surface (this lane):
  `Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean`
  (fixed-N w-law, exact drift map, no-phantom bound on the monotone branch,
  strict phantom-to-capacity-loss readback, CPL projection; axioms limited
  to propext, Classical.choice, Quot.sound).

## Proposed row content (verbatim-ready cells)

The cells below follow the field names of
`claims/frozen_prediction_register.json`. Cells marked OWNER SLOT are fixed
by the owner at freeze time and are absent here by design.

### id

FZ-13

### owning_issue

742

### milestone

Pre-DESI-DR3 freeze: the row is eligible only if frozen and anchored before
the DESI DR3 cosmology release.

### status

OWNER SLOT (set only at freeze; the current pending index row remains
`registered_pending_freeze` and carries no verdict weight).

### content (proposed, verbatim-ready)

Conditional fixed-capacity dark-energy model-discrimination stance, not an
independent derivation of fixed capacity: conditional on premises P1
through P4 of
`proof/epic_wins/dk01_wlaw/DK01_FIXED_N_WLAW.md` (fixed total record
capacity N over the scored range, the de Sitter capacity relation
Lambda l_P^2 = 3 pi / N, constant l_P, no repair-rate drift and a closed
dark sector). Fixed branch: (w0, wa) = (-1, 0) exactly in the CPL
parameterization, with zero continuous freedom within that conditional
branch; the capacity value, the
open closure equation N = F(N), and the capacity-vs-display gap all cancel
out of w(a). This result is the algebraic consequence of identifying
rho_DE = kappa/N and defining effective w through the continuity equation;
the Lean theorem does not derive those identifications or the fixed-N
premise. Monotone branch (capacity nondecreasing in the scale factor):
w(a) >= -1 at every scored scale factor; no phantom epoch. Exposure map
(exact, both directions): w(a) = -1 + (1/3) d ln N / d ln a, so any
established w(a0) < -1 forces capacity loss dN/da < 0 at a0, and any
established w(a0) > -1 forces local capacity growth, conditional on the
same density and closed-sector premises. Machine-checked in
`Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean`.

### comparison_protocol (proposed)

Post-freeze published DESI DR3 and Euclid w0waCDM posteriors only. The
estimator is the published posterior in the (w0, wa) plane from the
collaboration chains of the frozen combination; no likelihood re-analysis
or target-dependent nuisance change is permitted. Such a posterior is a
model- and likelihood-dependent CPL projection of the expansion history,
not a direct pointwise measurement of the capacity law. The frozen
combination
(BAO release, CMB likelihood, SNe compilation, and the named fallback if
the combination is unavailable at decision time) is a declared choice fixed
by the owner in the anchored row before the DR3 cosmology release. DESI DR1
and DR2, and every contour published before the freeze, are seen data:
eligible for postdiction bookkeeping only and excluded from any prospective
comparison this row scores. Historical event and release names appearing in
the sources of record are provenance references to already-seen data.

### kill_band (proposed structure; numeric cells are OWNER SLOTs)

Scope premise: premises P1 through P4 (fixed branch) or P2 through P4 plus
monotone capacity (monotone branch) as named in the theorem file. Which
premise absorbs a failure is a declared exit named in the theorem file.

- FZ13-R01 (FAIL, fixed branch): the point (-1, 0) lies outside the
  two-dimensional credible region of the frozen combination's (w0, wa)
  posterior at a declared level (OWNER SLOT, credible-set construction and
  numeric threshold fixed at freeze).
  The direction of the posterior mean is not an extra gate: a sufficiently
  established exclusion in any direction conflicts with the fixed branch.
- FZ13-R02 (FAIL, monotone branch): the monotone-compatible CPL set
  `{(w0,wa) | w0 + wa(1-a) >= -1 for every a in the frozen scored range}`
  is excluded in the declared statistical sense at the declared level
  (OWNER SLOT, scored range, statistic, and numeric threshold fixed at
  freeze). This is the CPL image of the pointwise machine-checked bound and
  avoids turning a posterior-mean quadrant label into a physical crossing.
- Systematics rule (must be frozen before exposure): every nuisance,
  fallback, audit, and no-verdict condition is direction-neutral. The
  already-seen display-bias mock may motivate a nuisance check, but a
  freezing-direction posterior is not automatically discounted and no
  post-exposure direction-specific rider may veto FZ13-R01 or FZ13-R02.
- Compatibility and no-verdict bands (proposed structure): (-1, 0) inside a
  declared inner credible region is COMPATIBLE with the fixed branch but
  supplies no confirmation credit for fixed capacity or OPH, because the
  same point is the standard LambdaCDM null. Intermediate or protocol-failed
  outcomes carry the exposure forward with no verdict. Numeric cells are
  OWNER SLOTs.

Per the corpus kill-condition protocol, the audit, re-audit, and repair
ladder runs before any failure is banked as a framework kill.

### Eligible freeze boundary

An owner may legitimately freeze the conditional target, its named premise
exits, the post-freeze data release and likelihood combination, the CPL
projection, a direction-neutral nuisance protocol, and complete decision
thresholds before accessing the comparison. Freezing cannot turn P1
(fixed N), P2 through P4, or the monotone-capacity alternative into derived
facts, and compatibility with the LambdaCDM point cannot be scored as
positive evidence for those premises. The row is not anchorable while its
combination, scored scale-factor range, statistic, thresholds, fallbacks, or
systematics rules remain OWNER SLOTs.

### custody

OWNER SLOT (source commit, custody commit, and stamped package path are
fixed at anchoring).

### attestation

OWNER SLOT (OpenTimestamps or equivalent anchoring is an owner action).

### content_sha256

OWNER SLOT (computed over the adopted content bytes at anchoring).

### frozen_utc

OWNER SLOT (set at anchoring).

## What this document is

This document proposes row content, a comparison protocol, and the
structure of the kill bands for a pre-DESI-DR3 freeze of the
fixed-capacity dark-energy w-law, with the machine-checked Lean surface as
the proof artifact and the DK-01 theorem file as the source of record. It
is a proposal pending the owner's freeze. Every threshold left numeric is
an OWNER SLOT; every choice already made here (domain, estimator class,
seen-data exclusion, band structure) is a declared proposal and is
renegotiable by the owner until anchoring. Nothing in this document or in
the Lean surface examines, fetches, or numerically uses comparison data.
