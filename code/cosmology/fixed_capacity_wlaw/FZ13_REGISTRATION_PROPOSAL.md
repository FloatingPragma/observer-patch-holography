# FZ-13 registration proposal: fixed-capacity dark-energy w-law, pre-DESI-DR3

Status: PROPOSAL. This document is not a registration. It claims no
registered, frozen, or scored prediction. Anchoring, adoption into
`claims/frozen_prediction_register.json`, the numeric kill-band cells, the
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

OWNER SLOT (set at adoption; this document is a proposal and carries no
status in the register).

### content (proposed, verbatim-ready)

Fixed-capacity dark-energy stance, conditional on premises P1 through P4 of
`proof/epic_wins/dk01_wlaw/DK01_FIXED_N_WLAW.md` (fixed total record
capacity N over the scored range, the de Sitter capacity relation
Lambda l_P^2 = 3 pi / N, constant l_P, no repair-rate drift and a closed
dark sector). Fixed branch: (w0, wa) = (-1, 0) exactly in the CPL
parameterization, with zero continuous freedom; the capacity value, the
open closure equation N = F(N), and the capacity-vs-display gap all cancel
out of w(a). Monotone branch (capacity nondecreasing in the scale factor):
w(a) >= -1 at every scored scale factor; no phantom epoch. Exposure map
(exact, both directions): w(a) = -1 + (1/3) d ln N / d ln a, so any
established w(a0) < -1 forces capacity loss dN/da < 0 at a0, and any
established thawing signal forces capacity growth, contradicting the fixed
branch. Machine-checked in
`Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean`.

### comparison_protocol (proposed)

Post-freeze published DESI DR3 and Euclid w0waCDM posteriors only. The
estimator is the published posterior in the (w0, wa) plane from the
collaboration chains of the frozen combination; no re-analysis, no
re-weighting, no in-house likelihood evaluation. The frozen combination
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
  posterior at a declared sigma level (OWNER SLOT, numeric cell fixed at
  freeze), with the posterior mean in the thawing quadrant
  (w0 > -1, wa < 0). Kills the fixed-capacity branch of this registration.
- FZ13-R02 (FAIL, monotone branch): an established phantom crossing, that
  is w(a) < -1 established on the frozen combination in a declared
  statistical sense at a declared level (OWNER SLOT, numeric cell and
  statistic fixed at freeze). Kills the monotone-capacity branch as well,
  through the machine-checked readback that w < -1 forces capacity loss.
- Proposed rider (sign check, inherited from the DK-01 mock record): a
  posterior mean in the freezing direction (w0 < -1, wa > 0) is scored
  against the display-bias mechanism recorded in the theorem file before
  any physical reading, at whatever level it appears.
- Pass and no-verdict bands (proposed structure): (-1, 0) inside a declared
  inner credible region scores a forward pass for the fixed branch;
  intermediate outcomes carry the exposure forward with no verdict. Numeric
  cells are OWNER SLOTs.

Per the corpus kill-condition protocol, the audit, re-audit, and repair
ladder runs before any failure is banked as a framework kill.

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
