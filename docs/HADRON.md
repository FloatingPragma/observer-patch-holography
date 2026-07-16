# Hadron Data Policy

The operating policy for particle rows whose final numeric values depend on
low-energy QCD.

## Standing

The source chain fixes everything on its side of the frontier: each declared
P map has an interval-certified unique fixed point on its stated domain, the
charged-lepton transport term is a closed one-loop kernel, and the gauge-width
map lands 2.5×10⁻⁶ relative from the measured α⁻¹ with the hadronic transport
term open. That open term is a frontier of the entire field, a property of
the observable: closing CL-1/CL-2 asks for 4×10⁻⁹ relative precision on the
hadronic moment, beyond every method on Earth. The best data-driven
dispersive determinations carry relative uncertainties near 4×10⁻³, the
leading lattice programs reach a comparable order after campaigns at the 10⁷
core-hour scale, the region below 2 GeV is nonperturbative, and a source-only
payload additionally forbids measured hadronic input even for scale setting.
See "Why The Hadronic Test Is Hard" in
[OPH_FALSIFICATION_PROGRAM.md](OPH_FALSIFICATION_PROGRAM.md). Existing grid
work is exploratory.

## Closure-Program Note

The hadronic lane is generator G1 of the closure program
([CONSISTENCY_STACK.md](CONSISTENCY_STACK.md), ledger dependency table): a
genuinely target-blind Ward-projected transport computation under a valid
detached contract is the object that could close or move CL-1 and CL-2. The
target-contract provenance (externally timestamped v1/v2, v2's invalid scalar
scoring algebra, the permanently inactive post-target-access v3 erratum
scaffold, the exploratory non-blind V1 grid) is recorded in the reading rules
of [CLOSURE_LEDGER.md](CLOSURE_LEDGER.md). A valid experiment needs a
detached successor whose complete source method is frozen before genuinely
withheld data, or an audited clean-room producer with no target access. The
payload computation must not read the target file; the complete transitive
dependency cone and chronology require independent verification. No current
target declares an active pass/fail outcome.

## Decision

The particle program uses two public numerical surfaces.

1. Source-only OPH values stay primary for sectors that do not require a
   nonperturbative hadronic spectral payload.
2. Rows controlled by confined quarks and gluons use a clearly marked empirical
   hadron closure when a final display value is needed.

The empirical closure uses measured $e^+e^-\to\mathrm{hadrons}$ data through
the standard dispersion relation. This is the strongest available route in the
absence of a production OPH hadron backend, because the optical theorem turns
the measured hadronic cross section into the same electromagnetic spectral
function that enters the Ward-projected Thomson transport.

## Terminology

The hard barrier is nonperturbative QCD. Perturbative QCD is useful
for high-energy tails, matching windows, and consistency checks. It cannot
replace the low-energy resonance-region spectrum. That region contains the
$\rho$, $\omega$, $\phi$, charm threshold, and other hadronic structures
that dominate the dispersion integral.

## Fine-Structure Constant

The OPH fixed-point equation is

```text
P = phi + sqrt(pi) / A_T(P)
alpha(0) = 1 / A_T(P)
```

with

```text
A_T(P) = a0(P) + Delta_lep(P) + Delta_had(P) + Delta_EW(P).
```

The D10 source map emits `a0(P)`, and the charged-lepton kernel is implemented
as a closed one-loop transport term. The missing source-only piece is the
Ward-projected hadronic spectral measure plus a same-scheme endpoint
remainder.

For the empirical closure, the hadronic term is computed from $R(s)$, the
ratio of the bare cross section for $e^+e^-\to\mathrm{hadrons}$ to the
pointlike $e^+e^-\to\mu^+\mu^-$ cross section. The dispersion form is

```text
Delta alpha_had(q^2)
  = -(alpha q^2)/(3 pi) PV int ds R(s) / (s (s - q^2)).
```

For the OPH Thomson endpoint package this is converted into the same
subtraction convention used by `a0(P)` and the source kernel:

```text
Delta_had(P)
  = mZ(P)^2/(3 pi) int ds rho_Q(s;P) / (s (s + mZ(P)^2)).
```

The empirical spectral payload is not a source-only OPH theorem. It is the
best data-driven completion of the low-energy hadronic transport term.

## Current Fine-Structure Values

The public display value is

```text
alpha^-1(0) = 137.035999177(21)
alpha(0)    = 0.007297352564331425030245795264691683...
P           = 1.630968209403959324879279847782648941...
```

At the public endpoint pixel, the source audit records

```text
a0(P_C)               = 128.307965473286248209961108741756716187...
Delta_impl_exact(P_C) = 8.686567842708528400985442542885969768...
Delta_required(P_C)   = 8.728033703713751790038891258243283813...
Delta_source_residual = 0.041465861005223389053448715357314044...
```

The empirical $e^+e^-\to\mathrm{hadrons}$ closure should compute that
transport gap as a dispersion integral, with the uncertainty and compilation
choice shown in the final row.

## Data Sources

Primary ingest candidates:

- PDG hadronic cross-section and $R$-ratio data files:
  `https://pdg.lbl.gov/current/xsect/`
- HEPData records for inclusive $R(s)$ and exclusive hadronic channels:
  `https://www.hepdata.net/`
- alphaQED by F. Jegerlehner, including hadronic running of
  $\alpha_{\mathrm{QED}}$:
  `https://www-com.physik.hu-berlin.de/~fjeger/software.html`
- KNT19, `arXiv:1911.00367`, for a modern data-driven HVP compilation.
- DHMZ and HVPTools analyses for an independent data-combination path.
- Muon $g-2$ Theory Initiative 2025 white paper, `arXiv:2505.21476`, for
  current HVP context and lattice/data-driven tension policy.

The production pipeline should not rely on a single black-box number. It should
store the source compilation, correction convention, covariance treatment, and
integration kernel.

## Required Empirical Payload

The empirical payload should live beside the source-only payload, with a
different artifact name:

```text
oph_empirical_ee_hadronic_spectral_measure
```

Minimum fields:

```text
artifact
format_version
source_compilation
data_release
energy_grid
R_values
statistical_covariance
systematic_covariance
radiative_corrections
vacuum_polarization_undressing
final_state_radiation_policy
exclusive_channel_sum_policy
inclusive_region_policy
pQCD_tail_policy
kernel
integral_value
integral_uncertainty
guards
```

The guards must state:

```text
source_only = false
empirical_hadron_closure = true
external_cross_section_data_used = true
promotable_as_OPH_source_theorem = false
usable_for_public_final_values = true
```

## Output Classes

Every particle row should carry one of these classes.

```text
source_only_oph
oph_plus_empirical_hadron_closure
compare_only
work_in_progress
```

Use `source_only_oph` for rows whose forward value is emitted without importing
PDG, NIST, or hadronic cross-section data.

Use `oph_plus_empirical_hadron_closure` for rows where OPH fixes the source
frame and the remaining QCD spectral term is supplied by measured
$e^+e^-\to\mathrm{hadrons}$ data.

Use `compare_only` for rows that are useful diagnostics but do not define the
public prediction.

Use `work_in_progress` for in-scope rows whose computation is not part of a
public final value.

## Pipeline Architecture

The production pipeline, in its target form, comprises:

- a row-class field in generated particle status manifests and final
  prediction tables;
- public output split into `source_only_value`,
  `empirical_hadron_closure_value`, and `reference_value`;
- `code/particles/hadron/empirical_ee_hadrons_sources.yaml` as the source
  registry;
- a schema for `oph_empirical_ee_hadronic_spectral_measure`;
- an ingest script that reads PDG/HEPData tables without silently applying
  physics corrections;
- a normalization script that converts cross sections or $R(s)$ into the
  Ward-projected spectral convention used by `code/P_derivation`;
- a dispersion integrator with covariance propagation and compilation labels;
- `P_derivation/empirical_thomson_endpoint.py`, combining the empirical
  hadronic term with the OPH D10 source anchor and the lepton kernel;
- a final fixed-point solver that emits both the empirical display row and
  the source-only audit row;
- papers, README files, learning material, and figures whose public prose
  states the two-surface policy once and keeps ledger details in the status
  section.

The unbuilt pieces are work in progress.

## Claim Policy

The public claim is:

```text
OPH fixes the screen-cell fixed-point equation and the electroweak source
frame. The low-energy hadronic spectral term is supplied by the best available
empirical $e^+e^-\to\mathrm{hadrons}$ data until a source-only OPH hadron
backend emits the same spectral measure directly.
```

The public claim is not:

```text
OPH derives the full hadronic spectral measure without QCD computation.
```

The empirical closure can produce final display values. It cannot be used to
claim an arbitrary-precision source-only derivation.

## Learning-Material Summary

Use this wording for the particles section:

```text
OPH defines a screen-cell fixed-point map and computes certified roots of two
declared incomplete maps. Neither root is a physical Thomson-limit
fine-structure prediction. The charged-lepton transport term is explicit, but
the confined-hadron term requires the electromagnetic spectral function of
QCD. Since the OPH hadron backend does not emit that spectral function
directly, the public empirical-closure row uses the standard data-driven
dispersion method based on measured e+e- -> hadrons cross sections. That row is
an external-data validation surface, not a source-only prediction. Source-only
OPH and OPH plus empirical hadron closure are shown as separate output classes.
```
