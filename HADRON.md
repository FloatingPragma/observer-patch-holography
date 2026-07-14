# Hadron Data Policy

This file is the operating plan for particle rows whose final numeric values
depend on low-energy QCD.

## Closure-Program Note (2026-07-14)

The hadronic lane is generator G1 of the closure program
([CONSISTENCY_STACK.md](CONSISTENCY_STACK.md), ledger dependency table): the
blind Ward-projected transport computation against the frozen target is the
object that closes or moves CL-1 and CL-2. The blind target is frozen and
externally timestamped (OpenTimestamps,
`../falsification/frozen_targets/hadronic_closure_target_2026-07-14.json`).
The payload computation must not read the target file; the dependency-cone
audit checks this. Pass/fail outcomes are declared in the target file.

## Decision

The particle program uses two public numerical surfaces.

1. Source-only OPH values stay primary for sectors that do not require a
   nonperturbative hadronic spectral payload.
2. Rows controlled by confined quarks and gluons use a clearly marked empirical
   hadron closure when a final display value is needed.

The empirical closure uses measured \(e^+e^-\to\mathrm{hadrons}\) data through
the standard dispersion relation. This is the strongest available route in the
absence of a production OPH hadron backend, because the optical theorem turns
the measured hadronic cross section into the same electromagnetic spectral
function that enters the Ward-projected Thomson transport.

## Terminology

The hard barrier is nonperturbative QCD. Perturbative QCD is useful
for high-energy tails, matching windows, and consistency checks. It cannot
replace the low-energy resonance-region spectrum. That region contains the
\(\rho\), \(\omega\), \(\phi\), charm threshold, and other hadronic structures
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

For the empirical closure, the hadronic term is computed from \(R(s)\), the
ratio of the bare cross section for \(e^+e^-\to\mathrm{hadrons}\) to the
pointlike \(e^+e^-\to\mu^+\mu^-\) cross section. The dispersion form is

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

The empirical \(e^+e^-\to\mathrm{hadrons}\) closure should compute that
transport gap as a dispersion integral, with the uncertainty and compilation
choice shown in the final row.

## Data Sources

Primary ingest candidates:

- PDG hadronic cross-section and \(R\)-ratio data files:
  `https://pdg.lbl.gov/current/xsect/`
- HEPData records for inclusive \(R(s)\) and exclusive hadronic channels:
  `https://www.hepdata.net/`
- alphaQED by F. Jegerlehner, including hadronic running of
  \(\alpha_{\mathrm{QED}}\):
  `https://www-com.physik.hu-berlin.de/~fjeger/software.html`
- KNT19, `arXiv:1911.00367`, for a modern data-driven HVP compilation.
- DHMZ and HVPTools analyses for an independent data-combination path.
- Muon \(g-2\) Theory Initiative 2025 white paper, `arXiv:2505.21476`, for
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
\(e^+e^-\to\mathrm{hadrons}\) data.

Use `compare_only` for rows that are useful diagnostics but do not define the
public prediction.

Use `work_in_progress` for in-scope rows whose computation is not part of a
public final value.

## Refactor Plan

1. Add the row-class field to generated particle status manifests and final
   prediction tables.
2. Split public output into:
   `source_only_value`, `empirical_hadron_closure_value`, and
   `reference_value`.
3. Add `code/particles/hadron/empirical_ee_hadrons_sources.yaml` as the source
   registry.
4. Add a schema for `oph_empirical_ee_hadronic_spectral_measure`.
5. Add an ingest script that can read PDG/HEPData tables without silently
   applying physics corrections.
6. Add a normalization script that converts cross sections or \(R(s)\) into the
   Ward-projected spectral convention used by `code/P_derivation`.
7. Add a dispersion integrator with covariance propagation and compilation
   labels.
8. Add `P_derivation/empirical_thomson_endpoint.py` to combine the empirical
   hadronic term with the OPH D10 source anchor and the lepton kernel.
9. Add a final fixed-point solver that emits both the empirical display row and
   the source-only audit row.
10. Update papers, README files, learning material, and figures so public prose
    states the two-surface policy once and keeps ledger details in the status
    section.

## Claim Policy

The public claim is:

```text
OPH fixes the screen-cell fixed-point equation and the electroweak source
frame. The low-energy hadronic spectral term is supplied by the best available
empirical \(e^+e^-\to\mathrm{hadrons}\) data until a source-only OPH hadron
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
The fine-structure constant is fixed by a screen-cell fixed point. OPH computes
the outer pixel equation and the electroweak source anchor. The charged-lepton
transport term is explicit. The confined-hadron term is the electromagnetic
spectral function of QCD. Since the OPH hadron backend does not emit that
spectral function directly, the public final value uses the standard
data-driven dispersion method based on measured e+e- -> hadrons cross sections.
That is the strongest available empirical closure: it measures the spectral
function that the theorem requires. Source-only OPH and OPH plus empirical
hadron closure are shown as separate output classes.
```
