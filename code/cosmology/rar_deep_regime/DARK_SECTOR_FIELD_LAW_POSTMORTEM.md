# Fixed-input QUMOND response: Cassini and certification boundary

Status: executable-evidence note, 2026-08-24. This note is scoped to the two
interpolation functions and fixed inputs recomputed by
`qumond_quadrupole_crosscheck.py`. It is not a no-go theorem for every local,
screened, environment-dependent, nonlocal, or dynamical dark-sector law.

## What was computed

An earlier OPH continuation used the unscreened QUMOND response

    laplacian Phi = div[nu(|grad Phi_N| / a_0) grad Phi_N],
    laplacian Phi_N = 4 pi G rho_b,

with either the simple or radial-acceleration interpolation function. Applied
unchanged to the Sun in the Galactic external field, the radial-acceleration
function at the declared Park et al. inputs gives
`Q_2 = 3.39e-26 s^-2`, compared with the Cassini estimate
`(1.6 +- 1.8)e-27 s^-2`. The published fixed-input response is in strong
tension with that benchmark.

The inner-multipole computation in this directory reproduces the
Blanchet--Novak simple-function value to about two percent and the Park et al.
radial-acceleration central value to `8e-5`. This agreement is useful
diagnostically, but the local integrations emit thousands of warnings and
have neither a singularity-aware split nor certified radial-tail bounds. The
local receipt therefore does not establish that its residual error is small
and cannot independently certify the published tension.

## What the local receipt does not decide

The historical OPH continuation applied the unscreened QUMOND response at both
galactic and Solar-System scales. Its phantom halo is anisotropic in an
external field, and published calculations report a quadrupole incompatible
with the declared Cassini comparison. That published result can motivate
setting the continuation aside, but the current local receipt is not an
independent numerical retirement certificate. Such a certificate needs a
singularity-aware formulation, retained cubature errors, analytic or rigorous
tail bounds, cutoff and precision refinement, and an independent
high-precision replay.

Even a certified exclusion would apply only to that fixed response grammar
and its declared parameters. It would not establish that every possible local
field law is excluded. Rapid-transition functions, derived screening or
environmental dependence, extra degrees of freedom, nonlocal responses, and
dynamical completions remain logically viable until individually specified
and tested. Such a route must earn its extra structure prospectively rather
than add a switch after seeing Cassini.

## What remains viable

The paper's current continuation instead treats anomalous modular energy as a
density carried on overlap collars. Near a compact source it conditionally
uses an ambient-density branch. An isotropic ambient density has zero
trace-free quadrupole; for the declared local density its traceful tidal scale is
`1.9e-31 s^-2` and the enclosed anomalous mass inside ten astronomical units
is about `5e-15` solar masses.

The Cassini `Q_2` result constrains the trace-free quadrupole channel, not
this isotropic term. Their magnitudes may be displayed as a dimensional scale
comparison, but Cassini is not an observational bound on the isotropic tidal
channel.

That arithmetic does not derive the density branch. The galactic-scale
carrier attachment, transition profile, full relativistic stress (including
pressure and anisotropic stress), lensing relation, normalization, and joint
likelihood remain open. The current matched SPARC diagnostic excludes
bulge-ambiguous galaxies and gives a paired interval containing zero after
changing the observable, radius, sample, and combination rule together. This
is consistent with a mixed-proxy explanation of the older displacement but
does not causally isolate it and neither validates nor falsifies the density
continuation. The separate penalized objective is an uncalibrated sensitivity
analysis, not a source-value rejection band.

## Receipts

- `qumond_quadrupole_crosscheck.py`
- `receipts/qumond_quadrupole_crosscheck.json`
- `rar_deep_regime_diagnostic.py`
- `receipts/sparc_deep_regime_diagnostic.json`
- Blanchet and Novak, MNRAS 412, 2530 (2011), arXiv:1010.1349
- Park et al., arXiv:2602.17884 (2026)
