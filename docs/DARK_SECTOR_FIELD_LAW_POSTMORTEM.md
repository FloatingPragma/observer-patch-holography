# Dark sector: why the field-law formulation failed Cassini

Status: tracking record, 2026-08-21. This document explains the result that
retired the first OPH dark-sector law. The dark-matter paper does not discuss
it; the paper states the density formulation that replaced it.

## What was computed

The first OPH dark-sector law was written as a local nonlinear Poisson
equation in the Newtonian field,

    div[ nu(|grad Phi_b| / a_0) grad Phi_b ] = 4 pi G rho_b,
    nu(x) = 1 / (1 - exp(-lambda sqrt(x))),   lambda = 1 - P/24,

with the acceleration constant `a_0 = 1.029e-10 m s^-2` taken from the
capacity de Sitter branch. On the SPARC radial-acceleration relation it
reached 0.133 dex scatter, the same as the empirical interpolating function it
coincides with at `lambda = 1`. Applied unchanged to the Sun inside the
Galactic external field it predicts a Solar-System quadrupole
`Q_2 = 3.62e-26 s^-2` against the Cassini value `(1.6 +- 1.8)e-27 s^-2`, a
19 sigma exclusion at fixed inputs and about 11 sigma after the Gaia
uncertainty on the Galactic acceleration.

## What was right

The arithmetic. The retired audit
(`oph-physics-sim/oph_fpe/cosmology/cassini_external_field.py`) implements
Park et al. (2026) Eq. 9b and reproduces their published benchmark to
`7e-5`. An independent first-principles computation, the inner multipole
expansion of the phantom potential with one integration by parts
(`code/cosmology/rar_deep_regime/qumond_quadrupole_crosscheck.py`),
reproduces the Blanchet and Novak (2011) value `4.1e-26 s^-2` for the simple
function to 2 percent, agrees with Eq. 9b to four digits on every case run,
and reproduces the Park benchmark to `8e-5`. The exponential tail of the
interpolating function does not reduce `Q_2`: the quadrupole is generated in
the shell where the field is of order `a_0`, a few thousand astronomical
units from the Sun, and the radial-acceleration function and the simple
function agree there to 0.2 percent. Park et al. report the same tension
for modified-gravity MOND in general, at 3 to 15 sigma depending on the
galaxy sample.

## What the error was

The formulation. A local field law in `grad Phi_b` gives every mass its own
phantom halo, and in an external field that halo is anisotropic. This is the
external field effect of modified-gravity MOND, and the Solar System excludes
it. Nothing in OPH requires that form. The OPH source is a density, the
anomalous modular energy carried on overlap collars, and the collars that
carry it are a property of the screen's record state on galactic scales.
Near a compact source the anomalous density is the ambient value. An
isotropic ambient density has zero quadrupole; its tidal scale
`4 pi G rho / 3` with the local dark density `6.8e-22 kg m^-3` is
`1.9e-31 s^-2`, four orders of magnitude below the Cassini uncertainty, and
the enclosed anomalous mass inside ten astronomical units is `5e-15` solar
masses. The galactic phenomenology is carried by the profile of the density
sourced by the galaxy as a whole, which in the deep regime is
`rho_A = sqrt(M_b a_0 / G) / (4 pi r^2)`.

The field law was adopted because it was the shortest route from the collar
cut count to a computable equation. It imported the external field effect as
a side effect. The applicability predicate added afterwards tried to switch
the law off in the Solar System and was correctly judged a post hoc fit; the
rotor completion that followed posited an action and kept the same failure
until a screening branch was supplied.

## What carried over

- The generic dark-sector theorems: under universal coupling the geometry
  responds to total modular charge, so the luminous-only Einstein relation
  holds exactly when the non-luminous charge vanishes.
- The modular-anomaly source with its bound by the collar recovery defect.
- The deep-regime law `a_A = sqrt(a_b a_0)`, `v^4 = G M_b a_0`, as the
  consequence of the linear enclosed-mass profile.
- The SPARC deep-regime comparison with one fitted constant.
- The Cassini bound, as a bound the density formulation passes by four orders.

## What was discarded

- The interpolating function `nu` and every statement that depends on its
  form at intermediate accelerations.
- The settled-domain applicability predicate.
- The rotor action.
- The `a_0` built from `H_0` and `Lambda`; `a_0` is a fitted comparison value
  until the magnitude target closes.

## Lesson

Do not write an OPH response as a local nonlinear law in the Newtonian field.
Derive the density from the source, and compare the Solar System with the
density footprint. The open question is the one the paper states: whether the
anomalous charge is generated on galactic-scale collars and saturates near
compact sources. That is a theorem target on the collar recovery defect, not
a fit.

## Receipts

- `code/cosmology/rar_deep_regime/receipts/qumond_quadrupole_crosscheck.json`
- `code/cosmology/rar_deep_regime/receipts/sparc_deep_regime_diagnostic.json`
- `oph-physics-sim/docs/BEST_OF_PUBLIC_DATA_COMPARISONS.md` (retired audit)
- Blanchet and Novak, MNRAS 412, 2530 (2011), arXiv:1010.1349
- Park et al., arXiv:2602.17884 (2026)
