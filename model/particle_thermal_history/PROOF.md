# Conditional particle thermodynamics and finite-source controls

The particle input is the complete calibrated quadratic spectrum in
`model/conditional_particle_spectrum/results.json`. Its source bindings
are authenticated before use. It contains 17 types, with 90 fermionic
and 28 bosonic degrees of freedom per momentum. A Dirac quark has twelve
degrees including colour, spin and antiparticles; a charged lepton has
four. Each Majorana neutrino mass eigenfield has two, with no additional
antiparticle or sterile copy. The W entry includes both charges and has
six polarizations; Z, Higgs, photon and gluon entries have 3, 1, 2 and 16.
Their masses and conventions are retained without changing the inputs.

## The actual finite source

The source has two opposite nonzero momenta, 180 fermionic canonical modes
and 56 bosonic modes. Its global boson occupation cap is two. Enumerating
the vacuum, all single occupations and all two-mode multisets gives
`1+56+56*57/2=1653` orthogonal states. Their energy-class multiplicities
are checked against every projector rank of the original finite instrument.
For a supplied positive source temperature the complete Gibbs probabilities
are `exp(-E_j/T)/Z`, where `Z` sums all 1,653 weights. Every occupation
record is retained, including positive weights below binary64 range;
their finite log probabilities are stored explicitly.

The fermionic partition function is the exact finite product
`Z_F=product_(a=1)^180(1+exp(-omega_a/T))`. It accounts for all `2^180`
occupation states. Energies are quasiparticle excitation energies above
the selected quadratic vacuum; the additive vacuum constant is subtracted.
Both per-mode log probabilities are retained. The
Gibbs entropy identity `S=log Z+<E>/T` holds for both sectors. These states
are supplied equilibrium preparations of the existing quadratic action;
the native repair law is not claimed to prepare them.

The global boson cap prevents independent geometric occupations. At
infinite supplied temperature its mean total number is exactly
`(56+2*1596)/1653=3248/1653`, less than two. The uncapped occupation product
on the same 56 modes has unbounded mean as temperature increases. The
executed finite controls quantify this distinction. In particular the
continuum gas below is a declared mathematical extension, not the
finite-cap source's own abundance law.

## Declared continuum gas and its identities

Extend the fixed quadratic dispersions to independent free occupations
on three-dimensional continuous momentum, with a common supplied
temperature and zero chemical potentials. The volume limit excludes a
separately populated zero-momentum condensate. The energy density excludes
the zero-temperature vacuum contribution. Let `y=m/T`,
`epsilon=sqrt(x²+y²)`, and `f=1/(exp(epsilon)+1)` for fermions or
`f=1/(exp(epsilon)-1)` for bosons. Each physical degree contributes

```
n/T³   = (1/(2 pi²)) integral_0^infinity x² f dx,
rho/T⁴ = (1/(2 pi²)) integral_0^infinity x² epsilon f dx,
P/T⁴   = (1/(6 pi²)) integral_0^infinity x⁴ f/epsilon dx,
s/T³   = (rho+P)/T⁴.
```

The last identity follows at zero chemical potential from the Gibbs
partition function. Differentiating at fixed mass, or integrating by
parts, gives `dP/dT=s` and

```
(d rho/dT)/T³ = (1/(2 pi²)) integral x² epsilon² f(1-/+f) dx,
(d P/dT)/T³   = (1/(6 pi²)) integral x⁴ f(1-/+f) dx.
```

The minus sign is fermionic, the plus bosonic. All integrands are
nonnegative, and `x²<=epsilon²`; hence the heat capacity is positive and
the equilibrium sound-speed ratio satisfies `0<c_s²=(dP/dT)/(d rho/dT)<=1/3`.
The exact trace identity is
`rho-3P = g T⁴ y²/(2 pi²) integral x² f/epsilon dx >=0`.
These statements hold for the declared free model at every positive
temperature. Number, energy and entropy fractions use their own
normalizations and are not interchangeable.

The massless integrals give number weights one for a boson and `3/4`
for a fermion, and energy/entropy weights one and `7/8`. Thus this menu
has exact effective counts `g_number=191/2` and `g_energy=g_entropy=427/4`.
Its quark/lepton fractions are `108/191` and `27/191` for number, and
`252/427` and `63/427` for energy. The fermion-only quark fraction is
`4/5`. These are ideal equilibrium comparisons implied by the supplied
menu, not measurements of native defects or an independent cosmological
prediction.

## Numerical controls

Adaptive real quadrature calculates the massive integrals. Its error
estimates are numerical diagnostics, not rigorous interval certificates
or physical model errors. Independent fugacity-series controls use

```
n/T³   = y²/(2 pi²) sum_l sign_l K2(l y)/l,
P/T⁴   = y²/(2 pi²) sum_l sign_l K2(l y)/l²,
rho/T⁴ = y²/(2 pi²) sum_l sign_l [3 K2(l y)/l² + y K1(l y)/l].
```

For fermions `sign_l=(-1)^(l-1)`, and the alternating truncation error
is at most the first omitted positive term. For bosons all signs are
positive. The integral representation of `K_nu` gives
`K_nu((l+1)y)<=exp(-y)K_nu(ly)`, so the omitted sum is bounded by its
first term divided by `1-exp(-y)`. Controls at `y=0.1,1,10` check both
statistics; special-function roundoff is not given an interval guarantee.

For `y>=200`, exponentially suppressed contributions are omitted with
the displayed analytic per-degree dimensionless bound

`2(768+192y+16y²) exp(-y/2)/(2 pi²)`.

Indeed `epsilon>=(x+y)/2`, `epsilon<=x+y`, and both occupation and its
temperature-derivative factor are at most `2 exp(-epsilon)` in this
domain. Integrating the corresponding polynomial against `exp(-x/2)`
gives the bound. Its logarithm is retained even when binary evaluation
underflows. This omission never asserts an exactly zero massive density.

## Conditional expansion and physical admission

Under the additional homogeneous adiabatic assumption `a³s=constant`,
the calculated gas gives `a(T)/a(Tref)=[s(Tref)/s(T)]^(1/3)`. If a flat
Friedmann density law is also supplied, the gravitational constant cancels
from `H(T)/H(Tref)=sqrt(rho(T)/rho(Tref))`. The derivative identity gives
`-d log a/d log T=(d rho/dT)/(3s)`. The receipt uses `Tref=10 GeV`.
These are consequences of explicitly supplied cosmological equations;
they derive no native expansion, cosmic seconds or freezeout event.

Fixed running quark and charged-lepton coefficients at the parent MZ
scale are used as free action masses, without conversion to thermal
pole masses or temperature-dependent renormalization-group evolution.
The low-temperature free coloured gas omits confinement and hadron
replacement; its massless gluons alone prevent a nucleosynthesis or
late-neutrino interpretation. The fixed Higgs background omits the
electroweak thermal crossover. Common neutrino temperature and zero
chemical potentials omit decoupling and conserved asymmetries. Pions
are not added on top of free quarks without a matching construction.
No temperature in the computed curve is certified as a physical
early-Universe population of the native observer.

The separate interacting-equation-of-state comparison uses selected
rows from the primary table of Saikawa and Shirai,
[JCAP 05 (2018) 035](https://arxiv.org/abs/1803.01038), at exactly the
published GeV temperatures. The table includes interaction and crossover
effects absent from this ideal gas. Its quoted energy/entropy errors
and the free-model differences are retained separately. No parameter
is fitted to that table and no statistical agreement claim is made.
`comparison_inputs.json` records the selected rows, original data hash
and selection rule.
