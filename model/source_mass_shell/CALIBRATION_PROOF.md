# Physical units and the source clock

The source evolution has the convention
`U(tau)=exp(-i H_model tau)`. Write `T>0` for the number of SI seconds per
unit of this same `tau`. If `v>0` is the continuum signal speed in source
coordinates, identifying that signal with physical light fixes

\[
 \ell=\frac{c_{\rm SI}T}{v},\qquad
 E_{\rm GeV}=\frac{\hbar_{\rm GeV\,s}}{T}E_{\rm model}.
\]

Here `ell` is metres per source coordinate unit. The convention `v=1`
therefore gives `ell=c_SI*T`. This matching is an explicit physical
identification. It supplies no value of `T`.

## Exact SI coefficients

The [SI definitions](https://www.nist.gov/pml/special-publication-330/sp-330-section-2)
fix `c_SI=299792458 m/s`, `h=6.62607015e-34 J s`,
`e=1.602176634e-19 C` and the unperturbed cesium-133 ground-state hyperfine
frequency `nu_Cs=9192631770 Hz`. Thus `1 GeV=1.602176634e-10 J` exactly.
Write `J_GeV=1.602176634e-10` for the numerical conversion from GeV to joules.
The elementary charge here converts electronvolts into joules; no electron
rest mass enters the conversion. The code evaluates the resulting rational
coefficients exactly. The only irrational coefficient is
`hbar=h/(2*pi)`, enclosed by the exact identity
`pi=16 atan(1/5)-4 atan(1/239)` and alternating rational series bounds.

Suppose a positive interval `[epsilon_lo,epsilon_hi]` encloses the source
cesium transition gap of **this same generator** `H_model`, and a separate
physical-process identification relates this source transition to the
unperturbed atomic transition. Its phase is `epsilon*tau`, while the SI phase
is `2*pi*nu_Cs*T*tau`. Therefore

\[
 T=\frac{\epsilon}{2\pi\nu_{\rm Cs}},\qquad
 E_{\rm GeV}=\frac{h\nu_{\rm Cs}}{J_{\rm GeV}}\frac{g}{\epsilon}.
\]

For `0<=g_lo<=g<=g_hi`, monotonicity gives the exact enclosure

\[
 E_{\rm GeV}\in
 \left[\frac{h\nu_{\rm Cs}}{J_{\rm GeV}}\frac{g_{\rm lo}}{\epsilon_{\rm hi}},
       \frac{h\nu_{\rm Cs}}{J_{\rm GeV}}\frac{g_{\rm hi}}{\epsilon_{\rm lo}}\right].
\]

The factors `2*pi` cancel in the energy ratio. The clock gap is an angular
generator gap, and the physical frequency is in cycles per second. Using
`hbar*nu_Cs` in the ratio would omit a factor `2*pi`.

`CsClockInput` requires separate provenance for the source gap, the physical
process identification, and the shared time normalization. Mismatched
generator identities, absent provenance, nonpositive clock gaps and implicit
binary floating inputs are rejected. These structural checks do not prove
the supplied physical premises. Every conversion with such an input is
explicitly conditional. In particular, a modular clock on a subalgebra, a
different co-driven apparatus, or a host execution timer does not acquire
the target Hamiltonian's time normalization through this arithmetic.

The package supplies no physical cesium clock input for the radial source.
Its absolute energy and mass entries are therefore null. It retains the
exact interval coefficients `hbar*g` in the free-scale expression
`E_GeV=(hbar*g)/T`, and `h*nu_Cs*g/J_GeV` in the conditional expression
`E_GeV=[h*nu_Cs*g/J_GeV]/epsilon_Cs_model`.

## Absolute scale is not determined by matching light speed

Fix the dimensionless source Hamiltonian, state, observer law, scheduler,
readout records and source spacetime coordinates. For every `lambda>0`,
replace the SI interpretation by

\[
 T' = \lambda T,\qquad \ell'=\lambda\ell,\qquad
 E'_{\rm unit}=E_{\rm unit}/\lambda.
\]

The physical speed remains `v*ell'/T'=v*ell/T=c_SI`. Every dimensionless
source evolution and readout probability remains equal. The physical phase
is also equal because `E'_SI*t'_SI/hbar=E_SI*t_SI/hbar`. Nonzero energies in
GeV vary by `1/lambda`. Consequently the source specification and light-speed
matching alone do not select an absolute energy scale. This statement
concerns the map from source time to SI time; it does not rescale the
dimensionless scheduler while keeping its transition rates fixed.

A source result proving an identified clock frequency, or an independently
measured frequency with a proved map to this generator, breaks this freedom.
The SI definition of the cesium second alone gives no value of the source
dimensionless cesium gap. Reading the simulator's dimensionless clock and
repair laws supplies their specified rates, without supplying a physical
frequency absent from those specifications.

For a relativistic continuum branch, the energy and momentum conversions
give `E_SI^2-c_SI^2*p_SI^2=(E_unit*M)^2`, where
`p_SI=hbar*p_model/ell` and `M` is the model rest energy. Then the rest mass
is `(E_unit*M)/c_SI^2`. A finite localized charged gap has no dispersion
premise merely because its energy can be divided by `c_SI^2`. The units
module consequently reports energy equivalents, with physical rest-mass
identification false.

The retained radial exterior one-body witness has continuum rest energy
`M=64`, inherited from its supplied mass coefficient. The finite interacting
27-site charged gap `g` has its separate full-Fock interval. Both coefficients
are evaluated by the same unit-conversion formulas. Their source domains
are recorded separately: a one-body exterior continuum limit and a finite
interacting bound-state response. No equation identifies `g` with the
relativistic rest energy `M` or gives `g` a moving-pole dispersion.
