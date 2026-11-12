# Native baryon weak-current form factors

## Source and current

The source is exactly the positive native action in
`model/openwave_calibrated_hadrons`, on its complete 2024-dimensional
three-quark occupation sector. Its colour-singlet embedding has 120 columns.
The weak raising matrix is the literal native exterior-representation
generator, restricted to its quark module. For each cell, second quantize
`sigma_mu tensor I_plus` and `sigma_mu tensor I_3`, where
`sigma_mu=(I,sigma_x,sigma_y,sigma_z)`. The full occupation currents
intertwine the colour embedding. No composite transition matrix defines
the source current.

Prepare the lowest charge-zero, spin-up, even-cell-parity state in the
spin-half/isospin-half representation. Applying the actual spin-lowering
and weak-raising currents produces its neutron/proton spin quartet. The
isospin-breaking mass term is a scalar on each fixed-charge sector, so it
changes the energies without changing these representation vectors.
The total vector current has the exact two-dimensional matrix element
`p† I_plus n=I_2`. Rotational covariance of the internal spin gives
`p† J_plus^a n=g_A sigma_a`. This concerns internal spin on the declared
graph; a relativistic spacetime interpretation requires an additional map.

## Exact local calculation and finite-hopping bound

In the co-located colour singlet, start with an integer epsilon-colour
vector of charge zero and spin projection one half. Apply the exact
polynomials `15-4S²` and `15-4I²`, then add its cell-reflected vector.
Their omitted scalar denominators cancel in normalization. The stored
integer occupation vector gives exactly

`g_V=1`, `g_A=5/3`, `||(J_plus^z-(5/3)I_plus)n||²=8/9`.

The latter is actual current leakage into other baryon states. Projecting
the current onto the nucleon quartet and treating that quartet as an
invariant weak sector would discard a nonzero source term.

In fixed `S=I=1/2`, positive cell parity and fixed charge/spin projection,
the co-located space is one-dimensional. Its orthogonal complement has
split-site quarks. Its colour Casimir cost is at least `8U/3` and its
local baryon-number-square cost is at least `5A/9`. Positive kinetic and
local spin/weak terms can be dropped for a lower bound. The global weak
term is `3 eta/4` on the whole sector. The co-located trial energy is

`E_loc=3mbar+delta/2+A+3(Js+Jw+eta)/4+3t`.

Thus the complement exceeds that trial by at least

`gamma=8U/3-4A/9-3(Js+Jw)/4-3t > 0`.

The ground energy is at most the trial energy, and the off-diagonal
coupling norm is at most `3t`. The complement component of the normalized
ground vector is therefore bounded by `3t/gamma`. On three quarks,
`||I_minus||<=2` and `||J_plus^z||<=3`. The pure-state expectation bound
for `I_minus J_plus^z` gives

`|g_A-5/3| <= 12(3t/gamma)`.

All constants in this bound are exact rational numbers, including the
stored calibrated coefficients. The result is the enclosure
`1.665734201267 < g_A < 1.667599132066`. Floating evaluation gives
approximately `1.666666666668`; its extra displayed digits are diagnostic,
not a tighter certified interval.

## Physical comparison

There is also an exact restriction on scalar repairs. The co-located
colour-singlet three-quark space is
`Sym³(C²_spin tensor C²_weak)=(S,I)=(1/2,1/2) direct-sum (3/2,3/2)`.
The Cauchy decomposition of the symmetric cube gives these irreducible
representations with multiplicity one, of dimensions four and sixteen.
By Schur's lemma, any scalar deformation preserving local colour, spin
and weak symmetry acts as a scalar on the nucleon representation. It
cannot change its states or the fixed native-current reduced matrix
element `5/3`. An axial repair therefore requires different state mixing
(for example off-site, orbital or sea degrees of freedom) or a justified
change of the effective current. Merely refitting scalar co-located
energies cannot remove this obstruction.

The PDG 2025 neutron listing gives `lambda=g_A/g_V=-1.2753 +/- 0.0013`.
The comparison uses its magnitude because the native Pauli-current
coefficient uses the positive sign convention. This datum is held out of
the hadron mass calibration. The certified interval excludes the quoted
magnitude by much more than its uncertainty. The calibrated action with
this bare-current identification fails the physical axial-ratio test.
It supplies no renormalized QCD current, recoil corrections or observed
weak lifetime. This failure concerns this specified action and dictionary.

The computation preserves all excited current channels. Its use of a
finite projector invokes the existing complete Lüders-instrument theorem
once the projector identities hold; it adds no Lean theorem count.
