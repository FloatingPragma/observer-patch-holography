# Boosts of the common native record field

This is a theorem about `MFR-common-native-field-v1`, using its sealed
all-moment continuum reconstruction and actually executed common cubic
member. It is independent of the positive bank/publication-order route.
Neither route requires Poisson event times or primitive unitary repairs.
The finite source graph and attempt order remain disclosed regulator data.
No exact Lorentz symmetry of that finite graph is asserted.

## Source and retained finite theorems

The common source derives the massive scalar covariance
`(p_E²+lambda)^(-1)` and physical Maxwell curl covariance from occupation
records. The masses squared are the native internal eigenvalues
`5-sqrt(5), 6, 5+sqrt(5)`. Its full Gaussian hierarchy, reflection quotient,
positive energy and causal commutators are proved in the sealed continuum
packet. Finite block, contact, volume, spacing and first-hit clock errors are
separate in the fixed common-member specification.

The existing RER `Geometry/FlatLorentzVolume.lean` proves
`boost_preserves_square`, `boostMatrix_det_one` and `boost_future_iff`.
The canonical Hermitian source module and golden causal/count theorems
remain valid on their stated interfaces. This proof uses the same Lorentz
quadratic form through the independently reconstructed field covariance;
it does not identify the two native candidates by their common limit.
The older `model/openwave_boosted_records/PROOF.md` supplies a useful
dispersion-error comparison but installs a different matter source. Its
particle/preparation claims are not imported.

## Invariant effective action and observer frames

For each massive band put `m=sqrt(lambda)>0` and
`dmu_m(p)=d³p/[(2pi)³ 2 sqrt(m²+|p|²)]`. This measure is the positive mass
shell restriction of `theta(p0) delta(p²-m²) d⁴p`; the determinant and cone
identities above prove its invariance. Therefore

    (U(Lambda) h)(p)=h(Lambda^(-1) p)

is unitary on the same source-derived one-particle space `L²(dmu_m)`.
Translations act by their mass-shell characters. Their tensor powers give
the Gaussian field representation. The Maxwell two-form covariance similarly
transforms with its two tensor indices and preserves the physical null
quotient. These are representations of limiting observer records, not new
updates of the native loads.

Inverting the measured Gaussian covariance gives the effective Euclidean
quadratic form. Its continuation is

    S = integral d⁴x [sum_b (partial_mu phi_b partial^mu phi_b
                            - lambda_b phi_b²)/2
                       - F_mu_nu F^mu_nu/4].

The displayed action summarizes the measured inverse covariance. It is not
executed instead of the native source. Tensor contraction and invariant
four-volume prove its Lorentz invariance. The vacuum hierarchy and the
state/observable transformation rule have no preferred inertial frame in
this limit. Finite-spacing anisotropy remains a measured error, bounded in
the sealed covariance and dispersion theorems; it is not discarded.

The scalar retarded response and Maxwell commutator have support inside the
common cone `|x|<=|t|`. Consequently a compact local source intervention has
zero linear response at spacelike-separated local observables in the limiting
field. Massive group speeds are `|p|/sqrt(m²+|p|²)<1`; photon speed is1.
This is a causal response statement, stronger than a group-speed bound
alone. The primitive attempt counter is not substituted for this clock.

## A state with a measurable contracting envelope

Length contraction of an arbitrary freely spreading packet is not an exact
rigid-body identity. The following controlled rest-packet family supplies
an explicit statement that can be tested without that false assumption.

Choose a nonzero real radial `f` in `C_c²(R³)`, supported in `|q|<=Q`.
Take the source one-particle state with mass-shell amplitude
`h_sigma(p)=C_sigma 2E(p) f(p/sigma)`, where `C_sigma` normalizes its
`L²(dmu_m)` norm. The field matrix element is

    psi_sigma(t,x)=C_sigma integral d³p/(2pi)³
                         f(p/sigma) exp(-iE(p)t+ip.x).

Its normal-ordered scalar intensity is
`<h_sigma|:phi(x)²:|h_sigma>=2|psi_sigma(x)|²`. This positive finite signal
has the vacuum contribution subtracted by the same source reference. It is
an intensity observable, not a claim that bare primitive occupation is a
relativistic position probability.

Apply the unitary boost giving velocity `v` along axis1, `|v|<1`, with
`gamma=(1-v²)^(-1/2)`. Scalar covariance gives exactly

    psi_sigma,v(t,x)=psi_sigma(gamma(t-vx1),
                              gamma(x1-vt),x2,x3).

At lab time0 remove the carrier phase `exp(i m gamma v x1)` and use
`y=(sigma gamma x1,sigma x2,sigma x3)`. Apart from an overall constant the
remaining amplitude is the Fourier transform in y of

    eta1(q)=q1+v(E_sigma(q)-m)/sigma,   eta_perp=q_perp,
    E_sigma(q)=sqrt(m²+sigma²|q|²),
    D(q)=det(d eta/dq)=1+v sigma q1/E_sigma(q),
    g_sigma(eta(q))=f(q)/D(q).                              (1)

This is a diffeomorphism since `D>=1-|v|>0`. Both its map and first two
derivatives converge uniformly to those of the identity on the compact
support as `sigma/m ->0`. Thus `g_sigma -> f` in H¹. This last assertion
also follows directly from the derivative formulas below and dominated
convergence, including the vanishing boundary traces of f.

For a real Fourier amplitude the intensity centroid is zero. Parseval gives
the exact longitudinal and transverse width ratios relative to the rest
intensity at time0:

    W1(v)/W1(0)=(1/gamma) sqrt[(I1/N)/(I0/N0)],
    W2(v)/W2(0)=          sqrt[(I2/N)/(I0/N0)],              (2)

where `N0=integral f²`, `I0=integral (partial_1 f)²`,
`N=integral f²/D`, and `Ii=integral |partial_eta_i g|² D dq`.
Radial f makes the rest widths equal. Hence

    W1(v)/W1(0) -> sqrt(1-v²),   W2(v)/W2(0) ->1.          (3)

These are normalized root-mean-square intensity widths of actual boosted
states in the source-derived field representation. They are not dimensions
of an assumed rigid classical object. Uniformity holds for velocities in
any fixed compact subset of `(-1,1)`.

For reproducible finite checks define `a_j=v sigma q_j/E_sigma` and
`D_j=partial_j D`. The exact integrands are

    partial_eta1 g = f_1/D² - f D_1/D³,
    partial_etaj g = f_j/D - a_j f_1/D² - f D_j/D²
                             + a_j f D_1/D³,   j=2,3,
    D_j=v sigma [delta_1j/E_sigma
                         - sigma² q1 qj/E_sigma³].         (4)

All denominators are bounded below. Compact integration with certified
derivative bounds therefore computes each finite ratio to any desired
tolerance. Equations(1)-(4) retain the spreading correction instead of
declaring the finite-width state to contract exactly.

## Source state functional and finite approximation

A positive-Euclidean-time scalar insertion at tau>0 with spatial transform
`Jhat(p)=exp(tau E(p)) h_sigma(p)` reconstructs precisely this one-particle
vector. The same formula with the boosted mass-shell amplitude reconstructs
the boosted vector. Both supports are compact in momentum, so these smearings
have finite source norms. Their inner products and intensity matrix elements
are two- and four-point functions of the same native record limit.

Signed or complex insertion scores are not classical acceptance
probabilities. Positivity is the proved reflection/state positivity. This
specifies a reconstructed preparation and observable functional; an executed
physical detector or a primitive coherent preparation is not claimed here.
Bounded spatial/time smearings and finite resolution approximate the displayed
intensity. At a fixed finite packet and detector menu, source moment and
refinement errors can be made arbitrarily small by the common-member rule.
On a compact energy band, real-time phases are uniformly approximable by
polynomials in the measured transfer generator; finite spectral arithmetic
and its error are retained. The narrow-packet limit follows after those
fixed-menu source limits, or by a diagonal choice of their tolerances.

## Coulomb comparison and scope

The zero-frequency inverse of the same Maxwell covariance is `1/|k|²`.
For a declared external conserved source its static Green response is
`1/(4pi r)`, with its electric field and boosted magnetic response related
by the same two-form transformation. This is the external-source Green
identity. It does not identify a primitive relative charge with electric
charge, prove an interacting charged-matter action, or make the earlier
row20 candidate the same source.

The theorem supplies invariant free action, frame covariance, causal
response and controlled boosted-state contraction on one common native
record family. The charged source/response dictionary, executable moving
state/detector apparatus, owner pin, in-platform test and other-group review
retain their separate status. No numerical check or same-group review
constitutes OpenWave acceptance.
