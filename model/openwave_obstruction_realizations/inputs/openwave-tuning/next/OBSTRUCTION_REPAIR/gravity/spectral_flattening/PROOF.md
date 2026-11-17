# A positive relativistic spectrum with finite equal-time entropy

## Flat continuous source and complete generators

On doubled three-dimensional spinors write A = alpha dot p / |p|, beta A = -A beta, and phi = arctan(R^3 |p|^3). The unitary

V(p) = cos(phi/2) + beta A sin(phi/2)

satisfies V (alpha dot p) V† = H_flat, where

H_flat = (alpha dot p + beta R^3 |p|^4) / sqrt(1+R^6 |p|^6).

Both energies are exactly plus/minus |p|, each twice. This does not change the momentum variable or delete any band. Its negative projector is precisely the quartic source projector whose actual regional entropy was proved finite in `../continuum_source/`. Multiplication by V is unitary on the entire one-particle space. Conjugating the full standard Dirac Poincare representation gives a full unitary group representation, not just formal pairwise commutators. Its generators have the transported self-adjoint domains; a common invariant smooth core can be obtained by transporting a Dirac invariant core. The apparent angular dependence at p=0 is harmless for defining the unitary multiplier, with V(0)=I; the displayed ordinary derivative formulas are used only away from zero.

Spatial momentum is unchanged. V is a rotational scalar, so rotations retain their ordinary spin-plus-orbital form. The boost is

K_i = V [i(alpha dot p) partial_i + i alpha_i/2] V†
    = i H_flat partial_i + i V(alpha dot p)(partial_i V†) + i V alpha_i V†/2.

The second term is essential. The code checks analytic first/second derivatives, [H,K_i]=-i p_i, [K_i,K_j]=-i epsilon_ijk J_k including the spin zero-order part, and rotational covariance. The remaining translation commutators follow immediately from the derivative coefficient. `ALGEBRA.json` records roundoff-level identities and independent centered-difference derivative controls.

This one-particle representation preserves the vacuum projector and hence acts on the source particle/hole Fock representation. It does not assert unitary equivalence of the standard continuum Dirac *vacuum Fock representation*: their projectors differ at order zero in the ultraviolet. The finite-entropy source representation is the one constructed relative to the local beta filling in the preceding continuum packet.

## A constant mass must change the dressing

For a constant nonnegative mass m commuting with D, use

H_flat,m = [D + beta(m + R^3 D^4)] sqrt(D^2+m^2) / sqrt[D^2+(m+R^3 D^4)^2].

It has exactly the ordinary massive energies and retains a projector differing from the local beta filling by O(|p|^-3). The rotation angle is the difference arctan((m+R^3 r^4)/r)-arctan(m/r), with V_m(0)=I for m>0. By contrast, applying the massless V to D+beta m leaves a generally O(m/|p|) projector tail. The latter is not the same finite-entropy architecture. The numerical comparison retains both.

A commuting finite mass register can be treated blockwise, but this algebraic extension does not establish an interacting mass pole, covariant gauge-current implementation, native proper clock, or transfer of the raw432 species coefficient. The conditional gravitational normalization receipt remains in `../clock_join/`.

The flattened vacuum energy relative to the local beta product has ultraviolet difference per two-level block -|p|(1-b), with b=R^3 |p|^3/sqrt(1+R^6 |p|^6), hence O(|p|^-5). Its compact three-dimensional sum is ultraviolet convergent, unlike the quartic time-generator's O(|p|^-2) difference. This is a reference-energy component; no gravitational action or selection equation follows from it.

## Actual tilted planes are a separate test

For a boost speed u along x and proper observer momentum q on t=u x, both real laboratory roots are retained:

k_(s),x = gamma(q_x + s u |q|), k_(s),perp=q_perp,
E_(s)=gamma(|q|+s u q_x), s=plus/minus.

Their positive Jacobians are j_s=gamma(1-u s k_(s),x/E_(s)). If P_s(k) is the full rank-two source band projector, the sampled-plane CAR Gram and occupied covariance are

G(q)=sum_s P_s(k_s)/j_s, N(q)=P_-(k_-)/j_-.

Canonical positive whitening gives C_u=G^-1/2 N G^-1/2. It is an exact rank-two projector: the two rank-two band frames make an invertible four-column complete frame, whose polar factor is unitary. There are no extra quartic-dispersion real roots or discarded branch register. However G depends on q and generally differs from the ordinary constant Dirac Gram gamma(I+u alpha_x). Therefore this canonical observer map is spatially nonlocal. A different momentum-dependent unitary convention for the canonical apparatus is a distinct observer architecture, not an equality of this region map.

`observer_area.py` forms actual interval covariance matrices from C_u for each transverse momentum and integrates their Gaussian entropy. The area is in proper observer coordinates; each interval has two interfaces. Full four-component matrices and all bands are retained. The measured coefficient depends strongly on u and cut orientation. Resolution/control receipts accompany the verdict; this map fails the universal area gate despite the exact Poincare generator algebra.

## Region algebras must not be conflated

The original equal-time regional algebra uses chi_B, and it has the finite entropy proved earlier. The fully transported Dirac net V chi_B V†, evaluated in V P_D V†, has exactly the standard Dirac regional spectrum. It does not retain the original chi_B finite-entropy result. Neither exact momentum-space Poincare covariance nor static chi_A <= chi_B proves geometric locality of the original finite-entropy net.

The original source propagation has nonzero spacelike tails. `leakage.py` retains complete finite Fourier systems and operator norms as well as average norms; its sharp spectral edge is an explicitly separate approximation error. `continuum_leakage.py` independently computes the exact transverse-zero continuum kernel outside the light cone, subtracting its known Abel-regularized high-momentum term before convergent quadrature. Nonzero continuum tails rule out strict microscopic locality for this original regional net. Approximate locality at separations large compared with R needs quantitative domain/error bounds; finite averaged leakage alone does not establish it.

An alternative boost-orbit apparatus U_flat chi_B U_flat† will be tested additively after freezing this packet. It starts from the finite-entropy chi_B, and is therefore distinct from V chi_B V†. Its vacuum entropy is invariant by the actual unitary boost. The remaining task is an operational geometric interpretation, including cross-observer overlap and complete finite implementations, not another algebra-only assertion.
