# A finite spectral left-Weyl source and its free field limit

The declared program is `coherent-native-spectral-weyl-v1`. Its matter space
at each cell is the full exterior algebra of `C2_spin tensor R15`, in
spin-major order, where `R15=Lambda^2(C5)+Lambda^4(C5)` has the actual native
current attachment of the [charged-module proof](../field_identity_audit/PROOF.md).
The [sitewise source encoding](../chiral_source/PROOF.md) implements this
occupation algebra and its currents on literal flux/cloud carrier slots.
The antisymmetric spin factor and the fermionic grading are declared
architecture. No native binding energy or choice of a Higgs state follows
from the encoding. The dark native module has `K=0`.

## Regulator choice and complete finite domain

Choose an even side length N, physical torus length L and spacing `a=L/N`.
Use the antiperiodic spin structure in each spatial direction. Its complete
one-axis momentum menu is

`p_n=2 pi (n+1/2)/L`, `n=-N/2,...,N/2-1`.

The three-dimensional menu is its Cartesian cube. It contains N cubed
momenta, is closed under `p -> -p`, and has no zero mode. Every momentum,
spin and internal component is retained, giving `30 N^3` independent
fermion modes and Fock dimension `2^(30 N^3)`. There is no spectral
postselection and no deletion of an unwanted chirality from another source.
This finite menu defines a different kinetic law from the local sine,
Wilson and overlap sources.

In physical momentum coordinates the exact source symbol is

`h_L(p)=-sum_i sigma_i p_i tensor I15`.

Use the same continuous source time as this Hamiltonian. The Pauli identity
gives `h_L(p)^2=|p|^2 I`; there is one positive and one negative spin band
per native internal state. The positive excitation has helicity `-1/2`.
The antiperiodic spin structure lifts finite zero-mode degeneracy without
introducing a constant mass term. The first energy is `sqrt(3) pi/L`,
which vanishes as the physical box grows.

For a one-axis normalized Fourier matrix
`F_xn=N^(-1/2) exp[2 pi i(n+1/2)x/N]`, the position operator is
`P=F diag(p_n) F*`. Its entries are

`P_xx=0`,

`P_xy=-(i pi/L)(-1)^(x-y)/sin[pi(x-y)/N]`, for `x!=y`.

This follows by differentiating
`sum_n exp[2 pi i(n+1/2)r/N]=sin(pi r)/sin(pi r/N)` at nonzero integer r.
It is Hermitian. The cubic position Hamiltonian is the sum of these three
axis operators tensored with the same Pauli matrices. The receipt checks
the complete Fourier identities for N equal to 2, 4 and 6 and every column
of the side-two cubic operator. It exhibits a non-neighbor matrix element
`P_(0,3)=-i/2` at `N=6,L=2pi`.

Thus finite hopping has arbitrary range along coordinate axes. There is
no exact finite causal cone and no claimed uniform exponential locality.
The Fourier symbol is a truncated linear branch with a discontinuity at
the periodic momentum boundary. This is the explicit price of the spectral
regulator, rather than a counterexample to the local doubling obstruction.
The single physical low-momentum Weyl cone uses the entire retained menu;
all cutoff states remain part of every complete measurement.

## A positive reference on the same source

Let `P_±(p)=(I±h_L(p)/|p|)/2`. Fill the negative spin band at every momentum
and every internal row. The resulting finite Slater determinant is the
unique ground because the menu has no zero energy. It has `15V` occupied
source modes, with `V=N^3`, and energy

`E_sea=-15 sum_p |p|`.

Relative to this same reference, the complete occupation Hamiltonian is

`H-E_sea=sum_(p,A) |p| (b_(p,A)* b_(p,A)+d_(p,A)* d_(p,A)) >= 0`.

A particle adds a positive-band mode. A hole removes a negative-band mode;
its physical momentum is minus that removed mode's momentum. Both have
positive excitation energy. The full gap is `sqrt(3) pi/L`; the nonvacuum
relative-charge-zero threshold is twice that value. Neither is a nonzero
infinite-volume rest mass. The massless Weyl continuum has `E=|p|`.

Every traceless native subgroup current commutes with the flat free symbol.
The filled reference is invariant under their global action because each
complete internal multiplet is filled together. A hole transforms with
current `-T_A^T` and spin `-J_i^T`; at momentum `-p` its helicity is `+1/2`.
It is the antiparticle of a left-Weyl field, not an independent right-Weyl
field with the same native charges. The exact nonzero fifth hypercharge
trace changes sign under charge conjugation and separates those two cases.
The fifth trace is a representation diagnostic, not a four-dimensional
anomaly. The anomaly checks are the separate retained cubic and linear
current identities.

This is the free flat-link reference. It is not asserted to satisfy the
continuous local gauge Gauss law of a fluctuating quantum-link theory.
The interior source-code Gauss constraint and the full dynamical gauge
constraint are distinct. Replacing a quantum link by the identity is not
proved to be a physical invariant sector.

## Complete free records

At one momentum and internal row, the two spin modes have a four-dimensional
Fock space. Its reference has one negative-band occupation. Couple one
source mode to a fermionic reservoir, recording occupied/empty preparation
settings with probability one half and both final reservoir flags. The
complete Kraus set is

`a*/sqrt2`, `-n/sqrt2`, `a/sqrt2`, `(I-n)/sqrt2`.

The instrument is the finite reservoir operation of the
[mass-readout construction](../mass_instrument/PROOF.md), restricted to two
source modes. The code directly checks Kraus completeness and its full
four-state energy/charge effects. The joint detector has four nonzero
operator projectors: the vacuum of excitations, one particle, one hole and
a particle-hole pair. Combining them with four preparation records gives
sixteen records. Zero-probability records are retained.

For the actual retained momentum `p=(1/2,1/2,1/2)` at `L=2pi,N=2`,
`E=sqrt(3)/2` and `w=<e0,P_- e0>=(3+sqrt(3))/6`. The six populated original
weights are `(3-sqrt(3))/12,(2+sqrt(3))/12,1/12,(3+sqrt(3))/12,
(2-sqrt(3))/12,1/12`, up to detector ordering. Addition/removal success
gives energy `sqrt(3)/2` and relative charge plus/minus one. Failed
preparation has zero-energy and energy-`sqrt(3)` neutral components.
No probability is renormalized after a success flag. This ideal finite
apparatus does not establish scattering or an empirical energy standard.

## Controlled free-field refinement

Fix a physical torus and its antiperiodic one-particle Hilbert space. Let
`Pi_N` project onto the complete retained Fourier cube. The continuum
left-Weyl symbol and its negative-band covariance restrict exactly to the
finite ones. Both commute with `Pi_N`. Consequently, on every retained
momentum window, the free time evolution and two-point functions agree
exactly at every time. There is no kinetic or covariance approximation
inside that window.

For a smooth spinor f, the omitted modes satisfy
`|p|>=Lambda_N=pi(N+1)/L`. For every nonnegative integer s,

`||(I-Pi_N)f|| <= (1+Lambda_N^2)^(-s/2) ||f||_(H^s)`.

The anticommutation-relation norm is `||c(f)||=||f||`. It therefore bounds
the field truncation by that tail, uniformly in time. The same covariance
projection bounds two-point errors. Replacing each field in a product in
turn gives the sum of its tail times the norms of all other smearings.
This controls complete finite words, rather than only plane-wave energies.
The negative-band filled states are the restrictions of the same free
quasifree state, so no reference state is changed during the cutoff limit.

The limiting Weyl current obeys
`partial_t rho+div j=0`, `rho=psi*psi`, `j=-psi* sigma psi`.
The Pauli norm bound gives `|j.n|<=rho`, and integration outside an expanding
ball proves propagation speed at most one. Hence continuum fields
anticommute at spacelike separation and even observables commute. At finite
cutoff those anticommutators need not vanish. For spacelike-separated
smearings, before a torus wrapping can occur, their residual is bounded by
`||(I-Pi_N)f|| ||(I-Pi_N)g||`, since the retained and omitted projections
commute with the free evolution. Its vanishing follows from the displayed
tail bound, without claiming a finite source light cone.

For compact momentum smearings away from zero, increasing L gives ordinary
Riemann-sum convergence to the R3 covariance and dynamics. The norm bound
extends this convergence to all L2 smearings by density. This gives the
massless left-Weyl field with its positive particle/hole vacuum. Its
Lorentz covariance and helicity use the actual two-component spin/current
identities in the [Weyl-frame calculation](../../standard_model_source_map/KINETIC_ATTACHMENT.md).
Boosts act on fields and Cauchy surfaces together; they are not equal-time
unitaries on a fixed finite momentum menu.

## Gauge interaction and alternative regulators

A finite background covariant hop can multiply each long-range derivative
matrix element by a path transporter. Its endpoints transform with the
actual native representation, giving an exact background gauge-covariance
identity. The [source encoding](../chiral_source/PROOF.md) provides that
interface. This algebra does not bound the long Wilson paths, gauge
fluctuations, electric energy, or interaction-induced cutoff errors.
It supplies no nontrivial interacting gauge continuum.

This caution is intrinsic to the spectral choice. Drell, Weinstein and
Yankielowicz describe a lattice fermion construction retaining chiral
symmetry and the state count, with controlled lower-dimensional examples
([1976 paper](https://doi.org/10.1103/PhysRevD.14.1627)). Gauge-coupled
spectral derivatives need their own continuum analysis; finite free success
cannot supply that analysis. Rabin's
[SLAC gauge-theory study](https://www.slac.stanford.edu/pubs/slacpubs/2500/slac-pub-2676.pdf)
is a direct warning against transferring the free result to gauge dynamics.

The repository's [overlap route](../../common_source/chiral_candidate/REGULATOR_ROUTES.md)
and [spatial overlap proof](../../interacting_extension/chiral_lattice/PROOF.md)
provide alternative regulated interfaces. The first retains a gauge-dependent
Euclidean Weyl projector and an admissible-kernel locality bound. The second
retains vectorlike Hamiltonian matter and a modified chiral charge that
need not commute with electric motion, consistent with
[Creutz, Horvath and Neuberger](https://arxiv.org/abs/hep-lat/0110009).
Neither is silently substituted for the spectral source. Lüscher's
[nonperturbative abelian construction](https://arxiv.org/abs/hep-lat/9811032)
and [perturbative nonabelian construction](https://arxiv.org/abs/hep-lat/0006014)
have their respective hypotheses; anomaly-free finite currents alone do
not instantiate those constructions.

The finite theorem `Lean/Screen/WeylYukawaConventions.lean` fixes the
all-left charge-conjugation convention and allowed charge sums. The
spectral polynomial conservation theorem in `Lean/QFT/PeriodicMomentum.lean`
applies once the actual translation and clock commute, as they do in this
free source. Neither theorem selects a gauge vacuum, Higgs expectation,
Yukawa coefficient or GeV scale. Those inputs and a controlled interacting
chiral limit are outside the free regulator theorem.
