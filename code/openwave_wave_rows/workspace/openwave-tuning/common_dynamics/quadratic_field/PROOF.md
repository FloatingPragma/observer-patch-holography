# Fixed-coupling quadratic field on the actual carrier: bounded proposal

This package proposes `native-local12-dirac-noncompact-quadratic-v4-proposed`.
It replaces the **compact cosine field architecture** with a **noncompact,
quadratic field architecture**. It does not prove these architectures
equivalent. The native code services, actual carrier incidence and face maps,
matter routes, Gauss convention and operational phase-clock ancestry are
retained and pinned. No producer source, acceptance configuration or official
OpenWave row is changed.

The constructive result is a free quantum field on the actual complex, with
an exact charge-sector energy, shared-field current translations, a controlled
two-polarization photon frequency limit, and a bounded native free-photon
spectral-code/clock/readout implementation. A complete interacting matter
implementation and an official platform validation remain separate work.

## 1. Exact source and state space

Let `B` be the actual 96 by 288 incidence matrix and `A` the actual 288 by
208 face-boundary matrix from the frozen DYN08 size-two complex. Their integer
identities are

```
B A = 0, rank B = 95, rank A = 190,
ker B = im A direct-sum H, dim H = 3.
```

The frozen providers prove the ranks and the three integral homology
generators. This package rechecks the numerical projectors and builds its
rigorous finite spectrum certificate directly from the integer matrix
`K=A A^T`. A numerical eigensolver is a proposal mechanism for rational
integers, not a proof oracle.

The Gauss convention is exactly `n+B E=G`. Put `q=G-n`; periodic neutrality
`sum q=0` is required, with an explicit compensator if needed. For
`L=B B^T`, define

```
G_L = L^+ = (L + 11^T/96)^(-1) - 11^T/96,
P_q = B^T G_L q,
P_L = B^T G_L B,       P_T = I-P_L.
```

The physical field Hilbert space in a fixed matter occupation sector is
`L2(ker B)`, after quotienting longitudinal gauge coordinates. The electric
operator is `P_q+E_T`; its canonical coordinate is `theta_T`. In particular,
`[theta_e,E_f]=i(P_T)_ef`. This is a **real continuous field**, not an integer
electric tape relabeled as a continuous wave function. The embedding of the
matter sectors is the direct sum of these fibers, all using the same
transverse field.

The supplied, fixed microscopic coefficients are `a=1/2`, `b=3`. The source
Hamiltonian is the positive self-adjoint quadratic operator

```
H_q = (a/2)||P_q||^2 + (a/2)||E_T||^2
      + (b/2)||A^T theta_T||^2.
```

The coefficient choice is architecture. With the separately frozen actual
harmonic normalization `rho=2/3`, its canonical low-mode speed is
`sqrt(rho*a*b)=1` and its dimensionless charge coupling squared is
`rho*a=1/3`. This does not predict an observed fine-structure constant.

Let `U` be an orthonormal coexact frame diagonalizing `K`, with positive
eigenvalues `lambda_j`, and let `H` now denote an orthonormal frame of the
three-dimensional harmonic space. Write `omega_j=sqrt(a*b*lambda_j)`.
The specified normalizable field reference is the product of all 190
oscillator vacua and three **free Gaussian packets**. At time zero it has

```
mean E = P_q,                         mean theta = 0,
Cov(theta) = U diag(a/(2 omega)) U^T + H H^T/(4 v),
Cov(E)     = U diag(omega/(2 a)) U^T + v H H^T,
Cov_sym(theta,E) = 0,                 v > 0.
```

Thus `Cov(theta) Cov(E)=P_T/4`, the pure Gaussian minimum-uncertainty
identity on the physical canonical space. All 190 optical and long-wave
oscillators contribute; none is removed through a low-frequency label.
The packet choice in the finite receipt is `v=1/2`.

The three free modes have no normalizable stationary ground. Their electric
variance stays `v`; their coordinate variance is
`1/(4v)+a^2*t^2*v`, and their symmetrized cross covariance is `a*t*v`.
Their energy is `3*a*v/2`. They remain in the source and its evolution;
the native cap-one word below does not discretize them.

For each positive mode the exact equations are

```
theta(t) = cos(omega t) theta(0) + a sin(omega t) E(0)/omega,
E(t)     = cos(omega t) E(0) - omega sin(omega t) theta(0)/a.
```

For a harmonic pair, `theta(t)=theta(0)+a*t*E(0)` and `E(t)=E(0)`.
`source.propagate_mean` implements precisely these means, including the
harmonic drift. A mean trajectory is not by itself the quantum state;
the covariance, phase and matter amplitudes must also be retained.

## 2. Same-source charge energy and routed current

Since `A^T P_q=0` and `P_q` is orthogonal to `ker B`, every electric field
in the fiber has

```
||P_q+E_T||^2 = q^T G_L q + ||E_T||^2.
```

Therefore the charge-dependent part of the **quantum** reference energy
is exactly `(a/2) q^T G_L q`. The transverse oscillator zero-point energy
`sum omega_j/2` and retained free-packet energy are charge independent.
There is no integer-electric minimum being substituted for a quantum
vacuum. The finite static examples use opposite unit charges at vertices
0 and 1, 12, 36, or 84 and retain every oscillator and harmonic factor.

A matter move with actual integer routed current `J` satisfies
`n'=n-BJ`, hence `q'=q+BJ`. Its full electric update is `E'=E+J`.
Because

```
P_(q')-P_q = P_L J,
E_T'      = E_T + P_T J,
```

these are complementary pieces of the **same** current. On a transverse
fiber a momentum translation by `j=P_T J` is the Weyl operator
`T_j=exp(i j dot theta)`. Indeed `T_j^* E T_j=E+j`. If the branch order is
reversed, its overlap uses `T_j^*`; the executable characteristic function
is explicitly the latter convention,

```
<T_j^*> = exp(-j^T Cov(theta) j/2 - i j dot mean(theta)).
```

Changing this sign convention consistently has no effect on the real
time-zero interference examples. It must not be silently changed when
composing moving Gaussian branches or complex amplitudes.

At size two, the positive and negative DYN10 routes on an axis have the
same endpoint, but distinct integer currents. Their difference is closed.
Recombining their matter path register gives complete probabilities
`(1+Re overlap)/2,(1-Re overlap)/2`. The three positive outcomes are
approximately 0.897725, 0.881261, and 0.713056. Both outcomes are retained.
Independent-field substitution would give 0.5; replacing the field by a
deterministic character gives 1. These controls are separated. These are
specified logical interference measurements on the continuous source,
not yet native routed-matter measurement words.

The exact finite static energy supplies a force route: differences of
energies for actual charge configurations and the phase differences they
produce in the same clock. A controlled large-volume Green-function
asymptotic, preparation of moving charged states, and native force
measurements are not proved by the four finite energies. In particular,
this package does not claim a certified measured Coulomb coefficient or
a completed force row.

There is nevertheless an exact static normalization check on this same
graph. Let `C` be the 36 by 3 matrix placing a unit on each of the two
inter-cell edges in an axis, and `P_cycle=I-B0^+ B0`. For the actual
harmonic cycle frame `W` used below, the executable verifies
`C^T W=I`, `A0^T P_cycle C=0`, and
`C^T P_cycle C=(3/2)I`. Minimizing the cell potential corrector therefore
gives conductivity `3/2`, dual to the harmonic cycle Gram `2/3`.
The acoustic scalar Laplacian eigenvalue is `|k|^2/8+O(|k|^3)` on the
12-vertex normalized constant vector, so the root-to-root Green Fourier
pole is `(2/3)/|k|^2`. Multiplying by electric coefficient `a=1/2` gives
the same canonical charge-coupling squared `1/3` as the photon units.
This leading-symbol agreement is not a substitute for a quantitative
real-space remainder, a force instrument or an empirical strength test.

## 3. Actual-face photon frequency theorem at fixed coupling

`photons.py` uses the actual 36 by 26 Bloch face matrix `A(k)`, reconstructed
from the frozen oriented faces and their integer cell shifts. At zero
momentum `rank A(0)=22`. Its Moore-Penrose inverse is rational, with exact

```
||A(0)^+||_F^2 = 3383/240 < 16,
||A(0)|| <= ||A(0)||_F = sqrt(88) < 10.
```

Consequently its positive singular values are at least `1/4`. Set
`P=I-A0 A0^+` and `Q=I-A0^+ A0`. Write `A(k)=A0+E1+R`, where `E1` is
the exact linear Taylor term. Direct integer sums of the actual face
shifts give, for `r=|k|`,

```
||E|| <= 11 r = eta,        ||R|| <= (17/2) r^2.
```

The six coefficients of the quadratic matrix polynomial
`(P E1 Q)(P E1 Q)^*` are checked exactly. They equal the harmonic-space
form `W (|k|^2 I-k k^T) W^T`, where `W^T W=(2/3)I`.
The large left null space introduces no gauge singular values: the code
also verifies that `P E1 Q` lands in the three-dimensional harmonic
subspace. Thus its two positive singular values are exactly
`sqrt(2/3)|k|`; its other singular values vanish. This is an all-direction
polynomial identity, not an inference from sampled axes.

Here is the quantitative block elimination. Put

```
S = P E A0^+ - (P E A0^+)^*,
T = A0^+ E Q - (A0^+ E Q)^*.
```

These are skew-adjoint, with `||S||,||T|| <= 4 eta` because their two
off-diagonal supports are orthogonal. Expanding
`exp(-S)(A0+E)exp(-T)` to first order removes both high/low off-diagonal
blocks of `E`. The error from `S E+E T` is at most `8 eta^2`; the second
derivative integral remainder is at most `32 eta^2(10+eta)`. Replacing
the remaining low block `P E Q` by `P E1 Q` costs `(17/2)r^2`. Therefore
the whole singular-value perturbation is bounded by

```
D(r) = 8 eta^2 + 32 eta^2(10+eta) + (17/2)r^2,
eta = 11 r.
```

The high block has positive singular values at least `1/4-eta`.
The explicit admission inequalities are

```
eta < 1/8,
r+D(r) < 1/4-eta,
D(r) < r/2.
```

The actual incidence identity `B(k) A(k)=0` and rank `B(k)=12` for
sufficiently small nonzero `k` bound rank `A(k)` by 24. The block estimate
proves at least 22 separated high singular values and two positive low
ones. Hence there are exactly two acoustic positive singular values,
each within `D(r)` of `sqrt(2/3)r`. The remaining two singular values
are exactly zero. The incidence rank follows directly from connected
in-cell edges: a left null vector would be constant within the cell,
and a nonzero component of `k` forces that constant to vanish through
the actual inter-cell edges.

With `a*b=3/2`, the two **quantum oscillator** frequencies consequently
satisfy

```
|omega(k)-|k|| <= (5/4) D(|k|).
```

At cell spacing `h`, the declared common physical time is `t=n*tau0*h`;
the cell Hamiltonian generates physical frequency `omega(hp)/(tau0*h)`.
In units `tau0=1`, the frequency error is at most
`(5/4)D(h|p|)/h=O(h|p|^2)`, with fixed `a,b`, not a vanishing-coupling
limit. The receipt evaluates rigorous bounds at `h=2^-18,2^-20,2^-22`.
These are analytic bounds, not claimed dense native executions on a
millions-of-cells register.
The symbol theorem holds on real Bloch fibers. A finite periodic box
must additionally use its allowed discrete momenta `k=2*pi*m/L`; the
displayed `|p|=1` examples evaluate the uniform Bloch bound, rather than
claiming that momentum belongs to every finite periodic box.

For finite photon number at most `R`, bounded nonzero momentum, and
bounded time `T`, the excitation unitary error under the specified
spectral dictionary is at most `T*R*max frequency_error`. The reference
comparison **matches the same source vacuum scalar phase on both
sides**, and carries every other coexact vacuum factor and all three
free harmonic factors. It does not establish absolute continuum
vacuum-energy convergence. No scalar phase has been silently deleted
from the finite native implementation below.

The prepared one-photon Fock mode has unit vacuum-to-one-photon matrix
element for its normalized creation operator. This is nonzero quantum
spectral weight. A local-edge electric or magnetic pointer has the
additional actual mode function and oscillator normalization; native
implementation of that local observable is a separate encoder/readout
obligation. The present native Ramsey read is a spectral-code probe.

## 4. Finite Fock code, exact boundary and rational geometry certificate

For an oscillator occupation cutoff `N<=R`, use the truncated ladders
with their literal top-sector boundary:

```
[a_i,a_j^*] = delta_ij I
             - P_(N=R) (delta_ij I+a_j^* a_i).
```

The free Hamiltonian is declared as `H_cap=E0 I+dGamma(omega)`, where
`E0=sum omega/2`. Equivalently, it is the quadratic expression built
from the truncated canonical coordinates **plus**

```
(1/2) P_(N=R) (dGamma(omega)+sum omega).
```

This top-cap term is part of the finite architecture. It gives exact
free evolution of the vacuum and all retained photon states. It does
not restore the canonical commutator or guarantee Weyl-current accuracy.
For cap one, the commutator defect on a one-mode top state has norm 2,
which the executable checks. The code has 191 physical states (vacuum
plus all 190 one-photon states) and 65 retained padding states in eight
qubits. Padding has the declared Hamiltonian zero.

The geometry certificate begins with rational matrices `Q` and diagonal
`Lambda`, denominator `2^44`. Their integer numerators are selected from
floating eigensolver output and then checked through exact integer
products. Twenty-two-bit limbs keep every intermediate integer matrix
product within signed int64; recombination and residual multiplication
use arbitrary-precision integers. Explicit overflow guards are executed.
No floating residual is used as a proof bound.

Let `epsilon` bound `||Q^T Q-I||`, `r0` bound `||KQ-Q Lambda||`, and
`m=min Lambda`. Integer outward square roots of Frobenius squared norms
supply both bounds. Put `s=r0/m` and
`gamma=epsilon+2(1+epsilon)s+s^2`. Define

```
Y=K Q Lambda^(-1),          U0=Y(Y^T Y)^(-1/2).
```

`Y` lies exactly in `im K=im A`, `||Y-Q||<=s`, and
`||Y^T Y-I||<=gamma<1/2`. Thus `U0` is an exact orthonormal frame of
the full coexact space and `||U0-Q||<=s+gamma`. In particular, exact
Gauss support is defined by this polar frame, not claimed for rounded
floating eigenvectors.

Using the frozen all-refinement curvature bound `||K||<=1504/17`,

```
||U0^T K U0-Lambda|| <= eta_K
  = r0+(1504/17+max Lambda)(s+gamma).
```

The exact positive-frequency Fock space may use the frame `U0`; its
one-photon Hamiltonian is
`Omega=sqrt((3/2) U0^T K U0)`. Its vacuum is the actual Gaussian vacuum,
independent of a choice of orthonormal one-photon basis. The approximate
native code uses `Omega_hat=sqrt((3/2)Lambda)` in that **same** basis.
It does not identify rounded Euclidean coordinate basis vectors with
exact curvature eigenvectors. The positive matrix square-root
perturbation bound yields

```
||Omega-Omega_hat|| <= (5/4) eta_K/sqrt_lower(m) = delta_omega.
```

On cap one, the whole Hamiltonian error, including the vacuum scalar,
is at most `(190/2+1)delta_omega=96 delta_omega`. The executed rational
certificate is below `2.279e-5`. The vacuum-only part is not omitted.
`SPECTRAL_FRAME.json` retains the exact rational integer frame.

## 5. Native source, common clock and complete records

`native_photon.photon_phase` emits the actual M6 source compiler's
diagonal word for all 256 spectral-code basis states. Physical energies
are `E0,E0+omega_1,...,E0+omega_190`; padding energies are zero. All
energies use integer-outward square-root intervals. The Walsh coefficient
of every one of the 256 Pauli strings is checked by the existing native
template consumer, and every actual rounded control generator is charged
through the pinned source certificate. Scalar `I` pulses and zero-point
phases are retained.

The tested preparation is `(vacuum+one photon in frame mode 1)/sqrt(2)`,
implemented by a native Hadamard on bit zero. The free phase word is
followed by a native Hadamard analysis. The same operational clock uses
a separately prepared ninth bit, a native `Z` pulse of angle `-t/16`,
and configured X or Y analysis. Its relative phase is `t/8`, precisely
the frozen common-clock convention. The clock operations commute with
the free photon phase, so their placement in the bounded emitted word
does not change the joint experiment.

`recovery.execute_complete` executes the actual nine-bit word and
retains every first export-complement branch through its frozen literal
prefix matrices. Both native pointer laws C and E are executed. The
receipt stores all 512 coarse joint outcomes and every one-bit fine
effect; the complete nine-bit fine law is their tensor contraction with
the unnormalized successful source state. Export and reset failures
are additional retained outcomes. There is no success normalization.

Nine blank resets and nine fresh pointer-ready resets are charged by
the minimum success effect over all native carrier states. The frozen
`RESET_DOMAIN` rank-one Kraus proof permits this bound even for coherent
or correlated pre-reset carriers; no product-state independence is
assumed. Reusing fresh pointer auxiliaries after each completed read is
an explicitly reset apparatus operation.

The complete record total-variation bound is the sum of initialization
failure, native word/pointer bounds, and `t*96 delta_omega`. This is an
instrument bound, not a renormalized fit to the successful numerical
vector. Deleting the tested photon frequency changes the retained
photon bit by more than 0.05, well above its error bound. The finite
native packet does not implement the three harmonic factors, a local
edge-coordinate encoder, or a moving matter/current interaction.

## 6. Precise next obligations and interfaces

The shared source is now explicit enough for a separate interacting
extension. A finite matter route history can retain a tuple of matter
occupation/spin amplitude, `q=G-n`, complex phase, transverse Gaussian
mean and common Gaussian covariance. Each actual route uses `P_L J`
for the sector and `P_T J` for the same field displacement. Quadratic
free propagation preserves Gaussian form; distinct histories remain
coherent and must be combined through their full Gaussian overlaps.
The longitudinal energy contributes the sector-dependent phase.
No independent field copy may be assigned to each matter particle.

The following are remaining requirements, not conclusions of this package:

1. Select this named amendment through the shared source-contract process;
   the compact cosine proposal is not silently changed by this directory.
2. Compile the actual Dirac/Fock current-conditioned Weyl operations, with
   CAR signs, real Gauss embedding, common clock and retained amplitudes.
3. Bound photon occupation tails and the **intermediate** cap boundary
   under those operations. Exact cap-one free evolution cannot supply this.
4. Discretize or otherwise implement all three harmonic packets with their
   translation, drift and retained tails; or explicitly propose a different
   global sector choice as another architecture amendment.
5. Supply source-native local electric/current/density instruments and the
   mode-to-edge dictionary, with complete pointer and failure records.
6. Control branch growth, spatial/volume refinement, the preparation and
   observable transfer maps, clock/word-precision accumulation and the
   interacting state domain. The free frequency theorem supplies only one
   of those terms.
7. Derive and measure charged-particle dispersion, mass/vacuum stability,
   forces and empirical constants for the declared dynamics. A prescribed
   charge spectrum, Clifford frame or coupling is not an empirical derivation.
8. Run the official OpenWave criteria and obtain in-platform validation.
   No claim in this package promotes any of the 31 rows.

Stable construction interfaces are `source.geometry`, `source.reference`,
`source.propagate_mean`, `source.transverse_current`,
`source.gaussian_translation_overlap`, and `source.two_route_records`.
The actual-face bound is `photons.frequency_bound`; the finite exact
geometry certificate is `fock.spectral_certificate`; the native free
source and clock are `native_photon.photon_phase` and `native_photon.word`.
Numerical Gaussian reference APIs are explicitly distinct from the
integer-outward finite native certificate.
