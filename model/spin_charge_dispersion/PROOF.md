# Bounded source charges in a joint dispersive band

## 1. The source results consumed

Use the ordered ground of the spin-one-half cubic antiferromagnet, its
actual scaled positive generator L, and its same-history spatial translations.
The spin order obeys m>=m0>0. The
[canonical field construction](../spin_canonical_field/README.md),
[actual field energy bounds](../spin_canonical_energy/README.md), and
[source charge construction](../protected_defect_catalog/PROOF.md) supply
the following facts on one common bounded-history representation.

For a real unit-L2 test F with smooth Fourier support in
`r/2<=|k|<=r`, put

    x=sqrt(r) Q3(F),  y=sqrt(r) Q2(F),  p=P2(F)/sqrt(r),
    ell=L/r.

The configuration fields x,y strongly commute. Their circle action rotates
the pair, fixes Omega, and commutes with ell and spatial translations.
The original canonical Weyl group of p shifts x by m times its parameter
and fixes y; this last identity is exact before scaling because both
source y and p use spin component2. No second conjugate momentum field
is appended. The source state obeys

    ||x Omega||²<=6,   ||y Omega||²<=6,   ||p Omega||²<=1/4,
    integral E dmu_p(E)<=1/4,
    integral E^(-1) dmu_p(E)<=1/12.                       (1)

Here mu_p is the ell spectral measure of p Omega and includes all energies.
In particular p Omega has no zero-energy component. Its spatial spectral
support lies in the stated annulus. The source conjugation C sends
`x->x, y->-y`, fixes Omega and commutes with L and spatial translations.
The transverse source energy-kick bound gives
`q_ell(sin(u[x cos(theta)+y sin(theta)])Omega)<=3u²`.
All these are proved source inputs; neither Gaussianity nor a chosen
dispersion law enters them.

The external ordered-spin input used by these parents is precisely the
reflection-positive field-selected theorem of
[Koma, Theorem2.1 and equations2.11-2.16](https://arxiv.org/html/2105.04970v2).
Its excitation energy is a filtered mean. It does not identify a spectral
atom or prove a vanishing relative width. The argument below uses its
matched source inequalities through (1).

### A certified stronger order constant

The existing infrared producer in `spin_goldstone_sector/infrared.py`
bounds the same normalized three-dimensional integral I by exact rational
Darboux sums, with the origin cell retained and outward cosine/square-root
bounds. At256 divisions per direction its exact output satisfies
`I<0.356515121`. An outward enclosure of sqrt(24) therefore proves

    b_* =1/12-I_upper/sqrt(24)_lower >0.010559989>1/100.    (1a)

Here b_* is a lower bound for the liminf of the full-volume invariant
ground expectation b_N=<M²>, where M=N^(-1)sum epsilon_x Sx^1.
It is not a magnetization measurement or an inferred structure-measure
atom without volume control. The parent's finite-volume infrared estimate
and its singular Riemann-sum bound establish precisely this liminf.

The same variational ground argument admits any fixed rational bias t>0.
For a purification G of the invariant finite-volume ground projector,
symmetry gives <M>=<M³>=0 and ||M||<=1/2. The normalized vector

    psi_t=(I+tM)G/sqrt(1+t²b_N)

has norm1 and magnetization `2t b_N/(1+t²b_N)`. The exact source double
commutator bounds <M K M> by3J/N, so its excitation expectation is at
most `3J t²/[N(1+t²b_N)]`. The bias changes a trial state, not the source
Hamiltonian or the selected-state principle.

Project to energies `[0,J/sqrt(N)]` and normalize. The removed probability
is at most `3t²/sqrt(N)`, and the normalized-vector error tends to zero.
The projected vector has excitation-generator norm at most J/sqrt(N).
For every local operator A this makes its ground-state inequality defect
at most `(J/sqrt(N))||A||²`; every weak local limit is a ground state.
The invariant ground density, M and the spectral projection preserve the
source's even translations and odd-translation/axis-reversal symmetry.
Consequently the averaged magnetization equals the local staggered mean
and survives this limit. The actual field-selected state's maximal-order
theorem then bounds its magnetization by

    m >= 2t b_* /(1+t² b_*).

This expression is increasing in b_*. With t=10 and (1a) it is larger
than0.102723683, hence the exact selected source has

    m>m0=1/10,       kappa>1/200.                         (1b)

The second inequality consumes the inherited exact kappa>=m²/2 relation;
neither m nor kappa is assigned an approximate numerical value. As a
lower-cost independent control, the unchanged64-mesh integral and t=14
certify m>7/100. These are rational inequality certificates on the same
source. No empirical magnetization enters, and the biased finite vectors
are not identified as exact selected ground states.

## 2. Two bounded, exact-charge creation operators

Take the fixed scalar function

    f(t)=t/(1+t²)=integral_0^infinity exp(-u) sin(ut) du,
    B_+=(1/2pi) integral_0^(2pi) exp(-i theta)
             f(x cos(theta)+y sin(theta)) dtheta,
    B_-=B_+*.                                             (2)

These are bounded operators with norms at most1/2. The angular integrals
are strong integrals. Rotation shifts theta, proving
`R(phi) B_+ R(phi)*=exp(i phi)B_+` and the opposite relation for B_-.
They are exact charge±1 operators of the source configuration subnet.

Joint spectral calculus evaluates the same integral explicitly:

    B_+=(x-i y) h(x²+y²),
    h(v)=[1-(1+v)^(-1/2)]/v,  h(0)=1/2.                  (3)

For a scalar point `(x,y)=rho(cos(phi),sin(phi))`, rotate the angular
variable and use
`(1/2pi) integral (1+rho² cos²(theta))^(-1) dtheta
 =1/sqrt(1+rho²)`.
The sine part of the numerator integrates to zero. This proves (3),
including its continuous value at the origin.

The complete preparation has Kraus operators

    B_+, B_-, B_0=sqrt(I-B_+* B_+-B_-* B_-).              (4)

Both first operators have norm at most1/2, so the positive operator under
the square root is at least I/2. This construction also works before the
configuration fields commute: define the finite source B_+ by the angular
charge projection of f(x), and B_- by its adjoint. The same norm bounds
prove completeness. No branch is discarded or normalized in the complete
instrument. The isometry sum_q B_q psi tensor |q> preserves arbitrary
untouched references.

## 3. A positive original probability, with no amplitude-existence oracle

Let A=Re(B_+)=(B_++B_-)/2. The canonical Weyl shift of x and invariance
of y imply the weak commutator identity

    |< [p,A] >| = m < partial_x [x h(x²+y²)] >.

The scalar function in (3) has bounded continuous first derivatives.
Thus its difference quotient under the Weyl shift is uniformly bounded
and converges strongly by joint spectral calculus. Differentiating the
matrix coefficient on Omega in Dom p proves this weak identity without
postulating an unbounded product domain. Equivalently, the bounded
commutator criterion shows that A preserves Dom p.

Rotation invariance of the configuration spectral measure gives

    < partial_x [x h(R²)] >
      =<h(R²)+R² h'(R²)>
      =(1/2)<(1+R²)^(-3/2)>,      R²=x²+y².

The function `(1+v)^(-3/2)` is convex on v>=0. Since <R²><=12,
Jensen's inequality yields

    |< [p,A] >| >= d0=m0/(2*13^(3/2)),
    d0²=m0²/8788.                                        (5)

This lower bound pertains to the explicit bounded operators (2). It does
not rely on finding an unspecified sufficiently small sine amplitude.

Set

    a=3 d0²/2=3m0²/17576,
    b=2/d0²=17576/m0²,
    Pi=1_(r/2<=|P|<=r) 1_(a<=ell<=b).                    (6)

The first factor is the same operational field-translation spectral
measure as in the inherited position construction. By (1),

    ||(I-Pi)p Omega||² <= a/12+1/(4b)=d0²/4.             (7)

Since ||A||<=1/2, the part of the weak commutator (5) outside Pi has
magnitude at most d0/2. Its inside part has magnitude at most
`2||Pi p Omega|| ||Pi A Omega||<=||Pi A Omega||`.
Therefore `||Pi A Omega||>=d0/2`.

Pi commutes with the source circle action and conjugation. Its two
vectors Pi B_±Omega have orthogonal charges and equal norms. Consequently

    ||Pi B_+Omega||² = ||Pi B_-Omega||²
       =2||Pi A Omega||² >= m0²/17576.                   (8)

Equation (8) is the joint probability of a specified creation record and
an energy-momentum band in the original complete experiment. It is
uniform as r decreases to zero. The selected band's energy lies between
ar and br on the actual L clock, while its momentum lies between r/2
and r. The nonlinear B_±Omega need not have all their momentum in this
annulus; every momentum, energy and preparation complement is retained.
The positive lower bound is derived with that complement present.
Using the certified m0=1/10, the original weight is at least1/1757600,
and the energy interval is `[3r/1757600,1757600r]`. These deliberately
coarse band endpoints are not measured energies or a fitted dispersion.

The source energy form also bounds the full created state. From the
Laplace representation in (2), the source sine energy bound and the
triangle inequality in the energy-form norm,

    sqrt(q_ell(f(x cos(theta)+y sin(theta))Omega))
       <=sqrt(3) integral_0^infinity u exp(-u) du=sqrt(3).

The angular integral preserves this bound. Thus each B_±Omega has
unnormalized first L energy moment at most3r. This is a controlled
energy bound for the complete bounded preparation, without identifying
its mean with a rest mass.

## 4. Finite construction and the same-history limit

For any supplied finite source box and rational test coefficients, x is
a finite selfadjoint matrix. The rational function f(x) is computable
because I+x² has spectral lower bound1. Diagonalize the known total spin
charge, or use its complete native Fourier pointer. If P_q are its charge
projections, the exact finite operators are

    B_+=sum_q P_(q+1) f(x) P_q,  B_-=B_+*.                (9)

For n spins, all operator charge differences lie in [-n,n]. A discrete
Fourier average with power-of-two N>n+1 implements (9) exactly: no other
difference congruent to±1 occurs. Both (9) and the positive complement
in (4) are computable finite matrices. The existing
[finite instrument compiler](../spin_operational_position/PROOF.md)
therefore supplies a terminating source pulse program with its full
reference-safe error and resource bound. It consumes the declared coherent
source service; bare protected repair has no internal-spin records.

The finite certificate evaluates (9) on the actual selected eight-spin
H_B ground, using `x=sum epsilon_x Z_x/4`. It computes f(x) as exact
rationals, performs the charge projection through the actual integer
Walsh matrix, and checks `[Q,B_±]=±B_±` exactly. The cube's charge-reversing
automorphism exchanges the two matrices up to sign and preserves H_B.
Thus their energy measures are paired. All256 spin rows and three
preparation records are retained. The inherited ground-vector certificate
gives original probability error at most epsilon/2 because
`||B_±*B_±||<=1/4`. This finite example does not substitute its state or
selection field for the infinite ordered ground.

The continuum identification of (2) uses bounded source operations. First
truncate the Laplace integral; its operator tail is at most exp(-U).
Approximate the integral and circle average by finite quadratures on a
fixed compact parameter set. The inherited full-history amplitude and
test continuity make the errors tend to zero on every finite history
vector. The same common subsequence therefore carries these bounded
preparations and their adjoints into (2), without identifying an added
unbounded auxiliary vector as a field state. Completeness transfers
through the positive square root in (4).

This continuum quadrature statement is strong convergence on those
vectors, not uniform operator-norm continuity of an unbounded field.
At a fixed finite source, a Laplace quadrature with mesh Delta_u and exact
positive exponential interval weights has operator error at most
`||x|| Delta_u`, in addition to the tail. The resolvent representation
`f(x)=((x-iI)^(-1)+(x+iI)^(-1))/2` proves its operator-Lipschitz constant1.
An angular mesh Delta_theta therefore contributes at most
`Delta_theta(1/2+||x||+||y||)`. These finite norms are read from that source
and may grow under scaling. They give explicit finite compilation costs
and reference-safe errors without asserting a uniform continuum cost.

Pi in (6) is a joint spectral event. For a finite-time operational
approximation choose a smooth energy-momentum detector effect equal to1
on that closed inner band and supported in the wider band
`ar/2<E<2br`, `r/4<|P|<2r`. Its positive and complementary Kraus operators
are both retained. Smooth Fourier approximations use the same time and
translation histories; the finite instruments are covered by the same
native compiler. This effect has probability at least (8). A complete
unhalved trace error delta lowers that probability bound by at most
delta/2 before any conditioning. These are ordered qualitative state
and continuum limits. No effective thermodynamic ground-preparation
modulus or finite-source error to that ground is invented.

## 5. Reuse of the finite theorem stack and the shell boundary

The source premises in (1) are directly matched: one ordered ground,
one canonical pair, the added commuting configuration field, source
charge covariance, positive energy, same translations, and their exact
first/inverse moment bounds. This eliminates a separate assumed particle
wavefunction or a guessed energy operator for (8).

The following finite interfaces have different hypotheses and cannot
replace that source match.

| Existing theorem | Exact conclusion and premise relevant here |
| --- | --- |
| `Geometry/MassShellKinematics.lorentzQ_fourMomentum` | A declared mass times a Lorentz frame has that mass shell. The spin spectral measure is not identified with this declared vector. |
| `Screen/LightSignalAdequacySurface.physicalFrequency_eq_mode_branch` | A supplied oscillator lift and frequency identification give the carrier branch. The interacting H_J action is not that supplied lift. |
| `Screen/FieldSectorEnergyInnerProduct.mode_unitary_group_explicit_generator` | The declared linear carrier-mode flow has its exact generator. It does not linearize the source spin dynamics. |
| `QFT/SourceHistoryGNSDynamics.sourceHistoryHamiltonian_eq_changes` | An eight-history finite adapter fixes diagonal repair counts. Its state and clock differ from the ordered spin field. |

The [exact current identity](../spin_current_spectrum/README.md)
`E² nu_S(k,dE)=2e(k)nu_j(k,dE)` and inverse-current sum rule are stronger
matched dynamical inputs than those unrelated finite adapters. They
preserve all the weight used in (8), but do not force its frequency width
to vanish. The existing [shell criterion](../spin_shell_diagnostic/README.md)
identifies the required source-specific estimate: a common-speed filtered
current residual of squared norm o(p³), equivalently the stated force
residual o(p⁵), would concentrate the whole broad response. A distinct
reducing branch with independently bounded nonzero original weight could
instead carry a pole alongside a continuum. No arbitrary narrowing and
renormalization of an energy filter supplies that weight.

Equation (8) constructs bounded, measurable charged excitations with
energies proportional to their selected momentum scale. It does not
prove a sharp dispersion, one common Lorentz speed, a nonzero connected
scattering amplitude, a laboratory particle identity or a mass value.
