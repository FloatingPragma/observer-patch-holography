# A covariant position measurement on the same dynamical response space

Use the source, selected state and continuum clock of
`spin_canonical_field`, and the source Gram isometry in
`spin_covariant_localization/CONTINUUM_POSITION.md`. Denote the original GNS
space by `H`, vacuum by `Omega`, spatial translations by `U(a)`, and time
translations by `D(t)=exp(-itL)`, with `L>=0`. The two groups commute. Write
`V:L2(R3,dnu)->H_Q` for the existing source polar isometry, where
`dnu=d^3k/(2pi)^3`. In particular

    U(a)V=V M_a,       (M_a f)(k)=exp(-ik.a)f(k),
    V*Q(F)Omega=sqrt(s_Q(k)) Fhat(k),
    m^2/|k| <= s_Q(k) <= 3/|k|.                              (1)

The state and generator are not changed by the construction below. The
position coordinate has the source envelope units, and `hbar` is the chosen
action unit in the momentum generator. No laboratory calibration follows
from choosing it.

## 1. Complete the field response under the actual clock

Define

    H_resp=closure span {D(t) V f : t real, f in L2(R3,dnu)}.   (2)

This space is separable: take rational times and a countable dense family
of `f`, using strong continuity. It reduces both spatial and time
translations. Indeed time evolution shifts `t`, and spatial evolution
commutes with `D(t)` and multiplies `f` as in (1); the inverse group elements
preserve the same span. Consequently the orthogonal complement also reduces
both groups. The restricted time generator is the same positive `L`.

Let `F_P(A)` be the joint spatial-momentum spectral projection on a set
`A` of Lebesgue measure zero. Equation (1) gives `F_P(A)V=0`. Since it
commutes with `D(t)`, it vanishes on (2). Thus the entire dynamical response
has absolutely continuous spatial spectrum. This follows from the source
field lower and upper bounds and clock commutation; a separate particle
shell is not an input. The vector `Omega` is orthogonal to (2), because it
is orthogonal to `H_Q` and is fixed by time evolution.

## 2. Construct the spectral multiplicity from actual time histories

The following construction avoids postulating a new Hilbert space for the
source. Enumerate the rationals `t_0=0,t_1,t_2,...` and set

    T_n=D(t_n)V.

Each `T_n` is an isometry, and `U(a)T_n=T_nM_a`. Therefore `T_n*T_l`
commutes with all `M_a`, and is multiplication by a measurable bounded
function `g_nl(k)`. These functions are determined by the source correlation
forms

    <T_n f,T_l h>=integral conjugate(f(k)) g_nl(k) h(k) dnu(k).

They obey `g_nn=1`, `g_ln=conjugate(g_nl)`, and every finite matrix
`[g_nl(k)]` is positive semidefinite almost everywhere. To obtain a common
exceptional null set, first use the integrated positive quadratic forms
on a countable family of rational complex coefficient vectors and their
finite supports, then use continuity in those coefficients. The pointwise
Cauchy--Schwarz bound is `|g_nl|<=1`.

At each `k`, give the countable formal vectors `v_n(k)` this Gram matrix,
quotient its null space, and complete. Apply Gram--Schmidt in the fixed
order, skipping zero residuals. All residual inner products are finite
measurable combinations of the `g_nl`; the sets where a residual is zero
are measurable. Choosing the next nonzero residual by the first index
keeps the resulting orthonormal frame measurable. Its first vector is
`v_0`, since `g_00=1`. Write `b_n(k)` for the coordinates of `v_n` in this
frame, padded with zeros in `ell2(N)`. Then

    <b_n(k),b_l(k)>=g_nl(k),       b_0(k)=e_0.              (3)

For finite sums define

    J sum_n T_n f_n = sum_n b_n(k) f_n(k)
        in L2(R3,dnu;ell2(N)).                             (4)

Equations (3) and the source Gram identity prove that (4) is well-defined
and isometric. It extends to all of `H_resp`, with closed range `R`. It
intertwines the existing spatial translations:

    J U(a)=M_a J,               J V f=f e_0.               (5)

This explicitly constructs the needed spatial spectral multiplicity out of
same-clock source histories. No particular energy measure, smooth fiber
frame, constant multiplicity, Gaussian state or dispersion is assumed.
The ordering of rational histories fixes a frame convention on the extra
channels. Different such conventions can change the position observable
on those channels. All of them retain the existing `H_Q` position specified
by the positive source Gram factor. Uniqueness of a physical detector does
not follow from (2)--(5).

## 3. A complete position POVM with a time-invariant detection sector

On the ambient vector-valued `L2` let `E_0(B)` be the ordinary Fourier
position projection, acting identically on every multiplicity channel.
Define on the actual response space

    E_resp(B)=J* E_0(B) J.                                  (6)

This is a normalized positive-operator-valued measure: positivity and
countable weak additivity pass through the bounded isometry, and
`E_resp(R3)=J*J=I_resp`. It need not be a projection-valued measure, since
the source range `R` can depend measurably on momentum and need not reduce
`E_0`. Equation (5) proves exact covariance

    U(a) E_resp(B) U(a)*=E_resp(B+a).                        (7)

It also proves agreement with the earlier field position as a sharp
measure on that subspace. Since every `f e_0` lies in `R` and `E_0`
acts componentwise,

    E_resp(B)V=V E_scalar(B).                              (8)

Thus `H_Q` reduces the position effects even though it need not reduce the
time dynamics. For any fixed finite measurable partition `{B_j}` of space,
choose the complete Lueders instrument with Kraus operators
`sqrt(E_resp(B_j))`. Functional calculus in this reducing subspace gives
`sqrt(E_resp(B_j))V=V E_scalar(B_j)`, so that specified finite-partition
instrument has the same restricted back-action. No instrument for an
ideal uncountable sharp-position collapse is asserted. Evolving a state out of `H_Q` does not leave the complete
response detection sector: every normalized `D(t)Vf` has total spatial
detection probability one under (6).

For the full source space, include its orthogonal complement as a distinct
outcome:

    E(B)=inclusion E_resp(B) projection,
    E_dark=I-Pi_resp.                                      (9)

The resulting measure is complete. Both translations and the actual time
clock preserve `E_dark`. Its probability is exactly time independent, and
it equals one on the vacuum. A normalized covariant position measure on
`R3` without this outcome cannot cover the invariant vacuum: it would
induce a translation-invariant probability measure on Euclidean space.
This obstruction is independent of any particle mass or field equation.

## 4. Uncertainty applies to these position outcomes

Let `P_j` be the existing spatial generator on `H_resp`, in action units
`hbar`. Strong intertwining gives, including spectral projections and
operator domains,

    J P_j = (hbar k_j) J.                                  (10)

For a normalized source vector `psi`, its position outcome measure in
(6) is exactly the ambient Fourier-position measure of `J psi`. If
`J psi` has finite Fourier-position variance and `psi` has finite momentum
variance, the canonical Weyl relation on the ambient vector-valued `L2`
space gives

    Delta_(E_resp) x_j  Delta_psi P_j >= hbar/2.            (11)

One proof puts the two first derivatives of the Weyl relation on separate
vectors, as in the preceding position construction; only the separate
position and momentum domains are needed. Cauchy--Schwarz applied to the
centered vectors gives (11). Equivalently one can approximate in the joint
first-derivative and multiplication graph norms by smooth cutoffs and
apply the usual Fourier inequality componentwise. This is an inequality
for the complete measured position distribution, not for the variance of
a different sharp observable on a replacement state.

If the compressed first moment `A_j=J*X_0,j J` is defined, its variance can
be smaller than the outcome variance. The difference is the nonnegative
measurement-noise term

    ||(I-JJ*) X_0,j J psi||^2.

It is therefore the outcome variance in (11) that is justified without
additional sharpness assumptions. No smoothness of the source fiber frame
has been proved, so finite position second moments and their persistence
under time evolution require their own bounds. In particular an invariant
response sector does not by itself supply a dense common finite-moment
core on every multiplicity channel. The initial vectors `V f`, for
Schwartz `f`, do have finite position and momentum moments, by (5).
For full-space states containing dark weight, (11) is a statement about
normalized detected-sector states; it is not a constant lower bound on an
arbitrary assignment of a real coordinate to the dark outcome.

## 5. Scientific and operational boundary

The response sector in (2) is invariant under the actual source time
dynamics, and (6) is a covariant position measurement using its actual
spatial momentum. This discharges the mathematical time-invariant response
space requirement for this observable. A response sector can contain many
energy branches or a continuous energy measure. No sharp particle pole,
Lorentz covariance, surviving scattering interaction, observed species or
mass follows from the construction.

The operators in (6) lie in the ambient bounded-operator algebra. The
source correlation Gram defines them, but does not by itself place them in
a specified observer's measurement algebra or make their whitening local.
An operational realization requires admitted complete source instruments
approximating these statistics on the same histories, with all dark
outcomes retained. To identify the coordinate with a literal carrier
record requires the actual record map to converge to these effects.
Second-moment tail bounds are needed to transfer (11) to finite detector
data. Knowledge of the source specification does not replace these
instrument and convergence statements.
