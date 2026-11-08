# Relaxed spin response from the source ground

The ordered-spin continuum has a positive, finite energy cost for a prescribed
linear response in either canonical field. That cost is determined by the
ground's inverse-energy susceptibility. A regularized source resolvent constructs
energy-minimizing response vectors and gives a complete finite readout. The
construction uses the same spin interaction, reference state and clock as the
[canonical continuum](../spin_canonical_field/README.md). It needs no equality
between a microscopic twist energy and the energy retained in that continuum.

This is a constrained tangent response. It does not identify a helicity modulus,
propagation speed, sharp particle shell or interacting relativistic theory.

## The actual continuum vectors

The source is the cubic spin-one-half antiferromagnet

    H_J = J sum_<xy> (Sx.Sy-1/4), J>0,

in the zero-field ground selected by an infinitesimal staggered field along
internal axis one. The scaling is X=p*x and T=J*p*t. Write K for its positive
continuum time generator, Omega for its ground, and Q(F), P(G) for the actual
Stone fields of the bounded-history construction. The
[energy theorem](../spin_canonical_energy/README.md) places Omega in their
domains, proves that their vectors have no zero-energy component, and supplies
their finite first and inverse energy moments. These are the vectors in the
same continuum Hilbert space, rather than unprojected auxiliary limits.
Q and P are canonical field quadratures. They are not a particle's position
and momentum, and their canonical relation does not select particle masses.

Take a real test with unit L2 norm and Fourier support in
r/2<=|k|<=r. For either A=P(G) or A=Q(F), put

    a=A Omega,       n=||a||^2,
    h=||K^(1/2)a||^2,       chi=||K^(-1/2)a||^2.

Let 0<m0<=m, where m is the source's staggered order, and let
m^2/2<=kappa<=1/2 be its exact transverse bond coefficient. The inherited
source inequalities are

    P: n >= m0^2 r/24,   h <= kappa r^2/2,   chi <= 1/12;
    Q: n >= m0^2/r,      h <= 6 kappa,       chi <= 4/r^2.

Cauchy--Schwarz in the positive energy measure gives n^2<=h chi.
Consequently the same actual fields obey

    m0^4/(288 kappa) <= chi_P <= 1/12,
    m0^4/(6 kappa r^2) <= chi_Q <= 4/r^2.                 (1)

The positive lower bounds follow from the nonzero canonical commutator and
the source moment bounds. They do not assume a Gaussian state or saturation
of those bounds. They are uniform in the admitted envelopes and every selected
continuum subsequence. The code can replace kappa by its upper bound one half;
it does not supply its unknown numerical ground-state value.

## The exact energy-minimizing response

Define the tangent response of a finite-energy vector v perpendicular to
ker K by the real functional b(v)=2 Re<a,v>. For a prescribed real b,

    inf { ||K^(1/2)v||^2 : 2 Re<a,v>=b }
        = b^2/(4 chi).                                    (2)

This definition requires only a=A Omega. It does not assert that v lies in
the domain of every unbounded field product. On vectors where the expectation
is differentiable it is its usual first derivative along
(Omega+epsilon v)/||Omega+epsilon v||. Bounded sine functions of the Stone
field approximate that derivative on Omega.

For the lower bound in (2), write
<a,v>=<K^(-1/2)a,K^(1/2)v> and apply Cauchy--Schwarz. For eta>0 set

    chi_eta=<a,(K+eta)^(-1)a>,
    psi_eta=<a,(K+eta)^(-2)a>,
    v_eta=b (K+eta)^(-1)a/(2 chi_eta).

This vector has precisely the requested tangent response and finite norm.
Its energy is exactly

    q(v_eta)=b^2 (chi_eta-eta psi_eta)/(4 chi_eta^2).       (3)

Monotone convergence gives chi_eta->chi as eta decreases to zero. Moreover
eta/(E+eta)^2<=1/E and this integrand tends to zero for E>0, so dominated
convergence gives eta psi_eta->0. Equations (2)--(3) follow. For b!=0 the
unregularized infimum is attained precisely when K^(-1)a is a Hilbert-space
vector. No inverse-square moment is supplied by the inherited theorem, so
attainment is not presumed. With the positive form K+eta the minimizer
v_eta is attained for every eta>0.

Combining (1) and (2) gives finite, nonzero response-cost intervals:

    P: 3 b^2 <= inf q(v) <= 72 kappa b^2/m0^4;
    Q: b^2 r^2/16 <= inf q(v) <= 3 kappa b^2 r^2/(2m0^4). (4)

These inequalities identify the response's scaling in the two canonical
channels. They do not identify the full Hamiltonian with a quadratic free-field
Hamiltonian. The state can have other fields and spectral continua.

## An exact residual certificate

The regularized inverse response also has the dual characterization

    chi_eta = sup_v { 2 Re<a,v>-<v,(K+eta)v> }.

For any trial v for which the residual r=a-(K+eta)v is defined, completing
the positive quadratic form gives

    chi_eta-[2 Re<a,v>-<v,(K+eta)v>]
       = <r,(K+eta)^(-1)r> <= ||r||^2/eta.                (5)

Thus an approximate finite linear solve gives a two-sided certificate from
its actual residual, without a spectral-gap premise or an eigensolver's
convergence declaration. Complex trial vectors satisfy the same identity.
For positive spectral mass n and first moment h,

    chi_eta >= n^2/(h+eta n),
    psi_eta >= chi_eta^2/n.                               (6)

These follow from Cauchy--Schwarz with E+eta and then with the constant
function. They give a strictly positive regularized response and filtered
weight without postselecting an artificially narrow energy interval.

For a finite source whose actual positive gap is certified to exceed g>0,

    chi_eta <= chi <= (1+eta/g) chi_eta.                   (7)

The finite certificate uses the sharper identity

    chi=chi_eta+eta psi_eta
        +eta^2 integral dnu(E)/[E(E+eta)^2],

whose last term is between zero and eta^2 psi_eta/g. This is the exact
dual residual (5) applied on the positive-energy subspace with its gap.
The continuum has gapless spectrum; (7) is used only for the finite
positive-field source. There is no insertion of a finite-box gap into (1).

## Reading the actual continuum coefficient through bounded source probes

The susceptibility in (1) belongs to the actual continuum Stone vector.
It need not equal a limit formed directly from the unbounded source vectors
A_p Omega: that inference would ignore the projection issue in the inherited
energy construction. The following ordered limit identifies an accessible
response without making that inference.

Use the same selected bounded-history subsequence p_j and its actual fields.
For a fixed real s!=0, set

    A_(p,s)=sin(s A_p)/s,
    A_s=sin(s A)/s,
    L_p=K_source/(Jp).

These operators have norm at most 1/|s|. Their two-point time histories are
finite linear combinations of the same Weyl histories used to construct
the continuum. At fixed eta>0 the reciprocal on the positive spectrum has
the integrable Fourier representation w_eta(t)/eta from (8). Dominated
convergence therefore gives

    lim_j <A_(p_j,s)Omega,(L_(p_j)+eta)^(-1)A_(p_j,s)Omega>
       = <A_s Omega,(K+eta)^(-1)A_s Omega>.                (7a)

The dominating integrand is 1/(eta s^2) times a probability density; all
times use the same source scaling. Positivity of the source and continuum
generators is essential to replace |E| by E in this formula.

Since Omega belongs to Dom A, the spectral theorem and
|sin(sx)/s|<=|x| imply A_s Omega->A Omega in norm as s tends to zero.
Writing d_s=||(A_s-A)Omega||, the right side of (7a) differs from chi_eta
by at most 2||A Omega|| d_s/eta. Finally eta decreases to zero by monotone
convergence. Hence the actual relaxed susceptibility has the source formula

    chi = lim_(eta downarrow0) lim_(s->0) lim_j
       <A_(p_j,s)Omega,(L_(p_j)+eta)^(-1)A_(p_j,s)Omega>.   (7b)

Every innermost experiment is bounded. No exchange of these three limits
or first-energy uniform integrability is required. This uses the original
bounded-history representation throughout, so no discarded auxiliary vector
is declared to be the actual field vector.

The complete two-outcome probe has Kraus maps sin(s A_p) and cos(s A_p).
They are linear combinations of the admitted opposite-angle Weyl controls,
and their effects sum to the identity on the entire source, including
untouched references. The first branch prepares s A_(p,s)Omega; the second
branch is retained. Combining it with (10), the unnormalized score is
eta s^2 times the quantity inside (7b). An absolute score error epsilon
therefore gives inverse-response error epsilon/(eta s^2). This is division
by fixed declared probe parameters, without dividing by observed success.

There is no effective convergence modulus in j, s or the finite-box/
selection-field realization of the specified ordered state. Formula (7b)
is a proved common-source identification with an ordered limit; it is not
a finite runtime bound for a prescribed continuum susceptibility accuracy.
It also provides no shell concentration from susceptibility alone.

## A source-only complete filter and susceptibility readout

For a bounded finite source probe A with ||A||<=Lambda and the ground of its
own Hamiltonian, the successful target operator is a Heisenberg average.
For native implementation require Lambda to normalize an admitted complete
controlled block encoding or linear combination of unitaries for A.
An operator-norm bound alone supplies no such observer operation. The finite
Z probes below have explicit Pauli sums with normalization one; the bounded
sine probes use the specified controlled Weyl operations.

    B_eta = Lambda^(-1) integral w_eta(t) alpha_t(A) dt,
    w_eta(t)=integral_0^infinity eta exp(-eta s)
                        s/[pi(s^2+t^2)] ds.

This positive integrable density has total mass one. The Cauchy characteristic
function and the elementary Laplace integral give

    integral w_eta(t) exp(itE) dt = eta/(|E|+eta),
    B_eta Omega = eta (K+eta)^(-1)a/Lambda.                (8)

Here finite-source K=H-E0 and eta use that source's energy units. To implement
a fixed continuum regularizer eta_c for K/(Jp), use eta=Jp eta_c in (8).
This is the same clock rescaling, without a separately fitted energy scale.

In particular ||B_eta||<=1, its success weight is eta^2 psi_eta/Lambda^2,
and every complementary outcome has the remaining total weight. An additive
scalar in the source Hamiltonian cancels inside alpha_t(A); the filter needs
no absolute ground-energy phase. This does not remove the energy-origin
requirements of a subsequent independently implemented energy detector.

There is an explicit finite approximation. Truncate |t|<=T, partition that
interval into cells of width at most delta, and replace alpha_t(A) by its
midpoint value. Since

    integral_(|t|>T) w_eta(t) dt <= 2/(pi eta T),

the normalized operator error is at most

    2/(pi eta T) + delta ||[H,A]||/(2 Lambda)
                  + total coefficient error.            (9)

The tail follows by integrating the Cauchy tail bound 2s/(pi T) against
eta exp(-eta s). The midpoint estimate follows by differentiating the
Heisenberg operation. The bin weights are computable integrals: truncate
s in [epsilon,S], bound the omitted pieces by eta epsilon+exp(-eta S),
then integrate the smooth arctangent expression on that compact interval
with certified elementary-function intervals. Choose nonnegative rational
lower approximants to the weights. Their sum is at most one; two opposite-phase
identity terms fill the missing normalization. Coherent preparation, SELECT
and unpreparation implement this finite sum using the same services as
[the complete spin packet readout](../spin_packet_readout/README.md).
Every nonblank label is retained. Source evolution, control precision,
preparation and detector errors require their stated separate budgets.

The regularized susceptibility itself has a bounded complete score. First
apply the complete probe operation with successful vector a/Lambda. After
its successful branch, use the source's complete excitation-energy readout
and score each excitation by eta/(E+eta); score every probe or detector
failure zero. The expected unnormalized score is eta chi_eta/Lambda^2.
Use the bounded score eta/(max(E_label,0)+eta) on numerical energy labels.
Require a certified joint spectral-label law whose probability of
|E_label-E|>d_E is at most ell, uniformly on the source spectrum. For
unhalved complete output error d_out and joint input trace error d_in,
the score error is bounded by

    d_E/eta + ell + d_out + d_in.                         (10)

On the near labels use the score's Lipschitz constant 1/eta; on the far
labels its change is at most one. Thus a finite Fejer clock needs its
far-bin leakage and alias bounds in ell. A bin spacing or nominal resolution
alone does not provide this premise. A hard spectral interval service can
use ell=0 only if it proves that support property. The score is in [0,1];
no division by a success probability enters
the expectation. Independent repetitions need their sampling error added.

## Certified finite-source response

The executable consumes the complete 256-vector integer witness for the
unchanged open eight-site source, J=1 and staggered B=1/4. Its true unique
ground, complete polar spectral residual, gap and vector error are verified
before any response is evaluated. Global spin flip commutes with this source
Hamiltonian, and the three Z probes are odd under that flip, so their exact
ground means vanish. The positive selection field and its torque are retained;
this state is not substituted for the infinite zero-field ordered ground.

All 256 spectral weights, including exact or approximate zeros, enter the
regularized susceptibility and squared-resolvent response. Resolvent identities
propagate the full operator certificate. If the excitation-operator error is d,
the stored-vector polar error is g and the ground-vector error is e, then for
the norm-one probes the susceptibility error is at most

    d/eta^2 + [g(2+g)+2e]/eta.

The squared-resolvent and energy-numerator errors are respectively bounded by
2d/eta^3+[g(2+g)+2e]/eta^2 and
3d/eta^2+[g(2+g)+2e]/eta. Equation (3) and the sharpened form of (7) produce certified
response-cost intervals. The receipt retains the complete successful filter
weights and complementary rejections. It reports a finite target calculation;
the large native resolver tape is not claimed executed.

The construction supplies no effective box-size/selection-field convergence
modulus for the specified infinite ground. It does not evaluate the four-spin
inverse-current tail in [the exact current sum rule](../spin_current_spectrum/README.md),
prove its vanishing, or establish the sharp subbranch condition in
[the source shell test](../spin_shell_diagnostic/README.md). These are separate
questions from the positive relaxed response in (1)--(4).

```sh
python -W error -m pytest -q model/spin_relaxed_response
python -W error -m model.spin_relaxed_response --write
python -W error -m model.spin_relaxed_response --check
```
