# Positive finite-amplitude entropy response of the natural Weyl shells

FA106 extends FA101 on the **exact U=0 source**, with the same state,
regional algebras, normalized exterior preparations and eight retained
lattice valleys. It computes a finite, strictly positive second variation
and controls the actual finite-amplitude relative entropy in the continuum
limit. It supplies no area dynamics or generalized-entropy stationarity.

## 1. The same source states, with their actual normalization

Let P be the full filled-band projection, E the ball restriction,
A=EPE on EH, and g a unit exterior vector, Eg=0. Put f=EPg and m=||f||^2.
The local normal states are

    omega_theta(X)=omega((1+theta[n_g-omega(n_g)])X), |theta|<1.

The exterior number operator is even and commutes with the regional CAR.
The multiplier is positive and the states are normalized. The covariance
is exactly A_theta=A-theta ff*, and the states are gauge-invariant
quasifree: Wick's determinant update has only one rank-one column. This
is a local affine family, not an affine global family through a pure
vacuum. For each finite lattice ball use the full P_a and

    g_a=M_b J_a g / sqrt(r_a),  r_a=||M_b J_a g||^2,
    f_a=E_a P_a M_b J_a g / sqrt(r_a).

No continuum valley is removed from P_a. The chosen preparation approaches
one of the eight summands. The exterior occupation is not r_a. Cell
supports must remain exterior; the separated FA91/101 shells satisfy this
for the mesh below. Each finite density dominates (1-|theta|)rho_0.

The projection block B=EP(1-E) obeys BB*=A(1-A). Thus

    f=[A(1-A)]^(1/2)h, ||h||<=1, m<=1/4.

For t>0 define q(t)=<f,(A+t)^(-1)f> and
r(t)=<f,(1-A+t)^(-1)f>. Their endpoint limits satisfy

    q(0)+r(0)<=1,  q(t),r(t)<=m/t.                       (1)

Indeed their sum at zero is the squared norm of the component of h on
0<A<1. Unused endpoint subspaces cause no problem. This controls both
endpoints without a uniform covariance gap.

## 2. Exact Hessian and finite-amplitude remainder

Write D(theta)=D(rho_theta||rho_0). At finite dimension the Fock density
formula gives the one-body relative entropy

    Tr[A_theta(log A_theta-log A)
       +(1-A_theta)(log(1-A_theta)-log(1-A))].

This follows by diagonalizing each covariance in the Fock logarithm;
A and A_theta need not commute. Differentiating the trace, using the
resolvent derivative of log and the rank-one inverse identity, gives

    D(0)=D'(0)=0,
    chi(theta)=D''(theta)
      = integral_0^infty {q(t)^2/(1-theta q(t))^2
                         +r(t)^2/(1+theta r(t))^2} dt.   (2)

The denominators have opposite signs because the particle covariance
changes by -theta ff* and its complement by +theta ff*. The tested
finite example has [A,ff*] != 0, so (2) is not a commuting shortcut.

For |theta|<=b<1, put s=1-b and chi0=chi(0). Then

    4m^2 <= chi0 <= m <= 1/4,
    chi0/(1+b)^2 <= chi(theta) <= chi0/s^2,
    |D(theta)-theta^2 chi0/2| <= chi0 |theta|^3/(3s^3),
    2m^2 theta^2 <= D(theta) <= chi0 theta^2/(2s^2).      (3)

For the upper bound on chi0, Cauchy-Schwarz gives
q(t)^2 <= q(0) integral lambda/(lambda+t)^2 dmu_f(lambda).
Integrating gives integral q^2 <= m q(0); apply the complementary
version and (1). For the lower bound let a=<f,Af>/m when m>0. Jensen
bounds q>=m/(a+t), r>=m/(1-a+t), hence chi0>=m^2/[a(1-a)]>=4m^2.
The same Jensen argument applied directly to A_theta proves
chi(theta)>=4m^2. If m=0 all quantities vanish. Differentiating (2)
and using q,r<=1 gives |chi'(theta)|<=2chi0/s^3. Taylor's integral
remainder proves (3). These bounds also prove smoothness on |theta|<1;
all differentiations are dominated uniformly on compact theta intervals.
In particular this is a nonzero entropy Hessian for every nonzero shell
restriction, not just an upper bound or a formal first law.

## 3. Identification with continuum Araki relative entropy

The continuum states are normal on the same regional von Neumann algebra:
their positive bounded exterior multipliers give vector representatives.
Choose increasing finite-dimensional regional one-particle spaces H_N
with dense union. The embedded finite CAR(H_N) are unital matrix von
Neumann subalgebras whose union generates the regional algebra.
Their restricted covariances are C_N=P_N C P_N and f_N=P_N f.
The natural-vector inequality survives compression:

    f_N f_N* <= P_N C(1-C)P_N <= C_N(1-C_N).

Consequently every bound (1)-(3) is uniform in N. Strong resolvent
convergence at t>0 gives q_N(t)->q(t), r_N(t)->r(t). For |theta|<=b,
the low interval (0,tau) contributes at most tau/s^2; the high interval
(M,infinity) contributes at most 1/(8Ms^2). Dominated convergence on
[tau,M] proves uniform convergence of chi_N(theta), hence of

    D_N(theta)=theta^2 integral_0^1 (1-u)chi_N(u theta)du.

To identify this limit we import precisely the martingale continuity of
relative entropy for increasing unital von Neumann subalgebras with
weakly dense union: Hiai, *Quantum f-divergences in von Neumann algebras I*,
Theorem 4.1(v), with f(x)=x log x
([primary proof, pp. 12-13](https://arxiv.org/pdf/1805.02050)). The proof
uses the variational formula on the weakly dense union. Equivalently,
Kaplansky approximation replaces its finitely many operator values by
uniformly bounded strong approximants from that union, which then fit
in one finite member. The original full-algebra values need not belong
to any finite member. The theorem's hypotheses are verified above. Therefore the limit D(theta) is the actual Araki relative entropy
of omega_theta relative to omega, and (2)-(3) hold for that quantity.
We do not subtract two divergent continuum entropies.

## 4. A quantitative joint lattice bound

Rescale the ball to unit radius, eta=a/R. Let ||g||H1<=B, B>=1,
eta<=1/(6B), with exterior cell supports. Use FA101's proved estimates
for T_eta=E_eta M_b J_eta and C=chi_B Pi chi_B. For j<=N,

    ||f_eta-T_eta f||<=20B^2 eta^(1/6),
    ||(T_eta*T_eta-I)f||<=7K_0 B eta^(1/72),
    ||(A_eta T_eta-T_eta C)C^j f||,
    ||(A_eta T_eta-T_eta C)(1-C)^j f||
        <=12K_N B eta^alpha_N,
    K_N=2048*8192^N*12^(N(N+1)/2),
    alpha_N=1/(6*12^(N+1)).                               (4)

The first includes the actual cell normalization, using
1-r_eta<=3eta^2 B^2 and r_eta>=11/12. The two power families in (4)
have the same fractional Sobolev loss as FA101's (1-2C)^j family:
each extra sharp restriction costs at most 2048/s_j and the enlarged
FA101 K_j bounds both C and 1-C. No H1 preservation at a sharp boundary
is assumed.

For 0<tau<=1 and t>=tau, expand (A+t)^(-1) in powers of (1-A)/(1+t),
and the complementary resolvent in powers of A/(1+t). Truncate before
power N. The tail norm is <=(1+t)^(-N)/t. Power telescoping uses
sum j/(1+t)^(j+1)=1/t^2; vector replacement, (4), and the T_eta
inner-product error give, for both scalar resolvents,

    |q_eta(t)-q(t)|, |r_eta(t)-r(t)|
       <=64B^2 K_N eta^alpha_N/tau^2 +1/(N tau^2).        (5)

Here (1+tau)^(-N)<=1/(N tau). This comparison takes place through T_eta,
not by subtracting operators on unrelated spaces. It retains the full
finite lattice covariance, including spectator valleys.

The derivative with respect to q of q^2/(1-theta q)^2 is bounded by
2/s^3 on 0<=q<=1. Split the integral at tau and M>=1 to obtain

    sup_|theta|<=b |chi_eta(theta)-chi(theta)|
      <=[2tau+1/(4M)]/s^2
        +(4M/s^3)[64B^2 K_N eta^alpha_N/tau^2
                         +1/(N tau^2)].                 (6)

For 0<epsilon<=1 choose

    tau=epsilon s^2/16, M=2/(epsilon s^2),
    N=ceil(32M/(epsilon s^3 tau^2)),
    eta<=min(1/(6B),
      [epsilon s^3 tau^2/(2048 M B^2 K_N)]^(6*12^(N+1))). (7)

The two endpoint errors total epsilon/4 and each of the two middle
errors is at most epsilon/8. Thus (6) is <=epsilon/2<epsilon.
The positive mesh is retained symbolically; the implementation never
materializes these powers or substitutes floating-point zero.
It follows that |D_eta(theta)-D(theta)|<=theta^2 epsilon/2,
uniformly on |theta|<=b. Combining this with FA101's first-derivative
convergence also controls the finite entropy changes

    S(rho_eta(theta))-S(rho_eta(0))
       -> theta ell - D(theta),  ell=-<f,k(C)f>.

This last expression is a finite entropy variation, not an absolute
continuum von Neumann entropy. The limits allow a common shrinking
R/L and a/R schedule; a fixed nonzero theta is permitted for this
absolute relative-entropy comparison. This does not assert an o(R^4)
nonlinear correction at fixed theta. The leading normalized matter
first-variation statement continues to use FA101's shrinking-amplitude
schedule; (8) below is a lower bound, not a sharp asymptotic.

## 5. A strictly positive response for the actual shells

For every unboosted directional FA91 shell at scale L, unit normalization
gives |Pi g(0)|^2=c_gamma/(pi^3 L^3),
c_gamma=4249388/150088295. Its first spatial derivatives have norm at
most 1024 L^(-5/2). Along the segment to |x|<=R the change is at most
2048(R/L)L^(-3/2). Since c_gamma>1/64 and pi<4, the central norm is
greater than L^(-3/2)/64. For rho=R/L<=2^(-18), the change is at most
half that norm. Integrating over the ball yields

    m >= c_gamma rho^3/(3pi^2) >= c_gamma rho^3/48,
    chi0 >= c_gamma^2 rho^6/576,
    D(theta) >= c_gamma^2 theta^2 rho^6/1152.              (8)

The nine original directional preparations all satisfy (8); we do not
claim this unboosted lower constant for an arbitrary transformed shell.
The same normalization argument as (4) gives
|m_eta-m|<=24B^2 K_0 eta^(1/72). Adding the positive mesh condition
eta<=[m_lower/(48B^2 K_0)]^72 makes m_eta>=m_lower/2, so the finite
lattice Hessian and relative entropy are strictly positive as well.
These are lower bounds, not a claimed exact rho^6 asymptotic.

## 6. Attachment to the existing Einstein stack

The existing G2/G3/G4 matter-modular inputs are used, not replaced.
This closes a further sector-relative G5 input: the same prepared source
states have a controlled nonzero entropy Hessian and finite-amplitude
relative-entropy limit, including a joint shrinking-ball schedule.
It neither supplies G1's responsive geometry nor the edge/area and
stationarity inputs in EinsteinRegisterData. For example, when the
geometric readout is held equal to the flat ball area 4pi R^2, its state
derivative is zero, whereas the established shell entropy first derivative
is positive. No choice of positive constant Newton coefficient makes
that fixed-area variation stationary. A compensating, source-derived
area response remains necessary. We do not impose it as a new law.

This package does not derive a nonzero-U entropy limit, an interacting
relativistic sector, sharply boosted lattice algebras, a native entropy
circuit, or a physical particle identification. The finite controls below
verify algebraic identities and noncommuting density comparisons; they do
not numerically execute the astronomical refinement schedule.
