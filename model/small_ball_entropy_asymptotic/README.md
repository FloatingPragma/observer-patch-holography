# Fixed-amplitude small-ball entropy and its exact leading coefficient

FA107 continues FA106 on the exact U=0 filled-band source. Use all eight Weyl valleys;
one prepared normalized exterior mode occupies a selected summand, with
all spectators retained. The existing FA87/96 identification gives the
regional covariance C on the unit ball and its selfadjoint modular
one-body operator k=log((1-C)/C). We never assume a uniform covariance gap.

For a normalized compact shell g_1 supported d<=|y|<=D, d>0, let
 g_L(y)=L^(-3/2)g_1(y/L), rho=R/L,
 g_rho(y)=rho^(3/2)g_1(rho y), f_rho=1_B P_-g_rho.
Assume rho<=d/4. This includes the original nine shells (d=1,D=2) and
all81 observer-frame preparations (d=1/2,D=4). They remain preparations
on ordinary rest balls, not sharply boosted lattice region algebras.

## 1. An actual modular exponential domain estimate

Use [La Piana--Morsella, Theorem 3.2, equation 3.9](https://arxiv.org/html/2405.07888v2), on the ordinary signed-CAR
Cauchy space with the FA87 convention bridge already established. The
regional real flow is exp(-is k). At initial time zero its local formula
is the conformal spinor factor times the Weyl solution evaluated at

 tau_s(x)=cosh(pi s)^2-|x|^2 sinh(pi s)^2,
 t_s(x)=sinh(2pi s)(1-|x|^2)/(2 tau_s(x)),
 x_s(x)=x/tau_s(x),
 M_s(x)=tau_s(x)^(-2)[cosh(pi s)I+sinh(pi s)sigma.x].

The opposite chirality has the corresponding constant Pauli-frame
change; every estimate below is unchanged. At s=i a, |a|<=1/3,

 tau=cos(pi a)^2+|x|^2 sin(pi a)^2 >=1/4,
 |x_s|<=2, |Im t_s|<=2,
 ||M_s(x)||=tau^(-3/2)<=8.

The exact norm identity uses
(cos I-i sin sigma.x)(cos I+i sin sigma.x)=tau I.
In particular it is not a bound on an arbitrary complex spinor matrix.

The already derived spacelike negative-band Weyl kernel is
 K_-(t,z)=-i(t I+sigma.z)/(2pi^2(|z|^2-t^2)^2).
At t=i a and real z, its norm is exactly
 1/[2pi^2(|z|^2+a^2)^(3/2)].
We can use the weaker bound 1/[pi^2|z|^3]. Since y is in the shell,
 |x_s-y| >=d/rho-2 >=d/(2rho).
Thus the kernel is nonsingular throughout the vertical parameter segment.
Its denominator remains nonzero in a connected open complex neighborhood
of that compact segment, uniformly for x in the closed ball and y in the
compact shell. The kernel integral therefore defines an L2(B)-valued
holomorphic function psi(s), including at both endpoints s=+/-i/3.

For real s near zero, finite propagation makes the evolved truncated
initial spinor 1_B P_-g equal to the evolved full P_-g throughout the
interior diamond. Approximate the truncated L2 data by smooth data and
use the unitary real modular group to apply the imported local formula.
Therefore psi(s)=exp(-is k)f for those real s. The complex formula is a
local analytic continuation; no exponentially growing global Fourier
integral on the wrong energy half-plane is being used.

Let E_N=1_[-N,N](k). The identity theorem gives
 E_N psi(s)=exp(-is k)E_N f
on that connected neighborhood: both sides are holomorphic and agree on
the real interval. Hence at s=+/-i/3,
 ||exp(+/-k/3)E_N f||<=||psi(+/-i/3)||.
Spectral monotone convergence proves actual membership in both
Dom exp(k/3) and Dom exp(-k/3), and equality with the continued vectors.
This is the required operator-domain argument, not formal continuation.

By ||g_rho||_1<=sqrt(4pi D^3/3)rho^(-3/2), the preceding bounds give

 ||exp(+/-k/3)f_rho||^2
 <=65536 D^3/(9pi^2 d^6) rho^3
 <=K rho^3,                 K=8192 D^3/d^6.              (A)

The displayed K deliberately enlarges the bound. For all81 preparations,
K=2^25 is sufficient. Cauchy-Schwarz between the two exponential weights
also gives m=||f_rho||^2<=K rho^3.

## 2. An O(rho^6) bound at fixed amplitude

Since (1+e^x)^(2/3)<=1+e^(2x/3), (A) implies
 <f,C^(-2/3)f>, <f,(1-C)^(-2/3)f> <=2K rho^3.
For t>0, lambda+t>=lambda^(2/3)t^(1/3), and hence the same q,r as FA106
obey
 q(t),r(t)<=2K rho^3 t^(-1/3),
 q(t),r(t)<=m/t<=K rho^3/t.
Use the former on0<t<1 and latter on t>1 in FA106's exact Hessian:

 chi_rho(0)=integral_0^infty(q^2+r^2)dt <=26 K^2 rho^6.   (B)

For |theta|<=b<1, FA106 then gives the actual regional Araki entropy

 0<=D_rho(theta)<=13K^2 theta^2 rho^6/(1-b)^2.            (C)

Combining (B) with FA106's original-shell positive lower bound gives
 c_gamma^2 rho^6/576 <=chi_rho(0)<=26K^2 rho^6
for rho<=2^-18. Similarly D is bounded above and below by positive
constants times theta^2 rho^6. This is two-sided order, without yet a
claim about its exact coefficient. (C) is o(R^4) for fixed L as R->0;
for varying L one must require R^2/L^6->0. It is not an unconditional
statement for every joint L,R limit.

## 3. Controlled finite-amplitude lattice entropy response

Retain FA106's actual cell norm, full finite regional covariance and all
valleys. Choose its Hessian modulus and FA101's linear-entropy modulus
with epsilon=rho^6, B=6 for all81 shells, taking the minimum of the two
positive symbolic meshes and the support-separation mesh. No numerical
zero is substituted for this mesh. The source entropy identity is exact:

 S(rho_a(theta))-S(rho_a(0))=theta ell_a-D_a(theta).

The two proved errors are |ell_a-ell|<=rho^6 and
|D_a(theta)-D_rho(theta)|<=theta^2 rho^6/2. With
 N_R=n0 rho^4, n0=8c_gamma/(15pi),
and fixed nonzero |theta|<=b, the normalized finite entropy difference
approaches the same nine-direction response D_ij from FA101 with

 |[S(rho_a(theta))-S(rho_a(0))]/(theta N_R)-D_ij|
 <=[C_B+(1/n0)(1+b/2+13K^2 b/(1-b)^2)]rho^2.             (D)

This bound is uniform for0<|theta|<=b; at theta=0 use the derivative.
C_B is the already established FA91 Taylor constant. Thus amplitude no
longer has to shrink with rho. The same finite-matrix Neumann condition
9*(68/5)*C_total*rho^2<1/4 keeps the finite response matrix invertible
with inverse norm<=272/15. These are the same rest-ball algebras and
observer-frame preparations. No area, Newton, stationarity or Einstein
premise is manufactured.

## 4. Exact leading Hessian coefficient (additional derivation)

For a constant spinor u on the unit ball, define f0=1_B u. Its real
modular orbit inside B is M_s(x)u, since the constant Weyl solution is
constant throughout the diamond. Angular integration kills sigma.x:

 A_u(s)=<f0,exp(-is k)f0>
 =2pi|u|^2[cosh(t)/sinh(t)^2-t/sinh(t)^3], t=pi s,
 A_u(0)=4pi|u|^2/3.

Let p0(k)=1/[4cosh(k/2)^2]. The elementary beta integral gives
 integral exp(-isk)p0(k)dk=pi s/sinh(pi s).
Differentiating shows that the exact spectral measure of f0 is

 dmu_u(k)=|u|^2(k^2+pi^2)/[4pi cosh(k/2)^2] dk.            (E)

In the covariance coordinate lambda=1/(1+e^k), this is
 |u|^2/pi[(log((1-lambda)/lambda))^2+pi^2]d lambda.
The Hessian kernel is
 L(lambda(k),lambda(l))+L(1-lambda(k),1-lambda(l))
 =2(k-l)cosh(k/2)cosh(l/2)/sinh((k-l)/2),
where L(x,y)=(log x-log y)/(x-y), continuously at x=y.
Set h(x)=x/sinh(x/2), a(x)=(x^2+pi^2)/cosh(x/2).
Their Fourier transforms are
 hhat(s)=2pi^2 sech(pi s)^2,
 ahat(s)=4pi^3 sech(pi s)^3.
Parseval and integral sech(t)^8 dt=32/35 give

 chi[f0]=(1/(8pi^2))|u|^4 integral h(k-l)a(k)a(l)dk dl
        =64pi^4|u|^4/35.                                (F)

The separated-kernel Taylor expansion also holds at s=+/-i/3:
 rho^(-3/2)f_rho -> f0, u=(P_-g_1)(0),
in the graph norms of exp(+/-k/3). Indeed after source-coordinate change
 psi_rho(s,x)/rho^(3/2)=M_s(x)integral
 K_-(rho t_s, rho x_s-y)g_1(y)dy,
which converges uniformly on the compact imaginary segment and closed
ball to M_s(x)u; the kernel and its first derivatives are uniformly
bounded away from its singular set. The same argument gives an O(rho)
error for the fixed shell profile.

Weighted resolvent Cauchy-Schwarz gives pointwise convergence of
q_rho/rho^3 and r_rho/rho^3. Their squared integrands are dominated by
an integrable constant times min(t^-2/3,t^-2). Therefore

 chi_rho(0)/rho^6 ->64pi^4|u|^4/35,
 D_rho(theta)/rho^6 ->32pi^4 theta^2|u|^4/35,              (G)

uniformly for |theta|<=b<1. For the last conclusion note q_rho(t),r_rho(t)
->0 at every fixed t>0; the denominators in FA106 are bounded below by
1-b. Dominated convergence is uniform on the compact amplitude interval.
The finite lattice result has the same leading coefficient if its
Hessian approximation error is o(rho^6), e.g epsilon=rho^7, and the
corresponding existing positive symbolic mesh is used.

For each original normalized FA91 shell, |u|^2=c_gamma/pi^3. Thus
 chi_rho(0)~64c_gamma^2 rho^6/(35pi^2),
 D_rho(theta)~32c_gamma^2 theta^2 rho^6/(35pi^2)
for fixed nonzero theta. This computes the previously missing leading
small-ball nonlinear coefficient; it remains a free-matter result.

An independent Stieltjes calculation verifies the constant. For
L=log((1+t)/t), put r=t/(1+t)=exp(-L). Substituting
x=lambda/(1-lambda) and using the elementary beta integral gives

 I(a)=integral_0^1 [lambda/(1-lambda)]^a/(lambda+t)d lambda
     =integral_0^infty x^a[1/(x+r)-1/(x+1)]dx
     =pi(1-exp(-aL))/sin(pi a),                 |a|<1.

The formula first follows directly for -1<a<0 and extends analytically
across a=0 by dominated convergence. Therefore I(0)=L and
I''(0)=(L^3+pi^2 L)/3. Equation(E) yields

 q_u(t)=|u|^2[L^3+4pi^2L]/(3pi),  r_u(t)=q_u(t).

Since t=1/(exp L-1) and
integral_0^infty L^n exp(L)/(exp(L)-1)^2 dL=n! zeta(n) for n>1,
expanding q_u^2 gives precisely64pi^4|u|^4/35. The zeta identity follows
by integrating the positive series sum_{j>=1}j exp(-jL) termwise;
zeta(2)=pi^2/6, zeta(4)=pi^4/90, zeta(6)=pi^6/945.

For completeness the beta Fourier transform used above is
integral exp(-isk)/(4cosh(k/2)^2)dk
=Gamma(1+is)Gamma(1-is)=pi s/sinh(pi s), using x=exp(k).
Similarly integral exp(-isk)/cosh(k/2)dk=2pi sech(pi s).
Differentiation and Fourier inversion give both transforms used in(F).
All kernels in Parseval are integrable, so no distributional subtraction
or principal-value prescription is required.

The limiting constant spinor f0 is a modular analytic vector used to
compute the Hessian coefficient. It is not asserted to generate a
finite-amplitude positive rank-one family by itself: its endpoint
inverse-covariance norms can diverge. The actual finite-amplitude states
are always the original natural exterior-shell states.

## 5. Attachment to the existing Einstein stack

No change to the already proved G1-G6 implication or tensor completion.
The result strengthens only the source-realized G5 matter estimates.
It supplies no responsive area, central edge identification, Newton
normalization, generalized entropy stationarity, nonzero-U entropy law,
sharp boosted lattice algebra, measured particle species or complete
interacting Einstein/QFT sector.


## Validation and reproducibility

Run `python -m pytest -q model/small_ball_entropy_asymptotic` and
`python -m model.small_ball_entropy_asymptotic --write`, followed by a
fresh `--check`. Exact checks cover the imaginary spinor prefactor,
Weyl-kernel norm, beta/Stieltjes expansion, real-orbit characteristic
function, both exact Hessian calculations, rational error budgets and
symbolic refinement. Independent floating integrals check the same
spectral density and a two-variable Hessian without supplying the proof.
Their raw floating tails are not serialized as exact certificate data.
A one-mode natural-vector counterexample demonstrates why the new
modular exponential estimate is essential. The source manifest pins
current mathematical inputs narrowly; it claims no full-repository or
old-receipt replay.
