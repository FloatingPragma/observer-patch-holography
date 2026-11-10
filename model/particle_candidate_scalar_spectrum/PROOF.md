# A resolved positive radial-scalar transfer on the primary source

## Source and observable

Fix an authenticated member of the [primary contract](../particle_common_contract/PROOF.md).
Its physical Hilbert space, full interacting Hamiltonian H, selected polynomial
density rho, canonical quantum reference, and time parameter are unchanged.
Write a for its spacing, V for its site count, [L,U] for its authenticated
energy enclosure, D=U-L, and R=max(|L|,|U|). On the physical Hilbert space,

    rho >= c I,             c=2^(-N)>0.                         (1)

N is the certified density-floor exponent, not a capacity in nats.
The selected polynomial gives rho=p(H) with p decreasing. No Gibbs-state
replacement, compressed Hamiltonian, or filled-state preparation is used.
The [radial-scalar calculation](../particle_candidate_scalar_gauge_response/PROOF.md)
supplies the native one-site contraction

    A_x=(1/2)sum_b Phi_(x,b)* Phi_(x,b),
    Phi_b=(l_b0 e_1-l_b1 e_0)/sqrt(2),       0<=A_x<=I.       (2)

It preserves Gauss support and each native occupation. The actual full
six-mode radial spectrum is reused, including every occupation sector.
The apparatus has successful rows K_x=A_x/sqrt(V) and the retained failure
sqrt(I-sum_x A_x^2/V). All sites, failures, energy records and quantum
reference outputs belong to the complete experiment.

## Full-source commutator at every site

The normalized physical vector Omega fills negative-beta matter modes and
uses constant links. It is an operator witness, never a replacement for rho.
The radial calculation executes the native Wilson matrices and CAR signs
on a corner and proves

    ||[H,A_corner]Omega||^2=9/(16a^2).                       (3)

The same local identity applies to all sites. A site's z incident bonds
contribute mutually orthogonal one-link words. For each weak representation
of dimension2, a negative-mode hole changes A by1/4. For the singlet, it
changes A by1/2. Each representation contributes d/(2a^2) times the square
of this change. Thus one incident bond contributes

    (1/a^2)[(2/2)(1/4)^2+(1/2)(1/2)^2]=3/(16a^2).         (4)

For an incoming oriented bond the Wilson hopping matrix is adjointed.
Its beta-changing spin entries retain modulus1/(2a); Haar orthogonality
also retains 1/d for an adjointed link representation. Hence (4) includes
incoming and outgoing bonds. Different link coordinates, representations
and matter occupations give the same orthogonality as in(3). Every first
word lies inside the actual cutoff; projection and Gauss commute with A.
Onsite, electric and magnetic terms commute with A and contribute zero.

Every site of the actual open cubic source has z>=3, so

    ||[H,A_x]|| >= ||[H,A_x]Omega||
                 =sqrt(3z)/(4a) >= d0:=3/(4a).            (5)

For side length s, the degree3+j has multiplicity
binomial(3,j)2^(3-j)(s-2)^j. The executable checks both its V=s^3 total
and its degree sum against twice the actual source edge count. There is
no periodic boundary substitution.

## A frequency cut without a Hilbert-dimension factor

For lambda>0 define the real even integrable kernel

    f_lambda(t)=lambda sinc^2(lambda t/2)/(2pi)
               =(1-cos(lambda t))/(pi lambda t^2).        (6)

Its integral is1 and its Fourier transform, with frequency measured in
the unchanged H time, is (1-|omega|/lambda)_+. One verification integrates
the triangular function on [-lambda,lambda] in the inverse Fourier
formula; the value at t=0 follows by continuity.

For g(s)=(1-cos s)/(pi s^2), the same inverse transform on [-1,1] gives
|g'(s)|<=1/(6pi). Its integral on |s|<=1 is at most1/(3pi)<=1/9.
For |s|>=1 direct differentiation gives

    |g'(s)| <= 1/(pi s^2)+4/(pi |s|^3).

The two tails integrate to at most6/pi<=2. Consequently
||g'||_1<=19/9<3 and ||f_lambda'||_1<=3lambda.

Let alpha_t(A)=exp(itH)A exp(-itH), and define the analytical operator

    T=integral [2f_(2lambda)(t)-f_lambda(t)]alpha_t(A)dt.   (7)

Its multiplier on an energy difference omega is

    m(omega)=2(1-|omega|/(2lambda))_+
                 -(1-|omega|/lambda)_+.                 (8)

It equals1 on |omega|<=lambda, decreases linearly from1 to0 on
lambda<|omega|<2lambda, and vanishes on |omega|>=2lambda.
In particular 0<=m<=1. The kernel in(7) has L1 norm at most3 and
derivative L1 norm at most15lambda. Its boundary values vanish.
Integration by parts, with ||A||<=1, proves

    ||[H,T]|| <= 15lambda.                                (9)

The integrals exist in operator norm because the physical Hilbert space
is finite. They are proof devices; no infinite-time filter is executed.

Put B=A-T. Choose 0<lambda<=1/(40a). Equations(5) and(9) yield

    ||[H,B]|| >= d0-15lambda >= d0/2.

Subtracting the midpoint of [L,U] only in a norm estimate gives
||[H,B]||<=D||B||. This does not change the source energy zero or clock.
Thus ||B||>=(d0-15lambda)/D. In an energy basis B has multiplier
b=1-m, with 0<=b<=1 and b=0 on |omega|<=lambda. If mu_A is the original
two-energy success measure, then

    mu_A(|omega|>lambda) >= integral b(omega)^2 dmu_A
      =Tr(rho B*B) >= c Tr(B*B) >= c ||B||^2
      >= c (d0-15lambda)^2/D^2.                           (10)

All traces and norms are on the full physical source. The matrix B is
Hermitian, but the trace inequality also holds without that property.
Neither a count of energy levels nor an assumed small effective dimension
enters (10). The only population premise is the actual lower bound(1).

For E_j>E_i, Hermiticity gives |A_ji|^2=|A_ij|^2, while decreasing p gives
p_i>=p_j. Positive transfers therefore carry at least half of the paired
weight. Applied to(10), this proves at every site

    mu_(A_x)([lambda,D]) >= c q,
    q=d0^2/(8D^2)=9/(128a^2D^2).                          (11)

The sharper coefficient (d0-15lambda)^2/(2D^2) is also retained in the
receipt. Summing the physical menu rows K_x=A_x/sqrt(V) averages the
V bounds, leaving the same cq, with no factor1/V loss.

The default lambda=1/40 is admitted for every primary member, since
a<=1. It is a fixed positive source-clock cutoff. Its upper endpoint D
grows with capacity, and cq decreases with capacity. Equation(11)
therefore proves no surviving weight in a fixed compact energy band,
isolated pole, vacuum mass gap, relativistic scalar field or GeV mass.

## A large total radial response can be elastic

For each Hermitian A, its exactly zero-transfer measure contains all
energy-diagonal terms and hence obeys

    mu_A({0}) >= sum_i p_i |A_ii|^2 >= |Tr(rho A)|^2.     (12)

Degenerate energy blocks only add nonnegative contributions. Averaging
the site menu and applying Jensen once more gives

    mu_menu({0}) >= [(1/V)sum_x Tr(rho A_x)]^2 >= b^2,
    b=1-delta/(4V) >=59/400.                              (13)

Here delta is the authenticated weak-plus-singlet defect estimate reused
by the radial producer. Thus the same capacity-uniform bound3481/160000
certifies elastic weight. At member0 the exact b^2 exceeds7/10. These
lower bounds coexist with(11); total detection alone cannot separate
elastic response from a propagating particle contribution.

The inherited deleted weak-plus-singlet hopping comparison does not
replace c in(11): at member0 its source mean upper bound minus the
deleted-source lower bound is positive. It gives no strictly negative
hopping expectation from which to infer a macroscopic radial response.
That inequality is a limitation of this estimate, not a proof that
macroscopic inelastic response is absent.

## Complete finite clocks and exact resource accounting

Use the [complete joint spectroscopy](../particle_primary_joint_spectroscopy/PROOF.md)
energy clock before and after the radial instrument. Its J exact rows are

    F_j(H)=(1/J)sum_(r=0)^(J-1) exp(2pi i jr/J)exp(-ir tau H),
    tau=1/(4R),              eta=lambda/8.                (14)

The signed centers are 2pi(j if 2j<J else j-J)/(J tau).
The full source phase radius is at most1/4. The admitted tau eta<=1
gives no alias probability. With d=tau eta/7, the inherited complete
tail is at most min(1,1/[2(Jd)^2]+1/[2Jd]); Jd beta>=1 suffices to
bound it by beta<1. The r-substep Cayley compiler supplies one complete
clock error at most (J-1)/(768r^2), uniformly with a quantum reference.
These are approximations to the original H evolution, not a substitute
clock or dispersion.

Set beta=cq/16 and choose integer exponents

    k=max(0,ceil(log2(16/(dq)))),       b_clock=N+k,
    h=max(0,ceil(log2(1/(48q)))),
    s_clock=ceil((b_clock+N+h)/2),
    J=2^b_clock,                      r=2^s_clock.         (15)

The executable uses exact rational comparisons for the ceiling logarithms.
Then Jd beta>=1 and (J-1)/(768r^2)<=beta. The two tails and two complete
clock implementation errors consume at most4beta=cq/4. Outside those
tails a true transfer in [lambda,D] reports a difference in

    [3lambda/4,D+lambda/4],                               (16)

because the difference of two clock errors is at most2eta=lambda/4.
Equations(11) and(15) therefore prove that the original probability of
the reported event(16) is at least3cq/4. A record in(16) outside the tails
has true transfer at least lambda/2. No conditioned or renormalized
probability is substituted for this original event weight.

The complete apparatus retains (V+1)J^2 fine records, two b_clock-bit
energy registers, the address or failure, the source, and its reference.
The initial unread energy clock fixes the source marginal rho because
its rows are functions of H. It can change the joint source-reference
state; preservation of that joint state is not claimed. The intervening
probe changes the state according to its declared instrument.

These resources are explicit finite integers, represented by exponents.
They are enormous because beta contains2^(-N). The receipt does not
allocate the registers, execute their histogram or satisfy a compact
native-observer memory budget. Apparatus implementation by a narrower
native program is a separate requirement. The maximum duration of each
ideal clock is (2^b_clock-1)/(4R); it is recorded in factored form. The
source floor, degree counts, commutator and resource inequalities are
replayed for members0,1,2,8, with their exact typed contract identities.

## Physical scope

The result is a positive, resolved, original inelastic energy-transfer
event of the finite interacting source. It supplies no spin assignment
beyond the operator's finite scalar transformation law and no comparison
to the observed Higgs. A rest-mass identification requires a populated
branch in a controlled compact window, the corresponding momentum and
dispersion readout, a species dictionary, and finite-capacity errors.
The mathematical proof is finite-dimensional analysis with executable
source and arithmetic checks; it is not a Lean formalization of the
frequency-cut lemma.
