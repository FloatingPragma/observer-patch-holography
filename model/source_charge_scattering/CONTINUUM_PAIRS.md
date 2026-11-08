# Separated source charge configurations

The [bounded source creator](../spin_charge_dispersion/PROOF.md) supplies
normal operators B_+=B and B_-=B* in the joint configuration-field algebra,
with charges±1 and norms at most1/2. They use one ordered source ground Ω,
positive energy L and commuting spatial translations U(x). At a momentum
scale r>0, the complete three-outcome instrument is

    B_+, B_-, B_0=sqrt(I-B_+*B_+-B_-*B_-).

The source order bound m>1/10 follows from the replayed infrared/biased
state certificate in `spin_charge_dispersion/order.py`. Its bounded
creator theorem gives each sign an original probability at least

    p_*=m0²/17576=1/1757600                                (1)

in the specified positive joint energy-momentum band. In particular the
full unfiltered creation norm squared is at least p_*. This is an
established input on the same source. The probability p_* is distinct
from the microscopic current coefficient denoted kappa elsewhere.

## 1. Every record in a finite sequence of creations

Fix an integer n>=1 and positions x_1,...,x_n. The complete preparation
has all3^n Kraus operators

    K_s=alpha_(x_n)(B_(s_n)) ... alpha_(x_1)(B_(s_1)),
    s_j in{-1,0,+1}.                                     (2)

Their adjoint products sum to I by iterating the complete three-branch
instrument. This identity does not require the operations to commute.
Every classical sequence and every native implementation failure is
retained. There is no conditioning on a selected success history.
The zero-branch operator has source charge zero, so a record s carries
total charge sum_j s_j.

Set A=B*B. The whole equal-time configuration algebra is commutative,
including all translates and adjoints. Normality also gives BB*=A.
Consequently every one of the2^n sign sequences using only±1 has the
same original probability

    w_n(x_1,...,x_n)=omega(product_(j=1)^n alpha_(x_j)(A))
                    =||K_s Ω||² >=0.                     (3)

Each factor is positive and at most I/4. Hence w_n<=4^-n. A tuple with
positive w_n simultaneously supports all these sign sequences, whose
total charges are -n,-n+2,...,n. This is a statement about source charge;
the number of creation operations is not identified with particle number.

## 2. A finite averaging inequality without clustering

Let C_R=[-R,R]^3, R>0, and let

    A_R=average_(x in C_R) alpha_x(A).

The bounded strong operator integral exists by source translation
continuity. It is positive, at most I/4, and omega(A_R)=omega(A)>=p_*
by invariance of the selected state. Independent averaging of (3) over
n positions gives exactly

    average_(x_1,...,x_n in C_R) w_n=omega(A_R^n).

Fubini applies to bounded matrix coefficients. The scalar spectral
measure of the positive operator A_R and convexity of t->t^n give

    omega(A_R^n)>=omega(A_R)^n>=p_*^n.                    (4)

Equation (4) holds for every R. It uses neither a mean-ergodic limit nor
a one-dimensional invariant subspace. No clustering or asymptotic state
factorization is assumed. The receipt includes exact controls in which
some translations give zero probability or an invariant mixture gives a
mean larger than the squared first moment, preventing stronger inferences
from an averaging inequality.

## 3. An explicit cube containing a separated rational configuration

Prescribe D>0. Independent points uniform in C_R have, for each pair,

    Pr(||x_i-x_j||_infinity<=D) <=(D/R)^3,

since the forbidden conditional volume is at most(2D)^3. A union bound
therefore bounds the near-tuple fraction by binomial(n,2)(D/R)^3.
Its contribution to the averaged original probability is at most

    4^-n binomial(n,2)(D/R)^3.                            (5)

Choose the least positive integer K satisfying

    K³ >= 4 binomial(n,2)/(4p_*)^n,       R=D K.           (6)

For n=1, K=1 and the pairwise constraint is empty. The integer K is
computed by exact doubling and bisection, with every comparison rational.
Equations(4)-(6) imply that the integral restricted to tuples with every
pairwise infinity distance strictly greater than D is at least3p_*^n/4.
Since its normalized volume is at most1, an interior allowed tuple has
w_n>p_*^n/2. Boundary sets have zero volume. Strong continuity then permits
a nearby rational tuple inside C_R, with the same strict separations and

    w_n>=p_*^n/4.                                        (7)

A common translation can place the first centre at zero; the other
relative centres then lie in [-2R,2R]^3, with every separation and
probability unchanged. The same tuple works for every±1 sign sequence. Since D is arbitrary,
this supplies configurations at arbitrarily large pairwise separation.
Equation(6) is an explicit finite search domain, not a computation of a
good tuple or a modulus for evaluating the infinite selected-state
probability. Success is not asserted for every separated tuple.

The positions are centres of translates of the source's Schwartz
smearing profile. An annular Fourier support is not compact spatial
support. Thus (7) does not assert disjoint compact operator supports or
independence of the smeared states. Their equal-time commutativity comes
from the configuration-field construction, not from discarding tails.

For n=2, equation(7) is the original charge-pair bound p_*²/4. Both
(+1,+1) and(-1,-1) records have that weight; the opposite-sign records
have it as well. All nine original preparation records are retained.
The finite source cube's local A_q(x) operators have a different finite
preparation law and are not substituted for these field B_q operators.

## 4. Energy of the entire unnormalized configuration

The source creator is the bounded form integral

    B=(1/2pi) integral exp(-i theta)
        integral_0^infinity exp(-u) sin(u X_theta) du dtheta,

with X_theta=sqrt(r)[Q3(F)cos(theta)+Q2(F)sin(theta)].
For ell=L/r, the inherited all-input energy-form kick estimate gives

    sqrt(q_ell(exp(±iu X_theta)psi))
        <=sqrt(q_ell(psi))+sqrt(3)|u| ||psi||.

Taking the sine difference and the energy-form triangle inequality
preserves that upper bound. The Laplace and angular integrals converge
in the closed form norm: integral exp(-u)du=1 and
integral u exp(-u)du=1. Therefore, for every form-domain input,

    sqrt(q_ell(B psi)) <=sqrt(q_ell(psi))+sqrt(3)||psi||.   (8)

The same estimate holds for B*, and for all spatial translates because
translations commute with L. In a unit-sign sequence, after j operations
the vector norm is at most2^-j. Starting from q_ell(Ω)=0, iteration of(8)
gives

    sqrt(q_ell(K_s Ω)) <=sqrt(3)sum_(j=0)^(n-1)2^-j,
    q_L(K_s Ω) <=12r(1-2^-n)².                            (9)

No intermediate state is normalized or assumed to be a ground state.
This energy bound concerns the±1 sequences in (3); it is not an assertion
that every other preparation record has the same energy.

For a separated rational tuple from (7), define

    E_max=96r(1-2^-n)²/p_*^n.

Spectral Markov bounds the original weight above E_max by p_*^n/8.
Thus the original joint record and energy event has probability at least
p_*^n/8 in

    0<=L<=96r(1-2^-n)²/p_*^n.                             (10)

At n=2, (9) gives27r/4 and (10) gives the pair cutoff54r/p_*².
The lower endpoint is zero. Neither result excludes a zero-energy
component or confines the nonlinear total momentum to the one-field
annulus. Every energy and momentum complement remains an outcome.

## 5. Nonzero integer-charge sectors and their energy infima

For any nonzero integer q choose n=|q| with all signs equal to its sign.
Equation(7) gives a nonzero source vector of exact charge q. The vacuum
supplies the charge-zero sector. Hence every integer source-charge sector
is nonzero. Exact circle covariance puts these vectors in the stated
charge eigenspaces; no unbounded product-domain assumption is needed.

The normalized vector from (7) has mean energy at most

    48r(1-2^-n)²/p_*^n.

The probability floor is independent of r. Letting r decrease to zero
proves that the energy infimum in each integer-charge sector is zero.
The actual L commutes with source charge, so each such sector is reducing.
This conclusion allows a zero-energy kernel. It does not prove that
positive spectrum accumulates at zero for every charge, exclude a gap
above zero eigenvectors, or identify a positive-frequency particle.
The single-charge positive-band theorem is a separate stronger statement
about its specified preparation and detector.

## 6. Operational attachment and scattering boundary

The parent bounded creators and their adjoints are common strong limits
of finite bounded source histories. Their finitely many translated
products and positive square-root complements belong to that same
history construction. Every finite approximation has the complete3^n
branch instrument. The existing finite-instrument compiler supplies its
arbitrary-input/reference error with all first-failure labels.
No effective infinite-ground preparation modulus or good-location search
oracle is supplied by the explicit domain(6).

A complete experiment with unhalved trace error delta reduces the stated
original event probability by at most delta/2 before conditioning.
Source-continuum approximation and finite native pulse error are separate
budgets. Existence of a limiting nonzero record is not an executed native
tape or a laboratory observation.

The construction supplies arbitrarily separated, finite-energy source
charge configurations. Scattering additionally requires identified
stable one-excitation channels, an asymptotic identification into this
same Hilbert space, quantitative incoming and outgoing wave-operator
tails, and the associated complete detectors. Charge sectors, repeated
creation operations and positive original probabilities do not alone
supply a particle count, Fock factorization, mass shell or scattering
matrix. The exact finite response and conditional Cook attachment are
proved in [PROOF.md](PROOF.md).
