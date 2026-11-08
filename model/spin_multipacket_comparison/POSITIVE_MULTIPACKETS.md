# Positive energy in separated internal charge sectors

## 1. Same-source preparations and a finite averaging argument

Fix r>0 and either one of the time-filtered source creators C=C_(r,+)
or C=C_(r,-) from LOCAL_CREATORS.md. Write A=C*C, 0<=A<=I/4,
and let p=9/(200000000 M²) be its original vacuum probability lower
bound. For independently chosen centres x_j in a cube of half-width R,
let A_R be the spatial average of the translated A. Translation
invariance and the spectral Jensen inequality give

    average_(x_1,...,x_n) <A_(x_1)...A_(x_n)>
       =<A_R^n> >= <A_R>^n >=p^n.                         (1)

This identity does not assume that the A_(x) commute. The average is
an operator product of the same positive average A_R. Individual
integrands can be complex; their absolute values are at most4^-n.
All integrals exist strongly, with the stated uniform norm bound.

Let D>0 and discard tuples with max-norm separation at most D between
some pair. Their fraction is at most binomial(n,2)(D/R)^3, so their
contribution in absolute value is at most4^-n times that fraction.
On the retained tuples assume the actual source locality estimates give

    ||[C_(x_i)^eta,C_(x_j)^zeta]||<=epsilon,
    eta,zeta in {identity,adjoint}, i!=j.                  (2)

The operational original preparation weight is

    w(x_1,...,x_n)=||C_(x_n)...C_(x_1) Omega||².

Moving each right-hand C_(x_i) beside its own adjoint costs at most
n(n-1) adjacent inter-site swaps. There are2n-2 remaining factors in
each error term, all with norm at most1/2. Thus

    |w-<A_(x_1)...A_(x_n)>|
       <=n(n-1)4^(1-n) epsilon.                           (3)

Both the C/C and C/C* commutator bounds are needed. Equal-time
commutation of the unfiltered configuration operators alone does not
supply (3) after energy filtering.

Choose epsilon so that the right side of (3) is at most p^n/4. Choose
an integer K with K³>=4 binomial(n,2)/(4p)^n and put R=DK. The bad-tuple
integral is then at most p^n/4. Taking real parts in (1) shows that the
integral of w over retained tuples is at least p^n/2. Hence a separated
tuple has weight at least p^n/2. Continuity supplies a rational tuple
with weight at least p^n/4 and strictly greater than D separation after
an arbitrarily small slack in these constants. Equivalently one can
make the two error budgets strict before selecting the tuple.

The producer uses conservative non-strict error budgets and records the
rational lower bound p^n/4. A boundary tuple has zero integration measure;
a positive integral over the open separated region contains an interior
point above this smaller floor. For n=1 no pair separation is imposed.
This gives a bounded search domain, not an algorithm for evaluating the
selected continuum state or finding a good point with certified cost.

The source tail bound (10) provides the required support distance
explicitly. If every centre pair is farther than D, compact profiles of
radius2R_profile/r have support distance at least D-4R_profile/r. This
subtraction is retained in the producer. Every preparation retains all
three records, so an n-step experiment has all3^n original records.
No success-conditioned normalization is inserted in (1)-(3).

## 2. Positive energy and charge

Each C raises energy by an amount in [ar/2,2br]. Repeated use of its
exact spectral-transfer identity proves

    C_(x_n)...C_(x_1) Omega
       belongs to E_L([n ar/2,2n br]) H.                  (4)

The vector has exact internal charge+n for C_+, and -n for C_-.
Source charge conjugation maps the two products into one another,
so the same centres and weights work for these two equal-sign strings.
No statement about an arbitrary mixed-sign string follows from (1).
After filtering, the opposite-charge creators need not have equal
C*C, and they are not each other's adjoints.

For every nonzero integer q, choose n=|q|. Equations(1)-(4) provide
nonzero positive-energy vectors in that internal charge sector at every
r>0. Since the upper energy endpoint tends to zero with r, positive
spectrum in each such sector approaches zero. This is stronger than
an infimum argument allowing only a zero-energy kernel. It does not
classify these vectors as elementary particles, assign a mass, or prove
that products have a tensor-product inner product. Charge zero is not
covered by this equal-sign construction.

The source spectrum, state and clock in (4) are exactly those used to
construct C. The finite-cube nilpotent local control in finite.py is
separate; its finite spectrum is not substituted into this scaling
statement.

## 3. An exact asymptotic derivative estimate

For differentiable bounded source packet operators A_i(t), let
M_i(t) bound their norms, epsilon_i(t)=||dot A_i(t) Omega||, and
c_ij(t)=||[dot A_i(t),A_j(t)]|| for i<j. Differentiate the ordered
product, then commute each derivative to the right until it acts on the
vacuum. This gives the exact estimate

    ||d/dt [A_1(t)...A_n(t) Omega]||
       <= sum_i epsilon_i product_(k!=i) M_k
        +sum_(i<j) c_ij product_(k!=i,j) M_k.              (5)

No factorization is used. An integrable right side on each of the two
time tails proves existence of the corresponding vector limits by the
Cauchy criterion. Nonzero norms and a comparison-space identification
require additional source information; an integral bound alone does
not prove scattering completeness.

For the actual time-filtered C, choose spatial packets so that on the
vacuum their smearing is fhat(P) exp(-it omega(P)) C Omega, and apply
the actual time evolution exp(itL). Here omega is a candidate function,
not an asserted dispersion. The one-packet defect is exactly

    epsilon(t)²=integral |E-omega(k)|² |fhat(k)|²
                         dmu_(C Omega)(E,k).             (6)

The commutation of the actual L and P justifies their joint measure.
The compact transfer band makes (6) finite for smooth bounded omega on
the packet support. Its right side is independent of t. A fixed-width
energy band therefore gives no integrable one-packet tail by itself.
A true spectral dispersion or another controlled comparison dynamics
must reduce this defect.

For two spatial packets, (8)-(9) of LOCAL_CREATORS.md bound c_ij by
integrating the actual almost-local commutator against |f_i f_j| and
against |dot f_i f_j| as appropriate. The packet norms are bounded by
||C_i|| ||f_i||1. All these L1 growth factors are retained in (5).
If a source-derived smooth dispersion and disjoint velocity supports
supply packet separation, integration by parts in the Fourier packet
phase gives inverse-power tails outside those supports. Choosing a
sufficiently smooth time filter then controls the cross terms in (5).
This implication does not supply the missing source dispersion or
replace it with a chosen relativistic energy function.

These preparations are even in the microscopic fermion grading. Their
internal U1 character, commuting equal-time configuration ancestors and
source comparison Gram do not determine physical spin-statistics,
Standard Model species, a Lorentz mass shell or CERN cross sections.
