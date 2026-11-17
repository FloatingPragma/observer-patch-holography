# Native-record free-field limit and its causal reconstruction

This theorem completes the free-field interpretation of the sealed native
occupation, slice-refinement and Maxwell components. The inferred real-time
field is the effective description. It is not required to be the primitive
endpoint event law. This is a written extension, with no new simulation or
platform status claim.

## 1. Source moments and the order of limits

Use the explicitly declared native product-graph family of the row30 slice
theorem and the four-channel row29 packet. At every finite graph, the shared
occupation theorem proves convergence of all fixed joint moments to the
Gaussian covariance L^+-Pi/(4W). Independent native zero-load coin records
give the optional exact L^+ completion. Either protocol has the same local
infinite-volume limit; the coin cost is retained if that protocol is used.

For a fixed lattice spacing a, first take the block length n to infinity on
a fixed finite periodic graph. Then take the graph volumes to infinity.
Finally refine a to zero. The source moment theorem makes the first limit
an actual limit of native record statistics. A diagonal sequence of finite
resources realizes any prescribed finite list of moments to any prescribed
accuracy. The following effective rule makes that diagonal explicit, without
asserting an efficient uniform central-limit rate or finite-tape precision.

At level r, specify a finite rational lattice-smearing menu, a maximal moment
degree and tolerance2^(-r), together with the finite graph and its rational
native transition table. For an additive scalar readout f, propagate

    m_(t+1)^(k)(y)=sum_x P(x,y) sum_(j=0)^k
                    binom(k,j) f(y)^(k-j) m_t^(j)(x).

The multinomial recurrence gives every requested mixed moment. Independent
auxiliary coin moments, if used, are computed exactly and convolved with
these moments. Normalization introduces only computable algebraic factors.
The finite Green matrix and its Wick targets are likewise computed from the
finite source. Starting at n=1, double n and use exact rational/algebraic
arithmetic, or certified intervals with shrinking widths, until every moment
error is certified below2^(-r). The fixed-graph moment theorem proves
termination. All previously tested n values are retained if this procedure
is executed. No efficiency or resource bound independent of r is claimed.
For interval arithmetic, use a fixed finite precision sufficient to bound
each error interval's width by2^(-r-2) at that n. If its upper bound is not
strictly below2^(-r), advance to the next n. Equality at the tolerance can
therefore never stall an inner refinement loop. Eventually the true error
is below2^(-r-1), so this finite-precision rule terminates as well.

Choose graph volume and spacing at the same level so that the covariance
bounds below and a certified finite Riemann-sum bound meet the tolerance.
Increasing the rational smearing and degree menus gives one computable
diagonal source family. Rational approximations of a dense smooth test menu
have errors controlled by the displayed covariance and pairing bounds. This
puts the limiting resource convention into the candidate, rather than
leaving an unspecified successful subsequence.

An explicit computable dense menu is the rational linear span of rational
translations and positive rational dilations of the Fourier bump
eta(p)=exp[-1/(1-|p|²)] on |p|<1 and zero outside, including its derivatives.
Take real/imaginary symmetric combinations for real fields. Approximate
identities followed by rational Riemann sums prove density in the smooth
compact-Fourier test functions, which are dense in Schwartz space. Every
derivative bound is computable: differentiating gives a finite polynomial
times powers of(1-|p|²)^(-1) times the exponential; putting
u=(1-|p|²)^(-1), each term is bounded using the explicit maximum of u^N e^-u.

Finite Green Riemann sums have a concrete certification rule. Exclude a
momentum ball of radius rho. The continuum inner integral is bounded by a
known constant times rho² sup|f_hat|². For mesh h, the nonzero points in the
shell ||k||_infinity=j number64j³+16j, so

    h^4 sum_(0<|hk|<=rho) 1/(c|hk|²) <=80 rho²/c,
    c=4/pi².

This bounds the discrete inner contribution uniformly in h. Outside the
ball, the menu's certified derivative bounds and the positive denominator
give a computable Lipschitz bound and hence a standard cell-by-cell
Riemann-sum error. First choose rho to bound both inner errors, then h to
bound the exterior error. Compact support bounds the number of cells; all
exponential/algebraic values admit certified rational intervals. Thus the
volume/refinement selection used above requires no convergence oracle.

After projection onto a normalized internal eigenvector with eigenvalue
lambda, the infinite-volume lattice covariance symbol is

    G_(a,lambda)(p)=1/[lambda+q_a(p)],
    q_a(p)=(4/a²)sum_(mu=0)^3 sin²(a p_mu/2).                 (1)

The finite torus zero mode causes no ambiguity in this limit. On the first
Brillouin zone, sin(x)>=2x/pi for0<=x<=pi/2 gives

    q_a(p)>=(4/pi²)|p|².                                    (2)

The singularity1/|p|² is locally integrable in four dimensions. To justify
the torus Riemann sum, split momenta at |p|=delta. The nonzero lattice points
in successive radial shells give a uniform O(delta²) bound for the inner
sum of1/|p|²; the exterior sum is an ordinary smooth Riemann sum. Sending
the volume to infinity and then delta to zero proves the stated Green limit.
For lambda>0 the integrand is bounded at zero. The contact term1/(4W)
vanishes as the finite graph grows. No finite constant mode is relabelled
as a propagating photon or massive particle.

## 2. Controlled continuum smeared covariances

One common proposed source contains the four independent native components
already used for Maxwell. Its internal constant projections supply the curl
algebra, and its nonconstant projections supply the massive scalar algebras.
Orthogonal internal eigenvectors have zero limiting cross-covariance and
hence independent Gaussian limits. Their rulers, spacing and reconstructed
time are identical. This is a single declared product family, not a
different model chosen for each scalar mass or photon. It remains an
unpinned N+ architecture, distinct from the two-event controller A.

Fix its common reflection explicitly. On the internal constant sector the
four channel labels transform as a vector: A_0 is odd and A_1,A_2,A_3 are
even under time reflection. On every nonconstant internal sector the chosen
massive scalar readouts are even. These internal spectral projectors commute
with the native covariance. Their limiting Gaussian independence makes the
physical F algebra and massive scalar algebras a product of their respective
positive reflection algebras. The indefinite constant-sector potential A_0
is excluded from that physical algebra. No single field is silently assigned
two incompatible reflection signs. The direct-sum tensor/scalar interpretation
is a disclosed part of this record dictionary.

Use normalized lattice fields a^(-2) times the projected occupation
fluctuation, paired with test functions using the cell measure a^4. For a
smooth band-limited test f with Fourier support |p|<=K and aK<pi, Poisson
summation has no alias in the Brillouin zone. Its covariance is the integral
of |f_hat(p)|² times (1), with measure d^4p/(2pi)^4.

The cosine remainder gives

    0 <= |p|²-q_a(p) <= (a²/12)sum_mu p_mu^4.

Writing c=4/pi², equation(2) implies, for every native lambda>=0,

    0 <= G_(a,lambda)(p)-1/(lambda+|p|²) <= a²/(12c).

Consequently the smeared covariance error is at most

    (pi² a²/48) integral |f_hat(p)|² d^4p/(2pi)^4.            (3)

Cross-covariances follow by polarization. This bound holds simultaneously
for the twelve native internal modes, whose eigenvalues are0,5-sqrt(5),6,
5+sqrt(5), with their already established multiplicities. General Schwartz
tests follow by smooth Fourier cutoffs and dominated convergence, using(2)
and rapid Fourier decay. These are distribution-valued fields; a finite
pointwise continuum variance is not asserted.

For the Maxwell readout, let h_mu,nu be an antisymmetric test tensor and
Phi(h)=(1/2)sum_mu,nu integral h_mu,nu F_mu,nu. Its covariance is that of
four scalar channels against v_nu=sum_mu i p_mu h_hat_mu,nu. On the lattice,
replace p_mu by kappa_mu=2sin(a p_mu/2)/a. On the same band-limited domain,

    |kappa-p| <= a² |p|³/24,
    ||v_a-v|| <= a² |p|³ ||h_hat||_F/24,
    ||v_a||,||v|| <= |p| ||h_hat||_F.

Splitting the covariance difference into its numerator and denominator
differences, and using(2), bounds it by

    (pi² a²/24) integral |p|² ||h_hat(p)||_F² d^4p/(2pi)^4.   (4)

The occupation-only contact contributes an additional bound
integral |p|²||h_hat||_F²/(4W), which vanishes with volume. The completed
native-coin protocol has no such term in its Gaussian block limit.

All limiting moments are Gaussian. If a finite list of pair covariances
differs by at most epsilon and each is bounded in magnitude by B, each
degree2r moment differs by at most

    (2r-1)!! r epsilon B^(r-1),                              (5)

by telescoping each Wick pairing. Odd moments vanish. Thus the limiting full
hierarchy, not only one pole or one mode, is fixed by the source moments.

## 3. The reconstructed scalar fields

The continuum scalar covariance is

    1/(p_0²+|p|²+lambda).

It is Euclidean invariant, symmetric and tempered. Reflection positivity
follows either from the native completed Dirichlet proof followed by the
limit, or directly from its time kernel

    exp(-E_lambda(p)|t|)/(2E_lambda(p)),
    E_lambda(p)=sqrt(|p|²+lambda).

For positive-time tests the reflection norm is the squared modulus of the
Laplace transform, weighted by1/(2E_lambda). The Gaussian polynomial lift
proves the full positive hierarchy. Spatial translations and positive time
translations then act on the one-particle quotient, the latter with
multiplier exp(-t E_lambda). Their symmetric tensor powers give the field
representation. The continued positive generator has the Wightman kernel

    W_lambda(t,x)=integral d^3p/(2pi)^3
                 exp(-i E_lambda t+i p.x)/(2E_lambda).       (6)

This is a source-derived kernel: its spectral mass lambda is the native
internal response, rather than a mass coin or independently selected
oscillator. Equation(6) is a positive-energy mass-shell integral and is
Lorentz covariant. Applying partial_t²-Laplacian+lambda to it gives zero.
Thus the inferred scalar field obeys the Klein-Gordon equation, with
omega²=|p|²+lambda, as a distribution.

The commutator is the difference of the two boundary values of(6). It has
support in |x|<=|t|. One direct proof uses the Klein-Gordon initial-value
problem with zero initial position and a compactly supported initial
velocity. Its positive energy has flux bounded by the energy density;
integrating outside the expanding light cone shows that initially zero
exterior energy remains zero. Approximating the delta initial velocity by
compact supports gives the distributional commutator support. Consequently
spacelike separated smeared fields commute. The construction therefore
includes positive energy, locality and the full Gaussian observation law.

## 4. The reconstructed electromagnetic field

The sealed row29 proof supplies the direct two-form reflected norm

    (|p|/2)||P_T e_tilde-(i[p cross]/|p|)b_tilde||²,

with rank2 at nonzero spatial momentum. The full vector-potential algebra
with vector reflection is not asserted positive. The physical algebra here
is the field-strength algebra, with its null vectors quotiented out.

Equations(4)–(5) give the convergence of its full Gaussian hierarchy to the
continuum tensor covariance. Its two quotient states have energy|p| and
mode amplitudes E=i|p|epsilon, B=i(p cross epsilon) for epsilon perpendicular
to p. Direct differentiation gives both divergence constraints and the
vacuum Faraday/Ampere equations. Its commutator is a tensor of second
derivatives of the massless scalar commutator, so it also vanishes at
spacelike separation. In3+1 dimensions the massless scalar retarded kernel
is theta(t)delta(t-|x|)/(4pi|x|); differentiation preserves its causal support.
The two transverse vacuum waves therefore propagate at the same speed1
as the characteristic cone of every reconstructed massive scalar band.

No dimensionful laboratory c is inserted. Space and reconstructed time use
the same declared isotropic lattice spacing, and the shared native symbol
fixes their relative slope. Absolute units are a different comparison.

## 5. What is proved and what is not

On the stated N+ source and record dictionary, the limiting free-field
descriptions satisfy the mathematical wave clauses of rows29 and30. Their
Lorentzian time is reconstructed from native full moments; a primitive
unitary event update is not an additional requirement. The first-hit clock
packet supplies a direct source-record realization of the same Euclidean
transfer, strengthening that dictionary.

This is not a platform award or proof of the full executable-model
hypothesis. The declared architecture, product channels, readout, limits
and pin ancestry remain visible. Interacting charged matter, physical
species assignments, orbital interactions, a common detector/control
theory for all rows and comparison with independent laboratory data remain
separate obligations. No uniqueness of the architecture is required or
claimed. Another-group review, owner pin and in-platform execution are still
necessary for any official cell.
