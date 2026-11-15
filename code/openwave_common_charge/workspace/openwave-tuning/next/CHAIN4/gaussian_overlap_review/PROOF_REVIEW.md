# Independent Gaussian overlap review

The additive estimate in native_field_opposing passes. This review is
independent of that estimate's author. It consumes the exact positive spectral
certificate from field25 and does not claim to reprove the whole field audit.

Write C and Lambda for the exact and implemented positive curvature forms,
with norm(C-Lambda) at most eta and both bounded below by lambda0 I.
Let X=sqrt(C)-sqrt(Lambda). The identity

    sqrt(C) X + X sqrt(Lambda) = C-Lambda

has solution obtained by integrating
exp(-t sqrt(C)) (C-Lambda) exp(-t sqrt(Lambda)) from zero to infinity.
Taking the operator norm gives eta/(2 sqrt(lambda0)); no commutation is
assumed. Multiplying by sqrt(b/a), then conjugating by W0^(-1/2), bounds
R-I by rho=eta/(2 lambda0), where W=sqrt((b/a) C),
W0=sqrt((b/a) Lambda), R=W0^(-1/2) W W0^(-1/2).
The lower bound for W0 used here is conservative and valid.

The normalized real Gaussian has wavefunction proportional to
det(W)^(1/4) exp(-x^T W x/2). Direct integration gives squared overlap
2^d sqrt(det(W)det(W0))/det(W+W0). Factoring out W0 gives the product
of 2 sqrt(r)/(1+r) over the positive eigenvalues r of R. Consequently
minus the logarithm of squared overlap is the sum of
log(cosh(log(r)/2)). For norm(R-I)<=rho<1,

    |log(r)| <= rho/(1-rho),
    log(cosh(x)) <= x^2/2.

The second inequality follows by integrating tanh(x)<=x for x>=0 and
using evenness. Thus minus log(overlap^2) is at most
d rho^2/[8(1-rho)^2]. Since both wavefunctions are real and positive,
their squared vector distance is 2(1-overlap), bounded by that logarithm.
This proves sqrt(d/8) rho/(1-rho). Positivity fixes the harmless overall
phase; the native controlled vacuum scalar remains a separate source term.

The exact arithmetic review checks rho from the frozen field spectral
certificate and squares the supplied rational upper bound to verify its
outward rounding. It verifies each new full source and edge budget equals
the old complete budget minus the old vacuum allowance plus the new one.
It changes no source, preparation, cutoff, clock, current term, native export
or reset probability. Adding spectator reference factors preserves the
initial vector bound; channel contraction transfers it to record distance.

Twenty-four independent rational two-dimensional widths with nonzero
commutators check the determinant formula and the bound using outward
integer-square-root intervals. Off-diagonal-only changes retain a strictly
positive distance, rejecting a diagonal-only shortcut. These finite controls
support the general proof; they do not replace it.

The reviewed vacuum bound is 9.70825949977e-7. The huge full field programs
remain unexecuted finite iterators. This is a second review of the additive
vacuum inequality and budget substitution, not a complete physical-row award.
