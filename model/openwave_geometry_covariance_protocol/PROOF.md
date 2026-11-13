# Complete covariance records and a common normalization

Let `Hg=(omega/2) sum_i(I-X_i)`, `Q=(q/N) sum_i Z_i`, with q=1/8,
omega=1 and a declared positive integer N. The geometry carrier contains all
2^N states. Its zero-source vacuum is `|+>^N`. The positive native matter
operators are the full Gaussian Fock operators transported by
`openwave_fock_field_transport`; their vacuum energy is subtracted once.
Declare a positive two-polarization free-photon Fock operator K_gamma with
dimensionless small-momentum slope one. No photon-number cutoff is imposed.
The common extension is

`H=Hg+exp(Q)K_s+exp(F Q)K_p+exp(f_gamma Q)K_gamma`.

Every factor uses the same source-time parameter and spatial unit. F and
f_gamma are specified profiles in [0,1]. The finite native matter and
geometry operators are bounded. In a geometry-coordinate basis,
`exp(f_gamma Q)K_gamma` is a finite direct sum of positive multiples of
K_gamma, with common domain D(K_gamma). Adding the bounded operators gives
a self-adjoint, lower-bounded Hamiltonian on that domain. Its positive
terms give the lower bound zero. This is a declared extension; no transfer
from the native twelve-port process or interacting photon field is proved.

## The complete instrument

For every k=0,...,N retain the full eigenspace projector P_k of Q, including
its binomial multiplicity. The effects `P_k tensor I_s tensor I_p tensor
I_gamma` sum to the complete identity. With N+1 degenerate memory states,
`U=sum_k P_k tensor T^k`, where T is the cyclic shift, is a unitary on
every memory input. A blank memory and its complete occupation measurement
implement the count-Lueders instrument. No bit label inside a count
eigenspace is measured or discarded as a failure.

Since X_i changes the count by one, `sum_k P_k X_i P_k=0`. All other terms
of H commute with Q. Thus the exact mean-work operator is

`sum_k P_k H P_k-H=(omega/2) sum_i X_i`.

The vacuum work is N omega/2 source units per copy. An independent pair
costs N omega u GeV, where u is the declared energy unit. The general work
can have either sign; its value depends on the input geometry energy.
Blank memories and independent ground preparations are supplied resources.
The instrument specifies its work; a zero-cost autonomous memory reset is
not inferred.

The vacuum count probabilities are `p_k=binomial(N,k)/2^N`, with
`x_k=q(2k-N)/N`. For two independent copies retain all (N+1)^2 records and
assign

`D(x,y)=(x-y)(exp(x)-exp(y))/2`, `Z(x,y)=(x-y)^2/2`.

Independence gives `E D=Cov(Q,exp Q)=C` and `E Z=Var Q=v=q^2/N`.
Equal-coordinate records have D=Z=0 and remain in the sample. Both
statistics use the same complete measurements.

## Sample and resolution bounds

The Rademacher moments give `E Q^4=(3-2/N)v^2`, hence
`E Z^2=(3-1/N)v^2`. The mean-value theorem on [-q,q] gives
`exp(-q) Z <= D <= exp(q) Z`. Consequently

`C>=exp(-q)v`, `E D^2<=3 exp(2q)v^2`,
`Var(D)/C^2<=3 exp(4q)`.

For a block of b independent pairs, Chebyshev's inequality bounds the
relative-error probability at epsilon by `3 exp(4q)/(b epsilon^2)`.
Choose `b>=12 exp(4q)/epsilon^2`; each block fails with probability at most
1/4. For an odd number K of independent blocks, their median fails only if
at least (K+1)/2 blocks fail. Independent Bernoulli variables with success
probabilities at most 1/4 are stochastically dominated by K independent
Bernoulli(1/4) variables. This proves the exact binomial tail stored in the
receipt. No normal approximation, fitted noise variance or sample rejection
is used. All recorded blocks are retained even when their mean is not the
median. The sample bound is independent of N; the work per pair grows as N
and the required coordinate resolution shrinks with sqrt(v).

Suppose each reported coordinate has deterministic error at most eta.
Write Delta=x-y and
`D=Delta^2 integral_0^1 exp(y+t Delta) dt/2`. The interpolated coordinate
error is at most eta. Using `E |Delta|<=sqrt(2v)` yields

`|E D_reported-C| <= exp(q+eta)(2 eta sqrt(2v)+2 eta^2)`
`                    +(exp(eta)-1) exp(q) v`.

The mean reported lapse is between exp(-eta) and exp(eta) times its exact
mean. These are apparatus-resolution bounds, not measured physical noise.
The receipt gives an illustrative eta proportional to sqrt(v) and retains
the ancestry of the count selected using observed G. The median theorem
above concerns exact coordinate records; this deterministic bias must be
included when applying it to rounded records. It does not assert the same
variance bound for an arbitrary noisy detector.

The same pair records give `A=(exp(x)+exp(y))/2`, with mean M(1).
The exact product identity gives
`Var(A)/M(1)^2=([1+tanh(q/N)^2]^N-1)/2<=v/[2(1-v)]`.
Since v<=q^2, the receipt chooses a common block size that also bounds
the relative mean-lapse failure by 1/4. Two median estimates then have
simultaneous coverage at least one minus twice the displayed binomial
tail. This uses the union bound; the two statistics on a pair need not be
independent. For estimates C_hat and M_hat, their simultaneous intervals
are `C_hat/[1+epsilon_C,1-epsilon_C]` and
`M_hat/[1+epsilon_M,1-epsilon_M]`. The executable propagates a synthetic
model-mean example through the full unit map. It supplies no observed
metric-noise data.

## Rest, light and length use the same source

Set `M(f)=<exp(fQ)>=cosh(fq/N)^N` and
`C=q M(1) tanh(q/N)`. On the chosen matter energy branches define

`E_N(x,y,F)=min spectrum[Hg+x exp(Q)+y exp(FQ)]`.

The isolated zero-source geometry ground has gap omega. Differentiating
its normalized eigenvalue equation at zero gives
`dE_N(x,0,0)/dx|0=M(1)` and `dE_N(0,y,F)/dy|0=M(F)`.
The same argument for one photon of free energy lambda(k) gives
`E_gamma(k)=M(f_gamma) lambda(k)+O(lambda(k)^2)`; its small-momentum speed
is M(f_gamma) in the common source units. These are ground-branch
derivatives. No finite-momentum dispersion or relativistic continuum is
inferred from them.

The full spectral-response theorem of `openwave_gravity_noise_normalization`
gives the far-profile coefficient
`g_N=S C/(2 pi omega)`, multiplying the bare labels x y. The physical weak
rest responses are `M(1)x` at the source and `M(0)y=y` at the far probe.
Therefore the rest-normalized coefficient is `g_N/M(1)`. Omitting this
factor would mix a bare source label with its measured rest response.

Let h denote hbar c in GeV fm, u the common energy unit in GeV and a the
single spatial unit in fm. Then

`c_gamma/c=u a M(f_gamma)/h`,
`G_weak=(g_N/M(1)) a/(u h)`.

Eliminating a gives

`G_weak/(c_gamma/c)=S C/[2 pi omega u^2 M(1) M(f_gamma)]`.

In the far photon reference, f_gamma=0, speed matching uses a=h/u and
the mean lapse cancels exactly:

`G_weak=S q tanh(q/N)/(2 pi omega u^2)`.

For a photon calibrated at profile one, there is an additional divisor
M(1). The two conventions cannot be interchanged. The far reference permits
the inherited electron energy anchor on the far native branch, whose
geometry coupling has profile zero. A finite source rest energy is
E_N(x,0,0), not M(1)x. Promoting the weak coefficient to a finite-source
Newton law requires the finite-strength and particle-identification
errors. The complete unit map supplies neither that promotion nor a
physical metric-noise identification.

`propagate_strength` encloses this common-unit expression for positive
intervals on C, both slopes, u and the measured speed ratio. Interval
monotonicity proves the endpoints; no statistical independence between
these intervals is assumed. Their simultaneous coverage must be supplied
by the observation protocol.

The strictly decreasing function tanh(q/N) gives distinct coefficients
at N,2N,... with the same far photon speed, clock and energy calibration.
An independent measurement of v would give N=q^2/v and an independently
measured C would constrain the response. No such physical datum is supplied.
The small-count comparisons are declared mathematical controls. The huge
count inherited from observed G is labeled as calibrated in every resource
estimate. None supplies a new validation against Nature.
