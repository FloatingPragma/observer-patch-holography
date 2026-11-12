# Certified component mixing of the declared bath

The source is exactly `native_thermal_completion`: 1,296 states, 26
communicating components, the same 17,352 accepted/rejected/idle records
and the original preparation's component weights. Temperature is the
supplied cost parameter `q=exp(-2 beta)`, not a GeV temperature. Attempts
are the bath's discrete clock, not cosmological time. Native descending
repair is a different law.

## Exact symmetric operator

On each component C the strictly positive Gibbs law is proportional to
`q^((E-Emin)/2)`. For the complete row-stochastic bath P, set
`S=diag(sqrt(pi)) P diag(1/sqrt(pi))`. Detailed balance makes S symmetric.
Its diagonal is the exact rational Pxx. For x different from y its entry is

`Sxy=(sxy/12)*sqrt(q^(abs(Ey-Ex)/2))`.

All square roots are enclosed by integer square-root bounds on a dyadic
grid; floating conversion error is then measured exactly as a rational.
The complete graph and detailed-balance identities are inherited and
rechecked by the unchanged record digest. P has diagonal at least 1/2,
so its entire spectrum is nonnegative. Connected positive support gives
a simple eigenvalue one in each nonsingleton component. The twelve
singletons have no nonconstant spectrum and no mixing error.

## Certified spectral enclosure

A binary64 symmetric eigensolver supplies only a candidate V and diagonal
D. Neither orthogonality nor its numerical residual is assumed exact.
For n-component dot products, gamma_k=k u/(1-k u), u=2^-53, bounds all
rounded products and additions. The implementation uses k=2n+2 for Gram
and absolute-Gram products, and k=2n+4 for V D V^T. Positive absolute
products provide the sums of absolute terms, enlarged by 1/(1-gamma_k).
The stated arithmetic domain excludes overflow and underflow of the
products; nonzero eigenbasis product magnitudes are checked.

The maximum enclosed Gram entry gives
`epsilon >= ||V^T V-I||_2`, using `||A||_2<=n max|Aij|`. The maximum
enclosed reconstruction entry and exact input intervals likewise give
`rho >= ||S-V D V^T||_2`. All subsequent calculations use rational
arithmetic. When epsilon<=1/2 the polar factor U of V is orthogonal,
`||V-U||<=epsilon` and `||V||+1<3`. Consequently

`||S-U D U^T|| <= rho+3||D|| epsilon = delta`.

Weyl's inequality encloses every ordered eigenvalue of S within delta
of the corresponding diagonal value of D. In particular the complete
nonconstant gap is enclosed by `1-D[n-2] +/- delta`. The receipt rejects
any case whose lower endpoint is nonpositive. This is a full-basis
certificate, not an isolated approximate eigenvector test.

## Mixing and cooling lag

Self-adjoint contraction on the mean-zero subspace gives, for every
probability law p supported on C,

`TV(p P^n,pi)<=sqrt(1/pi_min-1)*(1-gap)^n/2`.

The displayed bound is capped by one. The complete kernel retains all
components, so there is no global unique equilibrium. For prescribed
component weights w_C the global bound is the weighted sum of the
conditional bounds. A deterministic attempt budget for epsilon=10^-3
uses `b=ceil(log2(1/(4 epsilon^2 pi_min)))` and
`n=ceil(b/(2 gap_lower))`. Since log(2)<1 and 1-gap<=exp(-gap), this is
a rigorous sufficient count with no numerical logarithm assumption.

The cooling calculation consumes every original 1,296-component state
law at all six cuts of both the 32- and 256-attempt schedules. It compares
these with the exact rational Gibbs mixture at the same q and the same
component weights. If the producer's L1 arithmetic enclosure is e, the
true total variation differs from the rationally evaluated central TV
by at most e/2. All original rejection and idle branches retain their
record hash. These bounds concern terminal state laws after forgetting
the bath history; the full retained history need not lose its information
about initial conditions. Neither a cosmic clock nor physical abundances
are derived.
