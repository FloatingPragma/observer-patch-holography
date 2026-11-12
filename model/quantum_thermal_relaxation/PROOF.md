# Complete quantum thermal relaxation on the calibrated source

The Hilbert space is precisely the calibrated two-momentum source of
`particle_joint_mass_identity` and `thermal_particle_records`: all 180
canonical fermion modes and the 56 boson modes with global occupation
cap two. The source dimension is `2^180*1653`. Its fixed excitation
Hamiltonian, empirical charged-fermion coefficients, tree boson background
and optional Weinberg inputs are unchanged. Temperature and collision
services are supplied bath data. This is a different dynamical law from
native descending repair and from the separate 1,296-state classical bath.

## Whole-colour fermion collisions

The 84 source joint groups have rank one or three. The exterior-algebra
factorization by these orthogonal groups is unitary on every occupation.
For rank r and frequency E, its complete Gibbs state on the 2^r occupation
words is `tau_g(w)=q^|w|/(1+q)^r`, `q=exp(-E/T)`. One collision is

`Phi_g(X)=X/2+tau_g Tr_g(X)/2`.

Its public flag distinguishes hold and replacement. Hold has Kraus I/sqrt2.
The replacement has quantum bath coordinates `(a,b)` and Kraus
`sqrt(tau_g(a)/2)|a><b|`; all d^2 coordinates are retained. Its dilation
places the old source state, including every reference correlation, in
the bath register b and supplies the Gibbs purification on source and a.
It does not delete or copy an unknown quantum state. A whole-colour group
has a thermal density scalar on each occupation subspace, so this channel
and its covariant bath dilation respect colour rotations. The paired
purification carries the compensating charge. Measuring the individual
colour coordinates is not part of the public readout.

After n sweeps the probability that all 84 groups have been replaced at
least once is exactly `(1-2^-n)^84`. This gives a CP decomposition with
that coefficient of the full fermion Gibbs replacement channel. For any
untouched reference, the unhalved trace distance of the source/reference
marginal to `tau_F tensor rho_R` is therefore at most
`2*(1-(1-2^-n)^84)`. The full bath archive retains the old information.

## Capped boson instrument

The actual circular-polarization Fock rotation is inherited from the
same-source number construction. Let Pi_n denote total boson number n,
for n=0,1,2, and let a_j be the capped source annihilators. The supplied
number monitor and jump apparatus have operators

`K_(n,j,-)=sqrt(alpha) a_j Pi_n`,
`K_(n,j,+)=sqrt(alpha q_j) a_j^dagger Pi_n`,
`K_(n,0)=sqrt(I-alpha R) Pi_n`,

where `R=N+1_(N<2) sum_j q_j(N_j+1)`. Set `Q=sum_j q_j`. The executable
chooses a positive dyadic alpha no greater than `1/(Q+2)`, using proved
upper bounds on Q. Every eigenvalue of R is at most Q+2; the Kraus
adjoint products sum to identity on all 1,653 states. All 339 fine
quantum-coordinate operators are retained, including the 112 exact
forbidden zero operators. The public pointer records n, hold or the
creation/removal joint group: 171 records, with mode coordinates inside
each group retained coherently in the quantum bath. The globally capped
operators are used throughout; no uncapped Bose sampler is substituted.

Gibbs detailed balance holds on every occupation edge because its forward
and reverse probabilities have ratio q_j. No-jump diagonal terms complete
each row. Thus the exact finite capped Gibbs density is fixed. Within
each whole colour fiber the frequencies are equal, the jump family
transforms as its colour representation, and the bath mode coordinate
transforms contragrediently. An entangled one-gluon/reference state is
transferred coherently to this bath on removal. Dephasing its two occupied
colour coordinates instead changes that conditional density by unhalved
trace norm one, as executed independently.

## Complete reference-safe mixing bound

Consider two boson sweeps and trace the bath only to form the terminal
source channel. On the input n=0 sector, two holds give vacuum replacement
weight `(1-alpha Q)^2`. On n=1, a removal followed by a vacuum hold gives
weight `alpha(1-alpha Q)`. On n=2, two removals give weight `2 alpha^2`,
because `sum_ij a_i^dagger a_j^dagger a_j a_i=N(N-1)`. The first number
monitor separates these sectors. Their Kraus subfamilies therefore
dominate the completely positive map

`delta R_vac`, with `delta=2 alpha^2` and `R_vac(X)=|0><0|Tr X`.

The remaining map divided by 1-delta is a channel. For any reference,
replacement cancels the difference between rho_SR and tau_B tensor rho_R,
so the unhalved terminal distance is at most
`2*(1-delta)^floor(n/2)`. Combining the commuting fermion and boson
channels bounds the complete source/reference marginal by

`min(2, 2*(1-(1-2^-n)^84)+2*(1-delta)^floor(n/2))`.

For every positive T, Q<=56 gives the uniform minorization scale near
1/1682; the receipt uses the sharper actual temperature bounds and exact
dyadic alpha. The displayed sufficient count uses floor(n/2)>=20/delta.
A positive Taylor lower bound proves exp(20)>10^8; the fermion term
is below 10^-8 at this count, so the total is below 10^-6. This is a proved
resource count, distinct from the executed 32- and 256-sweep schedules.
No claim of input-independent complete archived histories follows.

## Executed thermal response and arithmetic

The initial finite Gibbs state is supplied at 100 GeV. The six bath
settings are 10,1,0.1,1,10,100 GeV, each held for 32 or 256 sweeps.
Every cut retains all 84 full fermion occupation factors and all 1,653
boson populations, together with the 435 boson group-count probabilities.
Every history is specified by the complete local Stinespring factors in
their original order; no source or bath coordinate is postselected away.
The executed diagonal populations are marginal laws, not an explicit
array materialization of the exponentially large quantum archive.

The group replacements generally produce non-Gaussian fermion factors:
their factorial second moment need not equal a binomial with the same
mean. The first sweep of each stage records this independent control.
Species shares are ratios of expected counts across all trials.

For exact fixed-H Gibbs-preserving channels, data processing gives
`Delta S-beta Q=D(rho||tau)-D(Phi(rho)||tau)>=0`, with heat into the source.
For an untouched reference use tau tensor rho_R; the same identity applies
to the joint entropy. The displayed heat/entropy values are numerical
diagnostics of this exact identity, not certified physical error bars.

Boltzmann coefficients have rational enclosures: the alternating series
for exp(-x) at x<=1 is bounded by consecutive partial sums, and repeated
squaring uses outward dyadic rounding. The exact source energy intervals
are retained. Initial Gibbs normalization and each transition coefficient
are bounded against these intervals. Positive sparse multiplication has
the standard gamma bound, plus an absolute allowance for subnormal
underflow across all 1,653 outputs. Repeated steps propagate the L1 error
without renormalization. Fermion updates have contraction 1/2; their
coefficient and rounding bounds are summed through that recurrence in
the checked error<=1/4 domain. The product of factor norm bounds and the
boson bound encloses the complete source-state L1 error. Hardware errors,
native preparation and a cosmological clock are separate premises.
