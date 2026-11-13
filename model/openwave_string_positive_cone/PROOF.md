# A uniform positive source beyond coefficientwise positivity

## Native carrier and scope

At each site retain both spins of the two native colourful quark flavours
and their conjugates: four fundamental and four conjugate colour triplets,
with 24 creation modes. The local Hilbert space is their full exterior
algebra, of dimension 2^24. The number operator N counts all these modes.
Colour generators have the literal native normalization C^2=4/3 on each
one-particle triplet. No particle-number, local colour or multiplicity
sector is discarded.

The spatial source is a finite open chain with full SU(3) link coordinates
and zero external colour flux. The exact tree Gauss isometry from
`openwave_colour_string` identifies its physical constraint kernel with
the complete global-singlet matter Fock space. Every local representation
and every allowed total particle number remains. Nonphysical preparation
failures remain flagged in the kinematic apparatus; this theorem does not
declare them physical zero-boundary-flux states.

The coefficients mbar, delta, U, A, t, Js, Jw, eta and lambda are precisely
the five-input spin-gradient source. Replace its local colour coefficient
by U-kappa/2 and its global colour penalty by the open-chain electric action
kappa sum E_x^2. The covariant hopping, local spin/weak squares, baryon-density
square and spin-gradient squares are unchanged and positive. The flavour
number term obeys (delta/2)D >= -|delta|N/2. The negative adjusted local
colour coefficient is retained in the actual source; positivity follows
from the complete operator inequality below.

## Exact local Fock bound

The native Cartan weights are four copies of
(1,0),(-1,1),(0,-1) and four copies of their negatives. The complete
occupation character is the finite polynomial

`prod_w (1 + z exp(w))`.

For each particle number n, multiply its weight character by the Weyl
denominator with positive roots (2,-1),(-1,2),(1,1). The coefficient at a
dominant weight (p,q) is its irreducible multiplicity. This follows directly
from the Weyl character formula: among the terms of each alternating Weyl
numerator, only its dominant highest-weight term lies in the dominant
chamber after the rho shift. The executable uses exact integer convolution.
Every coefficient is nonnegative, and the sum of multiplicity times
dimension `(p+1)(q+1)(p+q+2)/2` is exactly binomial(24,n), for every n.
Thus the census covers all 16,777,216 states in 547 number/irrep rows.

On (p,q), the native quadratic Casimir is
`[p^2+q^2+pq+3p+3q]/3`. Every retained row obeys C^2<=3n. Equality away
from the vacuum occurs only at n=8, (p,q)=(4,4), with Casimir 24 and
irrep dimension 125. Consequently the exact operator inequality is

`C_x^2 <= 3 N_x`.

This is a one-generation bound. Adding interacting coloured generations
changes the local exterior representation and requires another bound.
Leptons do not carry these currents, but an interacting spectator transfer
is not inferred merely from that observation.

## Uniform Gauss bound and number floor

On the physical kernel, C_x=E_x-E_(x-1), with both endpoint fluxes zero.
For every colour component and every vector, the parallelogram identity
gives the exact positive factorization

`4 sum ||E_x psi||^2 - sum ||C_x psi||^2
 = sum ||(E_x+E_(x-1)) psi||^2 >= 0`.

The code checks the integer incidence/plus-matrix factorization on chains
of length 2 through 17; the displayed telescoping identity proves all finite
lengths and does not require commuting colour components. Therefore

`kappa sum E_x^2 + (U-kappa/2) sum C_x^2
 >= (U-kappa/4) sum C_x^2`.

If U-kappa/4 is nonnegative this colour operator is nonnegative. Otherwise
the local Fock inequality bounds it below by
`-3(kappa/4-U) N`. Adding every unchanged positive source term proves

`H >= [mbar-|delta|/2-3 max(0,kappa/4-U)] N`.

The bracket is strictly positive whenever

`0 <= kappa < 4U + 4(mbar-|delta|/2)/3`.

For the actual rational five-input source the upper endpoint is about
1130.13842094. This is a sufficient full-operator positive interval, not a
claim that its boundary is the optimal instability threshold. It includes
negative local colour coefficients. The empty physical finite matter
vacuum is unique at zero energy, and the inequality is uniform in chain
length. It is not a construction of a four-dimensional continuum or of the
pure Yang-Mills theory. There are no independent transverse pure-gauge
degrees of freedom on this tree constraint kernel.

## Why the three-quark shortcut alone is insufficient

The full three-quark singlet has C_x^2=(2/3)n_x(3-n_x), and its electric
energy is (4kappa/3) times the span of its occupied sites. Hence its
kappa-dependent colour remainder is nonnegative for every configuration.
This permits any kappa on that sector, but does not prove full-Fock
positivity.

A complete native six-quark witness occupies three sites with counts 1,4,1
and local irreps 3,bar6,3. The middle bar6 is the symmetric combination of
two literal antisymmetric native quark pairs. Contracting the two outer
fundamentals with these pairs gives a global singlet. The executable checks
every native colour generator on its complete 81-coordinate carrier, its
exact integer norm 24, and the antisymmetric-pair CAR intertwiner. The local
Casimir sum is 6, whereas the two prefix Casimirs sum to 8/3. Its electric
minus half-local colour coefficient is therefore -1/3. The complete
number-plus-colour bound retains this sector and remains positive at the
selected coupling. No negative full-source energy is claimed from this
colour-only control.

## Same two-cell source and one common physical length

For any two-cell total-colour singlet, C0+C1 annihilates the state and
C0^2=C1^2 on it. Thus

`(U-kappa/2)(C0^2+C1^2)+kappa C0^2 = U(C0^2+C1^2)`

for every kappa, without a coefficientwise sign assumption. This identity
preserves the complete 120 baryon, 64 meson and 2520 nuclear singlet operators,
their states and every colour-preserving native current matrix element.
The executed native N3/N6 intertwiners are complete-sector controls of the
identity. All five empirical strong calibrations, including the native axial
value and nuclear binding, therefore transfer exactly. The identity is not
asserted on the retained coloured kinematic complement or on larger chains.

The nucleus supplies an outward graph-quadrupole interval q. In the
explicit same-length comparison, Q=a^2 q and sigma=(4kappa u/3)/a, where
u is the electron GeV unit. Use Q=0.285699 fm^2 and the illustrative target
sigma=1 GeV/fm with the ancestry in inputs.json. Integer square-root
enclosures give the required kappa interval approximately
[1109.457359467683,1109.457360715105]. Every point in that interval has a
strict positive number floor above 15.510795 source units.

The rational selection kappa=1109.45736 and exact length
`a=4kappa u/(3sigma)=0.7559087290603967... fm`
give an outward quadrupole interval approximately
[0.28569899963170337,0.2856990002741564] fm^2, containing the comparison
central value. The string tension is exactly the declared illustrative
target under that length map. These two quantities select a retrospective
physical-scale attachment; they are not independent precision predictions.
No extra empirical input changes the five strong-source coefficients.

## Complete records and work

The local number and quadratic/cubic colour Casimirs admit joint central
spectral projectors. The cubic eigenvalue is
`(p-q)(2p+q+3)(p+2q+3)/18`; the executable checks that the two central
values distinguish all irreps at fixed number. Their full ranks are the exact multiplicities times
irrep dimensions. On the maximally entangled local Fock/dual-reference
state, the complete 547-row record law has probabilities rank/2^24 and sums
to one. The dual reference transforms conjugately, so the preparation is
gauge invariant. Multiplicity spaces and post-measurement quantum states
are retained. These local central instruments commute with the local
source terms; a spatial probe's energy change is evaluated on its actual
Hamiltonian instead of inferred from that local fact.

In the complete N3 singlet on L sites, span projectors have exact ranks
20L for span 0 and

`(L-r)[binomial(4(r+1)+2,3)-2binomial(4r+2,3)+binomial(4(r-1)+2,3)]`

for span r>=1. This is inclusion-exclusion on all symmetric three-particle
orbital occupations. Their ranks sum to binomial(4L+2,3). On the uniform
full native three-quark carrier, colour admission has success probability
binomial(4L+2,3)/binomial(12L,3); its entire failure is retained. The full
Fock apparatus also includes the complementary number-sector effect;
its probability on this explicitly fixed-number input is zero. Gauss-dress
each span projector on the admitted branch. The full instrument is complete,
with string observable 4kappa r/3, and retains an idle reference. Its
nonselective mean work on the uniform input is exactly zero because the
projective map is unital. General state or branch work is the actual energy
expectation difference and is not asserted zero. Static string energies do
not include the rest-energy terms or constitute a decay calculation.
