# Particle production under a declared expansion

This consumer uses the calibrated `conditional_particle_spectrum` mass kernel,
its 180 canonical anticommutation-relation modes, and the exact named-field
and joint-group construction authenticated by `particle_joint_mass_identity`
and `thermal_particle_records`. Charged coefficients are the supplied running
pure MSbar parameters at 91.1876 GeV. The optional neutral sector uses the
supplied Weinberg prescription, normal ordering and zero lightest mass. No
coefficient is inferred from the expansion response. Only the fermion sector
is evolved; the capped boson factor can be an untouched spectator.

## Literal source and canonical phases

For positive scale a, replace the source momenta (0,0,+/-1) by (0,0,+/-1/a)
without changing its mass kernel K, canonical metric, or Weyl speed one.
The exact parent Takagi matrix W satisfies W^T K W=M and W^dagger W=I.
With T=diag(I4 tensor W, I4 tensor W*), its full 360-dimensional Nambu
symbol is unitarily equivalent to the same real mass symbol with the
kinetic diagonal multiplied by 1/a. This is an identity for every a>0.
The executable independently compares the literal full matrices at a=1,2,4.

Every positive kinetic coordinate i has one mass partner j. The exact
symbol entry is s*m, where s=+1 or -1 and m is its parent mass coefficient.
The particle-hole partner pair has i'=tau(j), j'=tau(i), and entry -s*m.
The consumer verifies these equalities using the exact symbolic kernel;
all 180 named modes occur in precisely 90 disjoint physical pairs.
For the zero mass, a unit coefficient selects only a pairing convention.
Its actual Hamiltonian coupling and production are exactly zero.

Put E=sqrt(p^2+m^2), c=sqrt((E+p)/(2E)) and
nu=m/sqrt(2E(E+p)). The positive named columns are
`u_a=c e_i+s nu e_j` and `u_b=c e_i'-s nu e_j'`.
The analytic parent construction uses these phases, not an independently
chosen eigensolver basis. Every column is checked against the actual named
frame. The negative column tau u_b* is `c e_j-s nu e_i`.
Furthermore T^dagger tau T*=tau exactly, so the complex flavour matrix W
introduces no omitted phase into this pairing formula.

In the ordered Nambu block (u_a,tau u_b*) the positive-mass evolution has
`S=[[alpha,beta],[-beta*,alpha*]]`. The literal signed block is
`S_s=diag(1,s) S diag(1,s)`. In the even Fock basis (vacuum,a^dagger
b^dagger vacuum) its unitary is
`G_s=diag(1,s) X S X diag(1,s)`, where X swaps the two coordinates.
The odd states have phase one in the Hamiltonian convention
`E(N_a+N_b-1)`. The pair Hamiltonian in the bare even basis is
`[[-p,s*m],[s*m,p]]`. These identities are independently replayed by literal
two-mode Fock matrix exponentials for both signs, including thermal
vacuum/pair coherence. Fix the permutation listing each stored pair in its
stored order, then concatenate all 90 pairs. The exterior-algebra reorder
multiplies an occupation word by (-1)^I, where I counts occupied inversions
of this permutation relative to the original 180-mode order. This fixes
all fermionic signs, including words with odd occupation in several pairs.
Disjoint even pair unitaries commute. Conjugating their tensor product by
this fixed reorder is therefore the complete
2^180-dimensional source unitary and extends as identity on any reference.

## Declared pulses and complete readout

The scale changes from 1 to 4. A sudden quench has no intervening time
unitary; its positive/negative frame overlap gives

`P_pair=(1-(m^2+1/4)/(sqrt(1+m^2)*sqrt(1/16+m^2)))/2`.

For each duration tau=1 or 8 in the declared GeV^-1 parameter, the executed
law has N=64,256,1024 constant segments. Segment k uses the exact rational
midpoint of `a(s)=1+3(3s^2-2s^3)`, s in [0,1]. Matrix exponentials are
evaluated in the initial energy frame: the Pauli coefficients are
`h_z=(m^2+p)/sqrt(1+m^2)` and `h_x=m(1-p)/sqrt(1+m^2)`.
This retains tiny neutrino couplings without interval amplification from
large arbitrary rotations. Every finite pulse and endpoint frame is bounded
with outward integer intervals. Square roots use integer square roots;
sine and cosine use the parent certified Taylor bounds, angle halving,
double-angle identities and a Lipschitz allowance for the angle interval.
No floating trigonometric result supplies a certificate.

To avoid rectangular interval wrapping, each pulse interval has an exact
dyadic matrix center A and an operator error epsilon. Twice the largest
entry's real-plus-imaginary radius bounds this error. The accumulated
product also has a dyadic center P and operator error e. Multiplying the
centers with directed integer rounding gives center C and rounding error r
by the same entry-norm bound. Because the ideal pulse U is unitary,
`||UP-C||` propagation is bounded by
`e_next <= e + epsilon*(1+e) + r`; here `||P||<=1+e` follows from the
previous enclosure and never assumes that P is exactly unitary. Only
centers, rather than their rectangular uncertainty boxes, enter subsequent
products. The final frame change uses the same recurrence. Its operator
error supplies the final matrix entry intervals. Floating midpoint
unitarity checks are additional diagnostics.

The two initial states are the instantaneous vacuum and the complete finite
fermion Gibbs Gaussian at supplied temperature 1 GeV. For a pair with initial
occupation f and c=|beta|^2, its final diagonal in order (00,10,01,11) is

`( (1-f)^2-(1-2f)c, f(1-f), f(1-f), f^2+(1-2f)c )`.

The retained even coherence is
`rho_(00,11)=s*(1-2f)*alpha* beta*`.
Thus all four diagonal entries and the complex off-diagonal entry specify
the full pair density. The tensor product of all 90 pair states is the full
Gaussian output, including every zero branch. Neither a selected-pair
normalization nor a one-excitation truncation is used.

At a=4 the joint group projectors are those of the final instantaneous
positive-energy fields. Their group lineage, electric charge, colour Casimir
and helicity agree with the source labels; physical momentum is +/-1/4 and
energy is sqrt(1/16+m^2). The receipt distinguishes these projectors from
the initial source projectors with momentum +/-1. Whole-colour number
groups preserve colour coherence. Their commuting spectral projectors form
a complete number instrument on every source occupation and on every
source/reference input. On the two executed preparations its exact joint
count generating polynomial is the product of all 90 factors

`p00+p10 z_g+p01 z_h+p11 z_g z_h`.

This factorization specifies all correlated joint outcomes, including zeros,
without enumerating an exponentially large array. Directed interval
convolution retains all 181 total-number outcomes, all 84 group marginals
and all 12 species number laws. Displayed point probabilities instead
use a single consistent f and squared matrix entry for each pair and
convolve these normalized pair laws. They are never defined as independent
midpoints of correlated count intervals or renormalized after the fact.
The executable checks nonnegativity, normalization, expectation in [0,180],
numerical agreement with the certified intervals, and interval widths
below 10^-14. The intervals enclose the exact ideal laws. The independent
binary64 point convolution is a diagnostic; its 2*10^-14 comparison
tolerance is not claimed as a proved point-rounding enclosure.
The coherent final source state and the
subsequent complete number measurement are distinct objects.

## Comparison with the continuous ramp

The finite segment law above is the executed construction. A controlled
comparison also follows for the continuous cubic ramp, with H=p(t) sigma_z
+m sigma_x. On a duration tau, the conservative derivative bounds are
`|p'|<=M1=9/(2 tau)` and `|p''|<=M2=117/(2 tau^2)`, since a>=1.
For one interval of width h=tau/N, use the interaction picture centered on
its midpoint Hamiltonian H_c. Taylor's theorem bounds the integral of
`p(t)-p(c)` by M2 h^3/24: the linear term integrates to zero.
Conjugation changes the integrand by at most
`2m |t-c| |p(t)-p(c)|`, so its integral is at most
`2m M1 h^3/12`. Its integrated norm is at most A=M1 h^2/4.
The Dyson terms of order at least two have norm at most
`exp(A)-1-A <= (3/2) A^2` for the checked A<1 regime.
Unitary telescoping therefore gives the single-pair operator bound

`delta_pair <= N*[h^3*(M2/24+2m M1/12)+(3/2)*(M1 h^2/4)^2]`.

Use the certified mass upper endpoint. For the full Fock operator, sum
this bound over all 90 pairs, retaining their multiplicities, and cap at
two. The unhalved output trace distance, including any reference, is at
most twice the uncapped sum, capped at two. A single-pair bound is never
identified with a full-source bound. There is a sharper mean-count bound without multiplying a full-source
trace bound by 180. Each local number observable satisfies 0<=N_pair<=2.
Its mean changes by at most min(2,2 delta_pair), since the pair half-trace
distance is at most delta_pair. Summing over the relevant pairs gives
rigorous total and species mean intervals around their exact finite-pulse
intervals. At 1024 segments, these intervals strictly order the vacuum
mean of the duration-eight continuous ramp below duration one, and
duration one below the quench. The massless production is exactly
zero for every declared pulse and for the continuous ramp.

The scale function and its GeV^-1 evolution parameter are supplied apparatus
data. No physical cosmological clock, Einstein backreaction, cosmological
volume, confinement crossover, interacting mass-pole matching, or observed
Universe abundance is derived. These statements concern complete
simulation-generated records on the declared finite quadratic source.
