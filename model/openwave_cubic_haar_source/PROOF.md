# Complete interacting Haar source on a spatial cube

The graph is the elementary three-dimensional cube, with eight vertices,
twelve links and six plaquettes. The source is the existing pure SU(2)
Hamiltonian `H=alpha sum C_e²+beta sum_f(1-W_f)`, with
`alpha=g²/(2a)`, `beta=2/(g²a)` and `W_f=Tr(U_f)/2`.
Its Hilbert space includes every Haar representation and all local Gauss
constraints. The finite matrices below are proof coordinates only.

## Complete Gauss census and exact multiplication

Every cube vertex is trivalent. A spin assignment has one invariant tensor
precisely when its three doubled spins satisfy the triangle inequalities
and even total parity. That invariant is unique. The normalized network
uses `sqrt(2j+1) D_j(U)` on each link and the unit 3j tensor at each vertex.
Peter--Weyl orthogonality makes these states an orthonormal complete basis.

The census enumerates every nonnegative doubled edge spin with total electric
energy at most the declared threshold, checking every vertex. The energy
is `sum k(k+2)/4`; an edge above the enumerated maximum alone exceeds the
threshold. Vertex-by-vertex recursion retains every admissible assignment.
Through energy six there are exactly 32 states: the vacuum, six four-cycles,
sixteen six-cycles, six eight-cycles, and three pairs of opposite four-cycles.
All nonzero edge spins there are one half. The complete census through
energy `13/2` adds twelve states at `13/2` and none at `25/4`.

An independent low-energy representation uses products of fundamental
characters on those disjoint cycles. Integrating the at-most-three matrix
occurrences on every edge gives the identity Gram matrix and all six
Wilson matrices exactly. Odd edge multiplicity vanishes under the SU(2)
central sign. The remaining integrals use the exact two-matrix delta or
epsilon pairing, including its signs. The 3j and character bases agree
after the displayed, fully determined diagonal phases.

For arbitrary spins, orient each edge from the even to the odd cube vertex.
Multiplication by a fundamental Wilson trace changes each of its four
edge spins by plus or minus one half. Matrix-coefficient Haar integration
at each edge gives two Clebsch--Gordan factors divided by the output
dimension. With normalized network edges this leaves
`sqrt(d_input/d_output)` per affected edge. At each vertex contract the
input and output 3j tensors with the two Clebsch--Gordan factors. An inverse
fundamental matrix is converted to its epsilon dual; summing its two
indices gives the signed local contraction in `vertex_contraction`.
The product of the four vertex contractions and edge factors, divided
by two for the normalized trace, is `cg_entry`. No undetermined recoupling
sign or fitted matrix entry remains. Complete successor expansions check
all fifteen Wilson commutators on every state through energy eight,
including every intermediate and final representation outside that set.

## Full-source energy and state certificates

For the complete projection through electric energy eight, the entire
complement has electric energy at least `33/4`; this conservative floor
uses the quarter-integer spectrum. Magnetic positivity gives
`QHQ >= (33/4) alpha Q`. The off-block norm is at most `6 beta`.
The complete 86-coordinate frame is revalidated by the existing exact
Gram/residual enclosure. The same full-complement Schur theorem as in
[`openwave_coupled_haar_spectrum`](../openwave_coupled_haar_spectrum/PROOF.md)
therefore bounds the first full excitation, without dropping any sector.

The ground trial is a fixed dyadic vector supported on all 779 Gauss
coordinates through electric energy fourteen and the six spin-three-halves
square loops at energy fifteen. These 785 trial coordinates need not form
a spectral cutoff: their role is solely to specify an explicit normalized
trial vector. The separate 86-state spectral frame is complete. Every Wilson successor of
every trial coordinate is emitted, including outside this support.
Exact algebraic entry intervals precede dyadic rounding. For each
rectangular action matrix, the square root of the maximum row-error sum
times the maximum column-error sum bounds its operator error. Thus the
native full-source residual is enclosed without a separate truncated-tail
approximation. Temple's inequality and the independently certified first
excitation give the full ground interval and projector distance.

The same complete successor vectors give all six means, all six variances
and all fifteen mixed Wilson moments. Their overlaps retain the entire
outside component. Algebraic rounding and the full ground-projector error
are propagated before computing covariance. Adjacent and opposite faces
are kept separately. All spatial permutations and sign reversals preserve
the Hamiltonian; its unique positive ground is invariant, so symmetry-
related exact moments agree. The interval calculations check their overlaps.
The adjacent-face covariance has a positive lower bound at declared
`g²=8` and `g²=16`. Other reported intervals retain their actual uncertainty.

## Complete six-face record and its work

The 64 effects `product_f(I+s_f W_f/2)/64`, with every sign choice, commute,
sum to identity, and lie between `1/4096` and `729/4096` times identity.
Their positive square roots define a complete full-Haar instrument, also
with an arbitrary idle external quantum reference. The reported probability
bounds hold for every input; the state-specific pair-record covariance is
one quarter of the certified Wilson covariance. No outcome is postselected
away. Summing the six single-loop gradient budgets gives nonselective
ground-state work between zero and `alpha/2`. Cross terms cancel because
the complete single-loop Kraus derivative satisfies `sum f_s f_s'=0`.

The graph is three-dimensional but finite. Its full Haar spectrum and
connected records do not supply a uniform-volume gap, a physical GeV mass,
an interacting four-dimensional continuum, or the Clay construction.
