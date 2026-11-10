# A quantum scalar cell with all native quark and lepton channels

## Full action and the specified preparation sector

The carrier is `F(C^2 tensor R15) tensor C^5`, of dimension `5*2^30`.
Every fermion occupation and every scalar-memory state exists in the
Hamiltonian domain. The twelve native currents, five irreducible matter
blocks and invariant pair tensors are exactly those of
[`physical_attachment/chiral.py`](../../physical_attachment/chiral.py).
The scalar matrices are the same five-state quantum fields used by
[`chiral_higgs`](../../chiral_higgs/PROOF.md). All three Yukawa coefficients
and the scalar number penalty are one. There is one cell and no spatial
hopping. These coefficients and geometry are declared architecture.

The Hamiltonian is `H=N_H+H_up+H_down+H_e`. The executable compiles its actual
pair words from `periodic_family.yukawa_words`; no quark occupation or
interaction is set to zero. Its integer CAR action is defined on every
full-carrier column. The all-occupation Ward identities and Hermiticity
follow from these same invariant pair tensors and their adjoints.

Start with the literal empty fermion and scalar state Omega. Enumerate
its complete connected computational-basis graph under H, continuing until
no outgoing matrix element exists. The resulting reducing space has
97878 states and 933620 nonzero matrix entries, including diagonal entries.
Its maximum graph distance from Omega is thirteen. Every outgoing column
is retained. This is a complete invariant preparation sector; it is not a
claim that this graph exhausts all conserved sectors of the full carrier.

## Exact ground certificate without a full numerical diagonalization

A diagonal matrix D of signs switches every nonzero off-diagonal H entry
to a negative real number. Exact traversal checks each edge, including
cycles. The switched matrix is irreducible because the enumerated graph
is connected. Therefore `I-DHD` is irreducible and nonnegative; its Perron
eigenvector is strictly positive and its largest eigenvalue is simple.
Equivalently H has a unique ground in this reducing space.

An exact equitable partition compresses this graph to 240 classes. For
class b, let the full vector B_b have the switching sign on every vertex
of that class and zero elsewhere. Every full row is checked against
`H B=B Q`; `B^dagger B=G` is the positive diagonal matrix of class sizes,
and `G Q=Q^dagger G`. The class containing Omega is a singleton. A stored
partition is only a proposed classification: replay verifies its complete
intertwining identity and source hash with integer arithmetic.

There are 64 classes with scalar number zero and 176 with scalar number
one. In this order, Q has block form `[[0,A],[B,I]]`; its squared-minus-
itself action on the blank-scalar block is `S=A B`. The exact Krylov
vectors `Omega,S Omega,...,S^17 Omega` are independent, while `S^18 Omega`
has an exact rational expression in them. The code checks this rank and
closure, producing the degree-eighteen minimal polynomial

`p(x)=(x-4)(x^2-19x+52)`

` * (x^3-313x^2+31738x-1041088)`

` * (x^3-229x^2+16786x-392664)`

` * (x^3-157x^2+7714x-117672)`

` * (x^3-97x^2+2794x-23152)`

` * (x^3-49x^2+586x-1104)`.

All eighteen roots are simple and positive, verified by exact rational
root isolation. For each root lambda, the actual H energies are
`E_±=(1±sqrt(1+4lambda))/2`. The resulting physical H-cyclic space has
dimension 36. Its basis vectors are polynomials of the full H applied to
Omega, so every vector satisfies every Gauss generator. This is an exact
invariant compression, not an omission of an H leakage term.

The connected-graph ground has nonzero overlap with Omega by strict
positivity. It must therefore occur in this cyclic spectrum. Its energy
comes from the largest root of the cubic with coefficient 313 and is
approximately `-11.14102996` source units. Exact rational root isolation
and outward square-root bounds enclose it in
`[-1531211498821/137438953472,-12249691990567/1099511627776]`.
No numerical eigensolver or rounded spectral gap enters this certificate.
The conclusion concerns the full connected preparation graph; other
conserved sectors of the complete carrier are not compared.

The ground is annihilated by every actual Gauss generator. Indeed Omega is
a gauge singlet and every spectral projection of the full H commutes with
Gauss. The ground projection of Omega is nonzero by strict positivity.
The reducing graph has a unique ground, so that projection is proportional
to this ground vector, even if another graph shares its eigenvalue.
This argument does not assert that arbitrary coordinate vectors of the
graph are physical singlets.

## Complete quantum energy probabilities

Let `m_k=<Omega|S^k|Omega>` and write `p(x)=sum_j p_j x^j`.
The vacuum resolvent has the exact form `q(z)/p(z)`, where

`q_k=sum_(j=k+1)^18 p_j m_(j-k-1)`.

The numerator has degree seventeen, leading coefficient one, and no
common factor with p. For a root lambda the S spectral weight is
`w_lambda=q(lambda)/p'(lambda)`. The H energy weights are

`w_±=w_lambda*(1∓1/sqrt(1+4lambda))/2`.

These follow by diagonalizing each two-dimensional scalar-number block;
the metric is retained in the moment construction. Exact interval
arithmetic certifies all 36 positive energy weights, with outward energy
grid `2^-40` and probability grid `2^-48`. Their sum is one by the finite
spectral resolution and the unit leading resolvent coefficient. The full
instrument also retains the complementary cyclic-code outcome, whose
probability on this input is exactly zero. The cyclic projector commutes
with H and every Gauss generator; all unnormalized output states exist on
the original carrier. No branch is renormalized.

The ground filter succeeds with probability in
`[67165484451/140737488355328,134330968903/281474976710656]`, approximately
`0.00047723947`. Its failure is the complement, including all other energy
outcomes. This finite spectral controller is a declared service. No
efficient local implementation or physical cooling rate is asserted.
The initial Omega preparation and the filtered interacting reference are
kept distinct.

## Gauge-compatible production records on the full Hamiltonian

For a separate experiment, evolve the original Omega with the same H for
source time `1/1000`. The initial state is not replaced by the ground.
A complete detector records irreducible-block occupation numbers and scalar
occupation. Its three distinguished outcomes are one scalar excitation and
respectively `(n_uc,n_Q,n_ec,n_L,n_dc)` equal to

- `(1,1,0,0,0)`, the up-pair channel;
- `(0,1,0,0,1)`, the down-pair channel;
- `(0,0,1,1,0)`, the charged-lepton-pair channel.

Each projector retains every color, weak and spin state in that occupation
sector. Their full ranks are 288, 288 and 32. They commute with every gauge
current: the current matrices preserve the native irreducible blocks, and
the scalar currents preserve scalar number. Individual color-coordinate
projectors would not supply this property. The fourth outcome is the
orthogonal complement of their sum, including the unconverted input.
The apparatus labels quantum-number channels; it identifies no measured
particle or lepton generation.

Exact integer applications give H^n Omega through n=5. The degree-four
Taylor vector has full-source remainder bounded by

`delta <= |t|^5 ||H^5 Omega||/5!`,

using the unit norm of the exact unitary in the integral remainder.
Here `||H^5 Omega||^2=409247068`; the square root is bounded outward with
integer arithmetic. For each effect the probability error is at most
`delta*(2+delta)`. The receipt encloses the three positive channel weights
and retains the complementary probability with their combined error.
No output is renormalized, and no graph-truncated evolution is substituted
for the full-source remainder. The leading channel weights are of order
`12t^2`, `12t^2` and `4t^2`; these are production records on this specified
input, without a rest-mass interpretation.

## Physical boundary

All quark and lepton Yukawa channels act together on the actual quantum
scalar memory. The complete preparation-sector ground and controlled
channel records are finite-source results. They do not identify the global
vacuum across other sectors, supply spatial propagation, establish a
relativistic field, derive three generations, or convert an excitation to
GeV. The source parameters have no empirical mass inputs.
