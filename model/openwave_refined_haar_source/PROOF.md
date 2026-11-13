# Complete interior-resolved Haar source for one spatial refinement

The source is the full SU(2) Haar Hilbert space on the54 links of the27-vertex
cubical graph with coordinates0,1,2 in each of three directions, restricted by
independent gauge rotations at every vertex. Its36 elementary Wilson faces are
normalized fundamental traces, so each is a bounded self-adjoint multiplication
operator of norm at most1. For declared positive alpha,beta, the Hamiltonian is
H=alpha E+beta(36-S), where E is the sum of link Casimirs and S is the sum of all
faces. E has its usual nonnegative self-adjoint Peter–Weyl realization, and the
magnetic perturbation is bounded and nonnegative. No representation is removed.

## Full arbitrary-valence Gauss basis

Orient every link from the even to the odd sublattice. A spin-network label
contains every link twice-spin k=0,1,... and every local left-associated fusion
path ending in spin0. At a vertex the invariant tensor is the product of
normalized Clebsch–Gordan coefficients along this path. Successive CG unitarity
proves that these tensors form a complete orthonormal basis of the invariant
subspace of the tensor product, at arbitrary valence. Haar orthogonality on each
link, with sqrt(2j+1) edge normalization, then proves orthonormal completeness of
the graph basis. Degree4 and6 fundamental tensors have respectively2 and5
singlets; the executable checks their entire exact Gram matrices and also
nontrivial degree5 assignments. A vertex degree is not an assumed intertwiner
multiplicity.

Multiplying a spin-j matrix by a fundamental matrix uses the exact two-factor
CG product identity. The output spin on each traversed edge is j±1/2. At each
traversed vertex, `vertex_contraction` contracts the old and new normalized
invariants with the two edge CG coefficients and the fundamental epsilon.
Summing both fundamental magnetic indices gives the recorded expression. The
edge normalization ratio is sqrt((2j+1)/(2j'+1)), and normalized trace contributes
1/2. The code enumerates every permitted output fusion path, including several
paths with the same output link spins. This is a full multiplication action,
not a projected source rule. The independent delta/epsilon Haar oracle checks
every compressed face column, including its sign relative to the trace basis.

## Complete201-dimensional proof frame and complement floor

For any occupied edge with j>=1, existence of an invariant at each endpoint
requires the sum of its other incident spins to be at least j. The two sets of
other incident edges are disjoint in this simple graph. For every half-integer
k>=0, k(k+1)>=3k/2. Consequently the edge and these neighbors contribute at least
j(j+1)+3j>=5 to E. This is a lower bound, not an assertion that energy5 occurs.

Below5 every occupied edge therefore has j=1/2. An invariant then requires even
support degree at each vertex. Each occupied edge costs3/4, hence there are at
most six. On a simple bipartite graph of girth4, an even-degree subgraph with at
most six edges is a single four- or six-cycle: two cycles would need at least
eight edges, and a higher even-degree vertex would require such a decomposition.
Every vertex on a single cycle has its unique two-leg invariant. There are36
four-cycles and164 six-cycles, enumerated by complete simple-path traversal with
canonical rotation/reversal. Together with the vacuum they give all201 states
through E=9/2. Every omitted state has E>=5. Thus the full complement obeys
QHQ>=5 alpha Q; this bound includes all higher representations and intertwiners.

All36 times201 face columns retain their entire output lists, including outputs
outside this frame. Exact character integrations provide independent matrix
elements; complete CG coefficients provide all residual channels. The second
and fourth Haar moments are <S²>=9 and <S⁴>=963/4. The executable full second
action gives ||P S² Omega||²=357/4 and ||Q S² Omega||²=303/2. These are exact
identities for the Haar-constant probe, not interacting vacuum assertions.

## Complete readout and work

For every36-bit sign word s, take E_s=product_f(1+s_f W_f/2)/2 and K_s=sqrt(E_s).
The faces commute as multiplication operators. Every effect is positive,
independent sign sums give sum_s E_s=I, and every probability lies between4^-36
and(3/4)^36. These formulas specify all2^36 outcomes without omitting any record.
The nonselective channel retains arbitrary external quantum references.
For a length-L Wilson loop the electric carré du champ is L(1-W²)/4. The
binary square-root instrument has sum_s Gamma(K_s)=eta² Gamma(W)/[4(1-eta²W²)].
At eta=1/2 this is bounded by L/48. Summing36 faces of length4 yields work at
most3 alpha above an interacting ground state. The magnetic term commutes with
all these effects. This is the existing finite readout/work theorem on the
larger source, with no fitted time or mass conversion.

The source supplies genuine spatial interior modes and arbitrary-valence
intertwiners. It does not establish a repeated-refinement limit, rotation/Lorentz
invariance, continuum Yang–Mills identification, a physical GeV scale, or the
Clay mass-gap conclusion. The target remains interacting four-dimensional
quantum Yang–Mills, not an exactly scale-invariant flat vacuum with a positive
mass gap; see the [official Jaffe–Witten specification](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).
