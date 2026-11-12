# Complete native low nuclear bands on finite graphs

## Actual local representation census

Use the identical positive current action and calibrated coefficients of
`model/openwave_calibrated_hadrons` on a finite graph. Its edge term is
`eta sum_(xy)(I_x+I_y)^2` and its kinetic term is the positive graph
Laplacian with hopping t. Graphs are explicit construction inputs. No
coefficient is refitted when the graph changes. Only colourful Q modes
are occupied in the N=6 experiment; every other native Fock factor remains
an invariant complement of the source.

Each cell has twelve Q coordinates: three colours, two spins and two weak
components. Exact native integer Casimir polynomials compute the complete
colour-singlet projectors in the220-word three-quark and924-word six-quark
spaces. On zero Cartan weight these use only64 and216 words. Exact projector
identities and magnetic-weight differences give the full local census:

* Three quarks: `(S,I)=(1/2,1/2),(3/2,3/2)`, each once, total dimension20.
* Six quarks: `(1,0),(0,1),(3,0),(0,3),(2,1),(1,2)`, each once, total dimension50.

For a character weight w(S,I), the highest-weight multiplicity is exactly
`w(S,I)-w(S+1,I)-w(S,I+1)+w(S+1,I+1)`. Thus these counts are derived from
the actual CAR representation, not supplied particle multiplicities.

## Complete locally colour-singlet subspace

Let P project onto states with zero colour current at every cell. At total
quark number six, local colour triality permits only a six-quark cell or two
three-quark cells. On a graph with L vertices, therefore,

`dim(P)=50L+400 binomial(L,2)`.

Every pair of local three-quark multiplets is coupled with the complete
finite spin and weak representations. The implementation enumerates all
allowed total S and I, all their magnetic multiplicities, and every cell
placement. The displayed dimensions exhaust P. In particular its S=1,I=0
space at fixed magnetic projections has dimension `L+2 binomial(L,2)=L^2`.
The S=0,I=0 space has dimension `2 binomial(L,2)`.

A single quark hop leaves P. Hence P H P contains the diagonal kinetic
degree cost and the actual local and edge current energies, with no cluster
hopping term inserted. A colocated six-quark multiplet at vertex x has
internal energy

`4A+Js S(S+1)+Jw I(I+1)+[eta I(I+1)+6t] degree(x)`.

For three-quark multiplets of local spin/isospin a and b at x<y, the
internal energy is

`2A+(Js+Jw)[a(a+1)+b(b+1)]+3t[degree(x)+degree(y)]`
`+eta[a(a+1)degree(x)+b(b+1)degree(y)`
`     +adjacent(x,y){I(I+1)-a(a+1)-b(b+1)}]`.

The bare mass adds `6mbar-delta Iz`. These formulas are actual compressions
of native current squares. Local Casimirs and the complete angular coupling
ensure that P H P is diagonal in the stated basis. They do not define a
projector Hamiltonian or discard the quark hopping complement.

## Full-source spectral certificate

Within the complete global-colour-singlet N=6 sector, any vector in Q=1-P
has at least two nontrivial local colour representations. Their Casimirs
sum to at least8/3. Every other interaction is positive. If r_L is the
minimum of `sum_x n_x^2` over nonnegative integer occupations summing to six,
then, after the common bare mass shift,

`Q H Q >= C=8U/3+A r_L/9`.

Only hopping joins P to Q. If z is the maximum graph degree,
`||P H Q||<=b=6tz`. Let a be the largest eigenvalue of P H P and suppose
C>a. With `epsilon=b^2/(C-a)`, completing a block square gives

`H >= (P H P-epsilon) direct_sum a I_Q`.

Min-max comparison with the P trial space then encloses every one of its
d=dim(P) low eigenvalues:

`lambda_j(P H P)-epsilon <= lambda_j(H) <= lambda_j(P H P)`.

The codimension-d Q trial space also gives `lambda_(d+1)(H)>=C`.
Thus the count and bounds cover the complete low spectrum below an explicit
coloured-complement threshold. They do not enumerate the high spectrum or
assert a result for a different global-colour sector.

For a three-quark source the analogous formulas use
`C_B=8U/3+A min(sum n_x^2)/9`, `b_B=3tz`. The full local20L band is retained,
including both spin/isospin multiplets. Its same-graph nucleon thresholds
are used for each nuclear comparison. The fitted up/down scalar shifts
cancel in a fixed-charge breakup. No threshold from a different graph or
Hamiltonian is substituted.

## Complete preparation, records and references

The local colour projectors are exact finite Casimir polynomials. Total
spin, weak and magnetic projectors can likewise be written as finite
polynomials in their native commuting generators. They therefore give an
explicit P for each stated sector without enumerating the full wavefunction.
Prepare the uniform full N=6 occupation input and apply this filter. For a
fixed magnetic sector of rank d its success is
`d/binomial(12L,6)`. The entire complementary state and every other Fock
sector have retained failure outputs. Opposite native charge is stored in
an initially correlated idle reference register, which is unchanged.

Let P_low be the exact d-dimensional low spectral projector. Its orthogonal
complement has the actual energy floor C. The off-block Sylvester equation
and the separation `C-a` give

`||P-P_low||<=b/(C-a)`.

The complete low/high spectral instrument on the prepared P density thus
has high probability at most `[b/(C-a)]^2`. The receipt gives original
success-low, success-high and preparation-failure bounds, without conditional
renormalization. This comparison is uniform in time and under arbitrary
external quantum references. Spectral readout commutes with H and does zero
mean work; nonselective preparation on the uniform input also does zero mean
work. Exact census and rational bounds are executed. Full large-carrier
propagation is not claimed.

## Finite graph results and physical boundaries

Paths of two, three and four cells, and cycles of three and four cells, use
one unchanged calibration. The spin-one/isospin-zero channel is bound on
each graph. Their complete excited bands and finite graph effects are
retained. The spin-zero/isospin-zero channel is also bound: its local
nucleon-pair diagonal energy equals the triplet's and its rigorous binding
margin stays positive. It is a different native state from the lifted
spin-zero/isospin-one exchange partner. The complete nuclear catalogue has
not been physically recovered. Without a physical orbital and spacetime
dictionary, this additional channel is neither a proved new particle nor a
proof that the source contains a forbidden physical nucleus.

The exact local three-quark nucleon representation has multiplicity one.
The colour, spin, weak and magnetic Casimir projectors isolate its exact
rank-one spin-up neutron line. Direct integer CAR multiplication gives
`Tr(P_neutron I_minus J_plus,z)=5/3`, retaining the actual native axial
quark current. This supplies the numerical reduced matrix element as well
as the representation multiplicity.
By Schur's lemma any colour-, spin- and weak-scalar deformation confined to
that local irrep acts as a scalar. It cannot change the native axial reduced
matrix element5/3 by altering its state. Spatial, sea, orbital or current
modifications require their own calculation. This representation statement
does not impose uniqueness of the overall architecture.

Graph variation supplies finite spectral error bounds and counterexamples.
It establishes no continuum limit, physical graph length, covariant nuclear
spin, observed excited-nucleus spectrum or universal physical clock.
