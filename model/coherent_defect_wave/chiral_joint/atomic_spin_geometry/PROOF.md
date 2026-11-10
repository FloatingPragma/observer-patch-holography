# An isolated interacting ground and its first collective excitations

## Source and the controlled comparison

The source is the complete two-cell-per-axis Peter–Weyl construction:
eight cells with thirty fermion modes and a five-state scalar each, and
twenty-four oriented 119-state links. All twelve local Gauss generators
are imposed. The link and scalar penalties and all three Yukawa
coefficients equal one. The Hamiltonian is

`H_lambda=H0+lambda T`,

where H0 is the complete atomic Hamiltonian and T is the actual hopping
operator. The hopping coefficient is an explicitly declared construction
choice within the interval proved below. No measured particle parameter
is used. The result does not transfer to unit hopping or another size.

The [atomic reference proof](../atomic_reference_control/PROOF.md) supplies
the full physical atomic ground P0, rank 256 and energy `8 E0`, a gap
lower bound `g>0.04117` against every other physical state, and the full
hopping bound `||T||<=tau=360`. It includes link-compensated local charge
sectors. The local ground is an odd spin doublet; the eight-cell ground
space has even total fermion parity. All source states outside P0 remain
in the complementary resolvent.

For `|lambda|tau<=g/5`, the canonical invariant-graph isometry identifies
the complete continued rank-256 band with P0. It uses the same Hamiltonian
and clock. Its exact effective generator K_lambda has second-order term

`K2=sum_e [-A_edge I+J_edge S_x dot R_axis S_y]`,

where R_axis is the pi rotation about the bond axis. The scalar shift is
retained. The actual full charged-state resolvent determines A_edge and
J_edge; all edges have the same coefficients. The integral expressions
and strict rational bounds are established in the atomic reference proof.
The [charged response calculation](../charged_ground_exchange/PROOF.md)
evaluates the full pair resolvent as an exact element of the local
ground's algebraic number field. It gives
`J_edge in [41579260299095/281474976710656,
5197407537387/35184372088832]`, approximately 0.147719206819003.
In particular `J_edge>=J_lower>0`. The source is not replaced by a fitted
spin Hamiltonian.

The [even-remainder proof](../atomic_even_remainder/PROOF.md) uses the
actual side-two bipartite fermion parity to prove

`||K_lambda-8 E0 I-lambda^2 K2|| <= 25 (|lambda|tau)^4/g^3`.

This is a full complementary-space operator bound. The parity argument
fails for a nonzero even-displacement path at side four, so this estimate
is not asserted for that member.

## Removing the bond-axis rotation

At cell `x=(x1,x2,x3)` set

`O_x=R_1^x1 R_2^x2 R_3^x3`.

The three pi rotations commute as SO(3) matrices. On every side-two edge,
`O_x^T R_axis O_y=I`. Their spin-half lifts may differ by signs; these
signs disappear from conjugation of spin operators. The implementation
constructs the actual 256-dimensional tensor-product unitary from the
Pauli matrices and checks both unitarity identities and the complete
matrix intertwiner. Consequently

`K2 = U[-24 A_edge I+2 J_edge H_cube]U^dagger`,

`H_cube=sum_12_undirected_cube_edges S_x dot S_y`.

The factor two counts the two distinct oriented source links between each
neighbor pair on this periodic member. Source translation is transported
with this change of coordinates; its observable is not redefined.

## The complete cube spectrum

The exact matrix `4 H_cube` has integer entries. In a spin configuration
its diagonal is the sum of products of endpoint Pauli-z values; each
unlike neighboring pair contributes an off-diagonal entry two when
exchanged. All nine magnetization blocks are included, of dimensions
`1,8,28,56,70,56,28,8,1`, whose sum is 256.

Exact characteristic-polynomial factorization and rational root isolation
give twenty-one distinct energies, with their complete multiplicities.
Every factor, root interval and magnetic multiplicity is retained in the
receipt. The lowest energy epsilon0 is the smallest root of

`epsilon^3+5 epsilon^2-epsilon-9=0`.

It is about -4.82009. It has multiplicity one, appears only in zero
magnetization, and is a spin singlet. The first excited energy is exactly
-4, with one state in each magnetic sector -1, 0 and +1. It is therefore
a spin triplet. The full commutators with total spin are checked on all
256 states. The cube gap `delta=-4-epsilon0` is strictly greater than
4/5. The distance from the triplet to the next energy is also greater
than 4/5. These are rotated atomic-spin labels, with no Lorentz-spin
identification.

## Actual source momentum

The [atomic momentum proof](../atomic_momentum_response/PROOF.md) checks
that the graded fermion permutation and antiperiodic seam signs cancel
on the even-side atomic ground space, giving ordinary permutations of
the local spin-doublet factors. This source translation is used here.

The code embeds the complete three-dimensional polynomial kernel
associated with the cubic ground factor, applies U and checks the actual
translation matrices. All three translations act as +1 on that kernel,
and hence on its unique lowest state. The ground momentum is zero.

The three first-excited eigenvectors are obtained independently from
their complete magnetic blocks, embedded into the full space and acted
on by U. Their complete joint translation projectors have rank one at
each of

`(pi,0,0)`, `(0,pi,0)`, `(0,0,pi)`.

The other five periodic momentum outcomes have rank zero within that
triplet. All outcomes are retained. A spin-coordinate rotation cannot
be used to relabel these source momenta as zero.

## Transfer to the full interacting Hamiltonian

Let delta_star be the smaller of the certified cube ground-to-triplet
gap and triplet-to-next-level gap. Let `Delta2=2 J_lower delta_star`.
Choose the largest positive dyadic coefficient lambda_star satisfying

`lambda_star tau<=g/5`,

`25 lambda_star^4 tau^4/g^3 <= lambda_star^2 Delta2/4`.

Exact rational arithmetic makes this choice without experimental inputs.
For every `0<lambda<=lambda_star`, write
`r4=25 lambda^4 tau^4/g^3`. The min-max principle and the full operator
remainder imply a ground gap at least

`2 lambda^2 J_lower delta_lower-2 r4>0`.

The first triplet is also strictly separated from every other level in
the rank-256 band. The band itself is separated from all other physical
states by at least `g-2|lambda|tau`. Thus the full physical H_lambda has
a unique ground and an isolated first-excited cluster of rank three.
The scalar `-24 lambda^2 A_edge` is present in both absolute energies
and cancels only when taking their difference. The certificate retains
upper and lower gap bounds using the actual interval for J_edge.

Every symmetry commuting with H_lambda, H0 and P0 commutes with the
unique canonical graph construction. It therefore acts on K_lambda with
the same source representation. Dividing its centered generator by
lambda squared gives a continuous family tending to K2. The strict
isolation bounds make the ranks of its momentum projectors constant
throughout the admitted positive interval. The exact interacting ground
has zero momentum, and the first cluster has one state at each of the
three X momenta listed above. The actual cubic source rotations permute
these three one-dimensional momentum sectors and commute with H_lambda.
Their three energies are consequently equal. No unproved conservation
of the rotated continuous spin symmetry at finite lambda is required.

## Reference, preparation and physical scope

The unique physical ground projector is a stationary reference for the
same interacting Hamiltonian. A formal joint Gauss/energy filter on the
full maximally mixed carrier succeeds with probability `1/D` and fails
with probability `1-1/D`, where

`D=(5*2^30)^8*119^24`.

The Gauss and energy projectors commute. This identity states the
complete probability law of the specified spectral filter; it does not
claim an efficient native implementation or an executed full-source
histogram. The source Hilbert-space entropy ceiling is `log(D)`; it
excludes apparatus and program storage.

The result is a controlled finite interacting ground and an isolated
collective excitation sector of one declared small-hopping source.
The lowest excitations occur at nonzero source momentum. Their gap is
therefore not a zero-momentum rest-mass measurement. The three momenta
are not identified with three fermion generations. No thermodynamic,
relativistic, Standard Model or GeV identification follows from this
finite certificate.
