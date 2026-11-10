# Complete Gauss kernels of the two-cell control

This calculation supplies the full physical gauge subspace of a specified
finite control. It does not replace a charged source by a source-only Gauss
state, restrict to a vacuum cyclic orbit, or use a numerical rank threshold.

## Carrier and exact currents

The active carrier has twelve lepton CAR modes, the six actual link states
and two five-state quantum scalar memories, of dimension
`2^12 * 6 * 5^2 = 614400`. A key is `(mask,link,hleft,hright)`.
The low six mask bits are the left cell, and the next six the right cell.
At each cell the spin-major internal order is native indices `(9,10,11)`:
e^c and the two L components, repeated for the second spin. The actual link
order is the blank, e^c link, and four weak-doublet matrix-unit states, using
the parent construction's `active_terms` and `link_labels` functions.

The full R15 source is retained as the owning carrier. This is its reducing
control with both quark Yukawa couplings zero, empty quark occupations and
unexcited colored links. Quark hopping cannot leave that spectator vacuum.
Both local charged-lepton Yukawa couplings are one, and both local quantum
scalar factors remain active. The Hamiltonian is constructed separately;
these gauge kernels contain every occupation in the admitted active control.

For each of the twelve native Hermitian generators, the local Gauss current
is its second-quantized matter current, its actual five-state scalar current
and the oriented endpoint current of the actual six-state link restriction.
The eight color currents vanish identically on this reducing control, which
the code checks on every factor. The weak currents are the actual matrices
at indices 8, 9 and 10; hypercharge is at index 11. Matter, scalar and both
link endpoint matrices satisfy the exact SU(2) commutator identities.

Both boundaries may carry fixed weak/color-singlet compensating charges
`(q_left,q_right)`. Their contributions to the respective hypercharge
currents are `q_site I`. The API accepts any exact integer pair, with an
integer q preserving the shorthand `(0,q)`. The charged-probe cases use
`q=0,-1,+1` and the reflected `(1,0)` sector.
Each charge is a retained apparatus label, not an additional particle species.
The two nonzero signs permit complementary charged source probes. Fixing
the boundary occupation contributes no further dimension within each
individual boundary sector.

## Complete diagonal sectors

The code enumerates all 64 local fermion occupations and all five scalar
alternatives, and joins every pair through every one of the six link states.
It selects local total `Y=0` and `T3=0`, with the boundary included at the
right endpoint. Multiplying the exact Y eigenvalues by six and T3 eigenvalues
by two makes the arithmetic integral; no rounded weights are used. Both
boundary hypercharges enter before the diagonal sectors are enumerated.

The conserved control labels are `L=N_L-N_ec` and `2Sz`, the sum of the two
source spin occupations with opposite signs. Local weak generators preserve
both. In the specified one-z-edge Hamiltonian, the diagonal Pauli hopping
and spin-singlet Yukawa terms preserve them as well. This is a finite
source spin projection; no relativistic spin classification is inferred.

| Boundary q | Complete active carrier | Diagonal Gauss dimension | Full physical dimension | Largest physical (L,2Sz) block |
|---|---:|---:|---:|---:|
| 0 | 614400 | 4140 | 1053 | 153 |
| -1 | 614400 | 2496 | 718 | 96 |
| +1 | 614400 | 3358 | 783 | 103 |

Here q denotes the right-only shorthand. The reflected `(q_left,q_right)=(1,0)`
case has diagonal dimension 2496 and physical dimension 718. Its physical
block dimensions equal those of `(0,-1)` under `(L,2Sz)->(4-L,-2Sz)`, checked
by the same exact kernels and independent characters. This equality of
dimensions alone does not prove a Hamiltonian symmetry; that requires the
separate operator intertwiner.

## Raising maps give the entire weak kernel

On each cell form the actual raising operator `Jplus=G8+i G9` and retain
its image in the adjacent `T3=+1` weight sector. For each fixed `(L,2Sz)`,
stack the left and right raising maps into a rectangular matrix M. A row
records both the endpoint and its full target key. Its codomain is not
compressed into the original zero-weight basis. Such compression would
erase every raising image and falsely accept nonsinglets.

On a zero-weight vector v, the identity `[Jplus,Jminus]=2Jz` gives
`||Jminus v||^2-||Jplus v||^2=2<v,Jz v>=0`.
Therefore `Jplus v=0` also gives `Jminus v=0`. Conversely, every weak singlet
is killed by both raising maps. Since both diagonal local generators were
already imposed, `ker M` is precisely the common kernel of all six local
weak currents. Combining it with the exact local hypercharge constraints
and the identically zero color currents proves that it is the entire
physical Gauss kernel of the stated active control and boundary sector.

All current actions use ordinary CAR bilinears, including their exact
occupation signs. Link and scalar actions are their actual matrix columns.
No singlet basis, Clebsch coefficient or target particle spectrum is supplied
to the kernel calculation.

## Exact rank, basis and coordinate map

Two domain columns belong to the same constraint component if a raising
row couples them. This bipartite decomposition partitions M into disjoint
rectangular blocks, after row and column permutations. In this carrier each
component has at most 36 domain columns. The algorithm computes rational
RREF of every complete component, including constrained one-column and
unconstrained zero-row components.

Every free variable supplies one rational null vector. Their assembled
columns form B in the lexicographic `diagonal_basis` order. The code verifies
`M B=0`, `rank(M)+columns(B)=columns(M)`, and that selecting the recorded
free-variable rows of B yields the identity. Thus B spans the complete
kernel and has independent columns. Its Gram matrix `G=B^dagger B` is
strictly positive definite. Selecting those same rows is an exact coordinate
left inverse on the physical image; it is not an orthogonal projector.

For an operator H preserving that image, its exact coordinate matrix can be
formed as `A=E H B`, where E selects the free-variable rows. Its physical
Hermiticity condition is `G A=A^dagger G`. Equivalently, use the generalized
Hermitian pair `(B^dagger H B,G)`. Treating the unnormalized B columns as
orthonormal would produce incorrect energies.

## Independent completeness count

There is a second count which does not use matrix ranks. In a finite unitary
SU(2) representation, the multiplicity of the singlet equals
`dim(T3=0)-dim(T3=1)`. Each nontrivial integer-spin irrep contributes one
state to both weights; half-integer irreps contribute to neither. Applying
this identity independently to the two commuting endpoint representations
gives

`m00=w(0,0)-w(1,0)-w(0,1)+w(1,1)`.

Each weight count retains the same local hypercharge, boundary, L and Sz
labels. The code independently enumerates all four weight spaces from the
actual current diagonals and verifies agreement with the exact nullity of
every block, for all three boundary charges. Tests additionally apply all
twelve actual local currents directly to every basis column in the largest
block of each boundary sector, including the lowering components and the
boundary sign. Removing a nonzero compensator fails the source-only
hypercharge constraint, as expected.

These results establish a complete finite Gauss domain for the subsequent
same-H reference and response calculation. They establish neither a ground
state nor a relativistic particle, mass, or continuum limit.
