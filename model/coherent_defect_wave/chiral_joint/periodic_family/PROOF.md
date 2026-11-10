# A full periodic three-axis quantum source

This declared family joins the actual native R15 fields, quantum links and
quantum scalar memories on one periodic cubic carrier. It supplies three
commuting source translations, one Hamiltonian clock and complete physical
translation-phase records. Its interacting ground and relativistic field
limit are separate requirements. No open-source or free-field ground is
imported into the construction.

## Complete carrier and source units

For even N at least two, let the cells be `(Z/NZ)^3`. Each cell has the full
`F(C2 tensor R15)` space of thirty fermion modes and the five-state quantum
Higgs memory. Every positively oriented coordinate edge `(x,x+e_i)` carries
the actual sixty-state quantum link. There are `N^3` cells and `3N^3` links,
with complete logical dimension `2^(30N^3) 60^(3N^3) 5^(N^3)`.
All occupations and link/scalar alternatives are retained. The
[sitewise encoding](../../chiral_source/PROOF.md) supplies their native
current attachment and graded support. Its binding and fermionic grading
remain declared architecture, without a claim that native idle selects them.

Choose a positive exact spacing a, length `L=N a`, and `hbar=1`. The default
volume family holds `a=pi` fixed. At `N=2`, `L=2pi`, the antiperiodic spectral
momenta are `-1/2,+1/2`; no factor of pi or time rescaling is discarded.
A cutoff sequence at fixed L is a different explicitly stated choice
`a=L/N`. The chosen real Yukawa coefficients and nonnegative scalar/link
penalties are uniform over the cells. Default values are all one.

## Actual periodic path Hamiltonian

Let P be the exact antiperiodic spectral momentum matrix from the
[spectral regulator](../../chiral_regulator/PROOF.md). For every anchor x,
axis i and forward displacement `r=1,...,N-1`, write `y=x+r e_i mod N`.
Use the complete forward oriented path of r edges, including the wrap edge
when present. In each actual native irreducible block define

`W_xy,ab=sum_(c1,...,c_(r-1)) U^1_(a,c1) ... U^r_(c_(r-1),b)`.

Every primitive factor is the actual matrix `U_ab=|a,b><blank|` on its own
quantum link. No factor is replaced by a number or normalized. The reverse
operator in the adjoint term is `W_xy^dagger`, not a separately assumed
unitary inverse. Define

`H_kin=(1/2) sum_(x,i,r,spins,a,b) c_x^dagger [-sigma_i P(x_i,y_i)] W_xy,ab c_y + h.c.`.

The half is necessary because all ordered anchors and displacements are
present. Both orientations are distinct path operators. Under a separate
classical identity-transporter replacement, the ordered terms plus their
adjoints reproduce exactly the full free spectral matrix. A quantum blank
link is not that identity background.

The coefficient before sigma_i in the displayed term is exactly

`(-1)^wrap [-i pi/(2L)] (-1)^r / sin(pi r/N)`,

where `wrap=floor((x_i+r)/N)`. The code checks this formula against the
existing spectral matrix. At N=2 and a=pi, the two oriented coefficients
are `+i/4` and `-i/4`. Their one-axis lepton restriction matches all 240
CAR-bilinear/two-link matrix-unit words of the
[frozen periodic control](../two_cell/SYMMETRIES.md). All 320 local
lepton/scalar columns also match its actual Yukawa operator. Together with
the identical number penalties, this proves the stated control equality
on every occupation, rather than assuming a one-particle truncation.

The full Hamiltonian adds all three native quantum Yukawa channels at every
cell, the scalar number penalties and the excited-link penalties. Its
Hermitian adjoints are included explicitly. The sparse `column` function
is defined on every complete carrier key. A key comprises a full `30N^3`
fermion mask, all `3N^3` link labels and all `N^3` scalar labels. Local pair
operators have even CAR parity, so external site-order signs cancel; global
hopping signs are evaluated with the full occupation mask.

## Local Gauss constraints

At each cell the local Gauss generator is its actual native matter current,
its quantum scalar current, the dual endpoint currents of its outgoing links
and the direct endpoint currents of its incoming links. Each oriented link
has commuting left and right actions. On a path, the adjacent endpoint
actions cancel at every internal vertex, and the remaining endpoint actions
cancel the matter creation/annihilation currents. This proof applies to wrap
edges exactly as to interior edges.

The [common-source path proof](../PROOF.md) checks the actual native
one- and two-edge tensors and proves the arbitrary-path cancellation.
The [quantum Higgs proof](../../chiral_higgs/PROOF.md) supplies the full local
Yukawa Ward identities. Scalar number and excited-link projectors commute
with their currents. Thus every local Gauss generator commutes with the
entire H on all occupations. The empty matter, blank links and scalar
vacuum give a nonzero physical state, without being asserted stationary.

## Three exact translations on the same source

For direction j, the one-body fermion translation shifts each cell forward
and supplies a minus sign exactly when its j coordinate wraps. Scalar
memories and oriented links undergo the corresponding periodic permutation.
The exterior lift includes every occupied-mode permutation sign. No source
occupation is removed to obtain the translation algebra.

The three one-body shifts commute, hence their exterior lifts commute.
Their link and scalar permutations commute as well. On the full source,

`T_i T_j=T_j T_i`, `T_i^N=(-1)^F`, `T_i^(2N)=I`.

Translation moves each forward path to the path with the translated anchor.
Its antiperiodic coefficient transforms by the product of its two fermion
endpoint signs. Each local Yukawa term has two such signs and is unchanged
apart from its cell label. Number penalties merely permute. Therefore
`T_i H=H T_i` for all three directions. The code checks every path word at
N=2 and N=4, with 72 and 1728 exact path/coefficient identities, and all
one-body translation commutators. The exterior-lift and word arguments
establish the full-occupation statement.

Translations permute the local Gauss generators, so they preserve the common
physical kernel. They need not commute with each individual local generator.
This distinction is relevant to their complete physical phase measurements.
The finite polynomial-conservation theorem in
[`PeriodicMomentum.lean`](../../../../Lean/QFT/PeriodicMomentum.lean)
uses precisely the supplied `Commute H T_i` hypotheses. The source also
supplies the pairwise commuting T_i needed for simultaneous polynomial
phase effects. These facts do not identify a particle momentum or supply
rotations, boosts or a relativistic dispersion.

## Executed physical three-dimensional control

At N=2 use a local charged-lepton Yukawa pulse on the empty physical source.
The exact 320-column local equality above gives `K^2 Omega=4 Omega` for
that local pair operator K. Thus `exp(-i pi K/4)` prepares its normalized
pair/scalar state, up to global phase, with four basis entries and norm one.
All eight cells nevertheless retain every native field and every Yukawa
channel during subsequent evolution. Acting with the full H gives 849 basis
entries and exact norm squared `3321/16`, including quark creation. This is
not an evolution projected into a lepton-pair subspace.

All 96 local Gauss generators annihilate both the prepared state and H times
that state. The three actual translated H-state columns equal the H-action
on the translated preparations. Exact reverse coefficients are checked on
all 852 transitions leaving the four prepared basis states. These controls
supplement the complete Hermitian and translation word proofs.

The three commuting finite translation-phase projectors are

`Q_(k1,k2,k3)=product_i [(1/(2N)) sum_(m=0)^(2N-1) exp(-i pi k_i m/N) T_i^m]`.

The finite cyclic-group identities prove orthogonal complete effects on the
entire source. All 64 N=2 phase triples are evaluated on the stated physical
preparation. Eight all-even triples have probability 1/8; the other 56 have
probability zero. Every nonzero branch is checked against all 96 actual
local Gauss currents and against its three translation eigenphases. The
branch vectors sum to the original state and the probabilities sum to one.
All phase alternatives, including zeros, remain in the receipt.

Because `T_i^N` is the same total parity for all axes, mixed-parity phase
triples have zero effects; this is an algebraic constraint, not discarded
data. The preparation here is even and nonstationary. These are translation
phases of the complete gauge-constrained source, without a claim that an
individual charged particle has been isolated. A charged apparatus would
also have to retain and transform its compensating reference.

## Reference, scale and remaining field requirements

The exact finite H defines one clock `exp(-itH)`. If P is its complete Gauss
projector, the physical Gibbs state `P exp(-beta H) P / Tr(P exp(-beta H))`
exists for beta positive. The parent construction's complete physical
preparation/filter argument applies to this finite H as well. A safe full
norm bound at default N=2 parameters is 1896. For a path of length r its
triangle contribution is `4 |coefficient| sum_d d^(r+1)`; local Yukawa and
number terms give the remaining finite sum. No global interacting ground
or Gibbs distribution is computed here.

Spectral paths extend up to N-1 edges. This source has no asserted finite
interaction range or exact light cone. In the separate identity-background
calculation its retained linear Weyl symbol has the same one-cone free
limit as the owning spectral proof; its periodic ultraviolet branch is
discontinuous. Replacing it by a local sine symbol is another regulator and
introduces the corresponding extra corner branches. Neither replacement
identifies the excitations of the actual quantum-link Hamiltonian.

Raw path coefficients retain all intermediate native indices. For a block
of dimension d, a fixed-endpoint path on blank links has `d^(r-1)` orthogonal
link words. This growth is not normalized away; bounding interactions under
refinement therefore needs more than the free spectral calculation. A
different primitive normalization would be a separately declared source.

The construction closes the finite three-translation and periodic-source
definition. It supplies no interacting Lorentz limit, uniform local field
bound, source-selected Higgs phase, isolated particle pole, measured species
identity, or GeV mass. Those requirements must use this same complete source
or a proved full-residual comparison to an explicitly changed construction.
