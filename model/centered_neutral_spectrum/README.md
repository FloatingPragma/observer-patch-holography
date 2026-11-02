# Complete finite neutral gap and spin/density spectra

This package resolves the centered cube's **entire neutral sector** into
exact reducing spaces and certifies its first excited energy. It also
resolves the local hidden-spin and density response, including every
nonzero weight and exact zero-weight lines. These probes access neutral
states absent from a charged linear-CAR measurement.

The declared source is the unchanged open eight-site model

```
K = dGamma(A tensor I_2)/16 + sum_x n_xup n_xdown - N/2,
A=A^T, A^2=3I,                  c=1/8, U=1.
```

Energies use this source clock with hbar=1. The reference is its unique
N=8, hidden-spin-zero, eta-zero ground. It is different from the original
uncentered N=4 ground. The [charged-spectrum proof](../centered_particle_spectrum/README.md)
provides its exact source, finite symmetry representation, singlet
selection and integer ground witness. The [centered-ground proof](../centered_source_ground/README.md)
provides the applicable Lieb uniqueness theorem and an all-other-charge
energy separation strictly greater than 1/8. These are explicit imported
mathematical inputs; the new packet separately reports historical receipt
custody instead of claiming to rerun the prior experiments.

## Complete neutral-sector coverage

The source commutes with hidden spin SU(2) and eta SU(2). On each site the
four-dimensional Fock representation is `(1/2,0) direct-sum (0,1/2)`.
Eight sites therefore have integer labels `S,E>=0`, `S+E<=4` in N=8.
Every neutral energy has a state with spin weight zero, so it appears in
the complete `(Nup,Ndown)=(4,4)` space of dimension 4900. Conversely,
there is no assumption that a locally visible state exhausts that space.

For each pair S,E choose spin weight S and eta weight -E. The occupation
numbers are

```
Nup=4-E+S,                 Ndown=4-E-S.
```

The simultaneous kernels of `S+` and `eta-` in this sector contain exactly
one vector per spin/eta irreducible copy. Raising eta E times maps this
space injectively to N=8 without changing energy. Its spin multiplet
then contributes `2S+1` neutral states. Each eta raising changes the
translation character by 111, so the table distinguishes the lowest-eta
character k from the actual N=8 character `k xor ((E mod 2)*111)`.
The character-orbit multiplicities are unchanged by this complement.
On spin weight S,
`S-S+` has eigenvalues `j(j+1)-S(S+1)` for `j=S,...,4-E`.
The analogous eta operator has eigenvalues
`j(j+1)-E(E+1)` for `j=E,...,4-S`. Their exact kernel projectors are products
of `1-L/r` over the nonzero listed roots. They commute with K and with
each other. The implementation checks the polynomial identities,
weighted self-adjointness, rank traces and source commutation with
integer matrices.

The signed cube translations have commutators and squares in the central
one-body signs. Their Fock lifts on every even-N sector commute and square
to one. Their eight characters are `chi_k(b)=(-1)^(k.b)`, k,b in `(Z2)^3`.
Axis permutations conjugate these characters by permuting their bits;
thus representatives k=0,1,3,7 cover all characters with multiplicities
1,3,3,1. We diagonalize every spin/eta/character block, including blocks
invisible to the proposed probes. No spatial irreducibility assumption
is needed for this complete-spectrum step.

The 60 blocks have exact disjoint signed-orbit bases, diagonal positive
orbit metrics, and exact spin/eta projectors. Selected projector columns
are proposed numerically but their independence is verified modulo 1009;
the exact projector trace proves that they span the full range. The
weighted sum of all block dimensions is **4900**. Multiplication by
`2S+1` gives **12870**, the full N=8 dimension. Thus the construction
covers all neutral states, rather than only the cyclic range of a probe.

## Complete residual certificate and the finite gap

For each complete block the archive contains integer coefficients for
every proposed eigenvector and eigenvalue. Let Y be the resulting full
column matrix with the exact orbit metric included, D the real diagonal
proposal, and R=KY-YD. Exact integer/Fraction arithmetic bounds

```
delta >= ||Y*Y-I||,       rho >= ||R||,
epsilon = rho/sqrt(1-delta),                 delta<1.
```

Frobenius norm bounds R and maximum absolute row sum bounds the Hermitian
Gram defect. Outward rational square roots use integer square roots.
The matrix Y is invertible onto the entire reducing block, and
`Y^-1 K Y=D+Y^-1 R`, with perturbation norm at most epsilon. Every true
eigenvalue lies in the union of the epsilon discs around D. This follows
also from `(D-lambda)Y*v=-R*v` for a normalized true eigenvector v and
`||Y*v||>=sqrt(1-delta)`. Along `D+tY^-1 R`, every connected disc component
retains its initial eigenvalue count by continuity. This justifies isolated
simple discs, without an ordered-Weyl argument for a nonnormal similarity.
Repeated higher eigenvalues are enclosed as clusters and never silently
assigned individual eigenvectors.

The lowest complete block is S=E=0,k=0. The next is S=1,E=0,k=7; its
first eigenvalue is isolated from every other block and every other
level in that block. Complete interval comparisons give

```
0.015021254176 <= Delta_neutral <= 0.015021254177.
```

Its multiplicity in N=8 is exactly **three**, the spin triplet; its
spatial character is k=111. Since every N unequal to 8 is already proved
to lie more than 1/8 above the ground, this same interval certifies the
**global spectral gap of the full 65536-state source**. This closes the
previous missing *finite* quantitative neutral-gap premise. It is not a
volume-independent or continuum mass-gap theorem.

For any bounded finite-source operators A,B of zero ground mean and
Euclidean time t>=0, the spectral theorem consequently gives
`|<A G, exp[-t(K-E0)] B G>| <= exp(-0.015021254176 t) ||A G|| ||B G||`.
This is a quantified finite connected-correlation bound on the same
reference, with no spatial or continuum limit presumed.

## Neutral probes and exact dark lines

Use the Hermitian neutral probes `Sz_x=(n_xup-n_xdown)/2` and
`rho_x=n_xup+n_xdown-1`. The ground singlet implies their elastic
expectations vanish. SU(2) selection gives

```
Sz_x G lies in (S,E)=(1,0),
rho_x G lies in (S,E)=(0,1).
```

For smaller exact matrices use the highest spin component
`S+_0=c_0up* c_0down` and lowest eta component `c_0down c_0up`.
The local spectral weights of Sz are one half the S+ weights; those of
rho are twice the pair-annihilation weights. These factors follow from
the spin-one ladder norm sqrt(2), not a fitted normalization. Eta raising
also multiplies a site's neutral density component by its bipartite sign.
The pair-sector energy itself is unchanged because K commutes with eta.

The origin stabilizer fixes both local probes. Exact signed stabilizer
orbits and the spin/eta projectors reduce each probe to 145 dimensions.
The total S+ channel (translation character 000) annihilates G, and the
total eta- channel (pair character 111) also annihilates G. Removing
these **proved zero** components leaves a complete 136-dimensional
candidate space for each local response. The exact operator images of
every vector in the retained ground singlet space are checked to lie in
these spaces; no amplitude cutoff is used.

Each candidate space has a simple eigenline at K=-2 with zero response
weight. Its integer vector z satisfies `(H16+32 Gram)z=0`. If T denotes
the exact source-probe transition to the ground-singlet coordinate space,
the archive provides an integer a and positive integer d satisfying

```
Gram_ground a=d T* z,
(H16_ground+32 Gram_ground)a=0.
```

Therefore the probe pullback of this eigenline lies in the ground-space
energy -2 eigenspace. It is orthogonal to G, whose certified energy is
near -4.124. This proves the weight is exactly zero. A tiny floating
amplitude is not used as evidence for its absence.

Every other local candidate line has a strictly positive interval weight:
**135 spin lines and 135 density lines**. For an isolated normalized
column with residual r and separation gamma from all other spectral
discs, its phase-aligned eigenvector error is at most `2r/gamma`.
The probe operator norm is at most one, so the amplitude error is bounded
by the sum of the ground and response vector errors. Exact normalization
intervals and squaring then enclose the physical weights, including the
smallest ones. All lines, zero lines and spatial residue entries are
retained in the receipt; the CSV retains every local energy and weight.

## Spatial residues, multiplets and current readouts

On even particle number the eight translation lifts form an ordinary
representation. For a normalized local candidate eigenvector psi define
`m_b=<psi,R_b psi>`. A simple eigenline in the complete stabilizer-fixed
space generates one irreducible cyclic cube representation: two nonzero
irreducible components would give two independent fixed vectors at that
same energy. The site-permutation representation has dimensions 1,3,3,1.
The eight translated vectors are an invariant tight frame, hence

```
sum_b |m_b|^2=8/dimension.
```

Exact overlap intervals uniquely identify dimensions 1 or 3. The spin
residue is `B[x,y]=weight*m_(x xor y)`. Density residues include the
additional bipartite factor `(-1)^parity(x xor y)` from eta raising.
These are dimensions of the probe-cyclic visible spaces, not assertions
that no other eigenvectors share a visible energy. Each component reaches
357 dimensions (24 rank-one and 111 rank-three lines). The three spin
components and density reach 1428 orthogonal neutral dimensions. Of the
12869 excited N=8 states, 11441 lie outside this specified cyclic range.
Current divergences add no new cyclic directions.

The same receipt gives the outgoing current-divergence readout. For
source hopping t_xy, define

```
j_xy=i t_xy sum_spin(c_xspin* c_yspin-c_yspin* c_xspin),
dot rho_x=i[K,rho_x]=-sum_y j_xy.
```

The onsite contact and number-centering terms commute with density.
The hidden-spin current uses the same expression with weights +1/2 and
-1/2 for the two spins, giving `dot Sz_x=-div j^S_x`. Thus each
current-divergence spectral weight is exactly E^2 times the corresponding
neutral weight, with the same complete support. This covers current
divergences, not independent transverse edge-current channels.

The first-moment identity `M1(rho)=4 M1(Sz)` follows from the
double commutators with K: opposite-spin onsite operators commute,
and the two squared spin coefficients are 1/4 of the density
coefficients. Its interval sum is checked independently of the
onsite equal-time identity.

The exact onsite identity
`Sz^2+rho^2/4=I/4` is checked against the summed interval weights.
No gap or spectral line is interpreted as a Lorentz mass. In particular,
a nonzero fixed-cube gap is compatible with the separately proved
massless free continuum branch: uniform volume/refinement control would
be necessary to infer a continuum gap. Measured particle identification,
a surviving interacting relativistic massive sector and Einstein dynamics
require their separate source-to-continuum proofs.

## Replay

Run `python -W error -m pytest -q model/centered_neutral_spectrum`, then
`python -W error -m model.centered_neutral_spectrum --write` and `--check`.
The archived proposal is untrusted; complete integer/rational checks are
rerun in a fresh interpreter. `--propose` creates a new proposal and must
pass the same verifier before replacing the archive. No old package,
receipt, source law or global ledger is modified by these commands.
