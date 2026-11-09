# Physical composite multiplets on the primary source

## 1. Domain and exact local kernels

Use member zero of the [primary contract](../particle_common_contract/PROOF.md),
with its original Wilson action, selected state, idle quantum reference,
Gauss projector, link cutoff and energy clock. The native internal blocks
have charges q=(-4,-3,1,2,6) and dimensions (3,2,6,3,1). Internal indices
are the ascending exterior-mask chart of `native_quantum_gauss`. No
lexicographic exterior-link matrix is identified with that chart by index.
Spin and charge-conjugation matrices in this construction act on Dirac
rows and the identity on those fifteen internal indices.

Let Ω fill both negative-beta Dirac rows and empty both positive rows at
each site, with constant link functions. At one site it fills thirty of
sixty CAR modes and is a physical even singlet. For a signed native
increment Δ, take all minimal-defect masks with positive particles for
Δ_r>0 and negative-beta holes for Δ_r<0, allowing both spin rows. Stack
all twelve actual native Gauss generators on that finite mask space and
solve their common nullspace. Exact exterior CAR calculation gives:

| Δ, with either overall sign | Ambient dimension | Gauss kernel | Dirac spin Casimir |
| --- | ---: | ---: | --- |
| (0,2,0,0,1) | 12 | 6 | 3/4 twice; 15/4 four times |
| (1,0,0,2,0) | 90 | 6 | 3/4 twice; 15/4 four times |
| (0,0,0,3,-1) | 40 | 8 | 2 three times; 6 five times |

These are the complete kernels on these stated minimal local occupation
spaces, not a catalogue of all physical source states. The first kernel
is the weak epsilon contraction of two doublets times a charge-six
singlet. The second is a color epsilon contraction of three color
antifundamentals. The third is a color determinant times a charge-six
hole. Total gauge charge is zero in each case. Gauss invariance follows
from all twelve matrices, not from this charge sum alone.

In the actual Dirac chart define J_a=dΓ(-i γ_b γ_c/2), with (a,b,c)
cyclic. The program computes their full CAR action on each mask basis,
checks invariance of the gauge kernel, [J_a,J_b]=iJ_c, and the exact
Casimir multiplicities above. The two identical weak-doublet factors
are spin symmetric after their antisymmetric gauge contraction. Their
spin triplet times the singlet field's spinor gives 3×2=4+2. The color
pair has the same symmetry. For the mixed even composite the spin
quartet times the hole spinor gives 4×2=5+3. These interpretations agree
with the computed matrices; they are not a dimension-only assignment.

## 2. Actual cube transformations and their physical spatial orbit

The [source cube action](../particle_primary_channel_labels/PROOF.md)
uses S(q)=wI-xγ_2γ_3-yγ_3γ_1-zγ_1γ_2 and rotates the actual cells and
links. Its central element is CAR parity. The code acts with S(q) on
every creation factor and its conjugate on every annihilation factor,
including all exterior signs. The negative-beta sea has determinant one
under each proper spin rotation. The resulting kernel matrices are
unitary and have these exact character decompositions over all48
binary-cube elements:

    odd three-defect kernels: G1 ⊕ H,
    even four-defect kernels: T1 ⊕ E ⊕ T2.

The local continuous spin coefficient algebra is not a continuous
symmetry of the spatial source Hamiltonian. A localized vector also
moves under a physical cube rotation. To obtain actual full-source
labels, write φ_(x,v) for local kernel vector v at one of the eight
member-zero sites, with Ω at the other sites and constant links. The
vectors at distinct sites are orthogonal by local occupation. All
spectator site states are even, so no fermion sign is omitted in

    U(q) φ_(x,v) = φ_(q x,D(q)v),
    Ψ_v = (1/√8) Σ_x φ_(x,v),
    U(q) Ψ_v = Ψ_(D(q)v).

The uniform orbital factor is A1g. Physical inversion contributes
(-1)^(number of holes) because beta acts as -1 on each hole and +1 on
each positive-beta particle; its determinant on the anchor is +1.
Consequently the positive odd multiplets have G1g⊕Hg, their opposite
increments have G1u⊕Hu, and both signs of the mixed even multiplet have
T1u⊕Eu⊕T2u. These are finite-source symmetry labels. They do not prove
Lorentz spin, chirality, an isolated particle pole or a Standard Model
species assignment.

The full-source central group projectors

    Π_λ = (dim λ/96) Σ_(g in binary cube × inversion) χ_λ(g)* U(g)

form a complete orthogonal instrument and commute with H and every
native total occupation N_r. Their joint measurement with all N_r and
the same source energy clock is therefore admitted without changing the
ideal energy marginal. Retain every group label and every clock outcome.
The local minimal kernel itself is not asserted invariant under H.

## 3. Actual composite operators, currents and finite statistics

For an orthonormal local kernel column v define F_v by replacing each
mask coefficient with its corresponding minimal quasiparticle monomial,
including the CAR sign that makes F_v Ω=v. Positive defects use c†;
holes use c. All terms have the same degree d=Σ|Δ_r| and increment Δ.
Thus on the full source Fock carrier

    [N_r,F_v]=Δ_r F_v,        (-1)^N F_v=(-1)^d F_v (-1)^N.

The native gauge commutator of F_v stays within the same minimal
monomial space. Acting on Ω gives zero by the actual Gauss-kernel
calculation. The map from those monomials to their distinct masks on Ω
is injective, so [G_a,F_v]=0 as an operator, for each native generator.
The connected compact native gauge action therefore fixes F_v. This
argument does not promote an arbitrary zero-charge monomial to a singlet.

If x≠y, every factor in F_(x,v) anticommutes with every factor in
F_(y,w) or its adjoint, and there are d_v d_w such interchanges. Hence

    F_(x,v) F_(y,w) = (-1)^(d_v d_w) F_(y,w) F_(x,v),
    F_(x,v) F_(y,w)† = (-1)^(d_v d_w) F_(y,w)† F_(x,v).

In particular the two trilinear families have exact odd disjoint-support
exchange. They need not satisfy the canonical free CAR at one site.
For a single trilinear product the local anticommutator is a nonconstant
occupation projector. This finite graded locality is compatible with
fermionic readout; it is not the relativistic spin-statistics theorem.

A column with s nonzero normalized monomials obeys ||F_v||≤√s by the
triangle and Cauchy-Schwarz inequalities. Thus K_v=F_v/√s is a bounded
complete-instrument success operator, with complement
(I-K_v†K_v)^(1/2). On Ω its success probability is exactly1/s because
F_v†F_v Ω=Ω. A success flag carrying the compensating -Δ sector, with
an even combined source-and-memory map, preserves total occupations
and parity. All complementary outcomes and references remain present.
These finite algebraic services do not prove efficient autonomous control.

## 4. Literal charge conjugation

Use the [source charge conjugation](../particle_primary_charge_conjugation/PROOF.md),
which preserves the full action, Gauss carrier, cutoff and selected state.
In the beta chart its Dirac matrix is an exact signed permutation. The
program applies c_i† to Σ_j M_ij* c_j with full sixty-mode CAR signs,
starting C|empty>=|full>. Locally CΩ=-Ω; on the eight-cell source the
anchor phase is +1. The induced map relative to that anchor takes Δ
to -Δ, is unitary between the computed kernels, and intertwines every
J_a. It also intertwines every proper binary-cube matrix.

For odd composites inversion changes sign under C; for even composites
it does not. The uniform-orbit labels in section2 therefore match the
literal identity C I C*=(-1)^N I. The same H gives equal conjugate
sector spectra and, for conjugate prepared inputs and covariant complete
apparatus, equal energy-record laws. This does not identify a named
particle or its measured antiparticle.

## 5. Populated complete symmetry-and-energy records

The [sector threshold](../particle_primary_sector_threshold/PROOF.md)
proves full physical lower bounds L_Δ>U_0, where E_0≤U_0 is its global
ground upper bound. These bounds hold on the entire conserved native
sector and therefore on every finite-cube subspace above.

For every normalized uniform-orbit vector Ψ_v, the exact source mean is

    M_Δ = <Ψ_v,H Ψ_v> = -720+3d.

The onsite expectation gives this value. Hopping and magnetic terms
have zero constant-link Haar expectation. Between distinct localized
vectors a single H term cannot move all three or four occupation defects;
therefore cross terms vanish too. This includes every vector in G1 as
well as H. A one-spin-row extremal witness alone would not establish G1
population. The sector proof's dressed trial Rayleigh value is not assigned
to every multiplet here.

Let L_0 be the global source lower bound. In each populated finite-cube
sector the bottom excitation obeys

    0 < L_Δ-U_0 ≤ inf spec(H|_(Δ,λ))-E_0 ≤ M_Δ-L_0.

This is a sector threshold interval, not an isolated line or rest mass.
For a prepared normalized vector in one selected cube component, Markov's
inequality applied to H-L_Δ≥0 gives original spectral weight at least1/2
in [L_Δ,2M_Δ-L_Δ]. Use the inherited complete energy clock with resolution
η=1/64 and total tail plus implementation error at most1/512. The joint
event consisting of that component and reported energy in

    [L_Δ-η, 2M_Δ-L_Δ+η]

has original weight at least255/512. Measuring the exact component first
and applying the clock makes this joint bound explicit; all other records
are retained. Arbitrary mixtures retain the same bound after summing their
populated labels. The clock error is charged once to that joint event.

Preparation is the inherited paired-sector Householder on two idle
physical memories, followed by graded source exchange and the required
parity correction. One memory keeps the whole original source state and
its correlations with the untouched reference; the other retains an
opposite-sector compensator. The source marginal changes to Ψ_v. The
combined occupations and parity are preserved. This is a controlled
intervention on the same H and clock, not a claim that the original
selected state has weight255/512 in an odd sector. The memory resource
bound and work interval are included in the bound receipt.

## 6. Scope of the finite theorem use

The primary finite Born/readout and retained-reference theorems apply to
these actual complete instruments; no quantum probability rule is added.
The source's finite group, gauge and CAR constructions provide explicit
premise witnesses for symmetry, conserved charges and graded locality.
No new Lean formalization of the sixty-mode matrices is claimed. Source
spectral thresholds and populated finite-cube records require no continuum
limit. Recovering physical masses requires a controlled particle branch,
its dispersion and a physical energy reference; a GeV mass is not inferred
from these finite rotation labels or bare onsite coefficients.
