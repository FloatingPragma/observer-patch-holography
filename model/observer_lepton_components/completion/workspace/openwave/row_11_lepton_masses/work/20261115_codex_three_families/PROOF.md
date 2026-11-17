# A geometric triplet and three relational observer channels

Platform row11 remains untested. This proof concerns the actual finite
source representations and a constructive observer-frame extension. It
makes no new numerical execution claim. No full Lorentz or continuum theorem
is added as a prerequisite for the static row.

## 1. The selected band is exactly the response-frame quotient

Use the common port labeling of `Screen/A5FamilyBand.lean` and
`Screen/PortFrameGram.lean`. Let A be the twelve-port adjacency, C the
antipode permutation, and A2 the distance-two matrix. The source tables give

    X=10(I-C),       Y=2(A-A2),
    P3=(X+sqrt(5)Y)/40,
    G=I-C+(A-A2)/sqrt(5).

Consequently, entry by entry,

    P3=G/4.                                                   (1)

Diagonal, adjacent, distance-two and antipodal entries of P3 are respectively
1/4, 1/(4sqrt(5)), -1/(4sqrt(5)), -1/4. The source proves P3 is a symmetric
idempotent of trace3. Thus V=range(P3) has complex dimension3, G is positive
semidefinite, and the Gram quotient is isometric to V by

    [x] in C^12/ker(G)  ->  2P3 x,
    x†Gx=4||P3 x||².                                         (2)

This is the same finite response-frame space, rather than another
three-dimensional multiplicity extracted from the coefficient table.
`observer_reconstruction/native_geometry.py` independently obtains the
candidate projector polynomial from complete native current records under
its disclosed twelve double-load preparations and uniform scheduler. Its
estimator reads no hidden output loads or source coordinates. The same code
states that physical position, scale and clock are not inferred merely by
recovering G. The current-quotient Lean files likewise leave those physical
attachments open. Equation (1) identifies their finite coefficient object;
it does not assert that a response vector is already a lepton family.

## 2. What multiplicity-free means here

Let U_g be the listed sixty proper port permutations. Their action is
transitive on ports. The stabilizer of a port cycles its five neighbors
and its five distance-two vertices, and fixes the antipode. Hence ordered
pairs have exactly four orbits, the four distances. A matrix B commutes
with all U_g exactly when B_(gi,gj)=B_(i,j), so

    End_G(C^12)=span{I,A,A2,C}.                               (3)

These four linearly independent matrices are also polynomials in A: the
source's exact four spectral projectors resolve the identity and have four
distinct adjacency eigenvalues 5, sqrt(5), -sqrt(5), -1. Their span is the
same four-dimensional algebra. In particular (3) is commutative.

For a finite unitary representation over C,

    H=direct_sum_lambda V_lambda tensor C^(m_lambda),
    End_G(H)=direct_sum_lambda I_(V_lambda) tensor Mat_(m_lambda).

Commutativity forces every m_lambda=1. Since the four source projectors are
nonzero orthogonal invariant projectors resolving I, and the commutant has
dimension4, each is one irreducible summand. Thus the source decomposition is

    C^12 = V1 + V3 + V3prime + V5,                            (4)

with multiplicity one for each inequivalent summand. In particular,
`End_G(V3)=C I3`. These facts are also explicitly certified in the existing
family-band and family_naturality code. They are reused here, not a new
source execution. The positive Laplacian costs are 5-sqrt(5), 6, 5+sqrt(5).
The declared faithful-band/window/cost task selects the first band. Its
realization and physical comparison clauses remain the source's named
interfaces; the unrestricted minimum is the constant zero-cost band.

There is a narrower readout rule that does not use the physical window:
on the zero-total current-response space the constant band is absent, and
the minimum positive Laplacian eigenvalue is uniquely 5-sqrt(5), with
projector P3. Thus a task that reads the slowest nonconstant response selects
a rank-three band from that source spectrum. This is a specified spectral
observation task. It does not yet prepare that band coherently or identify
its three response coordinates with charged particle species.

A rank-three residue is therefore one pole with three response components
for any source transfer that is a function of this Laplacian. It is not
three distinct energy poles in that transfer. An invariant self-adjoint
operator on this bare V3 is scalar. This statement has the bare-module and
unbroken-action hypotheses; it excludes no larger record-memory module or
relational observer construction.

## 3. The C3 register does not identify the full action

The lepton-register code explicitly allocates a neutral three-state register
F and applies its cyclic shift R on an oriented C-record. On a chosen face,
the native proper third-turn has eigenvalues 1,omega,omega² on V3, where
omega=exp(2pi i/3). This follows from its proper three-dimensional rotation
angle2pi/3, or the source frame action. Therefore the restriction of V3 to
that cyclic subgroup is unitarily equivalent to the regular C3 register.
An intertwiner can match the three normalized character vectors, with a
separate phase choice on each character. Matching this one subgroup does
not specify the action of the remaining port rotations or a physical
family interpretation.

The actual selected family amplitude is

    M(theta)=aI+b(exp(i theta)R+exp(-i theta)R†).

It commutes with R. For a phase where its effect has three distinct
eigenvalues, that effect cannot commute with every rotation of the bare
irreducible V3. The existing source already retains symmetry-breaking
attachments as possibilities. Rotating the detector/background together
can produce a covariant family of conjugate instruments with identical
spectra; invariance of one fixed apparatus is a stronger condition.

There is a finite spin check when, and only when, this same port action is
identified with the physical frame rotations. Let the binary lift of a
third-turn have spin-half matrix S with eigenvalues exp(±i pi/3), so Tr S=1
and its cube is -I. The other lift is -S, has trace-1 and cube I. The central
element acts as -I on either spin-half representation. On the bare triplet
and spin factor,

    Tr(V3(g) tensor S(gtilde))=0,
    Tr(I3 tensor S(gtilde))=3  (or -3 for the other lift).     (5)

Characters differ, so these six-dimensional representations cannot be
intertwined under that same lift. Equation (5) distinguishes a geometric
vector times a spinor from three internal spinor copies. It does not apply
when the family group is an independent internal action, or when a retained
observer frame participates. The next section supplies such a frame.

## 4. A constructive relational-frame remedy

Let G be the same order-sixty group and keep an actual observer frame label
g. On C[G] tensor V use diagonal frame transformations

    T_h=Left_h tensor U_h,       Left_h|g>=|hg>.

First, a coherent invariant-subspace construction is

    Theta v = (1/sqrt(60)) sum_g |g> tensor U_g v.             (6)

Orthogonality of frame labels and unitarity give Theta†Theta=I. Relabeling
g->hg proves T_h Theta=Theta. Conversely, invariance of
sum_g |g> v_g forces v_(hg)=U_h v_g, hence v_g=U_g v_identity.
The invariant subspace therefore has dimension dim V, exactly3 for V3.
Its projector is the group average of T_h. It contains three invariant
channels even though the bare V3 has no invariant vector.

For any complete instrument with Kraus operators K_(r,a) on V, define

    Ktilde_(r,a)=sum_g |g><g| tensor U_g K_(r,a) U_g†.         (7)

Then Ktilde Theta=Theta K and sum_(r,a) Ktilde† Ktilde=I.
Every failure and ancillary/reference output is carried through the same
formula; no successful probability is divided out. Equation (7) commutes
with T_h. Thus the invariant-sector instruments may have a full Mat3
algebra, including a nonscalar family effect. The bare Schur restriction
in section2 does not apply to this larger observer system.

Coherence between different frame labels is not needed for a reversible
relational record encoding. A smaller classical-frame version is

    E(rho)=(1/60) sum_g |g><g| tensor U_g rho U_g†,
    D(X)=sum_g U_g† <g|X|g> U_g.                             (8)

E is trace preserving and completely positive; D is a complete frame read,
conditional inverse and frame trace, also completely positive and trace
preserving. Direct substitution gives D E=id. The same identities hold with
identity on any untouched reference system, so entanglement within the
encoded family is retained. Measuring the frame label does not measure the
family character. E(rho) is invariant under T_h, and

    Phi_tilde,r E = E Phi_r                                  (9)

for the entire outcome-resolved instrument (7). Equivalently, the classical
frame algebra is the block-diagonal algebra of functions g->A_g. Its
invariant subalgebra consists exactly of

    A_g=U_g A_identity U_g†,

which is *-isomorphic to End(V). This is an explicit common observer algebra
with complete preparation, instrument and decoding maps.

The label is essential. If it is erased before applying its inverse,
the source band is twirled to I3/3 for any normalized input, by section2.
Freezing D while changing every implemented orientation to identity gives
the same depolarized result. A complete rank-one test on its prepared ray
has probability1 under the proper construction and1/3 under that changed
law, with the complementary outcome retained. These are analytical control
predictions, not execution results. The matched identity D E=id alone is
an identity of the specified services and earns no source or physical claim.

The two versions have different frame observations. Equation (6) uses a
coherent pure subspace. Equation (8) is a classical-quantum flagged encoding
with raw Hilbert-space dimension60 dim V and a logical Mat_(dim V) algebra;
it is not asserted to occupy a dim V pure subspace. They agree for the stated
block-diagonal frame-observation menu. Arbitrary coherent frame measurements
would distinguish them and belong to a different specified instrument.

If a spin factor is added and the binary group acts by
`Left_(pi(h)) tensor U_(pi(h)) tensor S_h`, equations (6)-(9) reduce its
logical action to `I_dimV tensor S_h`. The central lift fixes the frame and
port factors and acts as -I on spin. Thus the relational V3 sector has three
spin-half copies, with third-turn trace±3 on its coherent invariant sector.
The conditional obstruction (5) is removed by an explicit intertwiner.
An independently specified charge operator on a matter factor also commutes
with these maps. This transports an existing charge assignment; it does
not derive that assignment or join a different charge source by notation.

These are finite exact constructions. They require an enacted frame
register, conditional U_g and inverse readout, source realization of the
chosen V, and a complete source/reference map. The uniform classical frame
could be sampled by six fair native coins with rejection of labels60..63;
a cap of R trials has retained timeout16^(-R), and successful frames are
uniform. This is an analytical sampler specification, not a new event run.
Sampling a classical frame does not implement the coherent port rotations
or the family detector automatically. No new pulse certificate is claimed.

## 5. The full source still supplies all its bands

The relational construction works on every finite unitary module, not just
V3. Applying it to the full twelve-port space gives a logical Mat12 algebra,
or a coherent invariant subspace of dimension12. All four source bands
survive, with logical dimensions1,3,3,5. The inherited Laplacian remains
scalar at its four costs. Tensoring a common charged matter factor gives
that same charged representation to each supplied band unless a further
source/readout theorem distinguishes them.

Selecting P3 is a real operation or observation-domain choice. A complete
projective band instrument has outcomes P3 and I-P3, and its complement is
retained. The existing faithful-band/cost task provides a conditional,
non-empirical rule for favoring P3. It does not prove that every other band
is physically absent. A current-difference read removes the constant mode;
a frame-quotient read also annihilates the other bands. Darkness for that
menu does not prove darkness for all admitted observer interventions.
A physical claim may use a justified low-energy sector or operational
quotient, but it must state its energy/window and why omitted sectors do
not add the charged minima claimed by the row. There is no requirement to
show architecture uniqueness.

The same distinction applies to memory. A larger actual record module may
decompose as direct_sum V_lambda tensor M_lambda with dim M_lambda=3.
An invariant rest observable then acts as I tensor H_lambda on that
multiplicity space and may split it. This is a sufficient algebraic route
if the module, complete source instruments, charge readout and observer
energy are implemented. Merely tensoring a new C3 factor and declaring its
interpretation does not establish those source identifications.

## 6. The existing single-valley theorem is constructive but conditional

The configured field in standard_model_source_map has one-particle space
`C2_spin tensor M15 tensor C3_family`; its own fermion attachment explicitly
states that the family multiplicity is configured. Its full CAR extension,
internal current intertwiners, charged/colourless labels and all Fock
complements are genuine conditional results. They do not derive the C3
factor from the rank of the source Gram.

The existing openwave_single_valley_leptons source also proves the needed
static minimum count on its declared action. For each positive coefficient m,

    E(p)^2=[sum_j sin(x_j)^2+(a m+sum_j(1-cos(x_j)))²]/a²,
    dE²/dx_j=2 sin(x_j)[a m+1+sum_(k!=j)(1-cos(x_k))]/a².

The bracket is positive. Hence p=0 is the sole local and global minimum on
the Brillouin torus for that branch. Other corners remain as saddles/maxima.
Three configured colourless charged branches therefore have three minima;
removing the supplied Wilson term restores eight per branch. The complete
source retains other colourful charged channels, neutral modes and hole
coordinates. The theorem does not require a finished interacting continuum
before stating its finite static result.

Its common kinetic law, scalar background, phase, family attachment and
Wilson term are declared inputs. Their retrospective ancestry remains
explicit. No theorem reviewed here transfers that minimum count to native
current histories or to the conditional first-hit family source merely by
sharing the labels family, spin and charge.

## 7. A sufficient remaining join

The constructive target is now precise. Realize either the classical-frame
encoding (8) on the source-selected P3 response sector, or an actual native
record module with a three-dimensional internal multiplicity. Bind its
complete preparation and outcome maps to the same charged source. Then
show that a measured, common static energy/position dictionary gives one
minimum for each of its three nondegenerate charged responses, with all
other accessible sectors classified and all detector complements retained.
The phase and mass-map selection obligations are those of the separate
sealed packets; no value is supplied by this representation theorem.

The bare rank-three band alone does not identify three leptons. The explicit
relational-frame construction shows why that fact is not an obstruction to
an observer-first realization. This packet establishes the intertwining
mathematics and its exact scope, not the physical row or platform acceptance.
