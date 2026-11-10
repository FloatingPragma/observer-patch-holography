# Original composite spatial and final-energy records

## 1. Source, fields and observation domain

This calculation uses member zero of the
[common source contract](../particle_common_contract/PROOF.md): eight cells,
the full physical Wilson Hamiltonian H, the degree-one link carrier, its
selected faithful polynomial density rho, the retained purifying reference
and the declared Cayley clock. No compression of the physical Hamiltonian
to the trial vectors below defines an evolution law. All energies are source
energies. The source-only operator inequality is

    rho >= c P_G,       c=2^(-582407133).                         (1)

Its [owning construction](../particle_spatial_modular_response/PROOF.md)
does not assert the corresponding inequality for the joint pure reference
state. Positive complete-event probabilities on the source use (1); the
apparatus acts as the identity on the retained reference.

The [actual composite kernels](../particle_primary_composite_spin/PROOF.md)
and [native current tests](../particle_candidate_composites/PROOF.md)
supply three increments in native charge order (-4,-3,1,2,6):

    Delta_W=(0,2,0,0,1),   Delta_C=(1,0,0,2,0),
    Delta_E=(0,0,0,3,-1),                                  (2)

and their negatives. They are local gauge singlets. Each odd kernel has
six columns, with proper binary-cube representation G1+H. Each even kernel
has eight columns, with representation E+T1+T2. These are source rotation
representations, not Lorentz spins or identified particles.

The code forms the literal central character projectors in each kernel and
orthonormalizes their columns. Let F be one resulting local polynomial,
with T nonzero monomials and squared coefficient sum one. Each polynomial
is homogeneous in quasiparticle creators: positive-beta particle creation
or negative-beta hole creation, represented by a native annihilator.
On the filled-negative-beta, constant-link physical anchor Omega,

    phi_x=F_x Omega,  ||phi_x||=1,  F_x* phi_x=Omega,
    ||F|| <= sqrt(T),  A=F/sqrt(T),  ||A||<=1.              (3)

There are forty such components across the six signed families. Their
monomial count is a rigorous contraction normalization; no expectation in
phi_x substitutes for an expectation in rho.

## 2. Exact full-source first action and its energy

Let p be the CAR parity of F, let x and y be distinct cells, and define

    G_xy=F_x* [H,F_y] - (-1)^p [H,F_y] F_x*.               (4)

For nonneighbors this operator vanishes by finite-range locality and
graded CAR. For neighbors only their shared Wilson bond contributes.
The actual Wilson blocks in the
[native-chart calculation](../particle_primary_composite_transport/PROOF.md)
are -I/2 on positive beta and +I/2 on negative beta. The annihilator
commutator supplies the opposite sign for holes. Thus each term of
[H,F_y] replaces one quasiparticle at y by the corresponding one at x,
with coefficient -1/2 and a native link matrix entry, conjugated for a
hole. Cross-beta terms have no contraction with F_x* and vanish in (4).

The code evaluates this identity in the literal 120-mode two-cell CAR
chart on sparse witnesses, including every filled spectator mode. It uses the actual native
chart and separately verifies a non-diagonal link and all three Wilson
spin matrices. Haar orthogonality for an irreducible native block r gives
norm 1/d_r to each matrix entry. Within each family a given native block
appears with only one particle/hole orientation. Distinct blocks and
entries are orthogonal. The resulting vector G_xy phi_x has exactly
one defect at x, d-1 defects at y, and one native link word, where d is
three for the odd families and four for the even family.

Every such one-link word is admitted by the actual degree-one cutoff.
F commutes with that cutoff and with the Gauss projector, so the full
physical commutator equals the compression of the product-carrier
calculation. This is the same admission argument used in the full-source
[odd graded response](../particle_primary_odd_sum_rules/PROOF.md#3-literal-neighboring-graded-response).
The program checks every column's positive norm; no zero is replaced by
a guessed lower bound.

Put g_raw=||G_xy phi_x||^2 and let chi be this vector normalized to one.
Its expectation in the *full* H is exact. The onsite term is -720+3d.
Hopping has zero expectation because it changes the definite sublattice
fermion parity. A plaquette has zero expectation because at least one
of its other three links is constant and has zero native first Haar
moment. The electric mean is the weighted native Casimir:

    M_chi=-720+3d + (105/2)
        sum_(r,i,j) C_r ||v_(r,i,j)||^2/d_r / g_raw.       (5)

The native Casimirs are the ones used in the
[full-sector bounds](../particle_primary_sector_threshold/PROOF.md).
Both signs give the following exact rows; repeated entries denote separate
orthonormal source rotation components.

| Family | Component | T | g_raw | M_chi |
| --- | --- | --- | --- | --- |
| Weak | G1, each of 2 | 3 | 3/8 | -5513/8 |
| Weak | H, 0 and 3 | 1 | 3/4 | -690 |
| Weak | H, 1 and 2 | 3 | 1/4 | -690 |
| Color | G1, each of 2 | 9 | 5/36 | -6739/10 |
| Color | H, 0 and 3 | 3 | 1/4 | -2021/3 |
| Color | H, 1 and 2 | 9 | 1/12 | -2021/3 |
| Even | E, 0 / 1 | 2 / 6 | 1/4 | -1353/2 |
| Even | T1, 0 and 2 | 4 | 11/16 | -14841/22 |
| Even | T1, 1 | 6 | 1/4 | -1353/2 |
| Even | T2, 0 | 2 | 1/4 | -1353/2 |
| Even | T2, 1 and 2 | 4 | 7/16 | -1353/2 |

For the contractive A in (3), the double polynomial in (4) is divided by
T. Its squared first-action norm is g=g_raw/T^2. A common lower bound for
all forty components is

    gamma=min g=1/972.                                    (6)

The means (5) are neither eigenvalues nor mass estimates.

## 3. A complete energy-resolved spatial experiment

Let L_Delta be the proved lower bound for H on the full conserved native
occupation sector n_0+Delta. The anchor occupations are
n_0=(48,32,96,48,16). Equation (4) has net increment zero, so its first
action remains in that sector. The inherited global ground bounds L0,U0
satisfy L_Delta-1/64>U0 for all six sectors.

For a positive spectral measure with lower bound L_Delta and mean M_chi,
the Markov inequality applied to H-L_Delta gives

    <chi, 1_[L_Delta,3 M_chi-2 L_Delta](H) chi> >=2/3.     (7)

The spectral and completeness premises are supplied here by the finite
self-adjoint physical H and the explicitly calculated normalized chi.
No thermodynamic or continuum theorem is required for this finite
inequality. Use the actual complete
[energy clock](../particle_primary_joint_spectroscopy/PROOF.md), with
resolution eta=1/64, tail allowance1/1024 and complete implementation
allowance1/1024. Let E_J be its positive event effect for the reported
interval

    J=[L_Delta-eta, 3 M_chi-2 L_Delta+eta].                (8)

The summed tail and one-clock implementation loss is below1/6, hence
<chi,E_J chi> >=1/2. This is a full-source final-energy event. It is not
a signed transfer spectrum, an initial-energy selection or a spectral
calculation in a compressed Hamiltonian. The receipt also gives its
ground-relative enclosure using L0,U0.

Take R>=||H|| from the actual source preparation plan. Use one exact
Cayley step U_t=(I+i tH/2)^(-1)(I-i tH/2), and its inverse, in

    K_xy(t)=(A_x* U_t* A_y U_t
                -(-1)^p U_t* A_y U_t A_x*)/2.            (9)

Both products are contractions, so K is a contraction. The success K
and failure sqrt(I-K*K) form a complete instrument. K is gauge invariant,
even and has zero native increment; the positive-square-root failure has
the same properties. For example the finite Julia unitary dilation of
K implements this instrument with a retained two-state pointer and
commutes with each combined native charge. The inverse clock in (9) is
an explicitly declared coherent apparatus action, not an autonomous
reversal of source evolution. Its operators are functions of the same H.

Follow this complete instrument by a complete measurement of all five
native occupations and the complete final energy clock. Retain every
outcome, including the failure and all other sectors and energy bins.
The event of interest has success, occupation n_0+Delta and final label
in J. The native label is an absolute source occupation sector, not a
single-particle decoder; several excitations can have the same total.
There is no initial source preparation or initial charge measurement.

The short-time bounds proved in the
[same-clock source response](../particle_primary_odd_sum_rules/PROOF.md)
apply to any contractive A of either parity:

    ||K_exact(t)-i t G_A/2|| <=2 R^2 t^2,
    ||K_xy(t)-K_exact(t)|| <=R^3 |t|^3/6,                 (10)

where G_A=G/T. For |t|<=1/R the sum is at most3R^2t^2.
Set

    t=c gamma/(64 R^2).                                  (11)

The code stores this as a rational prefactor and a binary exponent; it
does not allocate an integer with 582407133 bits or the full source
matrix. On the physical phi_x, (7)-(10) and the contraction E_J^(1/2)
give

    ||E_J^(1/2) K_xy(t) phi_x||
       >= t sqrt(g)/(2 sqrt(2))-3R^2t^2
       >= t sqrt(g)/4.                                  (12)

For the second inequality, gamma<=min(1,g) gives
3R^2t<=3sqrt(g)/64 and 1/(2sqrt(2))>1/3.
Including the commuting native-sector projector does not change this
vector or its energy event.

Apply (1) to the positive complete-event effect, with the unmodified
original input rho. This gives the original, unnormalized probability

    p_neighbor(J,n_0+Delta,success) >= c t^2 gamma/16.     (13)

This is not a comparison between reference density operators. All
reference branches remain in the complete instrument's output, whose
source marginal generally differs from its input after intervention.

For distinct nonneighbors G_A=0, so the same event obeys

    p_nonneighbor(J,n_0+Delta,success) <=9R^4t^4,
    p_nonneighbor/(c t^2 gamma/16) <=9c gamma/256<1/2.     (14)

Thus the original joint spatial and final-energy probabilities differ
by at least c t^2 gamma/32>0 for every one of the forty components.
This is proved using the complete source's exact finite law; its huge
histogram is not explicitly sampled. A uniformly selected retained
neighbor/nonneighbor setting bit halves these joint weights. A uniformly
selected retained component label divides them by forty. Neither operation
renormalizes successful events or discards failures.

## 4. What the discriminator establishes

The six signed conserved-sector families have nonzero, energy-resolved
original spatial responses under this complete apparatus. Their recorded
sector, source rotation representation and CAR parity distinguish the
specified candidate channels. The source current and conditional baryon/
lepton bookkeeping remain exactly those of the
[candidate identity calculation](../particle_candidate_composites/PROOF.md).
No Standard Model field dictionary is inferred from those source labels.

The original odd neighboring first-shell moment m1 is a *signed* selected
expectation of a Wilson-dressed quartic transfer operator, as shown in
[the source sum rules](../particle_primary_odd_sum_rules/PROOF.md).
Equations (4)-(5) calculate a nonzero vector and a positive effect; they
do not determine that signed expectation. Faithfulness transfers positive
effects, not a sign for an arbitrary quartic. Likewise prepared third
moments do not determine the original two-point spectrum. The exact
calculation supplied here is the distinct four-point complete record
(13)-(14).

Its source-specific time and original probability floors are extremely
small. The result establishes neither a capacity-uniform finite-energy
residue, an isolated pole, a charged threshold, conserved momentum,
relativistic dispersion nor a GeV mass. No neutrino, neutron or antiparticle
of a measured species is identified. Signed source sectors and their
conjugation can be read without promoting them to that physical dictionary.
