# Localized fermion fields on the compact texture pair

The existing compact-pair source supplies a positive defect-to-mode map:
four spatial modes form an isolated band, its CAR algebra carries the actual
M15 currents, and the declared symmetric pairing restricts exactly to that
band. This is a map on a fixed designed texture. Its frequencies are not
the free homogeneous charged action's Dirac masses. No occupation operator
creates the background winding or supplies a moving-defect law.

## The actual source and an improved exact band enclosure

Use the literal source in `model/compact_defect_pair/{PROOF.md,spectral.py}`:
two cores (0,0,+1),(0,0,-1), the normalized compact Pauli texture
f=(x,y,z^2-1) inside radius2 and north outside, zero amplitude at the cores,
onsite magnitude M>=8 elsewhere, and centered Weyl hopping coefficient1.
The two surrounding shells have degrees +1 and -1, with constant exterior.
The topology concerns this declared texture. The spin-core fill, zero
fermion amplitude and exterior geometry retain their owning premises.

Let P select the four core coordinates, Q=1-P, h=T+V, B=QhP, D=QhQ.
The owning theorem gives exactly four finite- and infinite-volume eigenvalues
in [-2,2], no other spectrum in the gap, ||B||^2=7/4, ||T||<2, and

    E u=F(E)u,  F(E)=-B*(D-E)^(-1)B,
    |E|<=r=(7/4)/(M-4).                                  (1)

It also computes the exact leading matrix F0=-B*V_Q^(-1)B. Its diagonal
core blocks vanish and its offdiagonal block is -sigma3/(4M), so its
four eigenvalues are -alpha,-alpha,+alpha,+alpha, alpha=1/(4M).
The present calculation reuses those identities and the cluster-count
proof. It does not rerun or fit the retained 686-state eigensolver proposal. The
producer authenticates the compact owner and its public action fixtures;
internal matrices use that owner's existing action basis. The proved P15
dictionary in native_matter_attachment relates it to the literal native
current basis. An independent literal-current check is separate from the
bounded fixture replay.

There is a stronger source-specific remainder. Both cores have odd lattice
parity, hence B has range on the even sublattice. Removing the cores does
not change bipartiteness. V_Q-E and its inverse preserve parity, whereas
T_Q changes it. Every compressed odd-length term in the convergent Neumann
series for F(E) is therefore zero. Put q=2/(M-r)<1. Resolvent subtraction
for the zero-hop term and the sum of the even terms n>=2 give

    ||F(E)-F0|| <= eta(M),
    eta(M)=(7/4) r/[M(M-r)]
             +[(7/4)/(M-r)] q^2/(1-q^2).                (2)

For M=64 this is the exact rational number

    eta=2365097/60283957504 < 1/25000.                   (3)

Every one of the existing four roots lies within eta of +1/256 or -1/256.
The particle-hole symmetry hC=-C conjugate(h), C=epsilon_x i sigma2,
pairs opposite energies. Since eta<1/256, exactly two roots lie in each
sign interval. This is an enclosure of the full spectrum in that isolated
band, not a statement that its two positive roots are distinct. The same
proof works in every admitted exterior cube and on the infinite lattice.
The inherited exterior-volume error remains valid.

## Same-current CAR fields and exact time evolution

Let W be any complete four-column orthonormal frame for this band,
A=W* h W, S=W* C conjugate(W). The owning identities are

    hW=WA, C conjugate(W)=WS,
    S^T=-S, S*S=I, AS=-S conjugate(A).                  (4)

For internal/family index alpha in M15 tensor C3, define the fermion field
b_(j,alpha)=sum_x conjugate(W_xj) a_(x,alpha). These are canonical odd
CAR fields, not auxiliary even occupation-register operators. Every
spatially uniform native current intertwines through W tensor I45.
All current commutators follow on the full Fock algebra. The same spatial
band supports every internal M15 component; the texture does not select
the charged-lepton line or its three supplied family copies. This describes
180 localized CAR coordinates, not 180 independently identified particles.

For the already declared pairing Delta=C tensor conjugate(K), where K is
any symmetric 45-by45 internal/family kernel, the complete Nambu inclusion
E=diag(W tensor I45,conjugate(W) tensor I45) obeys

    B_full E=E B_band,
    B_band=[[A tensor I45, S tensor conjugate(K)],
            [(S tensor conjugate(K))*, -conjugate(A) tensor I45]]. (5)

It follows that the embedded CAR subalgebra has exact closed Heisenberg
Bogoliubov evolution, on the same clock, for every initial full-source state
and reference. Restricting that state supplies the band state; no initial
product factorization is needed. This is an algebra/dynamics statement.
It is not a Fock-vector embedding with a vacant exterior: complementary
pairing can create pairs from that proposed exterior vacuum.

Equation (5) includes A. Removing it has operator-norm defect at least
`1/256-eta` in this source. The discarded coefficient is neither rounding
noise nor an optional zero-mode convention.

### Complete feedback on the compact source

The spatial number-conserving Hamiltonian is independent of the family
coordinate. Extend that same Hamiltonian to every ancillary flag copy in
the [full-CAR analyzer](../charged_lepton_response/CAR_KOIDE.md). Its vacuum
ancillas reduce the Hamiltonian, and its uniform flag/family frame commutes
with the spatial operator. The common-Hamiltonian premises of that theorem
therefore hold for the declared staggered local pairing used here.
For a finite region with V sites, a conservative kick norm is
`6 V |v_H|`. Committed classical flags with all charged quantum complements
give the complete unhalved error bound `108 V^2 |v_H|^2 t^2/N`, plus
the summed compilation error. The corresponding coherent archive has the
separate weaker bound in the owning theorem.

This constructs the chosen paired action using the same spatial operator
and declared pulse clock. It supplies ancillary copies and feedback
services; it does not prepare the paired ground state, select the texture,
or identify a native attempt count with physical time.

## Complete charged frequency groups for the configured phase-zero action

Instantiate the actual charged-Higgs line with the parent's coherent Y at
phase0, neutral norm1. Its family eigenvalues are

    m_l=(15-10sqrt(2))/216, multiplicity2,
    m_h=(15+10sqrt(2))/108, multiplicity1.               (6)

The native rank-one directed charged tensor makes K have each family
singular value twice and 39 zeros. Its charged supports have electric
charges -1 and +1. The neutral pairing preserves that charge; all twelve
current identities require the associated Higgs covariance, not invariance
of a fixed neutral Higgs under broken generators.

The complete spectrum in (5) is E_j plus or minus each singular value of K.
One direct proof conjugates the full Nambu matrix by diag(I,C*) to obtain
identical h blocks and an offdiagonal internal K block; these commute.
The analogous band calculation uses S. Equivalently h Delta+Delta h^T=0:
the unpaired one-body Hamiltonian and this pairing commute on full Fock
space. This differs from the homogeneous Dirac anticommutation. Thus each
family contributes every
positive frequency

    m_family+E_j, j=1,2,3,4,                            (7)

once for each conjugate internal charge. Their negative Nambu counterparts
are retained and are not extra physical species. The complete positive
charged count is24, with156 further positive spectator Nambu rows in the
localized band. These counts include all spatial, family and charge rows.

Equations (2)-(6) rigorously separate four positive frequency groups:
`m_l-alpha`, `m_l+alpha`, `m_h-alpha`, `m_h+alpha`, with error eta and the
explicit rational sqrt(2) enclosure used by the code. Their respective
positive charged row counts are8,8,4,4. No member of any group is selected
or divided by its weight. The groups are disjoint; there are at least four
positive frequencies, not a certified three-species spectrum. Distinct
splittings inside a group are not inferred from numerical eigenvalues.

The rational bound sqrt(2)<707107/500000 gives

    m_l > 42893/10800000,
    m_l-(1/256+1/25000) > 1/40000.                     (8)

Hence this entire static paired Nambu source has a gap greater than1/40000:
charged band rows obey (8), spectator band rows have modulus at least
alpha-eta, and exterior rows have modulus greater than61 since |E|>=62
and ||K||<1. This is an exact gap in the chosen lattice clock. Neither M64,
the core separation, the coherent phase nor a frequency branch is fitted
to measured particle data.

## Original readout weights and the reference boundary

The existing per-eigenvector core probability is at least
p_core=1/[1+(7/4)/(62-r)^2]. Summing the four exterior-column norms gives
||Q_core P_band||^2<=4(1-p_core). Equal-rank principal-angle identities then
imply, for every normalized core coefficient f,

    ||P_band f||^2 >= 1-4(1-p_core) > 99/100.           (9)

This is an original full retarded CAR weight, with every band row and the
complement retained. It is not a separate addition or removal lower bound.
The CAR identity {a_f,a_f*}=I gives a complete addition/removal instrument.
For Pi=1_[0,1](H-E_ground), the four Kraus maps Pi a_f, Pi a_f*,
(I-Pi)a_f and (I-Pi)a_f* have total effect I on every Fock input.
This is a complete energy readout, not a renormalized band selection.
On the chosen paired quadratic ground, a complete energy projection after
each operation yields the positive excitation interpretation of (9). The
window [0,1] contains every band excitation and no exterior one-particle
excitation, since the latter have energy greater than61. All
other-energy outputs remain part of that instrument. On a general input,
the operator/CAR restriction remains exact but its records need not have
the ground-excitation interpretation. No native execution of this
spectral readout or preparation of that paired ground is asserted.

The unpaired filled sea cannot silently be called the paired ground.
On an actual two-mode pairing block with unpaired levels e,-e, the full
four-state Hamiltonian has energies -m,-e,+e,+m. Its unpaired filled-sea vector
has energy -e and is stationary, but is not a ground when m>e. The phase0
charged band satisfies precisely m>|E_j| by (8). The finite exact control
keeps all four occupation states and verifies electric-charge conservation.
The positive action and closed algebra therefore do not remove the state
preparation obligation.

## Which identification is excluded

The literal [hopping operator](../compact_defect_pair/spectral.py) uses the
same three Pauli matrices in its Weyl differences and its onsite texture
`M n(x).sigma`. This is a two-component source coordinate, with no physical
spin assignment required for the following algebra. For a Hermitian
`beta=b0 I+sum_j b_j sigma_j`,

    {beta,sigma_j}=2 b0 sigma_j+2 b_j I.

Anticommutation with all three kinetic matrices forces every coefficient
to vanish. The conclusion also holds for operator coefficients on unchanged
internal copies. Thus this onsite Pauli texture cannot be a nonzero
Lorentz-scalar Hamiltonian mass for that kinetic block. For uniform north
orientation the actual symbol is
`h(k)=-sum_j sigma_j sin(k_j)+M sigma3`; its square contains
`-2 M sin(k_3) I`, rather than a momentum-independent rest-mass square.
This does not invalidate its proved lattice gap.

Translation fails on the fixed compact texture itself. For the source's
unit core separation, its [direction function](../compact_defect_pair/topology.py)
gives `V(0)=-M sigma3` and `V(e_3)=0`. If `S_3|x>=|x+e_3>`, then on an
interior source column
`[h,S_3](|0> tensor v)=|e_3> tensor M sigma3 v != 0` for nonzero v.
The hopping commutator cancels exactly; this calculation uses no finite-box
boundary. Consequently these localized eigenvalues do not constitute joint
energy-momentum shells of the fixed Hamiltonian. A dynamical collective
position would require its own operator and kinetic law.

This compact source uses the staggered spatial pairing C=epsilon_x i sigma2.
The homogeneous free charged construction uses i sigma2. They are distinct
configured actions. Even replacing the former by the latter would not erase
the texture: the hopping obeys T J=J conjugate(T), while its nonzero onsite
Pauli texture obeys V J=-J conjugate(V). Consequently

    hJ-J conjugate(h)=2VJ != 0,                        (10)

so the homogeneous free Clifford-square cancellation fails on that texture.
The code checks (10) on the actual core-neighbor block, all three Pauli
axes, and a complex symmetric pairing control. The additive spectrum (7)
is not sqrt(p^2+m^2), and a static trap has no spatial momentum dispersion.
Translating the externally supplied texture gives a unitarily related
Hamiltonian, not a dynamical collective-position operator or its kinetic term.

The shell degrees remain properties of a declared background. Creating a
CAR excitation leaves that background unchanged. The exact current/charge
intertwiner does not attach an individual winding sign to an electric charge
or prove that a natural lepton is this texture. The map supplies localized
fermionic modes and a complete source frequency table; formation, motion,
stability, physical clocks and species matching remain separate. Gauge
charge endpoints with dressed CAR paths constitute a different source map,
not a reinterpretation of these two classical winding numbers.
