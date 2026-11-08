# Source-family compatibility for particle quantum numbers

## 1. Reusing the full-ground projection theorem

Use exactly the signed full-Fock source in `strong_coupling_refinement`
and `protected_spin_limit`, on a finite graph with E>=1 edges:

    K=T+UD-(U/2)N,   |t_xy|=c/2,   tau=cE,   U>=4cE>0.

Here D counts doubly occupied sites, N is total occupation and V is the
site count. The comparison operator H=K+UV/2 differs only by a scalar.
The one-particle-per-site space P has fixed N=V. Min-max on P gives
nonpositive low-band energies for H. Every sector with N different from
V has energy at least U/2-tau, so every full source ground has N=V and

    L=K-E_ground >= g Pi_(N!=V),  g=U/2-cE>0.                 (1)

This is the source theorem's exact charged-sector bound. The exchange coefficient on P is J_eff=c²/U. It is the same
coefficient used by the existing ordered spin construction, not a fitted
velocity or physical coupling.

Let R=U/c and let p be the spatial field scaling. Under one common positive
clock factor q, the ballistic spin time coefficient is b=q J_eff p. Then

    qg/b = g/(J_eff p)=R(R/2-E)/p
                         >=R²/(4p)>=4E²/p.                  (2)

Thus any sequence in this domain with p tending to zero and liminf b>0
has qg tending to infinity. This argument does not assume the particular
trajectory c=J E^4, U=J E^8. It applies to every c,U,E satisfying
the displayed full-ground hypotheses. Conversely, bounded charged energy
forces b to zero at least as rapidly as C p/[R(R/2-E)]. This says the
existing ballistic spin clock scale collapses; it does not by itself
assert convergence of another source dynamics when its approximation
errors have not been controlled.

The local dressed theorem in `bipartite_spin_limit/LOCAL_LIMIT.md` has a
different domain. It supplies volume-independent spin dynamics on dressed
inputs at fixed lattice spacing. It does not identify those inputs with
finite-c source grounds or transfer the bare repair instrument. Equation(2)
is not applied to that theorem outside the full-ground gap hypotheses.

## 2. Complete observable weight and fermion grading

For every normalized state rho with Tr(rho qL)<=C, positivity in (1) gives

    Tr(rho Pi_(N!=V)) <= min(1,C/(qg)).                       (3)

No conditioning or thermal reference enters this inequality. The source
fermion parity is Pi_F=(-1)^N. Every odd CAR polynomial changes the parity
relative to the N=V ground, so its output lies in the charged complement
of that ground sector. For an admitted odd creation contraction A with
unconditioned energy W=<A Omega,qL A Omega><=C,

    original creation probability ||A Omega||² <= C/(qg).    (4)

Complete preparation keeps A and its complementary outcome. Dividing by
the vanishing creation probability would change (4) and would not supply
a nonzero original particle residue. Every fixed finite-energy spectral
window is eventually disjoint from these source-odd sectors under (2).

On the bare single-occupancy inclusion, P A_odd P=0 exactly. Its spin
operators are even bilinears c_up* c_down and preserve N. The same charge
separation holds for the dressed low band, because the source conserves
N and its entire low band lies at N=V. Local even spin fields cannot be
renamed the original odd CAR fields. Their regular limiting CCR are
bosonic field algebra relations; neither those relations nor internal
spin prove a Lorentz-spin assignment.

This is not an exclusion of emergent nonlocal or topological fermion
sectors. Such a route needs an explicit sector, observable field operator
and exchange grading. Those data are not supplied by (1), and no claim
about all possible OPH sectors follows from this projection.

## 3. A chemical reference cannot retain both conjugate charges

Keep the same N=V reference and set K_nu=K-nu(N-V). For a charge-one
addition and removal their reference-relative energies transform exactly as

    Delta_+(nu)=Delta_+(0)-nu,
    Delta_-(nu)=Delta_-(0)+nu,
    q[Delta_+(nu)+Delta_-(nu)]
                       =q[Delta_+(0)+Delta_-(0)]>=2qg.       (5)

Consequently both conjugate gaps cannot stay bounded under (2), for any
choice of nu. Centering one branch can move the other in the opposite
direction. If a reference-relative excitation becomes negative, the
reference has ceased to be a ground against that excitation. It cannot
silently be replaced by a new vacuum while retaining this proof's state.
Equation(3) uses the positive ground-subtracted original K; it is not
applied to an indefinite shifted reference generator.

## 4. The finite source is not forbidden to have odd excitations

The producer reconstructs the existing signed two-site CAR source on all
sixteen Fock states, with c=1/8 and U=1. Its complete spectrum proves a
unique centered ground with N=2 and

    E0=-(4+sqrt(17))/8,
    Delta_+=Delta_-=(2sqrt(17)-1)/16,
    Delta_spin=(sqrt(17)-4)/8.

All four local modes have addition and removal probability1/2. Both
conjugate charge sectors and the spin triplet have finite positive
excitation energies on this one finite clock. The exact CAR and parity
identities, the whole spectrum and all creation weights are retained.
This control uses the original action; it changes no Hamiltonian to
produce a desired mass. It demonstrates why the exclusion is a statement
about the specified scaling and protected projection, not about the
existence of fermionic source readouts. These finite gaps are not physical
masses, and the dimer is not a relativistic vacuum.

## 5. What the finite gauge and statistics results actually attach

The existing constructions supply different exact attachments:

| Result | Supplied source attachment | Further physical requirement |
| --- | --- | --- |
| Finite compact gauge bracket and action composition | Classified bracket, chosen invariant metric and explicit coefficient bundle; `GaugeSectorBracketCompletion.lean` retains its parameter premises | The selected source action, state, clock and spacetime readouts must instantiate the same bundle |
| Native M15 matter currents | Additional auxiliary spin lift and nonvacuum even exterior register; fifteen exact internal weights | An interacting spacetime field with physical chirality, statistics and particle sectors |
| Complete gauge transport | Declared fifteen-component matter and sixty-state link, endpoint Gauss intertwiners and all hold outcomes | A source-selected operation schedule, collision law and common propagation clock |
| Native CAR source | Exact odd exchange grading and conserved occupation charge on full Fock space | Its finite-energy survival and relation to Lorentz spin in the chosen continuum |
| Ordered spin field | Actual neutral, source-even spin observables, regular canonical fields and complete history limit | Physical spin/statistics, nonabelian current attachment and common surviving matter sectors |
| Existing chiral trial | Declared all-left spectral kinetic symbol, Higgs tensors, finite CAR execution and clock records | A source-selected interacting chiral regulator, physical vacuum and common gauge dynamics |

The action-composition theorem explicitly keeps PR-35, PR-36, PR-46 and
PR-47 outside its conclusion. A finite gauge theorem is not a claim that
these field attachments are inhabited by the ordered spin source. This package changes no canonical premise or claim classification.

The M15 internal occupation masks all have degree two or four. Their
internal even parity does not establish the exchange statistics of a
spacetime matter field. The source CAR attachment and the auxiliary Spin
lift must be connected on the same operator/state domain. Likewise,
the literal raising link obeys the existing partial-isometry identities;
the separate rotor completion changes its coordinate dynamics. Native
Gauss covariance does not turn that link into a Wilson coordinate.

There is a concrete parameter countercontrol on a fixed extended
gauge source. The producer replays its color-conjugate block, declared dual
reference and complete target/hold detector. The same endpoint Gauss singlet
gives target probabilities16/25 and144/169 for the two existing admitted
rotations (3/5,4/5) and (5/13,12/13). All48 endpoint singlet checks pass,
and both outcomes are retained. Thus finite covariance and completeness
permit distinguishable operation parameters. Reading a specified simulator
reveals which program it runs; it does not make these finite covariance
premises select that program or a capacity-dependent coupling trajectory.

## 6. The remaining common-source route

The original full CAR action and complete field repair exist at
finite c,U. A candidate retaining those original odd fields must either
leave the full-ground projection regime used in (2), or change the target
field/sector with an independently justified grading. Leaving that regime
does not establish a useful continuum. It requires a source-selected
capacity-to-coupling and clock law, the corresponding common ground/state,
controlled response and a nonzero finite-energy sector.

The centered ballistic class with U=o(c^(2/7)) is an existing example that
retains odd free Weyl fields on its own state/clock. Its full hierarchy is
free and it supplies no surviving nonabelian interaction. The ordered
spin branch supplies different fields and a different scaling. Tensoring
these conclusions or relabeling their excitations does not derive a
common interacting action. The exact constraints above identify which
source-family and quantum-number premises must be supplied without
choosing a law to fit measured masses.

## 7. Positive spatial and exchange information on the spin branch

There is a genuine proper-cubic source action, independent of internal spin
rotations. On an even periodic cubic torus, each proper signed permutation
R of the three site coordinates preserves nearest-neighbor bonds and
epsilon_x. Site permutation therefore preserves H_J and its positive
staggered internal-axis-1 selection term. The unique selected finite
ground is invariant as a state. The stated ordered limits preserve this
action. The producer checks all24 rotations on the complete side-four
source graph: all64 sites, all192 bonds and every selection sign remain.
This is an exact graph/operator argument, not a diagonalization of2^64
states. The inherited positive-field uniqueness premise remains explicit.

Site rotations satisfy alpha_R beta_x alpha_R^(-1)=beta_(Rx), because
coordinate parity is preserved and the internal pi rotation in beta acts
independently of the site permutation. At finite p and in the inherited
field limits,

    alpha_R Q_a(F)=Q_a(F composed with R^(-1)),  a=2,3.

The same covariance holds for the bounded charged creator B_q(F) by its
functional calculus and angular integral. A rotation-invariant radial
profile consequently gives a proper-cubic singlet. Pulling this action
back along the binary-cubic double cover makes its central2pi element
the identity; no spinor factor is introduced by that site action.
The operators carry the separately proved **internal** charge.

This finite conclusion does not assign continuum spin0. The exact harmonic
polynomial

    h4(X)=X1^4+X2^4+X3^4-(3/5)|X|^4

has zero Laplacian and is invariant under all24 proper cubic rotations,
yet a continuous rotation changes it. It is a cubic singlet in the
spin4 harmonic representation. This concrete control prevents inferring
an SO3 spin solely from a finite point group. If the same charged field
obtains a strongly continuous **scalar-pullback** SO3 extension, with one
scalar intrinsic fiber, that fiber has spin0: every continuous
one-dimensional representation of connected SO3 is trivial. The existence
and physical attachment of that extension are separate premises.

The local source charge creators are also actual even CAR bilinears,
`A_+=S3-iS2`, `A_-=A_+*` in the physical spin frame. They preserve total
occupation and commute exactly on disjoint sites. The full two-site CAR
calculation verifies both assertions without compressing away empty or
doubly occupied states. The limiting configuration creators inherit their
equal-time commutative field algebra. These facts supply a source exchange
grading and spatial representation information. A relativistic
spin-statistics identification still needs the same continuum spacetime
action, positive spectrum, locality and field/sector attachment; it is
not proved by an onsite spin-one-half label.


## 8. Exact interfaces to the sector and Spin theorems

`QuantumMatterIntegration.FiniteQuantumMatterBridge` requires an actual
projective partition, nonzero component projectors, a component-to-sector
map, a compatible central-weight map, the declared block-pinching readout
and a chiral anomaly-free selection mask. Its coordinate-line construction
proves that this interface is mathematically inhabited. It does not
identify those lines with the selected spin-source excitation sectors.
The source U1 charge readouts and the M15 register are concrete inputs to
a possible bridge; matching their labels without the required operator
maps does not instantiate it.

`ExteriorComponentBridge` proves creation anticommutation in the supplied
five-mode exterior algebra. `MatterGrammarIndexBridge` matches its table
to the frozen Weyl table by an explicit permutation.
`FermionSectorAssembly` combines finite complex tuples and declared
coefficients into an invariant functional. These theorems retain the
physical Spin/chirality attachment as a premise. None identifies the
source-even spin creators with spacetime Weyl fermions.

`SL2F5PortCover` gives the canonical finite cover of the committed port
permutation group and its two-element central kernel. It explicitly does
not identify that cover with the executable auxiliary spin lift, construct
a complex spinor representation or attach it to physical rotations. The
proper-cubic site action in Section7 is a source representation with a
trivial central element under pullback; the finite port-cover theorem does
not replace it by a spinor representation.

`SectorInheritance` proves that a bijective endomorphism of the finite
full matrix algebra is inner and records a correspondingly restricted
localized equivalence. To obtain physical superselection sectors one
needs the selected vacuum representation, the actual local observable
net, transportable localized endomorphisms and their intertwiners. The
finite innerness theorem neither constructs this category nor rules out
sectors in a limiting net.

The compact-group reconstruction in the technical supplement requires a
transportable refinement-directed sector colimit, a rigid symmetric
C*-tensor category and a faithful bosonic fiber functor. It first yields
some compact group. The SU3, SU2 and connected U1 conclusion additionally
consumes the category/current identification and the separate local Lie
type and charge certificates. The ordered source has actual field,
positive-energy and spatial-action constructions. Attaching its sectors
to this category and its currents is the precise missing reuse step;
no absence claim about those existing source constructions is needed.
