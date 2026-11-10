# A finite open-cube field and native Dirac current

## Common source interface

The spatial complex is the boundary and interior of one open unit cube:
eight vertices, twelve edges and six oriented faces. The vertex order is
lexicographic on `{0,1}^3`; each edge goes from x to x+e_j. D has -1 at
the tail and +1 at the head. A face follows axis i, axis j, minus axis i,
minus axis j. `field.py` constructs these integer matrices and verifies
CD=0, rank D=7 and rank C=5. The oriented volume boundary is
`D2=(-1,1,-1,1,-1,1)` in the executed face order and satisfies D2 C=0.
Thus the magnetic field B=CQ has exactly zero magnetic divergence on
every retained photon state. This geometry is a declared spatial
architecture. No equality with the thirteen-vertex cone or the original
PW119 source Hamiltonian is asserted.

The matter provider supplies the complete native carrier S, its reducing
four-mode canonical code at each of the eight vertices, its full
Hamiltonian and the complete outside-code sectors. Let P_code project
onto the global code. Four native constituents per vertex give a full
outer port charge Q_native,x with integer range [-4,4]. The extended
free Bogoliubov unitary U_FW is identity on the global code complement,
is charge preserving and acts on the code by the verified finite Wilson
transformation. Define

```
rho_x=U_FW Q_native,x U_FW^dagger.
```

These commuting charges have full range [-4,4] and sum to the same
conserved total native port charge. On the code,
`rho_x=Psi_x^dagger Psi_x-2`, with range [-2,2]. The full complement
range determines the electric and reference capacities. The Gauss
embedding does not use the smaller selected-code range on other sectors.
The matter provider must establish these exact identities and identify
S and every native reference; the field provider does not infer them
from a Dirac label.

Let r=0 be the retained compensator position. There are seven independent
longitudinal labels `ell_x in {-4,...,4}`, x!=r, and a static reference
charge `R in {-32,...,32}`. Their dimensions are `9^7` and 65. The full
static reference Hamiltonian is `4|R|`, retaining both charge signs at
equal declared energy. The net-neutral experiment has R=0 throughout;
a net-charge-one source sector retains reference cost 4. The full
field photon carrier is

```
F_3 = direct_sum_(n=0)^3 Sym^n(C^5),  dim F_3=56.
```

Thus the kinematic carrier, before preparation counters and detector
registers, is `S tensor C^(9^7) tensor C^65 tensor F_3`. No electric label,
native complement or photon-cap boundary is removed from it.

Set `L=D^T D`, `G=(L+ones)^(-1)-ones/64` and
`PT=I-D G D^T`. On the full label carrier define

```
rho_ell=sum_(x!=r) ell_x (delta_x-delta_r),
Ebar_L=D G rho_ell,
Gamma_x=(D^T Ebar_L)_x-rho_x-R delta_(x,r).
```

The simultaneous Gamma kernel is exactly
`ell_x=rho_x` and `R=-sum_x rho_x`. These integer constraints are
independent of coupling g. They retain the physical carrier at g=0;
using the kernel of an identically zero unnormalized constraint there
would admit a different source. If Pi_rho are the joint matter-charge
projectors, the exact physical isometry is

```
Y psi=sum_rho Pi_rho psi tensor |rho_(x!=r)> tensor |-sum rho>.
```

The labels are unique for every joint charge tuple in the entire native
range. Hence `Y^dagger Y=I` and `YY^dagger=P_(all Gamma=0)`, including
all source charge degeneracies. On every retained matter sector,

```
D^T Ebar_L=rho+R delta_r,
C Ebar_L=0,
||Ebar_L||²=rho_ref^T G rho_ref,
rho_ref=rho-(sum rho)delta_r.
```

This is the outer integer-Gauss projector. Y also intertwines the
existing native Gauss constraints; it does not discard native
Gauss-violating kinematic states. Those states have zero weight in the
physical preparation and remain in the native code complement.

Every endpoint charge-changing matter matrix element is accompanied by
the same signed longitudinal shifts. The finite shifts are partial
isometries, zero at invalid capacity moves, and satisfy
`[ell,S_plus/minus]=plus/minus S_plus/minus` on the full label carrier.
No modulo wraparound is used. A tail-to-head matter hop changes rho by
the corresponding column of D and changes ell by its nonreference
entries. A net charge pump also changes R by the opposite net charge.
These give zero change of every Gamma on the complete kinematic carrier,
including invalid moves. `gauss.py` checks all orientations and 157464
two-endpoint boundary cases. Physical moves inside the actual source
charge spectrum stay inside the label range; a forbidden source charge
move has zero matrix element.

## Finite transverse field and its exact boundary

The transverse stiffness `K=C^T C` has spectral projectors

```
P4=(6 PT-K)/2,  rank P4=3,
P6=(K-4 PT)/2,  rank P6=2.
```

There is no harmonic one-form on this contractible complex:
`ker C=range D`. The constant vertex mode is a gauge redundancy and is
killed by G; it is not a photon oscillator. No periodic global-flux
sector is silently dropped. A periodic construction requires its own
retained harmonic variables and is a different member.

Choose any real orthonormal frame V for P4 and P6, and let
`nu=(2,2,2,sqrt6,sqrt6)`. The executed frame is exact algebraic. Define
literal total-cap ladders a_j and

```
Q_e=sum_j V_ej (a_j+a_j^dagger)/sqrt(2 nu_j),
P_T,e=i sum_j V_ej sqrt(nu_j/2)(a_j^dagger-a_j),
B_f=(C Q)_f.
```

Q is the dimensionless link coordinate `a A_phys`; P is the electric
flux `a² E_phys`. On the full photon carrier the chosen Hamiltonian is

```
Hph=[Tr(nu)/2+dGamma(nu)]/a
   =[sum_e P_T,e²+sum_f B_f²]/(2a)
      + P_(N=3)[dGamma(nu)+Tr(nu)]/(2a).
```

The positive completion is part of the declared Hamiltonian and retains
the top-cap states. The actual photon ground is the unique Fock vacuum,
with energy `(3+sqrt6)/a` and gap `2/a`. Literal finite matrices verify
the full Heisenberg equations for Q and P under this completed Hph.
This imports the finite-cap completion identity from
`model/affine_gauss_source`, with the actual cube modes substituted and
replayed. It does not import the cone spectrum.

The photon ladder algebra on total cutoff K is exactly

```
[a_i,a_j^dagger]=delta_ij I-P_(N=K)(delta_ij I+a_j^dagger a_i).
```

It is checked on all 56 states at K=3. Distinct Q_e generally do not
commute on the top boundary. `field.py` supplies a nonzero matrix element
of `[Q_0,Q_1]` and verifies that all its nonzero entries lie at N=3.
This boundary does not change the Abelian charge/Gauss algebra, whose
generators act on matter and longitudinal labels. It does exclude an
unrestricted canonical-commutator or commuting-Weyl-link claim. All
Peierls operators are literal finite exponentials `exp(i g Q_e)`.

Total occupation cutoff is invariant under U(5) rotations of the mode
space. A coordinate with vacuum variance c reduces to a single-mode
coordinate with coefficient sqrt(c), with spectator occupation fixed.
The active blocks have sizes 1,2,3,4. The largest squared eigenvalue of
`a+a^dagger` on these blocks is `3+sqrt6`. Thus

```
<0|Q_e Q_f|0>=(P4/4+P6/(2sqrt6))_ef,
c=<0|Q_e²|0>=1/16+1/(12sqrt6),
||Q_e||²=(3+sqrt6)c < 9/16.
```

The complete cutoff has no additional oscillator states. Vacuum words
of degree at most six match the unbounded oscillator comparison because
a returning path cannot pass level three. Longer words and finite Weyl
operators are evaluated on the finite carrier or bounded explicitly;
Gaussian identities are not transferred without such a bound.

## Coupling and orientation

One coupling g enters the electric flux, Coulomb term and Peierls phase:

```
P=g Ebar_L+P_T,
H_C=g² rho_ell^T G rho_ell/(2a),
W_e=exp(i g Q_e).
```

The matter source supplies the open Wilson matrix with its actual mass
m and a declared spacing a=1/m. For an edge in direction j,
`B_e=-(beta+i alpha_j)/(2a)`. The positive Peierls phase multiplies

```
K_e=Psi_head^dagger B_e^dagger Psi_tail W_e
```

and the negative phase multiplies K_e^dagger, with the corresponding
longitudinal shifts on the full kinematic carrier. Reversing the phase
relative to this D convention would violate the parent gauge
orientation. The exact matter charge flow is
`I_flow,e=-i(K_e-K_e^dagger)`, so `dot rho=D^T I_flow`. The derivative
current is `J_phase,e=+i(K_e-K_e^dagger)=-I_flow,e`; the first interaction
term is `g sum_e Q_e J_phase,e`. The scalar square-root
kinetic current is not used for this source. The finite Q commutator
boundary is retained in any electric-field equation; a continuum
Ampere equation is not asserted on the full cutoff.

On the selected physical code, Y identifies the same finite coupled H
with the Wilson matter matrix, Hph, H_C and the literal Peierls bonds.
The complete outside-code Hamiltonian and its current couplers are
specified by the native matter provider; all such sectors remain in the
source. A convenient declared choice makes the Dirac couplers zero on
the global code complement while retaining its native Hamiltonian and
charge. This is a global spectral-code architecture and does not prove
local tensor causality on the full native carrier.

The common source chooses the complete kinematic extension
`H_ext=Y H(g) Y^dagger direct_sum 0` on the outer Gauss-orthogonal
sectors. The unchanged native Hamiltonian on the logical native code
complement is distinct from this zero outer-constraint extension.
The full partial-shift controls prove physical intertwiners and boundary
consistency; they do not assert that an unguarded partial-shift
Hamiltonian has been adopted on invalid outer-Gauss sectors. Every
sector remains in the carrier and the final outer-constraint flag.

## Preparation and detector contract

The comparison input is the free positive-energy matter vacuum times
the finite photon vacuum, physically embedded by Y. It generally has
local rho fluctuations and entangled longitudinal labels despite total
charge zero. The actual coupled H(g) need not have this ground state.
The prepared input is explicitly nonstationary under H(g). No ground
placement or actual coupled-vacuum claim is imported from the free
matter or free field.

The entire native preparation instrument is lifted through Y on every
retained record, with identity completion outside its physical image.
The pre-code native Fock-vacuum seed is in the global code complement,
where U_FW is identity and all port charges vanish; with ell=R=0 it is
physical. The matter provider must verify this seed condition. The
lifted preparation preserves every outer Gauss constraint as well as
the native constraints. Its source copies, spin references, counters,
success probabilities and state error are retained. Exact preparation
and an efficient circuit are not inferred from a spectral formula.

For globally compressed couplers, an exact support argument can avoid
an initial ground filter. The deterministic native cube trial acts as
identity on every added local port reference R_loc. The matter provider
verifies that the intersection of their shared-vacuum support with the
global code is exactly the code vacuum `g^tensor32`. Therefore the
projection of the actual trial onto the code is
`g^tensor32 tensor chi`, with code probability at least `p0^32`, where
p0 is the certified single-source trial-ground overlap. Both Y and the
full coupled Hamiltonian reduce the code. The final code/count effects
are block diagonal, so their target probabilities are exactly the code
weight times the ideal code probabilities. The complete trial error
has its own final native-complement record. This reasoning uses the
actual source support; it neither prepares exact g on every branch nor
adds an unrecorded initial P0 measurement. A filtered preparation and
its fidelity bound form an independent alternative control.

On the global code the complete detector records
`(N_b,N_d,N_ph)`, with `N_b,N_d=0,...,16`, `N_ph=0,...,3`. The native
code complement is a separate retained outcome, with all four photon
counts. These physical count projectors give `(17²+1)*4=1160` labels.
One additional Gauss-orthogonal outcome completes their sum to the full
kinematic identity, so the retained pointer dimension is 1161. That
constraint-failure outcome has exact probability zero on the physical
preparation. The selected event `(1,1,1)` is production of one
particle, one antiparticle and one transverse quantum from the specified
free input. It is a charged-current/field interaction record in a net
neutral sector, not a weak or electromagnetic scattering cross section.

Let Pi_alpha be these detector projectors, including the retained
Gauss-orthogonal outcome, and rho the actual
prepared density, including retained references. All records are

```
preparation failure flags,
q tr(Pi_alpha exp(-iT H(g)) rho exp(iT H(g)) Pi_alpha), every alpha.
```

Their sum is one because every preparation flag is retained and the
Pi_alpha sum to the identity on the full source. A finite readout uses
`Uread=sum_alpha Pi_alpha tensor Shift^alpha` with a retained pointer,
a finite Hermitian logarithm L and the declared source-clock control
`V(s)=exp(-isH(g)) L exp(isH(g))` for one source unit. The total readout
is `exp(-iH(g)) Uread`; the record refers to its interval start.
The same autonomous H(g) acts on all no-operation/failure branches.
Controller work and synthesis workspace are separate resources.

A nonzero first current coefficient and the complete finite-time and
preparation error inequalities are supplied by the matter/current
calculation on these exact field matrices. Free Maxwell propagation or
a separate Coulomb energy identity alone does not certify that coupled
record. The source U(1) is named by its exact charge and current;
electromagnetic identification is a separate empirical comparison.
