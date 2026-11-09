# Cubic energy density and an exact boundary improvement

## Same physical action, state and clock

Use the [primary contract](../particle_common_contract/PROOF.md), with its
open Wilson action, all sixty matter modes per cell, product link cutoff,
Gauss polynomial density, idle conjugate reference and source clock.
The [physical cube action](../particle_primary_channel_labels/PROOF.md#the-physical-cube-action-and-the-open-boundary)
acts on the ordinary open interior graph. Every dangling link has a
one-link ghost Gauss constraint, which projects it to the constant
function. Its electric operator annihilates the physical carrier.
Deleting that zero operator from a physical identity changes neither
the Hamiltonian nor its original energy zero. It is not an ambient
operator identity before the ghost projection.

Write M_x for onsite matter, E_b for a complete interior electric term,
B_b for a complete Hermitian hopping, and P_p for the signed magnetic
plaquette term. For geometric endpoints V_b and all four vertices V_p set

    h'_x = M_x + (1/2) sum_(b: x in V_b) (E_b+B_b)
                 + (1/4) sum_(p: x in V_p) P_p.                 (1)

Every whole term is the same compressed parent operator as in the source.
Each has total share one; hence sum_x h'_x=H on the physical space.
No interaction, matter state, cooling law, clock or reference is replaced.
Equation (1) specifies a local energy convention, not a dynamical metric.

The finite binary-cube representation maps M_x to M_Rx, edges to edges
and real plaquettes to real plaquettes. Edge reversal uses U to U*;
its electric Casimir and complete Hermitian bond transform accordingly.
The product cutoff intertwines this action. Therefore

    R h'_x R* = h'_(R x),
    R J'_xy R* = J'_(R x,R y),    J'_xy=i[h'_x,h'_y].           (2)

The central binary element is CAR parity and acts trivially on these even
operators. Thus (2) is a scalar-energy/current covariance under the24
proper cube rotations. No continuous rotation or Lorentz tensor is inferred.

## Exact improvement of the owner allocation

Let h_x be the [owner allocation](../particle_primary_stress/PROOF.md):
electric energy belongs to the positive-link owner and plaquette energy
has three owner shares1/3. Put delta_x=h'_x-h_x on the physical carrier.
For each positively oriented interior edge u to v route E_b/2 from u
to v. For a plaquette with vertices a,b=a+e_i,c=a+e_j,d=a+e_i+e_j,
route P_p/24 along a to b and a to c, and P_p/8 along b to d and c
to d. Add the contributions to an antisymmetric operator I_uv=-I_vu.
Then exactly

    delta_x = sum_y I_yx,       sum_x delta_x=0.                (3)

The plaquette divergence is (-1/12,-1/12,-1/12,+1/4) at(a,b,c,d),
precisely the change from its three-owner to four-vertex allocation.
All improvement routes follow nearest edges. For any set R of cells,

    H'_R-H_R = sum_(u outside R,v in R) I_uv.                  (4)

Every interior route cancels. Thus the difference is an explicit boundary
operator, not a change in total energy or a discarded interaction.
The word "boundary" is an incidence statement here; its norm can grow
with cutoff and boundary size. Equation (4) does not prove a negligible
continuum improvement or an entropy area law.

The current change is also exact:

    J'_xy-J_xy = i[h_x,delta_y]+i[delta_x,h_y]+i[delta_x,delta_y]. (5)

Summing y gives i[delta_x,H]. Thus the two density/current pairs satisfy
the same source continuity equation with their respective regional
energies and full boundary fluxes. The operator-valued endpoint increment
and integrated antisymmetric current instantiate
[`CoarseEnergyWard`](../../Lean/QFT/CoarseEnergyWard.lean); its finite sum,
coarsening and time telescoping proofs apply without modification.
The actual primary matrices are not formalized as Lean values by this leaf.

## Supports, norms and complete records

Although all four plaquette vertices receive energy, its operator support
is its three positive-link owners. Every owner in h'_x lies within graph
distance two of x. Even disjoint operators commute, so J'_xy=0 at
distance greater than four. The full-carrier estimates give

    ||h'_x|| <= K_n=[180+(1827/2)ell_n^2+3/7]/a_n,
    ||J'_xy|| <= 2 K_n^2.                                    (6)

There are at most six incident electric terms, six bonds and twelve
plaquettes, with shares1/2,1/2,1/4. Source norms are90/a for matter,
609ell^2/(2a) for electric,30/a for a bond and1/(7a) for a plaquette.
The code evaluates boundary-specific supports and caps without allocating
the source Hilbert space. It independently enumerates source incidence,
all24 rotations and every improvement coefficient on small open boxes.
The proof establishes the formulas at arbitrary finite side.

For a Hermitian energy or current A with M>=||A||>0 use the complete
square-root instrument for effects (I+A/M)/2 and(I-A/M)/2, retaining
both labels +M,-M. Their mean is Tr(rho A). Zero operators admit an
identity outcome labeled zero. Use invariant caps such as(6) when
comparing rotated instruments. Functional calculus makes these
instruments transform by the same R as their observables. The selected
rho commutes with R and H, so complete rotated histories have identical
classical record distributions on the same source clock. Their retained
quantum outputs are related by R on the source and conjugate(R) in the
canonical-reference comparison; no reference gate is part of the experiment.
The source net regional flux has zero stationary mean. Individual pair
currents are not asserted to have zero mean. After interventions the
operator continuity law survives with the actual disturbed source state.

For a declared real scalar control f_x(t), put H_f=sum_x f_x h'_x.
It has the exact spatial transformation R H_f R*=H_(f composed with R^-1)
and H_1=H. The same finite source work theorem gives

    d<f_x h'_x>/dt + sum_y f_x f_y <J'_xy> = dot(f_x)<h'_x>,
    d<H_f>/dt = sum_x dot(f_x)<h'_x>.

These equations supply a cube-covariant energy coupling, including its
applied work, on the actual carrier. The prescribed f is an apparatus
control; identifying it with a dynamical proper lapse is a separate claim.

Owner-cell partial trace is a different operation. Rotation can move a
positively owned link into another owner cell. Therefore covariance of
H'_R does not by itself prove covariance of the owner's regional density
or assign an entanglement region to the geometric set R.

## Conjugate records and regional response

The [literal unitary charge conjugation](../particle_primary_charge_conjugation/PROOF.md)
fixes every complete source term individually. It is linear, so it
preserves i in the commutator. Consequently C h'_x C*=h'_x and
C J'_xy C*=J'_xy. All specified square-root instruments and source waits
commute with C. Conjugate input states therefore have identical complete
classical energy/current histories, with covariant retained quantum outputs.
Additional noncovariant detector back-action is outside this statement.

The site CAR conjugator is even because every site has60 modes; link
coordinate conjugation is local to its owner. Thus C factors across any
fixed ambient owner region and its complement, up to an irrelevant scalar.
On the selected state C rho C*=rho. Its supported regional marginal sigma_R
is C_R-invariant, as are its support and K_R=-log sigma_R on that support.
For a supported trace-zero Hermitian tangent X_R with
C_R X_R C_R*=-X_R, the finite regional first law gives

    dS(sigma_R+tX_R)/dt|_0 = Tr(X_R K_R)=0.                  (7)

Any C-odd admissible source tangent also has zero first-order expectation
of h'_x and J'_xy. This is a symmetry selection rule, not a vanishing
quadratic response or zero energy of a charged excitation. Physical
one-parameter paths must preserve positivity and the actual Gauss support.
No arbitrary odd matrix is asserted to be a preparable tangent.

The same-state first law, energy continuity and cube covariance are finite
premises for a stress reconstruction. They supply no proper geometry,
geometric diamond modular kernel, fixed-volume area response, Einstein
equation or refinement-uniform stress estimate.
