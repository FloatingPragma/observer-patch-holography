# Spatial Weyl hopping with the same fifteen-state gauge representation

The fifteen internal states and twelve currents are those of
[`native_matter_attachment`](../common_source/native_matter_attachment.py).
They are constructed from the declared native current, spin lift and selected
even exterior register. The separate physical canonical anticommutation
relations and spatial two-component spinor are explicit architecture.
The auxiliary even occupation parity does not supply physical fermion parity.

Three exact constructions have different field domains. The same-current
kinetic/Higgs intertwiner in `native_matter_attachment.kinetic_attachment`
handles arbitrary supplied complex Yukawa coefficients. The
[`native_gauge_matter_transport`](../common_source/native_gauge_matter_transport.py)
compiler implements background-linked physical CAR hopping, and
`native_gauge_continuum` gives its conditional constant-background limit with
two Dirac tastes and both chiralities. The
[`native_quantum_gauss`](../native_quantum_gauss/PROOF.md) construction has
operator-valued finite link registers and an exact local Gauss law. The
attachment below equips that last Hamiltonian with spatial Weyl hopping.
It preserves its configured link representation and all fifteen components.

## One representation and a full-Fock Gauss identity

Write the internal space as `M15=direct_sum_r V_r`. Its source-computed block
dimensions are `6,3,3,2,1`. Let `G_A`, `A=1,...,12`, be its anti-Hermitian
native current matrices. The executable compares these matrices with the
literal spin/exterior lift in the selected occupation-mask order. Physical
Hermitian current matrices are `i G_A`, as in the background compiler.
The spatial spin matrices `sigma_j` act on a separate `C2` factor and commute
with every `G_A`.

At each vertex use thirty physical CAR modes `c_x,s,a`, where `s=1,2` and
`a=1,...,15`. For an oriented edge `x -> y`, retain the owning link space

    L_xy=C|0> direct_sum_r(V_r,x* tensor V_r,y), dim L_xy=60,
    U_ab=|bar a,b><0|,                    a,b in one block,
    L_x(G)=0 direct_sum_r(-G_r^T tensor I),
    L_y(G)=0 direct_sum_r(I tensor G_r).                  (1)

Cross-block `U_ab` vanish. For real `c>0`, `epsilon>=0`, and a specified
edge direction `j`, define

    W_j=i c sigma_j/2,
    H_xy=sum_stab c_x,s,a^dag (W_j)_st U_ab c_y,t,b
             +adjoint+epsilon(I-|0><0|),
    Q_x(G)=sum_sab c_x,s,a^dag G_ab c_x,s,b,
    G_x(G)=Q_x(G)+sum_incident_edges L_x(G).             (2)

Every finite sum of these terms is Hermitian on the full finite CAR/link
carrier. The local generators satisfy the same Lie brackets as the native
currents, and generators at different vertices commute. For every vertex,
axis and native current,

    [G_x(G_A),H]=0.                                     (3)

The CAR identities give `[Q_x(G),c_x,s,a^dag]=sum_i G_ia c_x,s,i^dag` and
`[Q_y(G),c_y,t,b]=-sum_i G_bi c_y,t,i`. The literal link commutators are

    [L_x(G),U_ab]=-sum_i G_ai U_ib,
    [L_y(G),U_ab]= sum_i G_ib U_ai.                      (4)

Reindexing the sums in (2) cancels (4) against the two fermion charges;
`W_j` is unaffected. Adjoint terms obey the same identity. The occupied-link
projector commutes with both endpoint actions. This proves (3) on every
occupation sector, with no one-fermion hypothesis. Gauge-singlet onsite
interactions or loop contractions can be added if their invariance is proved;
their values and inclusion are additional action data.

Any supplied finite family factor tensors the same current with its identity
and repeats the hopping sum over that factor. The proof is unchanged; it
does not select a family count. The executed control uses one M15 copy.

`ward_certificate()` extracts the independent coefficients of
`c_p^dag c_q tensor |l><k|`. Their array has dimension `3600=60*60`, whereas
the two-vertex physical Fock space has dimension `2^60`. Its vanishing is
a faithful test of the operator identity: all `c_p^dag c_q` are linearly
independent, as their one-particle matrix elements show. It is not a
truncation of the represented operator to one particle. All 72 combinations
of twelve currents, two endpoints and three spatial axes are checked.
The tests also apply the full sixty-mode CAR action to occupied spectator
states. Reversing the link representation sign breaks (3).

Exponentiating the finite generators gives local gauge covariance. Every
spectral function of the same `H`, and every measurement Kraus operator
commuting with all `G_x`, preserves their joint kernel. This construction
does not identify local gauge redundancy with a primitive repair theorem.
Finite quantum-link registers have the context of
[Chandrasekharan and Wiese](https://arxiv.org/abs/hep-lat/9609042);
the identities here are verified for the displayed source register.

## A dressed spatial-spin transport on this Hamiltonian

Supply the owning static dual reference at `y`. For each block of dimension
`d`, and each spatial spin `s`, retain the full contracted states

    A_s=sum_b |y,s,b>|bar b>_reference|0>,
    B_s=sum_ab |x,s,a>|bar b>_reference|bar a,b>_link.    (5)

They have norms squared `d,d^2`, and every local Gauss generator annihilates
them. The four-dimensional space they span reduces the same Hamiltonian:

    H A_t=sum_s (W_j)_st B_s,
    H B_s=d sum_t (W_j^dag)_ts A_t+epsilon B_s.          (6)

In the basis `A_up/sqrt(d),A_down/sqrt(d),B_up/sqrt(d),B_down/sqrt(d)`,
the Gram matrix is `diag(1,1,d,d)` and

    H_r=[[0,d W_j^dag],[W_j,epsilon I2]],
    H_r^2-epsilon H_r=(d c^2/4)I4.                     (7)

The two finite levels, each of multiplicity two, are
`(epsilon +/- sqrt(epsilon^2+d c^2))/2`. The separation is a finite
matter/link/reference transport frequency. It is not a particle rest mass.
The code checks (6) by operating the owning sparse link Hamiltonian on all
sixty columns across five blocks and three axes, with 480 literal Gauss
state checks. Neither unwanted source blocks nor occupied link levels are
discarded to obtain this reduction.

Prepare the complete M15/reference singlet with spatial spin up and link
vacuum, then measure its internal block. The original block probability is
`d/15`. At `c=epsilon=1`, apply the exact same-Hamiltonian Cayley step at
`tau=1/8`, and measure position and spatial spin. All twenty one-fermion
outcomes per axis are retained, including zeros. They sum to one; each block
sums to `d/15`. There is also a complete-Fock complement: writing `P_1` for
total fermion number one, the Kraus maps are `P_x,s U P_r P_1` and
`U(I-P_1)`. Here `P_1` is the full one-fermion projector, with identity on
every link and reference level, not a projector onto (5). Their effects sum
to identity on the full carrier. The last
outcome has exactly zero original weight in this preparation and is retained
as the twenty-first record. All these maps commute with the local Gauss
generators, and `U` is the same Cayley clock in every occupation sector.
These ideal finite services are specified inputs. No native pulse
execution, prepared gauge vacuum, detector reservoir or continuum limit is
claimed by this calculation. The energy of `H` is preserved by the Cayley
step; a subsequent position/spin measurement need not preserve it.

## The entire local anomaly tensor on M15

For an all-left Weyl interpretation of one M15 representation, the linear
and symmetric cubic algebraic coefficients are

    ell_A=Tr_M15 G_A,
    D_ABC=Tr_M15 G_A(G_B G_C+G_C G_B).                  (8)

`anomaly_certificate()` verifies all twelve `ell_A=0` and all 364 unordered
triples `D_ABC=0`. Cyclicity of trace makes the cubic tensor fully symmetric,
so this checks every ordered component. Multiplying currents by `i` to
use the Hermitian convention leaves every zero unchanged. The individual
blocks have nonzero abelian cubic traces; the cancellation uses the entire
fifteen-state representation. There are four weak doublets counting color,
consistent with the even parity criterion for the
[SU(2) anomaly](https://doi.org/10.1016/0370-2693(82)90728-6).
The two spatial spinor components are not two independent chiral flavors.

The same selected-mask arithmetic is certified in
[`ExteriorSelection.lean`](../../Lean/Screen/ExteriorSelection.lean),
`parity_sectors_survive` and `witten_automatic`; the explicit matrix tensor in
(8) is an executable attachment, rather than a Lean evaluation of those
matrices. These are representation-level cancellations. A regulated
interacting chiral quantum measure, all global anomaly questions for a
chosen quotient gauge group, and a mirror-decoupling theorem do not follow
from (8). Finite-dimensional Gauss covariance itself holds for arbitrary
representations, including ones with a nonzero tensor (8).

## Separation from the relativistic continuum premise

The coefficient `i c sigma_j/2` is the literal spatial Weyl coefficient of
the cubic source. With supplied classical identity links it gives the free
symbol `-c sum_j sigma_j sin k_j` and its eight valleys. That statement
does not describe the spectrum of (2): here every `U_ab` is a nonzero
nilpotent raising operator with `U_ab^2=0`. It is not a unitary Wilson
link and admits no eigenvector with nonzero eigenvalue. In particular, its
vacuum cannot be set equal to an identity link in the Hamiltonian or state.

Equation (3) supplies a rigorous finite gauge interaction on the same
physical spinor/internal CAR representation. A relativistic interacting
limit also requires a compatible gauge-field state and scaling, control
of the quantum links, and a regulator with the desired chiral spectrum.
Neither the background two-taste construction nor the free eight-valley
construction establishes these properties for (2). The declared quantum
links, spatial graph, physical CAR functor, coupling and clock are explicit
premises. No numerical coupling, physical species or mass is selected by
this finite join.
