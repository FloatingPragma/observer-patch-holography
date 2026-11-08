# Complete two-charge response on the selected source

Use the actual eight-site cube, its unique selected ground Ω, Hamiltonian
H and clock from `native_ordered_ground`. The source charge and local
creation operators are

    Q=(1/2)sum_x epsilon_x X_x,
    A_q(x)=Z_x(I-q epsilon_x X_x)/2,       q=±1.             (1)

`protected_defect_catalog` proves [H,Q]=0, QΩ=0, the full-source repair
attachment, [Q,A_q(x)]=q A_q(x), and commutation of disjoint local
creations. These are internal source charges, without a particle species
or electric-charge assignment. Reading the source implementation determines
these operators; it does not require fitting an effective field action.

## 1. All two-charge preparation outcomes

For x≠y, the four Kraus operators A_q(x)A_r(y) form a complete instrument.
Apply the local X pointer and Z creation at x and y, retaining both local
outcomes. Its branch vector has charge q+r and original weight

    w_qr(x,y)=||A_q(x)A_r(y)Ω||².                            (2)

No fine all-site X measurement is inserted into this preparation. Such a
measurement would dephase the coherent state inside its charge sector.
For the separated cube corners x=0,y=7, the four weights are enclosed by

| First charge | Second charge | Original probability |
| --- | --- | --- |
| -1 | -1 | [0.120159438343,0.120159438345] |
| -1 | +1 | [0.715315632586,0.715315632587] |
| +1 | -1 | [0.044365490725,0.044365490726] |
| +1 | +1 | [0.120159438343,0.120159438345] |

The two zero-total-charge records are distinct preparation histories and
are retained separately. Conserved total charge does not turn them into
one coherent branch.

## 2. The pair preparations fill an actual reducing source sector

For q fixed, enumerate all28 pairs x<y and put

    B_q e_(x,y)=A_q(x)A_q(y)Ω,
    G_q=B_q*B_q,    M_q=B_q*H B_q.                           (3)

Each column has norm at most1 and charge2q. The exact Walsh/X basis of
the full source gives charge multiplicities binomial(8,4+2q)=28.
Its actual integer Hamiltonian block is read from W_q* H W_q/256;
the certificate verifies H W_q=W_q H_q and retains all source labels.
No eigenstate or spectral window is selected to create this block.

The complete source eigensystem certificate bounds the ground proposal
by ||g_hat-Ω||<=e_g. Project its pair columns onto the exact charge2q
space only for this error calculation. Their omitted proposal norms are
bounded by e_g²; the true vectors have zero omitted norm by (1).
Both true and proposal columns have norm at most1, and each column error
is at most e_g. Therefore every Gram entry changes by at most2e_g and

    ||G_q-Ghat_q|| <=56 e_g.                                (4)

The proposal Gram entries are exact rational numbers. Fraction-free
Sylvester elimination proves all28 leading principal minors of
Ghat_q-10^-9 I positive. Every integer division and every pivot sign is
checked. Consequently

    G_q >= gamma I,    gamma=10^-9-56e_g >0.                 (5)

The receipt gives gamma and all integer principal minors. This proves
that the actual pair vectors are independent. Their28-dimensional range
is the full Q=2q sector, hence is reducing for the same H.

There is therefore a closed source response, without a discarded
orthogonal memory sector:

    W_q=B_q G_q^-1/2,    h_q=G_q^-1/2 M_q G_q^-1/2,
    W_q*W_q=I,    H W_q=W_q h_q,
    K_q(t)=B_q*exp(-itH)B_q
          =G_q^1/2 exp(-it h_q)G_q^1/2.                   (6)

In nonorthogonal preparation coordinates its generator is G_q^-1 M_q.
The small Gram floor makes stable numerical inversion costly; no
conditioning problem is removed by renaming these states particles.
The receipt computes complete kernels by the exact reducing source block,
which avoids numerical inversion of G.

This dimension identity is specific to the eight-site source. On12
sites, the Q=2 sector has495 states while there are only66 local pair
preparations. The eight-site argument cannot assert closure of that larger
space, a two-particle Fock structure, or an asymptotic particle count.

## 3. Reusing the finite algebra and instrument theorems

The exact Hamiltonian graph of each charge2q block is connected. The
receipt provides a55-vertex walk through every one of its28 basis labels,
using only nonzero source matrix entries. Complete final local-X records
give every diagonal basis projector P_i. If H_ij≠0, then

    P_i H P_j = H_ij |i><j|.                                (7)

Division by that explicit nonzero scalar extracts the matrix unit.
Products along the source walk and adjoints yield every matrix unit.
This is an analytic application, with exact executable input witnesses,
of the generic Lean theorem
`QFT.SourceOperatorGeneration.adjoin_eq_top_of_walk`. No additional Lean
instantiation is claimed compiled by this package. The source H, every
edge and every separating projector are supplied; maximality of the
observable algebra is not an assumption. The resulting star-algebra is
M28 on each reducing charge sector.

Algebraic generation alone is distinct from a physical instrument.
`spin_operational_position/PROOF.md`, section3, supplies the needed
finite compiler: every computable finite instrument on the source spins
has a complete native control word with an explicit all-input/reference
error bound. It uses the declared addressed coherent controls, blank
records and pointer commitments of `native_spin_observer`. It retains
all source-export failures. The finite source ground algorithm and the
strict Gram bound (5) make the matrices here computable. Thus this
compiler applies without an added state-preparation or measurement oracle.
It does not provide efficient resources, bare-repair access or a continuum
scattering detector.

`Dynamics.ProtectedCharge` supplies the dual-fixed-observable criterion
for charge expectation conservation. The actual source commutators and
complete branch preservation discharge its conservation input here.
Its symmetry-only counterexample prohibits replacing those checks by a
claim of covariance alone. The specialized rich-history state in
`QFT.SourceHistoryGNSDynamics` is not substituted for Ω. The CAR and
Weyl/contact scattering packages use different source actions and reference
states; their particle interpretation and wave operators are not
transported to (1).

## 4. A complete pair-response detector

For one sign q use the computable effects

    E_j=(1/28)|B_q e_j><B_q e_j|,  j=1,...,28,
    E_dark=I-sum_j E_j.                                     (8)

Every column has norm<=1, so sum E_j<=I. Equation (5) also gives
sum E_j>=gamma P_(Q=2q)/28. Thus there is no undetectable vector in that
sector. Choose the finite Lüders instrument sqrt(E_j),sqrt(E_dark);
this specifies its backaction. The29 effects include the complement on
the whole256-dimensional source, not only a selected charge block.

For incoming pair i, retain all four outcomes in (2), evolve H and apply
(8). The original joint probability for the q,q branch and detector j is

    Pr(q,q,j;t)=|K_q(t)_(j,i)|²/28.                          (9)

The same-branch dark probability is w_qq-sum_j Pr(q,q,j;t).
The other three preparation branches have zero click probability and
retain their original weights in dark. There are116 coarse
preparation/detector records, with additional native first-failure labels
in a compiled implementation. The receipt encloses every one of the28
click probabilities and the same-branch dark outcome for the separated
incoming corners, at source times0 and1/4, for both signs.
A detector label j refers to overlap with a source preparation history;
it is not asserted to be the positions of two physical particles.

For the independent numerical enclosure, each exact source bond is minus
a unit projector and the selection field has norm<=1. Hence
-13I<=H<=I and ||H+6I||<=7. Degree-N exact integer Taylor evaluation has
operator remainder (7|t|)^(N+1)/(N+1)!. Each kernel entry has in addition
at most2e_g input/output error. Modulus squaring uses outward rational
bounds. The constant energy shift affects phases, not (9).
The full pair detector and evolution tape is not executed here; the
native compiler establishes its finite implementation with prescribed
complete error, on the same controller contract.

## 5. A nonzero source interaction and finite propagation response

Apply local Z kicks at x,y and compare their ground energy costs. The
mixed inclusion-exclusion operator is

    D_xy=Z_x Z_y H Z_y Z_x-Z_x H Z_x-Z_y H Z_y+H.

All source terms except a shared bond cancel exactly. On an actual edge,

    D_xy=-(X_x X_y+Y_x Y_y),

and on a nonedge D_xy=0. The receipt verifies these full operator
identities, including the selected field. For the actual edge(0,1),

    <Ω,D_01 Ω> in [-1.105415773065,-1.105415773063].          (10)

This is a nonzero source interaction observable, without a fitted
contact interaction. The bound ||D_01||<=2 supplies ground-proposal error
at most4e_g.

A complete finite response measures

    chi_xy(t)=-i<Ω,[Z_x(t),Z_y]Ω>
             =2 Im<Ω,exp(itH)Z_x exp(-itH)Z_y Ω>.           (11)

A native controlled word for the unitary inside (11), with a blank
interference pointer, has Y-plus probability1/2+chi_xy(t)/4.
All other pointer labels and failures are retained. Negative-angle source
controls implement the inverse in the word; passive time reversal is not
assumed. The compiler from section3 applies to this finite experiment.
Exact integer polynomial vectors enclose both unitaries. If the
single-unitary remainder is r, their scalar-product error is at most
2r+r²+2e_g, before multiplying by2 in chi.

At the same time t=1/8, the complete source gives

| Sites | Graph distance | Retarded response interval |
| --- | --- | --- |
| 0,1 | 1 | [0.137040588801,0.137040588803] |
| 0,3 | 2 | [-0.000948863729,-0.000948863727] |
| 0,7 | 3 | [0.000756454563,0.000756454565] |

For distance d, the exact nested commutators of orders0 through d-1
commute with the distant Z. This is the finite source support rule,
checked on the complete matrices. For an edge,
[[H,Z_x],Z_y]=-D_xy, so (10) fixes the leading connected-response
coefficient. These responses prove interaction and propagation in this
finite source. They do not constitute a continuum cross section or an
asymptotic scattering matrix.

## 6. Source locality and the asymptotic attachment

The existing ordered-spin construction supplies the volume-uniform bound

    ||[alpha_t(A_X),B_Y]||
       <=||A|| ||B|| |X| 2^-distance(X,Y)(exp(48J|t|)-1).    (12)

The selected onsite field can be put in the interaction picture, leaving
bond norms and supports unchanged. For a complete detector effect E,
a remote unitary kick Z_y changes its expectation by at most
||[alpha_t(E),Z_y]|| in every source state. For two disjoint unitary kicks,
the connected inclusion-exclusion of the four complete detector laws is
bounded by twice the smaller of their two bounds from (12): the two
conjugations commute and ||Ad_Z-I||<=2. This is a state-independent
causal control. It does not assume the selected ground factors across
separated regions. For a two-branch local charge instrument, summing
K_r*[A,K_r] instead gives an additional factor at most2 before any
conditioning. Keeping original weights avoids division by a rare branch.

There is a fixed-spacing infinite-volume construction for complete
local experiments. The ordered state is selected by the established
box-to-infinity passage followed by B decreasing to zero, or a justified
diagonal giving that same selected state. A family at fixed B=1/4 is not
identified with the zero-field ordered state. For that finite-ground
family omega_Lambda and its selected weak-limit state omega, choose the same
finitely supported creation Kraus operators K_a and local detector effects
E_r in every sufficiently large box. Their original probabilities are

    p_Lambda(a,r;t)=omega_Lambda(K_a* alpha_t^Lambda(E_r) K_a).

The inherited spin dynamics and (12) give norm convergence of the local
Heisenberg operator, uniformly on each compact time interval. Removing
the onsite selection field uses the same ordered passage; it is not an
exchange of the thermodynamic and large-time limits. Weak state
convergence then gives the actual source limit

    p(a,r;t)=omega(K_a* alpha_t(E_r) K_a).                  (13)

Positivity and sum_(a,r)p=1 follow from the complete instruments. Finite
sequences follow by composing their bounded completely positive maps;
norm continuity and telescoping preserve the same argument. Thus local
two-charge creation and complete finite-time collision readout extend to
the infinite ordered spin source. This is a consequence of the existing
fixed-spacing state and dynamics theorem, without an assumed particle
vacuum. It gives no effective selected-state/box convergence modulus.
The eight-cube numerical values need not equal their infinite-source
counterparts. In particular, the global ground-dependent rank-one frame
detectors (8) do not acquire a local infinite-volume limit from (12).

`spin_canonical_field` separately supplies the long-distance common
bounded field-history limits and a continuous-time local field net.
Only the controls and readouts actually included in that history
construction pass to its field scaling. Equations (10)-(13) do not prove
a nonzero connected field-scaled response, a transport of the finite
cube's ground to the field vacuum, or a scattering cross section.

For genuine scattering, a source-derived asymptotic Hilbert space and
identification J must be attached to the same selected-state generator L.
Given a self-adjoint comparison L0 and a common differentiable domain,
let D=LJ-JL0. The elementary identity

    d/dt [exp(itL) J exp(-itL0)psi]
       =i exp(itL) D exp(-itL0)psi

proves a quantitative wave-operator tail whenever
integral_T^infinity ||D exp(-itL0)psi|| dt is finite. For example an actual
source bound <=C(1+|t|)^(-1-delta), delta>0, gives tail
C/[delta(1+T)^delta]. Both incoming and outgoing tails are required to
compare a finite echo with a scattering operator. An isometric
identification and suitable range/completeness statement are separate
premises. The exact finite algebra in section3 discharges neither
asymptotic independence nor this integrability.

The supplied Weyl/contact action has such an attachment in
`scattering_packet_readout/PROOF.md`; its Hamiltonian, packet and reference
are explicitly different from (1). Its Cook machinery is reusable once
this source supplies L0,J and the estimates; its S matrix is not evidence
for those premises. `QFT.FiniteUnitaryScatteringNoGo` excludes the ordinary
limit of a single nontrivial finite unitary power and explicitly permits
relative comparisons. It is not a blanket no-go for scattering.
The finite construction above supplies the source operators, complete
pair response and detector needed to test a proposed asymptotic
attachment, with physical masses and cross sections unassigned.
