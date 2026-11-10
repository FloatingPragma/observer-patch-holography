# Exact equal work with an elastic law on the selected source

The actual primary source admits a complete parity apparatus whose mean
work equals the critical phase apparatus at every capacity. Its original
transfer measure converges to the elastic atom while its work has a
strictly positive lower bound and its second moment grows as inverse
spacing. The construction changes the apparatus, not the source, state,
reference, Hamiltonian or clock. It does not prove collapse of the
critical phase apparatus.

## 1. A physical local parity and its exact work

Use the [primary contract](../particle_common_contract/PROOF.md), with
`a=2^-n`, `R=1/a`, `L=2R^2` and `V=L^3`. Keep its full sixty-mode matter
carrier, quantum links, physical Gauss constraints, product cutoff,
selected polynomial state `rho`, idle canonical reference and H.
The source hopping is `H_hop=sum_b T_b`, with `||T_b||<=30/a`.
At cell x set `P_x=(-1)^N_x`. Literal CAR gives

    P_x c_(y,j) P_x=(-1)^(delta_xy)c_(y,j),
    P_x H P_x-H=-2 sum_(b incident to x) T_b.                (1)

Every other source term commutes with P_x. This unitary preserves all
native occupations, total parity, every Gauss constraint and every link
cutoff; no charged external memory or alternative preparation is needed.
Each physical open bond has two endpoints. Averaging all retained origins
therefore gives the exact selected work of a parity kick,

    W_parity=-(4/V) Tr(rho H_hop).                         (2)

No translation invariance of the finite-box state is used in (2).

For the original unit tent let kbar be the exact cosine sum in the
[critical-profile proof, equations4–5](../particle_primary_critical_refinement/PROOF.md#2-exact-selected-work-and-finite-spatial-errors).
Its selected work is

    W_critical=-(kbar/V) Tr(rho H_hop)=h_n kbar/a,
    h_n=-(a/V) Tr(rho H_hop) >= h_* >0.0424.                (3)

The expectation in (3) is the actual selected source marginal. Its
bound uses a variational inequality; the variational vector does not
replace that marginal. The directional expectation equality needed for
general rectangular profiles follows from the actual cube symmetry.

## 2. An exact source-defined coin and complete instrument

Choose the probability `p_n=kbar/4`, and use all2V Kraus operators

    K_(x,0)=sqrt((1-p_n)/V) I,
    K_(x,1)=sqrt(p_n/V) P_x.                              (4)

Their effects sum to I. The origin, inactive/active coin, source output,
idle reference and complete energy records are all retained. The
finite coin/isometry services are declared apparatus operations; no
large source histogram is asserted to have been executed. Their unitary
extension lies within each finite conserved block. The inactive outcome
is not discarded. In particular, this is an apparatus intervention and
not a change of the native cooling or repair law.

For the unit tent the three directional sums coincide. With
`m_R=1`, `m_u=2` for `1<=u<R`, the exact finite coin is

    p_n=(R/2) sum_(u,v=1..R) m_u m_v
                        [1-cos(pi*u*v/R^4)].              (5)

The2R edges of a one-dimensional tent have slope magnitude1/R.
The transverse tent values u/R have multiplicities m_u, proving (5),
including the wrapped coefficient list at R=1. The source bonds remain
open. At n=0, (5) gives exactly p_0=1. At n>=1, the exact tent sums and
`pi^2<10` give

    0<p_n <= (5/12)a^2 D2 <=(45/32)a<=45/64<1.             (6)

Here `a D2=2(2+a^2)^2/3<=27/8` for a<=1/2. Thus (4) is admitted at
every selected member. It uses algebraic roots of unity of order2R^4,
their real parts and positive square roots. No selected-state expectation,
mass datum or fitted parameter enters this coin. The code certifies
compact rational enclosures using Machin's identity for pi and finite
alternating cosine sums; those enclosures do not redefine the exact coin.

Combining (2)–(4) proves the original, unconditioned identity

                      W_control=W_critical.              (7)

The equality holds at every capacity, with the actual unknown h_n
cancelling in the coin. More generally a fixed admitted rectangular
profile and rational amplitude use p=kbar(b,s)/4. For any fixed finite
set of such settings, p=O(a) and hence p<=1 eventually. Their entire
spatial work form and coefficient-cancelling work ratios are reproduced
by this distinct apparatus. This is not a spatial-response or dynamics
equivalence between the apparatuses.

## 3. Original spectral collapse with surviving mean work

Initial exact energy recording preserves the source marginal because
rho commutes with H; it can change the source/reference correlations.
Let nu_n be the complete original transfer law of the origin-averaged
parity kick. The complete law of (4) is exactly

    mu_n=(1-p_n) delta_0+p_n nu_n,
    mu_n(omega !=0)<=p_n ->0.                            (8)

The total spectral mass remains one. Thus its transfer laws converge
in total variation to delta_0, despite the positive work (7).
The [selected work bound](../particle_primary_critical_refinement/PROOF.md#1-actual-state-clock-and-complete-profile-apparatus)
gives `W_control>=w_*>0`, with `w_*=(16/9)h_*` sufficient. Since all
nonzero transfer lies in the active coin, Cauchy–Schwarz gives

    M2 >= W_control^2/p_n.

The degree z is3 at n=0 and at most6 otherwise. Equation(1) implies
`||P_x H P_x-H||<=60z/a`. Source stationarity identifies its squared
expectation with the second transfer moment. Therefore

    W_control^2/p_n <= M2 <=3600 z^2 p_n/a^2,
    (32 w_*^2/45)/a <= M2 <=182250/a.                    (9)

The uniform lower bound also holds at n=0, since p_0=1 and32/45<1.
This is an actual selected-source example with `M2=Theta(1/a)`,
the order of the [critical apparatus bound](../particle_primary_weighted_spectral_tail/PROOF.md).
It is stronger than a separate-carrier counterexample: it uses the exact
same selected state and matches the exact work, without its evaluation.

For every fixed energy cutoff K>0, positive transfer below K contributes
at most Kp_n. Negative transfer only strengthens the following inequality:

    integral_(omega>K) omega dmu_n >=W_control-Kp_n,
    liminf_n integral_(omega>K) omega dmu_n >=w_*.          (10)

The control consequently fails uniform energy-weighted integrability.
Every positive fixed-energy band loses its original probability. No
conditioning on the active coin can remove that cost in (8)–(10).

## 4. The original reference and unchanged complete clocks

Compare (4) with a matching idle experiment using the same origin and
coin probabilities, but I in both branches. For an arbitrary input
including its quantum reference, their recorded classical-quantum output
distance is at most p_n in half trace norm: the inactive blocks coincide,
and each active block differs by at most its own probability. Subsequent
complete instruments contract that distance. A fixed history of k such
recorded interventions has distance at most kp_n by telescoping.
If the coin and analyzer are kept coherently unread, the difference of
the two isometries instead has norm at most2sqrt(p_n); the corresponding
reference-inclusive half trace distance tends to zero as well. These two
bounds retain their respective recorded and coherent domains.

Use precisely the critical apparatus's authenticated finite energy-clock
plan on the same H, with resolution eta_n and error allowances. Its
uniform [complete-clock estimates](../particle_primary_joint_spectroscopy/PROOF.md)
apply to each unitary in (4), and convexity preserves the budgets.
The signed mean error e_n is the same certified O(a^2) bound, so

    recorded_mean >=W_critical-e_n>0.                     (11)

For any fixed positive band with lower endpoint greater than2eta_n,
the idle branch enters only through the two-clock tail/implementation
allowance. Its reported probability is at most p_n plus that vanishing
allowance. All initial/final clock labels and zero-weight outcomes remain.
The receipts bind the identical clock digest and original2V clock-product
record count. They neither subtract an uncontrolled absolute error from
a rare event nor replace the original state by a conditional one.

## 5. The specific unresolved critical-probe calculation

Equations(7)–(11) show that the actual source, cooling law, complete
records, positive Dirichlet work form and an O(1/a) second-moment cap do
not determine nonelastic survival. Even exact knowledge of that work at
all fixed profile settings cannot choose between their spectra. This
does not exclude survival of the original critical phase experiment:
its many small coherent number rotations are not the rare parity kicks.

Its missing data are the selected current correlations, for example
`average_X Tr rho(C_X^2+S_X^2)` and the current double commutators in the
[selected spectral-tail proof, §5](../particle_primary_weighted_spectral_tail/PROOF.md#5-complete-original-finite-windows-and-the-remaining-limit-question).
The coin control evaluates neither quantity for that original apparatus.
A uniform second moment there would imply positive finite-energy weight
from the same positive-work argument, while its joint spatial attachment,
dispersion, charged or neutral threshold and physical mass identification
would require their own proofs. No mass or universal source-gap theorem
is asserted here.
