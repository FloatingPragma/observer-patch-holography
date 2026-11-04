# Local spin dynamics with a fixed positive native repair rate

## Source, state and clock

Use the complete two-flavor signed cubic source, in the onsite spin frame,
on periodic cubic tori whose side lengths are multiples of four and at
least four. All occupation sectors are present. The hopping signs are
`(-1)^(sum_(j>axis) x_j)`, so translations by two sites preserve the
source Hamiltonian. Write

```
G = nu N+T,   nu=U/2,   N=sum_x Q_x,   Q_x=(n_x-1)^2,
J=c^2/U>0,    U=c^2/J,   B=144c,       nu>=2^20 B.       (1)
```

Every edge has its own independent rate-Gamma Poisson clock, with fixed
`Gamma>0`. Each event applies the unchanged full-domain C repair or the
unchanged full-domain E repair; one law is fixed throughout the process.
The complete branch action is
[`full_domain.branch_action`](../defect_field_bridge/full_domain.py).
There are two coin values and sixteen raw slots, including collision
holds and zero Kraus branches. No branch is discarded or normalized
separately. The quantum coin is prepared in `|+>` in the source dilation.
Its classical readout is explicitly the dephasing in the coin/raw basis.

The initial source marginal is invariant under translations by two sites,
and has positive energy density bounded by a constant `e_plus`, independent
of c and volume. An untouched reference system can have arbitrary
correlations with the source. Symmetry is required of the unconditional
source marginal, not of each conditional trajectory. There is no
record-dependent feedback in the theorem. Uniform edge rates, the uniform
law and the symmetric unread coin preserve the stated symmetry. In
particular, with `q(t)=Tr(rho_t N)/V`, each site's defect probability is
at most `8q(t)`.

The physical comparison uses the source time of G, with fixed exchange J.
Only the lattice spacing is held fixed. This domain does not assert a
continuum field equation, a relativistic shell or a particle assignment.

## Exact conservation and the energy estimate

Every native repair Kraus operator commutes with N. For binary endpoint
occupancies the rule either holds or exchanges occupancies 01 and 10.
For a nonbinary endpoint it holds and records raw slot 8. Thus the number
of empty and doubly occupied sites is conserved by each branch, including
its fermionic sign. On `P_e=product_(x in e)(1-Q_x)`, every branch is
`delta_(raw,0) P_e/sqrt(2)`. The unread channel is unital. These statements
follow for arbitrary volume from the local branch rule; `algebra.py`
checks the complete finite matrices, including nonadjacent seams.

For a Hermitian local A with `P_Z A P_Z=0`, put
`p=Tr(rho(1-P_Z))`. Separating the P/Q blocks and applying the state
Cauchy-Schwarz inequality gives

```
|Tr(rho A)| <= ||A|| [p+2sqrt(p(1-p))] <=3||A||sqrt(p).  (2)
```

Each hopping bond has norm at most c and zero protected compression.
Its support has at most two sites and defect probability at most 16q.
There are three bonds per site, so

```
|Tr(rho T)|/V <=36c sqrt(q).                            (3)
```

Let `e(t)=Tr(rho_t G)/V`. Hamiltonian evolution conserves G, and repair
conserves N. A repair edge meets at most eleven hopping bonds. Each
increment `Phi_e*(T_f)-T_f` has norm at most `2c`, support at most three
sites and zero protected compression. Disjoint terms vanish because the
maps and interactions are even local operators. By (2), its expectation
has absolute value at most `6c sqrt(24q)`. Therefore

```
|e'(t)| <=3 Gamma*11*6c sqrt(24q(t))
        <=990 Gamma c sqrt(q(t)) <=1000 Gamma c sqrt(q(t)). (4)
```

Let `M_T=sup_(0<=t<=T) sqrt(q(t))`. Combining (3)--(4) with positivity
of `nu N` proves

```
nu M_T^2 <=e_plus+(36+1000 Gamma T)c M_T,
M_T <=m_T:=(36+1000 Gamma T)c/nu+sqrt(e_plus/nu).        (5)
```

For fixed J, Gamma and T, this is `q(t)=O(c^-2)` uniformly in volume.
Bare one-particle-per-site spin inputs have e_plus=0. Period-two averaged
actual source ground states also have e_plus=0: every bare protected
vector is a zero-energy variational trial. The full Gibbs state at fixed
inverse temperature beta has `e_plus=log(4)/beta`, by its variational
free energy and `S<=V log(4)`. No source spectral gap is assumed.

## Local rotation and the unital spin algebra

Use the exact five transformations Y of
[`LOCAL_LIMIT.md`](../bipartite_spin_limit/LOCAL_LIMIT.md). Its norm is
`||Phi||_k=sup_x sum_(Z containing x) exp(k|Z|)||Phi_Z||`, with connected
supports. The bounds, which also hold on these tori, are

```
Y G Y* =nu N+D5+V5,
||D5||_(1/2)<=2B,  ||V5||_(1/2)<=r:=B(4096B/nu)^5,
sum_(j=0)^4 ||Sj||_(7/12)<=2B/[nu(1-rho)], rho=4096B/nu,
||Y* A Y-A||<=8|supp A| ||A|| B/nu                    (6)
```

for a strictly local A. The square-root projection estimate is
essential here. For `rho'=Y rho Y*`, the Hilbert-Schmidt triangle
inequality applied to `Q_x Y sqrt(rho)` gives

```
sqrt(Tr(rho' Q_x)) <=sqrt(Tr(rho Q_x))+||[Q_x,Y]||,
Tr(rho'_t Q_x) <=p_T:=(sqrt(8)m_T+8B/nu)^2.            (7)
```

This is `O(c^-2)`. Bounding only the difference of the two projection
expectations would give an insufficient linear dressing error.

At each site identify the four-state carrier with charge times spin:
`|0,0>=|up>`, `|0,1>=|down>`, `|1,0>=|empty>`,
`|1,1>=|up down>`. Set `j(A)=I_charge tensor A`. This is an even unital
star homomorphism, commuting with every Q_x. On the charged block it can
change particle number by two; it is a comparison readout algebra without
a particle-number superselection restriction. The positive spin state
`omega_c(A)=Tr(rho'_0 j(A))` is an ordinary partial trace over charge.
No projection or conditioning of the input state defines this spin state.

The local-series proof underlying (6) gives, before compression,

```
D5=D0+K2+E,  D0=A(T),  K2=A([S0,V0])/2,
||D0||_(1/2)<=B/3,  ||K2||_(1/2)<=k2:=48B^2/nu,
||E||_(1/2)<=eta:=2^25 B^3/nu^2.                       (8)
```

Here A is averaging under the onsite integer penalty and `V0=T-A(T)`.
D0 has edge supports, commutes with their penalty and annihilates their
protected space. K2 has supports of at most three sites. The norm bound
for K2 follows from the interaction commutator bound with loss 1/12,
`||S0||<=2B/nu`, and averaging. Its protected compression is exactly
`Hspin=J sum_edges(S_x.S_y-1/4)`. In fact the first omitted step is
bounded by `2(288B/nu)^2 B`, and all later diagonal increments by
`B(4096B/nu)^2`; their sum is at most eta.

For each K2 support Z, subtract the unital extension of its compressed
spin term. The difference C_Z has norm at most twice that of K2_Z,
commutes with the local penalty and vanishes on P_Z. Thus
`C_Z=(1-P_Z)C_Z(1-P_Z)`. The same reduction holds for D0_Z. Since the
spin evolution and its algebra commute with every Q_x, the commutators
of these terms with every evolved spin readout retain this reduction.
Their expectations have a factor `|Z|p_T`, rather than its square root.

## Completely bounded repair maps

All local map norms in this section are completely bounded operator
norms. This retains arbitrary untouched reference systems. Let
`L_R=Gamma sum_e(Phi_e*-id)` be the original unread repair generator.
Its local norm at decay 7/12 is at most `48 Gamma`: six incident edges,
map norm at most two, and `exp(7/6)<4`. Each term annihilates the identity.
Local maps with disjoint source supports commute.

The connected-support commutator estimate also applies to these maps:
`||[M,L]||_(k')<=4||M||_k||L||_k/(k-k')`. The derivation of S_j has
map norm at most `2||S_j||`. Allocate decay loss `1/60` to each of the
five conjugations. The convergent Baker-Campbell-Hausdorff series has
ratios `q_j<=1440||S_j||_(7/12)`, and

```
sum_j q_j<=4320B/nu<1/2,
||Ad_Y L_R Ad_(Y*)-L_R||_(1/2)<=2^20 Gamma B/nu.        (9)
```

Indeed the product of the five majorants is at most
`1/(1-sum q_j)`, and its excess is at most `2sum q_j`.
Multiplying by `48 Gamma` proves the displayed constant.

For a spin readout, `(Phi_e*-id)(A)` commutes with N_e and vanishes
on P_e. Its expectation is bounded by `2p_T` times its map norm.
Thus the original repair contributes an effective map strength at most
`96 Gamma p_T` in the local comparison below. This uses (7) on the
exact evolved state, without assuming that this state is protected.

## Spin propagation and the quantum estimate

For a local identity-annihilating map M_Z of completely bounded norm
m_Z, Haar averaging over the local algebra on Z bounds
`||M_Z tau_spin^t(A)||` by m_Z times the distance of the evolved A from
the commutant of Z. The commutator-chain estimate in LOCAL_LIMIT.md
therefore gives the bound

```
2||A|| m_Z sum_(x in supp A,y in Z)
       exp(8d|t|-dist(x,y)/4),   d=48J,
```

together with the trivial cap `m_Z||A||`. Choosing
`L=ceil(64d|t|)` and using the cubic sphere count and exponential tail
from that proof bounds the sum over local maps by
`2^16 s||A||||M||_(1/2)(1+64dT)^3`, for `s=|supp A|` and `|t|<=T`.
The estimate holds for spin readouts tensored with an arbitrary reference
or record algebra: the averaging and norms are completely bounded.

Apply Duhamel's formula to the expectation of
`tau_spin^(t-u)(A)` in the exact dressed state at time u. The penalty
commutator is zero. For D0 and the charged complement of K2, multiply
the local commutator bound by `2p_T` and `3p_T`, respectively, using
the reduction following (8). Hamiltonian derivations have map norm
twice the interaction norm. A factor four in the following prefactor
dominates all these factors and the original unread repair contribution:

```
Delta_T=(2B/3+6k2+96 Gamma)p_T+eta+r+2^20 Gamma B/nu,
E_quantum(T)=2^18 s||A|| T(1+3072JT)^3 Delta_T.         (10)
```

This bounds the local expectation error between exact continuing native
repair and spin evolution on the initial positive pullback omega_c.
The comparison uses spin propagation with velocity controlled by J;
the large charged hopping coefficient c enters only the state-dependent
factor `Bp_T=O(c^-1)`. Equation (10) is `O(c^-1)` at fixed J, Gamma,
T and e_plus, independently of the total volume.

## Complete local records and time stamps

Fix a finite observed edge set F. Keep the complete coin/raw output and
time stamp of every event on these edges. Unobserved edges continue to
repair at rate Gamma. In a coherent tape, the ideal appended output is
`|+><+|_coin tensor |0><0|_raw`. For the classical record algebra both
actual and ideal records are dephased, giving the fair coin and raw zero.
These conventions agree with the exact protected source dilation.

Let I_e denote the actual recorded append, E_e the unread channel with
identity on the existing archive, and I0_e the ideal independent append
with identity source channel. At the event time, each append includes
that time stamp. The actual generator on the source and observed archive
is `L_G+L_R+Gamma sum_(e in F)(I_e-E_e)`. The ideal comparison has
generator `L_spin+Gamma sum_(e in F)(I0_e-id)`. This identity retains the
unread marginal on every edge and all observed events.

For an embedded spin times archive readout O put
`C_e(O)=(I_e-E_e-I0_e+id)(O)`. Every branch commutes with N_e, O
commutes with each Q_x, and the complete protected instrument is the
ideal append. Consequently C_e(O) commutes with N_e and has zero P_e
sandwich. Since P_e is the zero eigenspace of the positive N_e,

```
C_e(O)=(1-P_e)C_e(O)(1-P_e),
||C_e||_cb<=4,  |Tr(rho'_t C_e(O))|<=8p_T||O||.        (11)
```

This holds in the dual norm on the local source algebra on X tensored
with the full finite-cylinder archive algebra and an independent untouched
reference. Arbitrary bounded joint readouts include correlations with
the final spin readout, the complete archive and the reference. Only the
source marginal enters p_T. It does not infer conditional translation
symmetry or a trace norm for the entire source.

Conjugating a local completely contractive map by `exp(S_j)` changes
it by at most `8||S_j||_0` when its source support has two sites. This
follows by integrating its commutator with the derivation; only supports
meeting the edge contribute. Conjugations are completely isometric, so
the five changes telescope. Equation (6) bounds their sum by `32B/nu`
for each of I_e and E_e. The ideal append commutes with Y. Thus the
additional joint record error is

```
E_record(T)<=Gamma |F| T ||O|| (8p_T+64B/nu).           (12)
```

The ideal archive evolution is contractive, commutes with spin evolution,
and preserves its embedded algebra. Equations (10)--(12) therefore apply
to its backward-evolved readouts. Time-dependent append maps are covered
by the same Duhamel integral, or by the predictable Poisson compensator.
Finite event-count truncations converge by the Poisson tail, so the
estimate includes all record lengths on [0,T]. An infinite lattice archive
is described by consistent finite-region/time cylinder laws; no global
trace norm for an infinite archive is asserted.

## Bare readouts, faithful inputs and growing time

For a unit-norm local spin readout add `8sB/nu` to (10)+(12) to remove
the final local dressing. For a bare protected input, comparison of the
initial pullback with its original spin state costs at most

```
2^19 s(B/nu)(1+3072JT)^3.                              (13)
```

To prove (13), telescope the five isometric conjugations of the backward
spin readout. Each difference is the integral of its commutator with S_j.
Apply the spin propagation sum above and
`||ad(S_j)||<=2||S_j||_0` and `sum||S_j||_0<=3B/nu`. The constant
in (13) dominates the resulting `6*2^16`. The ideal archive factor is
a contraction and has no source
support of its own, so the same bound covers joint record readouts.

For any period-two spin state sigma, the full-source input

```
rho_0=(1-delta) sigma_P+delta I/4^V,
delta=(c/nu)^2,    e(0)=delta nu/2=J                   (14)
```

is faithful. Every occupation state has positive weight. Both Hamiltonian
and unread native repair are unital, so
`rho_t>=delta I/4^V` for every finite volume and time. For a correlated
reference replace the mixed branch by `I/4^V tensor sigma_reference`.
The source marginal remains faithful; the reference need not be faithful.
Compared with the spin evolution of sigma, add `2delta` to (13).
The coefficient delta is an input mixture weight, not a retained-outcome
probability. No input or repair event is postselected.

The choice of tracial sigma gives a faithful local stationary reference
for the interacting spin limit, with continuing independent Poisson/coin
records. This does not make the finite-c input (14) stationary or identify
its modular flow with source time. General full-source Gibbs and ground
inputs obey the energy-domain estimate (5) on their exact initial states.

At fixed J, Gamma and e_plus, let `T=O(c^alpha)` with `0<=alpha<1/6`.
Then `p_T=O(c^(2alpha-2))`, (10) is `O(c^(6alpha-1))`, (12) is
`O(c^(alpha-1))`, and (13) is `O(c^(3alpha-1))`. Every error tends to
zero. When alpha is positive the expected number of records on every
observed edge diverges. The rational example alpha=1/12 in `bounds.py`
certifies this expanding-time regime without executing an exponentially
large source.

All constants are uniform in torus volume. The source and spin generators
at fixed c are local and have thermodynamic dynamics; the connected
interaction bounds also define Y on local observables in this limit.
Taking local weak limits of the stated symmetric input families therefore
gives the same comparison on the infinite fixed-spacing lattice.

## Common-clock boundary

The source penalty and spin exchange have ratio
`nu/J=c^2/(2J^2)`, which diverges on (1) and is unchanged by a common
positive clock rescaling. A clock keeping the penalty finite sends the
exchange to zero; a clock keeping the exchange nonzero leaves the penalty
divergent. This is a ratio of coefficients in the source generator,
not an identification of its exact many-body charged gap with nu.
The theorem does not supply simultaneous finite massive charged and
massless sectors, a sharp Lorentz shell, a continuum interaction or
Einstein dynamics.
