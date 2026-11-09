# A prepared finite field and its complete occupation records

This construction composes the existing gauge preparation, fixed-spacing
locality and representation-cutoff theorems on one declared experiment.
It specifies a finite source, its prepared mixed state, one Hamiltonian
tick and a complete detector. Exact rational enclosures give a nonzero
occupation-transfer event without allocating its enormous Hilbert space.
The result concerns the chosen finite polynomial action. Its absolute
comparison to a periodic exact-sign parent does not resolve that rare
event in the parent, and supplies no interacting relativistic limit or
physical particle identification.

## 1. One periodic preparation and one action

Use the source representation `M15` of `G=S(U(3)×U(2))`, all sixty
Dirac/fermion modes at every cubic cell, and coordinate links `L²(G,Haar)`.
The [owning action](../particle_observer_outputs/FIELD_ROUTE.md) is

```
H_m = H_E + dGamma(h_m) - kappa sum_p Re Tr_M15(U_p)/15.
```

The electric coefficient is one half in the actual M15 trace metric.
The example fixes spacing `a=1`, polynomial order `m=0`, interior box
radius `R=5` and periodic side `s=21`. All sixty modes at the origin are
filled; all other fermion modes are empty. The determinant character of
the complete local representation is trivial, so this occupation state
is gauge invariant. The conserved total fermion number is sixty. The
Hilbert space retains every occupation, including sectors outside this
chosen input. There is no scalar pairing term in this experiment.

On the periodic gauge carrier prepare the normalized polynomial

```
Psi_L = P_Gauge product_edges [1-||R(U_e)-I||_F²/60]^L / norm,
L=n².
```

The [finite-word preparation](../particle_field_preparation/PROOF.md)
gives an explicit finite code and a deterministic unitary preparation.
Its [conserved-energy refinement](../particle_field_preparation/REFINEMENT.md)
chooses `L` and the fixed coefficient `kappa` for target `eta=2^(-20)`.
With `E=3s³`, initial electric bound `E0=(29/10) E L²`, kinetic bound
`B=60(2+epsilon_0)`, and `epsilon_0=240/17`, it chooses

```
kappa = 96000 (E0+2B)/eta,
<S>_0 <= eta/96000,
S = sum_p [1-Re Tr_M15(U_p)/15].
```

The represented-plaquette admissibility theorem gives
`1_{sigma_min(H_W)<1/2} <= 48000 S`. Conservation of this same Hamiltonian
therefore bounds the bad-gap probability by eta at every time. This is
a particular high-plaquette-cost architecture, with its coefficient
displayed in the certificate. All locality and cutoff estimates below
use that coefficient. The coefficient `1/100` of another control is
not substituted into any estimate.

Preparation and detection are declared controlled services with their
work accounts. They are separate from autonomous observation evolution.
A narrower repair implementation is governed by the
[finite readout contract](../particle_observer_outputs/READOUT_CONTRACT.md).
No numerical particle data chooses this action, state or clock.

## 2. The interior state, boundary charges and two reference roles

Attach each positive-direction link to its tail cell, as in the
[same-field locality proof](../particle_field_causality/PROOF.md).
Let `B_R` be the entire radius-R cube and retain all its outgoing link
registers. The open input is exactly

```
rho_B = Tr_periodic exterior |Psi_L,filled_origin><Psi_L,filled_origin|.
```

This is the whole-box marginal. Replacing it by independent Haar links,
an unrelated vacuum or only the origin's one-cell marginal would not
satisfy the comparison theorem. Links on deleted outgoing boundary edges
remain electric spectators. The smaller Wilson operator deletes hopping
through its boundary, so the onsite electric terms agree under the
tensor embedding.

Gauge averaging acts on each link by left/right translations and
preserves its matrix-word degree. The prepared ket lies in a tensor
product of per-link spaces of degree at most L. If Q projects onto that
local tensor product, `(Q tensor I)|Psi_L>=|Psi_L>`, hence
`Q rho_B Q=rho_B`. Tracing exterior factors cannot increase ket support.
The bra degree in a density matrix is not an additional ket degree.
With `E_B=3(2R+1)³`, total degree

```
K0 = E_B L                                                     (1)
```

is sufficient. The gauge support dimension is at most `451^(E_B L)`.
The Fock part is the specified pure occupied-origin state.

The marginal commutes with the induced gauge transformations, but need
not be supported only on box singlets. Boundary electric charges are
correlated with the exterior. A finite simulator prepares the same
mixed state with its canonical purification

```
|sqrt(rho_B)> = sum_ij (sqrt(rho_B))_ij |i> |j*>.
```

The reference transforms in the conjugate representation. Gauge
invariance of rho_B makes the joint vector invariant. Its coefficients
are computable from the finite source polynomial, exact compact-group
Haar inner products, the positive Gram quotient and a finite matrix
square root. A unitary taking a chosen invariant blank to this invariant
vector extends by identity on its complement. For example, after aligning
the overlap phase, a Householder reflection implements that map. It
retains all input columns and has no accepted/rejected preparation
branch. Zero Gram directions are coordinate redundancy, not discarded
physical outcomes.

There are two distinct uses of a reference. The periodic/open comparison
can retain any external system that is untouched by both experiments.
The original periodic exterior evolves and is traced out; its whole
output is not controlled by the local comparison. The simulator's
conjugate purifier is retained through its finite experiment, but need
not be the original periodic exterior. The comparison with the periodic
source concerns the local detector and a genuinely untouched external
reference, tracing this simulator purifier when the output spaces are
identified. No claim of global-state closeness is made.

## 3. Periodic locality, open dynamics and representation cutoff

The order of comparison is

```
periodic exact sign -> periodic polynomial
                   -> interior open polynomial -> finite word code.
```

The periodic exact-sign comparison uses the same prepared state, electric
term, plaquette coefficient and Hamiltonian duration T. The owning bound
is

```
e_sign <= min(1,60 T [epsilon_0+(2+epsilon_0)2^(-10)]).          (2)
```

The bad-gap estimate is applied to the periodic parent. An open Wilson
operator can have boundary modes; no periodic free-gap theorem is
applied to that operator.

For the polynomial, the [locality theorem](../particle_field_causality/PROOF.md#5-volume-control-for-complete-local-histories)
uses `d=2m+1`, `b=max(d,2)`, boundary separation `D=R-2b`, rooted kinetic
norm J_root and weighted interaction rate L_rate. Its boundary estimate
uses

```
C_R = 2 J_root [2(R+d)+1]³ + 3 kappa [2(R+2)+1]³,
e_box <= min(1,2 C_R T 2^(-D) [exp(2 L_rate T)-1]).             (3)
```

At m=0, `J_root=360/17` and `L_rate=2880/17+192 kappa`.
To apply the open-box proof within a periodic graph, cut all interaction
terms joining the box to its exterior and retain the exterior Hamiltonian.
Only roots in the same bounded boundary collar differ. Terms whose paths
leave and return through the boundary are included in C_R. Take the
conservative periodic side

```
s >= 2(R+2b+d)+1.                                             (4)
```

The entire relevant collar has a unique nonwrapping coordinate lift.
Every differing support is at periodic distance at least D from the
origin, including supports crossing the chosen torus coordinate seam.
The commutator-chain proof therefore gives (3) with the same constants.
All terms farther away agree and cancel in the Duhamel comparison.
Its operator bound holds for the actual whole-box marginal and untouched
references, irrespective of entanglement with the traced exterior.

Let `V_K` be the source's total matrix-word space of degree at most K,
choose `K=K0+12`, and define the finite Hamiltonian

```
H_K = P_K H_open P_K                                         (5)
```

on the full Fock carrier times `V_K`. Total-degree compression is not
exactly local, as the owning
[counterexample](../particle_field_causality/PROOF.md#6-why-the-total-degree-cutoff-needs-its-own-error)
shows. The parent locality estimate is used before compression.

Electric terms preserve word degree. At m=0, kinetic multiplication
raises total degree by at most one, magnetic plaquettes by at most four,
and the final number pointer preserves it. On the invariant N=60 sector,
put

```
B_open = 60(2+epsilon_0)+kappa E_B,
M_K = (29/10) K²+B_open >= ||H_K|| on N=60.                    (6)
```

The electric estimate follows from
`sum_links degree_link² <= K²` and the native trace-metric Casimir bound.
The [complete-record cutoff theorem](../particle_observer_outputs/FIELD_ROUTE.md#a-controlled-representation-cutoff-for-complete-observer-histories)
then gives

```
e_cut <= min(1,2 sum_(j>=4) (B_open T)^j/j!).                  (7)
```

The first omitted order is four because `K-K0=12` and the maximum word
step is four. This theorem applies to the mixed marginal by purification.
It keeps the full purifier and every detector outcome. The action remains
defined on every Fock occupation; N=60 is an exact invariant property
of this prepared input, with no filtering or conditional normalization.

## 4. A fixed finite tick and sixty-one detector outcomes

Choose

```
T = min(1, 1/(128 M_K²), 1/(4 L_rate)).                        (8)
```

One complete macro tick is `exp(-i T H_K)`. If repeated, the same fixed
unitary service defines the operational tick count. Preparation and
the final pointer have their own service boundaries. No physical seconds
or GeV are assigned. Equation (8) is a conservative choice for the
explicit probability enclosure; the high-kappa gap theorem itself holds
for every model time. This small-tick construction does not establish
a uniform observation interval along a continuum refinement.

The final detector measures the origin's total occupation `n=0,...,60`.
Its projectors have local Fock ranks `binomial(60,n)`, tensor identity
on every other field and reference register. Their sum is identity.
A six-bit pointer implements the complete controlled addition; labels
61,62,63 have zero ideal probability. For every input, the unitary
preserves all pointer columns. Implementation failure records, when
present, are retained separately. The original complete law is

```
p_n = Tr[(P_n tensor I) exp(-i T H_K) rho_B exp(+i T H_K)],
sum_(n=0)^60 p_n=1.                                         (9)
```

The code certifies a nonzero output of (9). The actual m=0 spin hopping
blocks obey

```
B_mu,± = -(4/17) beta (I±gamma_mu)/2,
B_mu,±† B_mu,± = (16/289)(I±gamma_mu)/2.
```

There are fifteen internal components and a rank-two spin projector.
All sixty origin modes start filled and all six distinct neighboring
cells start empty. The CAR makes different destination/hole labels
orthogonal, while unitarity of the native link representation cancels
its coordinate dependence. Thus, for every normalized prepared link
state and retained reference,

```
||P_59 H_K psi||² = 6*15*2*(16/289)=2880/289.                 (10)
```

All other terms preserve origin occupation. A single hopping insertion
cannot put it below59. The initial code is contained in V_K and one
action insertion has degree at most K0+4, so (10) is unchanged by the
compression. No flat-gauge configuration replaces the prepared state.

Let `r=T² M_K²/2`. Unitary Taylor remainder gives
`||exp(-iTH_K)psi-psi+iTH_K psi||<=r`, and (8) gives `r<=T/2`.
Since `3<sqrt(2880/289)<4`, the exact complete probabilities satisfy

```
(3T-r)² <= p_59 <= (4T+r)²,
sum_(n=0)^58 p_n <= r²,
1-(4T+r)²-r² <= p_60 <= 1-(3T-r)².                           (11)
```

Every individual n below59 is exported with interval `[0,r²]`; the
shared tail bound and exact total one accompany those intervals. This
is an enclosure of the original distribution, with all outcomes present.
It is not a list of normalized selected branches. The finite source
event in (11) has a strictly positive lower bound.

## 5. One error budget and its resolution

For `x=2 L_rate T<=1/2`,
`exp(x)-1<=x/(1-x)<=2x`. The code uses this rational envelope in (3).
For (7), it reuses the owning exact tail routine with `x=B_open T<=1/128`:
two series terms are evaluated and the tail starting at order six is
bounded by a geometric series with ratio `x/7`. These compact enclosures
avoid constructing enormous exponential numerators. Resource counts
stay factored; neither the
prepared large vector nor its propagator is materialized.

The final local-record total-variation bound, including untouched
external references in its half-trace form, is

```
e_total <= min(1,e_sign+e_box+e_cut).                         (12)
```

The fixed example certifies `e_total<10^(-6)`. Each individual parent
probability lies within e_total of its finite counterpart. Its p59
interval includes zero: the absolute comparison precision is coarser
than the rare transfer event in (11). The finite event is proved, while
positivity of that particular parent event is unresolved by this bound.
At m=0, epsilon_0 is large, so the polynomial is not presented as a
precise relativistic dispersion approximation.

If preparation, one macro tick and one detector use have complete
half-diamond errors `e_p,e_u,e_d`, with input half-trace error `e_i`,
the [existing composition contract](../particle_observer_outputs/READOUT_CONTRACT.md)
adds `e_i+e_p+e_u+e_d` to (12), capped at one. It compares full channels
with failure and unused labels retained. A supplied primitive has its
declared exact action; a compiled implementation needs the stated
accuracy and resource proof. No enormous compiled tape is asserted here.

The construction establishes a coherent finite experiment with an
explicit state, one action and clock, complete records and a composed
comparison budget. A physical sector requires correlation and
energy-momentum control along an appropriate interacting refinement,
including chiral dynamics and a physical detector interpretation.
The all-time gap estimate supplies one such premise at fixed parameters;
it does not supply those further premises by itself.
