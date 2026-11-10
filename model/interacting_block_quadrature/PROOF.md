# Relativistic neutral records from collective interacting blocks

## Declared physical source and exact collective algebra

Each constituent is the complete physical eight-cell source with native
hopping `2^-34`, unique ground g, and isolated first triplet e_a from the
immutable block-propagation parents. Its actual excitation energy is
Delta. Its native spin current supplies the exact dressed matrix unit

```
A_a=P_g O_a P_t/(60 sqrt(r))=|g><e_a|.
```

The normalization uses the actual positive source residue r. It is an
exact finite spectral service, with no efficient compilation asserted.
All higher physical block sectors and all native gauge constraints are
retained. Every A is gauge invariant, even, and preserves raw L. The
internal triplet index and internal cube momenta stay identifiable.

Place N identical constituent blocks at each site of a finite periodic
cubic graph of V sites and degree z=6. Let

```
A_xa^(N)=N^(-1/2) sum_(i=1)^N A_(x,i,a),
X_xa^(N)=A_xa^(N)+(A_xa^(N))^dag.
```

For dimensionless spacing `0<a<=1`, declare the full physical Hamiltonian

```
H_(a,N)=a^(-1) sum_(x,i)(h_(x,i)-E0)
         -t_a sum_<xy>,alpha X_xalpha^(N) X_yalpha^(N),
D_a=Delta/a,
t_a=Delta(1-a^2)/(12a).
```

Its clock is `exp(-it H_(a,N))`. The replication, exact spectral dressing,
pair creation and annihilation, and spacing-dependent source clock are
architectural choices. They differ from the number-conserving block
exchange. They use no empirical mass, speed, length or phase input.

The subspace with every constituent in g or its first triplet is exactly
reducing. In the subspace symmetric under the N permutations at each
site, normalized occupation vectors have `n_x=sum_alpha n_xalpha<=N`.
Direct counting of the constituent matrix units gives

```
(A_xalpha^(N))^dag=b_xalpha^dag sqrt(1-n_x/N),
A_xalpha^(N)=sqrt(1-n_x/N) b_xalpha.
```

These identities are exact, including the boundary `n_x=N`. No canonical
commutation relation is assigned to a single finite-capacity block.
The onsite energy on this subspace is exactly `D_a sum_x n_x`.

## The free quadratic reference

Replacing the square roots by one gives a precise comparison Hamiltonian,

```
H_free=D sum_x,alpha b_xalpha^dag b_xalpha
       -t sum_<xy>,alpha (b_xalpha+b_xalpha^dag)
                           (b_yalpha+b_yalpha^dag).
```

Here and below D=D_a and t=t_a. Write `b=(q+ip)/sqrt(2)`,
`Q=q/sqrt(D)`, `P=sqrt(D)p`. Up to its specified additive constant,

```
H_free=1/2 sum_x,alpha [P_xalpha^2+m^2 Q_xalpha^2]
       +D t sum_<xy>,alpha (Q_xalpha-Q_yalpha)^2,
m^2=D(D-12t)=Delta^2>0.
```

The finite Gaussian ground is unique and normalizable, including the
zero momentum mode. Its gap is m=Delta. Its complete normal-mode spectrum
and centered two-point function are

```
omega(k)^2=D^2-4Dt sum_i cos(k_i),
<Omega|Q_alpha(k,t) Q_beta(-k,0)|Omega>
 =delta_alpha,beta exp(-it omega(k))/(2 omega(k)).
```

The momentum menu is `k_i=2pi n_i/L`. The additive vacuum energy is
`E_free=(3/2) sum_k (omega(k)-D)`. It is nonpositive because the unsqueezed
vacuum has energy zero. Stability requires `12t<D`. Setting `12t=D`
produces a zero-frequency finite periodic mode and supplies no normalized
Gaussian ground. The positive mass in the declared scaling avoids that
failure.

## Quantitative finite-capacity ground comparison

The assertions in this section concern the actual finite hard-core source,
with all complementary physical sectors retained.

First, its global ground energy occurs in the all-low constituent sector.
Every pattern of higher constituent eigenstates is invariant and carries
a strictly positive onsite cost. Replacing those fixed states by g gives
a variational state in the all-low sector with identical expected
interactions and the higher onsite cost removed. Thus a higher-sector
ground cannot minimize the full energy.

In the all-low occupation basis, every off-diagonal Hamiltonian entry is
nonpositive. Replacing a vector by its entrywise absolute value cannot
increase its energy. Averaging a nonnegative global ground vector over
all constituent permutations gives a nonzero symmetric global ground.
Consequently the full global ground energy equals the minimum in the
exact symmetric occupation representation. No exclusion of nonsymmetric
sectors is assumed in this argument. Uniqueness or a uniform gap of the
full finite-capacity source is not claimed.

All off-diagonal entries of the symmetric Hamiltonian have magnitude at
most the corresponding free boson entry. The diagonal entries agree.
For a nonnegative minimizing vector psi_N, extended by zero to Fock space,

```
E_N=<psi_N,H_N psi_N>
    >=<psi_N,H_free psi_N> >=E_free.
```

This is a ground-energy comparison, not an operator ordering between the
two Hamiltonians on arbitrary signed vectors.

For an explicit upper bound, let Omega be the free Gaussian ground,
`n=sum_x n_x`, `M2=<Omega,n^2 Omega>`, and `P_K=1_(n<=K)`.
Take `K+2<=N`, `q=M2/(K+1)^2<1`, and use the normalized trial
`phi_K=P_K Omega/||P_K Omega||`. Every vector in its support belongs to
the symmetric finite-capacity source.

The quadrature inequality gives the free form bound

```
H_free <= (D+2zt)n+3ztV,
|E_free|<=3ztV.
```

Since `(H_free-E_free)Omega=0`, the projected trial energy error equals
the discarded component's positive energy divided by the retained norm.
The elementary tail bounds

```
Pr(n>K)<=M2/(K+1)^2,
<n 1_(n>K)><=M2/(K+1)
```

therefore give

```
epsilon_cut=[(D+2zt)M2/(K+1)
             +6ztV M2/(K+1)^2]/(1-q).
```

On the total-occupation subspace `n<=K`, the exact collective factors
obey `1-sqrt(1-u)<=u`. With the two intermediate occupation steps kept,
this yields

```
||(H_N-H_free)P_K|| <= epsilon_HP,
epsilon_HP=12 z V t (K+1)(K+2)/N.
```

The factor includes all three native polarizations and all `zV/2` bonds.
It is deliberately conservative and independent of a numerical fit.
Thus

```
0<=E_N-E_free<=epsilon_E=epsilon_cut+epsilon_HP,
1-|<Omega,psi_N>|^2<=epsilon_E/Delta.
```

The second bound uses the actual free gap, applied to the nonnegative
symmetric global ground chosen above. It supplies a stationary physical
reference of the full finite-capacity source with an explicit comparison
error. Finite spectral preparation of that reference is a declared
service; no efficient algorithm or observed vacuum interpretation follows.

An explicit dilute bound also holds for every such global ground. For the
collective lowering operators, counting occupations gives

```
sum_alpha (A_alpha^dag A_alpha+A_alpha A_alpha^dag)
 =2n_x+3-(2n_x^2+2n_x)/N <=2n_x+3.
```

This exact symmetric-sector identity suffices for the chosen global
ground. The inequality `(A+A^dag)^2<=2(A^dag A+AA^dag)` and the two-site
quadrature bound imply

```
H_N >=(D-2zt)n-3ztV,
<psi_N,n psi_N>/(NV) <=3zt/[N(D-2zt)]
                         =3(1-a^2)/(2Na^2).
```

The coefficient `D-2zt=Delta a` is positive. The declared N can increase
so this occupation per constituent vanishes together with epsilon_E.

## Complete low-density dynamical control

All norm comparisons use the same time parameter t in the actual source
and the free comparison. Let `P_N` denote the local-capacity Fock embedding
and let a normalized free input psi have

```
M4(T)=sup_(|s|<=T) <psi(s),(n+2)^4 psi(s)>,
psi(s)=exp(-isH_free)psi.
```

Finite Gaussian states and their finite normal-mode excitations have
finite M4(T), calculable by their exact Gaussian moments. For eigenstates
this quantity is time independent. In a real normal-mode basis each
vacuum mode has canonical occupation mean
`u=(D/omega+omega/D-2)/4`. Its first four number cumulants are

```
u, 2u(u+1), 4u(u+1)(2u+1), 8u(u+1)(6u^2+6u+1).
```

They follow by differentiating the exact generating function
`[(1-q)/(1-q z^2)]^(1/2)`, with `q=u/(1+u)`. Cumulants add over
the three polarizations and every real normal mode. A one-particle
normal-mode excitation replaces that mode's generating function by
`z[(1-q)/(1-q z^2)]^(3/2)`. These formulas supply M2 and M4 directly,
without a sampled tail or fitted truncation. For complex momentum modes,
cross terms of functions of total occupation vanish between distinct
real-mode one-excitation states, so the largest-mode bound applies.

Set `C=D+12 t |edges|`. The elementary
four creation/annihilation terms on each bond give

```
||H_free v||, ||H_N P_N v|| <= C ||(n+2)v||.
```

For `K+2<=N`, splitting the free state at total number K yields

```
||(H_N P_N-P_N H_free)psi(s)||
 <=epsilon_HP+2C sqrt(M4(T))/(K+3).
```

Duhamel's formula gives the explicit comparison

```
||exp(-itH_N)P_N psi-P_N exp(-itH_free)psi||
 <=|t|[epsilon_HP+2C sqrt(M4(T))/(K+3)].
```

Normalizing the projected input and removing the final projection adds
at most `2 sqrt(<psi,(1-P_N)psi>)` when the tail bound is uniform along
the free orbit. For a stationary or single-energy free input that follows
from `M4(T)/(N+3)^4`. If a separate actual-ground preparation replaces
the projected free reference, add its established state-vector error,
at most `sqrt(2 epsilon_E/Delta)` after phase alignment.

Every complete finite positive-operator-valued detector on this common
embedded carrier then has total-variation error at most the displayed
state-vector error. Complementary records and preparation failures are
retained. The result concerns bounded records and finite normal-mode
experiments. It does not replace an operator by a canonical boson at a
fixed nonzero occupation per constituent.

For a literal bounded field detector, take a finite-volume free Weyl
unitary W(f,g), the exponential of real smeared Q and P. Its finite source
representative is the contraction `C_N=P_N W(f,g) P_N` on the admitted
symmetric carrier, extended by zero on its full complement. The Halmos
matrix

```
U_C=[[C_N, sqrt(I-C_N C_N^dag)],
     [sqrt(I-C_N^dag C_N), -C_N^dag]]
```

is an exact unitary on the source and a retained two-state apparatus.
A controlled Hadamard test with the apparatus initially in its first
state uses a control phase phi and has four retained Kraus records,

```
K_(s,0)=[I+s exp(i phi) C_N]/2,
K_(s,1)=s exp(i phi) sqrt(I-C_N^dag C_N)/2,  s=+1,-1.
```

Their squared operators sum exactly to I. Their probabilities on a
normalized admitted vector psi are

```
p_(s,0)=[1+||C_N psi||^2+2s Re(exp(i phi)<psi,C_N psi>)]/4,
p_(s,1)=[1-||C_N psi||^2]/4.
```

Summing the internal apparatus flag gives
`p_s=[1+s Re(exp(i phi)<psi,C_N psi>)]/2`. The phase selects the real
or imaginary component. No internal flag or control outcome is discarded.
The compression is a declared finite spectral service and uses the same
admitted source, reference and clock as the evolution bound.

Let `tau=||(I-P_N)W psi||` for psi in the admitted carrier. The full
uncontrolled Halmos output differs from `W psi tensor |0>` by norm
`sqrt(2) tau`: both the success compression and retained leakage have
norm tau and occupy orthogonal flag sectors. The normalized controlled
Hadamard construction places half the squared amplitude in that branch,
so its full output difference is exactly tau. Recording the two registers
preserves this norm comparison and bounds their complete total variation.

For telescoping on arbitrary ideal free states, extend C_N by zero on
the Fock complement. The controlled-unitary difference is at most

```
||(I-P_N)W psi||+(1+sqrt(2))||(I-P_N)psi||.
```

This follows by splitting psi into admitted and complementary parts;
the uncontrolled unitary difference on the latter is at most twice its
norm, and the normalized control supplies `1/sqrt(2)`. Thus a common
tail bound tau_* before and after each Weyl operation gives a simple
per-instrument error `4 tau_*`. This conservative factor includes every
retained flag.

Here is an explicit moment bound for entire sequences. Let the canonical
displacement of kick j have complex amplitude eta_j. Put

```
R=a^(-1) sum_j ||eta_j||,
F(R)=1+8R+6R^2+4R^3+R^4.
```

Free evolution amplifies a canonical displacement by at most `1/a`:
in each real normal mode its phase-space matrix has norm at most
`max(D/omega,omega/D)<=1/a`. Every ideal coherent control history is
therefore a displacement of norm at most R acting on the stationary
Gaussian reference or an evolving single normal-mode packet. The
occupation identity
`W^dag n W=n+eta.b^dag+eta^*.b+||eta||^2`, expanded twice, gives

```
||(n+2)^2 W(eta) psi|| <=F(||eta||) ||(n+2)^2 psi||.
```

For example the cross terms `(n+2)L+L(n+2)`, with
`L=eta.b^dag+eta^*.b`, are bounded by
`4||eta|| ||(n+3)^(3/2)psi||<=8||eta|| ||(n+2)^2psi||`;
the other terms give the displayed coefficients. Consequently
`M4_sequence<=F(R)^2 M4_base`. The base one-particle bound takes the
largest canonical occupation among the actual normal modes.

Before each fresh control Hadamard, its alternative histories are
orthogonal in the retained controls. A source-only moment is the average
of their displaced-state moments. Subsequent control rotations preserve
it. This controls the coherent superpositions and all complete records
without conditioning on a rare outcome. Individual normalized
postselected records would require their own probability factors.

Set `tau_*=sqrt(M4_sequence)/(N+3)^2`. With m evolution segments of
total duration T and r Weyl tests, telescoping gives

```
complete vector error <=initial_vector_error
  +T[epsilon_HP+2C sqrt(M4_sequence)/(K+3)]
  +(2m+4r) tau_*.
```

The two tails per evolution segment follow by comparing the projected
Duhamel relation with the full ideal state. This bound applies to all
four records at every test, with their full conditional source outputs.
The executable supplies exact rational moment and resource bounds for
two kicks of canonical displacement norm at most `1/8` each. Unbounded
field moments are interpreted only on their stated finite-moment domain.

## Retained native translation and orientation reference

The two detector bits record the control outcome and Halmos flag. A Weyl
kick can also change the native internal translation character and cubic
component. The closed instrument therefore includes a separate retained
finite reference. Its records are joint relational records. Its source
Hamiltonian and time parameter are the ones specified above.

The concrete group is the finite cover

```
G=(Z/4Z)^3 semidirect 2O,     |G|=64*48=3072,
```

where 2O is the binary octahedral group and its proper signed-permutation
image acts on translation vectors modulo four. On one complete native
constituent, the actual antiperiodic translations satisfy
`T_i^2=(-1)^F`, `T_i^4=I`, and commute. Proper cube rotations have the
native two-spinor lift, whose central minus identity also acts as fermion
parity. Their exact relation is `R T_i R^dag=T_j^(sign)` when
`R e_i=sign e_j`. The representation of G can have a kernel; no faithful
3072-element action on the source is required. The smaller translation
image on the admitted even sector is `(Z/2Z)^3`. Using that image on the
entire complementary carrier would omit the odd-sector parity relation.

These are symmetries of the actual complete 119-state-link Hamiltonian.
`atomic_momentum_response/checks.py` proves the antiperiodic hopping-word
identities and rotation-translation relations. `oriented_links/checks.py`
checks all proper cubic spin, onsite singlet and hopping-word transformations,
including the antiperiodic seams. The full coordinate-current and
orientation identities in `peter_weyl_links/checks.py` apply those words
to every link column. Occupied-column controls also replay the actual
Hamiltonian compiler. `collective_spin_response/model.py` verifies the
native Fourier-current transformations `O_alpha` on that same source.
The local `reference.py` replays these checks, verifies all 576 spatial
rotation compositions and all spin-lift compositions up to the central
sign, and checks an exact complete three-component detector example.
Thus cubic covariance here does not rest on a postulated rotation of a
reduced triplet.

Apply the same native group element diagonally to every constituent of
`H_(a,N)`. The onsite Hamiltonians commute with this action. The isolated
spectral projectors commute with it, and the three dressed matrix units
transform by real signed permutations. The sum over their equal-component
quadrature products is invariant. This proves `[H_(a,N),U_g]=0` on the
full source. The action preserves the common Gauss kernel, raw L and
fermion parity. Its cubic part is discrete; no continuous spin symmetry
or Lorentz-vector representation is asserted. These internal constituent
translations differ from translations of the macrocell graph.

Let the retained reference be `R=C[G]`, with orthonormal basis `|g>`,
left regular action `L_h|g>=|hg>`, and reference Hamiltonian zero. Define

```
J psi=|G|^(-1/2) sum_g U_g psi tensor |g>,
C_tilde=sum_g U_g C_N U_g^dag tensor |g><g|.
```

Then `J^dag J=I`, `(U_h tensor L_h)J=J`,
`[C_tilde,U_h tensor L_h]=0`, and `C_tilde J=J C_N`.
The same identities hold for the adjoint and hence for both positive
defect square roots by finite functional calculus. The four Kraus
operators and the Halmos unitary therefore have complete lifts satisfying
`K_tilde J=J K`; the control and flag bits carry the trivial group action.
Their squared effects sum to the joint identity on the entire reference
and source carrier. They preserve covariance also on inputs outside
`range(J)`. All complement and failure outputs are retained.

Since `(H_(a,N) tensor I)J=J H_(a,N)`, every sequence of these lifted
operations has exactly the original record probabilities on `J psi`.
The entire final source-reference vector is J applied to the original
source vector in each retained record. J preserves norm errors, so all
finite-time, preparation, projection and total-variation bounds above
hold unchanged for these joint records. The ideal Gaussian comparison
uses the same signed-permutation action and the same J; its embedding
intertwines the group action.

For the chosen actual global ground `psi_N`, the lifted pure state
`J psi_N` is exactly stationary at the same energy and is jointly
invariant. Its source marginal is
`|G|^(-1) sum_g U_g |psi_N><psi_N| U_g^dag`, a stationary ground-supported
density. This need not equal the original source density. The displayed
Weyl expectations describe the joint relational instrument, rather than
unreferenced source-only Weyl measurements on that marginal. No uniqueness
of the full finite-capacity ground follows from the lift.

Preparing this joint state is part of the declared finite spectral
service. A complete preparation can be lifted by the same Kraus rule,
with its input represented by J and all outputs retained. Equivalently,
finite unitary state preparation within the invariant subspace `range(J)`
can prepare the selected ground from a declared invariant seed. A
controlled application of `U_g` to a uniform reference realizes J as an
isometry, but does not by itself prove a covariant preparation from an
arbitrary unreferenced input. There is no efficient native implementation
claim for this service.

One shared 3072-state register suffices for the entire finite sequence.
With r Weyl tests and both bits retained separately at each test, the
declared field apparatus dimension is `3072*4^r`. Its capacity is
`log(3072)+2r log(2)` nats, in addition to the full kinematic source
capacity. The resource rows report both counts. Other preparation
controllers or services require their own storage accounting.

## Massive relativistic refinement and its scope

Put `k=ap` for physical block momentum p. The exact free dispersion is

```
omega_a(p)^2=Delta^2 [1+(1-a^2)/(3a^2) sum_i (1-cos(ap_i))].
```

The limit is the massive Klein-Gordon dispersion

```
omega_0(p)^2=Delta^2+c^2|p|^2,
c^2=Delta^2/6.
```

Delta is the rest frequency in the declared source clock. Conversion to
an observed mass unit or to a measured speed requires an independent
physical calibration that is absent from this construction.

On any bounded physical momentum window the explicit error is

```
|omega_a(p)^2-omega_0(p)^2|
 <=Delta^2 a^2 [|p|^2/6+sum_i p_i^4/72],
|omega_a(p)-omega_0(p)|
 <=Delta a^2 [|p|^2/12+sum_i p_i^4/144].
```

The second bound uses both frequencies at least Delta. For a centered
Fourier-mode covariance at `|time|<=T`, the absolute error is at most

```
(1/(2Delta^2)+T/(2Delta)) |omega_a(p)-omega_0(p)|.
```

The canonical field normalization is `Q_cont(x)=a^(-3/2)Q_x`, with the
same normalization for P and the corresponding lattice smearing weights.
Finite momentum-window smearings therefore converge with these bounds.
There is also a direct volume error. For smooth Fourier smearings supported
inside `[-P,P]^3`, set
`F(p)=fhat(p)^* ghat(p) exp(-it omega_a(p))/(2 omega_a(p))`.
Let `h=2pi/ell` be the momentum spacing in a physical box of side ell.
The midpoint-cell estimate gives

```
|h^3 sum_p F(p)-integral F(p)dp|
 <=(2P+2h)^3 (3h/2) max_(i,p)|partial_i F(p)|.
```

The factor `(2pi)^(-3)` converts this to the field covariance convention.
The derivative is bounded uniformly in a because
`|partial_i omega_a|<=Delta |p_i|/6` and `omega_a>=Delta`.
The smearings vanish smoothly outside their support, and the momentum
window is contained in the Brillouin zone. Thus this error tends to zero
as the physical box grows, independently of the spacing error above.
Its exact limiting two-point distribution is

```
W_alpha,beta(t,x)=delta_alpha,beta
 integral d^3p/(2pi)^3 exp(-i sqrt(Delta^2+c^2|p|^2)t+ip.x)
                       /[2 sqrt(Delta^2+c^2|p|^2)].
```

Together with the canonical commutator and Gaussian higher moments, this
is the free massive scalar field covariance, with limiting speed c after
rescaling spatial units. Its mass-shell measure supplies the usual
Poincare covariance. The triplet labels are internal labels of these
three field components. No Lorentz-vector or observed particle-spin
identification is inferred from the native cube current.

At each finite box and spacing, K and then N can be chosen from the
explicit bounds so the finite-capacity record error is arbitrarily small.
A diagonal sequence choosing these integers before taking the spacing
and volume limits transfers this free relativistic record subsector to
finite systems composed of the actual interacting blocks. This is a
controlled dilute, increasing-capacity limit. Fixed N, the uniform
original lattice, a Lorentz-covariant interacting continuum, and measured
Standard Model masses receive no such theorem.

The finite-N square roots produce genuine occupation-dependent many-body
interactions. Their controlled contribution is bounded by epsilon_HP and
the explicit tails on every admitted experiment. The diagonal refinement
sets that contribution to zero. The constituent eight-cell interactions
retain their exact role in Delta, r, g and the dressed modes, but no
nonzero scattering interaction between the limiting Gaussian particles
is established.

With fixed D and t strictly below criticality, spacing refinement freezes
the kinetic coefficient. At criticality the periodic Gaussian zero mode
loses its normalizable ground. The spacing-dependent native clock,
pairing coefficient and increasing constituent capacity stated above are
essential hypotheses of this construction.

## Relation to spin-wave results

The exact collective square root is a multicomponent version of the
[Holstein-Primakoff representation](https://doi.org/10.1103/PhysRev.58.1098).
Its finite constituent matrix elements are derived and checked here.
Correggi, Giuliani and Seiringer's rigorous
[spin-wave free-energy result](https://arxiv.org/abs/1312.7873) concerns
the three-dimensional Heisenberg ferromagnet and a low-temperature free
energy. It supplies no direct theorem about this paired neutral source,
its selected reference, or a relativistic record limit. The inequalities
above prove the needed ground and finite-time bounds for the declared
Hamiltonian rather than transferring those different hypotheses.
