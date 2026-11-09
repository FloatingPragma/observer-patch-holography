# Local influence on the same gauge-coordinate field

## 1. Reused source and finite theorems

The Hamiltonian is the declared coordinate-link program in
[FIELD_ROUTE](../particle_observer_outputs/FIELD_ROUTE.md), with the same
native M15 representation, full Dirac Fock space, electric Casimirs,
polynomial `h_m` and Hamiltonian duration. Its coordinate completion is an
admitted architecture choice. It is distinct from the native raising-only
link program. No action parameter or time is fitted to a particle mass.

The following existing results fix the proof strategy and its boundaries.

| Existing result | Reused statement | Domain retained here |
| --- | --- | --- |
| [Full-Fock field, local dynamics](../full_fock_field/README.md#local-thermodynamic-dynamics-at-fixed-spacing) | Even-interaction CAR commutator chains and onsite interaction picture | Replace its two-mode bond constant by the actual Dirac/M15 polynomial constants; preserve odd-field grading |
| [Complete local matter records](../source_family_compatibility/LOCAL_MATTER_READOUT.md#a-local-complete-history-bound) | Defer instruments coherently and telescope the same local controls | Keep every record, feedback branch and untouched reference |
| [Finite causal observer net](../../Lean/QFT/FiniteCausalObserverNet.lean) | Declared local observable algebras, isotony and remote repair hypotheses | Equal-time even algebras commute; the present finite-time influence estimate has tails |
| [Refining lattice cones](../../Lean/Geometry/RefiningLatticeCausalCone.lean) | A supplied layered move law has inner/outer reachability bounds | Hamiltonian time evolution is not a finite-depth circuit or an exactly vanishing cone |
| [Native causal dataflow](../native_causal_dataflow/README.md) | Literal program event dependencies and scheduler scope | No redefinition of the Hamiltonian clock as FIFO depth |
| [Gauge-coordinate cutoff](../particle_observer_outputs/FIELD_ROUTE.md#a-controlled-representation-cutoff-for-complete-observer-histories) | Complete-record representation-cutoff and numerical bounds | Add those errors to the parent locality law; a total-degree cutoff is not a local tensor cutoff |

The parity-sensitive chain argument is the one used in
[Nachtergaele, Sims and Young, Theorem 3.1 and its proof](https://arxiv.org/html/1705.08553v3#S3).
Even interactions preserve parity. Disjoint operators commute if at least
one is even, while two odd operators anticommute. The proof below evaluates
the interaction weights for this source and treats the unbounded electric
generators by bounded onsite approximation. It does not infer an
interacting relativistic theory from the imported locality theorem.

## 2. Cells, local terms and uniform constants

Attach each outgoing positive-direction link `L²(G,Haar)` to its tail
cell, where `G=S(U(3)×U(2))`. Every cell contains all 60 Dirac/M15 CAR
modes. Distances are Manhattan distances in lattice units, on an open
cubic box or the corresponding periodic graph. Periodic distances must
include wraparound; the open-box volume comparison below uses no wraparound.
Local fermion algebras retain their CAR grading. Gauge-coordinate operators
are parity even. Physical gauge-invariant instruments are a subalgebra of
this ambient construction; Gauss restriction does not increase norms.

Write the source Hamiltonian as

```
H = H_on + sum_x Phi_x + sum_p Phi_p,
H_on = sum_cells (outgoing electric Casimirs + declared onsite terms),
Phi_p = -kappa Re Tr_M15(U_p)/15.
```

The onsite terms are coordinate independent, bounded and parity even.
The electric Casimirs are positive self-adjoint operators on their usual
compact-group domains. Each finite-volume `H_on` is self-adjoint, and its
bounded perturbation above is self-adjoint on the same domain. Couplings
need not be small. A declared onsite finite scalar and its even charged
pairing interaction can also be included in `H_on`; doing so does not
prove that it is the source's physical vacuum or particle sector.

For `d=2m+1`, let `P_x` project onto the 60 one-body modes at x and set

```
K_x = (P_x h_m + h_m P_x)/2,     Phi_x = dGamma(K_x).
```

The polynomial's paths starting or ending at x visit only cells in the
Manhattan ball `B_d(x)`. Link coefficients on every such path are supported
there. Consequently `Phi_x` has that support, including paths that return
to x. Also `sum_x K_x=h_m`; this is an exact decomposition of the source
Hamiltonian. For every coordinate configuration,

```
||K_x||_1 <= 60 ||h_m||,
||dGamma(K_x)|| <= ||K_x||_1 <= J_m,
J_m = 60 (2+epsilon_m)/a,     epsilon_m=16(15/17)^(m+1).
```

The first inequality uses the rank-60 projector and the trace-norm
triangle inequality. The second follows by summing occupied eigenvalues
of a Hermitian matrix on the complete exterior algebra. It holds on every
occupation sector, not only one-particle states. Taking the essential
supremum over coordinate configurations proves the same bound on the
full coordinate/Fock Hilbert space. This replaces the volume-dependent
global Fock norm by a local interaction norm.

The ball cardinality is
`N_d=(4d^3+6d^2+8d+3)/3`. A cell belongs to at most `N_d` rooted stars,
each containing at most `N_d` cells and diameter at most `2d`. A plaquette
may be assigned its four vertices, has diameter two and norm at most
`|kappa|`; each cell belongs to at most twelve plaquettes. Therefore

```
L_m := sup_x sum_(Z contains x) |Z| ||Phi_Z|| 2^diam(Z)
     <= J_m N_d² 4^d + 192 |kappa|.                    (1)
```

This constant is independent of the box volume, electric Casimir ceiling,
particle number, initial state and untouched reference dimension. It is
not uniform in m or inverse lattice spacing.

At m=0 the actual source projectors give a much smaller constant. With
`P_mu^±=(I±gamma_mu)/2`,

```
h_0 = (25/(17a)) beta - (4/(17a)) sum_mu beta W_mu.
```

Move the first term into `H_on`. Each nearest-neighbor block has thirty
singular values `4/(17a)` and thirty zeros. The full-Fock Hermitian bond
therefore has norm `J=120/(17a)`, the sum of the thirty positive one-body
eigenvalues. The code verifies the exact native spin projector and its
two-site spectrum before multiplying by the fifteen retained internal
components. This yields the sharper source constant

```
L_0 <= 24 J + 192 |kappa|.                              (2)
```

Six incident bonds, two cells per bond and distance weight two give the
factor24. No electric coefficient enters (1) or (2).

## 3. Electric domains and the graded influence estimate

The onsite interaction picture preserves each support, parity and operator
norm. Unbounded electric evolution need not be operator-norm continuous on
all local bounded link operators. It is insufficient to assert the finite
matrix theorem applies without handling this domain.

For each fixed finite box, replace every onsite self-adjoint generator by
its bounded spectral clipping. These clipped onsite generators preserve
supports, parity and the exact local gauge action. Their evolution is norm
continuous and the chain proof applies. The finite sum of clipped onsite
generators converges in strong resolvent sense to `H_on`; adding the fixed
bounded interaction preserves that convergence. Their unitaries and
Heisenberg evolutions converge strongly, uniformly on compact time
intervals on each vector. An operator-norm inequality uniform in the
clipping ceiling passes to this strong limit. Thus the estimate below holds
on the uncut coordinate Hilbert space. No electric-energy moment of the
initial state is required for these bounded-observable estimates.

For homogeneous `A_X,B_Y` on disjoint cells, define
`[A,B]_gr=AB-(-1)^(parity(A) parity(B))BA`. Put

```
M_xy = sum_(Z contains x and y) ||Phi_Z||.
```

The iterated even-interaction commutator argument bounds the contribution
of n connecting interaction supports by

```
2 ||A|| ||B|| (2|t|)^n/n! sum_(x in X,y in Y) (M^n)_xy.
```

Each chain is counted at least once by choosing an intermediate cell in
every adjacent support intersection. The bound therefore also covers
multiple-site interactions. Since
`sup_x sum_y M_xy 2^dist(x,y) <= L_m`, the triangle inequality for distance
gives the corresponding row bound `L_m^n` for `M^n`. The n=0 term vanishes
for disjoint supports. Summing the factorial series yields

```
||[alpha_t(A_X),B_Y]_gr||
 <= 2 ||A|| ||B|| |X| 2^-dist(X,Y) (exp(2 L_m |t|)-1).    (3)
```

The trivial normalized bound is two. Equation (3) is an ordinary
commutator estimate when either operator is even and an anticommutator
estimate when both are odd. Using the ordinary commutator for two odd CAR
fields would be false even at t=0. Local detector effects, currents and
gauge transformations used here are even.

Equation (3) gives exponentially small influence outside an approximate
cone with speed `2 L_m/log(2)` cells per model time. The tails are not
identically zero. Multiplication by a gives a speed in the declared length
units; neither this bound nor a spherical Lorentz-group parametrization
identifies it with the speed of light or constructs boost symmetry.

## 4. A complete detector consequence on the same state and clock

Use the source preparation with all60 modes filled at the origin, other
cells empty and constant Haar wavefunctions on the links. Its determinant
character is trivial, so it is gauge invariant. Let a remote operation be
the unitary phase of a gauge-invariant Wilson plaquette at support distance
D from the origin. Compare the two experiments that apply this operation
or identity, then evolve by this same H for the same time t and measure
the complete local occupation PVM. Every occupation record, including
zero-probability labels, is retained.

For any union of detector labels, its effect E is a positive contraction
on X. The change in its probability is at most `||[alpha_t(E),U_Y]||`.
Taking the supremum over unions bounds total variation of the entire
record law by

```
min(1, 2 |X| 2^-D (exp(2 L_m |t|)-1)).                   (4)
```

The same bound applies to the retained local quantum output and arbitrary
untouched references by effect duality. Reference dimensions introduce
no volume factor. Local complete apparatus must be supported on the
stated cells; a global energy projector is not a local detector by name.

For the exact m=0 program, a=1, kappa=1/100 and t=1/64, the code obtains
`L_0=72816/425`. At D=32, the complete-record TV bound is below `10^-7`.
It is a rational upper certificate on the large field experiment, not an
assertion that its exponential-size state or pulse tape was allocated.
This finite polynomial has a loose sign approximation; this calculation
is a locality control, not a relativistic accuracy test.

## 5. Volume control for complete local histories

Let X lie in the cube of radius r and contain s cells. Compare the actual
polynomial Hamiltonian on the open cube of radius R with any larger box,
on the same clock. The initial joint reduced state on the entire
radius-R simulation box and untouched reference must agree, for example
by restricting one common initial state. Agreement only on the observed
cells X is insufficient: a different neighboring preparation can affect
X at short times. Put
`b=max(d,2)` and assume `D=R-r-2b>=1`.
Keep link registers on deleted boundary edges as electric spectators;
their onsite electric terms are unchanged. They do not acquire an
unmentioned boundary interaction or an independent time convention.

It is necessary to compare the polynomials of the actual Wilson kernels.
A smaller-box polynomial drops paths that leave and return. Their
difference is supported in the width-d boundary collar. Every differing
rooted kinetic star is bounded by twice `J_m`. A valid upper bound on the
sum of boundary perturbation norms is

```
C_R = 2 J_m (2(R+d)+1)^3 + 3 |kappa| (2(R+2)+1)^3.       (5)
```

For the sharp m=0 decomposition use `J_m=3J`, grouping the three forward
bonds at each root. Formula (5) deliberately counts excess roots. Its
dependence is on the chosen collar, not the unbounded exterior volume.
The inside/outside split leaves each onsite electric term unchanged.
All boundary terms have distance at least D from X. Duhamel's identity
and (3), integrated for at most |t|, give

```
B_R(t)=min(2, 2 s C_R |t| 2^-D (exp(2 L_m |t|)-1))        (6)
||alpha_t^large(A)-alpha_t^R(A)|| <= ||A|| B_R(t).
```

For k local Stinespring interventions and a final local output, use the
same feedback policy and retain its coherent records. In the interaction
picture each evolved local intervention differs by at most B_R(T), where
T is the whole duration. Telescoping the product on both sides of an
arbitrary dual output effect and including the final local field effect
gives the conservative bound

```
half_trace(output_large,output_R) <= min(1,(2k+1) B_R(T)). (7)
```

This is a complete local output comparison, including arbitrary reference
correlations and all outcomes. The exterior field is not part of that
output. An initial joint box/reference half-trace error adds directly;
it is not included in (7). Apparatus and cutoff errors are also separate.
The rational planner finds an R for any positive requested error
at fixed m,a,kappa,T,r,s,k; a polynomial in R times `2^-R` tends to zero.

Equation (6) constructs fixed-spacing local dynamics by a norm-Cauchy
volume limit on each local observable, for a consistent local interaction
family. On all local `B(L²(G))`, this need not be a point-norm-continuous
time group because of electric evolution. No stronger continuity claim
is used. Finite normal-state detector probabilities have the stated
strong-unitary continuity. Selecting an infinite-volume interacting
vacuum or a uniform spatial-refinement state is a separate requirement.

## 6. Why the total-degree cutoff needs its own error

The source cutoff `V_K` bounds the total representation-word degree over
all links. Its projector is not a tensor product of local link projectors.
Consequently compressed distant coordinate operators need not commute.

An exact two-rotor witness uses the K=1 span
`|0,0>,|1,0>,|-1,0>,|0,1>,|0,-1>`. The parent shifts T1 and T2† commute.
After total-degree compression P,

```
[P T1 P, P T2† P] |0,1> = |1,0>,
||[P T1 P, P T2† P]|| = 1.                              (8)
```

The two links can be arbitrarily far apart. The symbolic code verifies
the complete compressed matrices. The actual nonabelian source has the
same nonzero witness: its native q=6 singlet is a nontrivial unit-modulus
character f of G, with Haar mean zero. On distinct links, `f1 f2` is
orthogonal to the whole total-degree-one space, while `conjugate(f2) f2=1`.
Thus for that actual projector P, `[P M_f1 P,P M_conjugate(f2) P]f2=f1`,
so its norm is at least one. The exact norm-one matrix above is the
two-rotor control; equality is not asserted for the larger M15 compression.
This excludes exact support locality
of this compression, not locality of the uncut source. One must apply
(3)-(7) to the parent, then add the separate complete cutoff errors.
For two experiments with cutoff TV errors eta1 and eta2, their compressed
record contrast is at most `min(1,parent_bound+eta1+eta2)`. The field-route
cutoff proof supplies such errors on its stated initial degree and bounded
interaction domain. Numerical Taylor errors enter in the same way.

A product of local Casimir cutoffs could retain exact support locality;
that is a different approximation specification and would need its own
convergence bound. It is unnecessary for the parent result above.

## 7. Physical boundary

The result supplies a finite, volume-uniform causal influence bound and
complete local record control for the same declared interacting field.
It supplies no exact light cone, Lorentz boosts, vacuum selection,
interacting chiral continuum, particle mass, or Einstein stress tensor.
The constants do not control the combined m→infinity and a→0 limit.
That limit needs a common prepared-state family and uniform local
correlation, gap/accuracy and field-renormalization estimates. Physical
particle readout additionally needs the same-source spectral and detector
dictionary. These are output requirements, not demands that the declared
architecture be unique.
