# Complete local records at fixed polynomial action and time

The [mixed gauge state](../particle_local_gauge_state/PROOF.md) has exact
matching local marginals across periodic volumes. Combining that fact
with the [fixed-spacing locality theorem](../particle_field_causality/PROOF.md)
constructs the limit of every specified finite local record experiment
at fixed polynomial order, spacing, coefficients and model duration.
Finite source-word compressions approximate those records with explicit
factored resource bounds. The time is not shortened as volume grows.
This is a thermodynamic local-record statement for the polynomial action,
not an exact-sign, chiral or relativistic continuum theorem.

## 1. Exact source marginal consistency

Fix the normalized one-link polynomial f_L, local maximally mixed
sixty-mode Fock density and site-gauge twirl of the owning state proof.
The same L is used in every periodic volume. A retained finite patch
consists of cells, all their outgoing links, and their full Fock modes.
The patch and endpoints of its retained links have a nonwrapping
embedding in every compared torus.

For any fixed gauge transformation, its unitary action factors over
link and Fock registers. Partial tracing the exterior removes the
unitary conjugations on exterior factors and their normalized product
densities. The remaining integral depends only on the independent Haar
variables at vertices incident to retained links, and on gauges at
retained fermion cells. Every other Haar variable integrates to one.
Its density is therefore exactly the same finite operator in every
such torus. No limiting estimate, normalization ratio or boundary
projection is needed. The local maximal Fock trace is also exactly
consistent. For nested patches, ordinary partial tracing gives the
smaller marginal.

This argument retains the shared gauge variables at vertices incident
to several links. Replacing the site twirl by independent endpoint
twirls on each link would change these correlations. Partial tracing
also retains boundary-charge mixtures: the patch is not projected onto
field-only singlets. Fermionic restriction is well defined on the
even local algebra used by the detector; the full tracial Fock state
has the usual compatible CAR restrictions.

The positive, normalized cylinder functionals consequently define one
state on the algebraic local source algebra and its norm completion.
This is a state of the declared mixed field family, without a vacuum,
Gibbs or Lorentz-invariance identification.

## 2. Reference systems and local preparation

The full source family is physically realized with its conjugate
purifier, as specified by the owning state proof. Purifiers in different
volumes have different spaces and need not have matching local joint
marginals. The default comparison traces these purifiers and compares
complete local source records. The dynamical exterior is also traced;
its full final state is not compared.

An independent fixed reference R can be added to all experiments.
More general correlations can be constructed by applying the same
fixed local completely positive preparation or instrument to a fixed
patch and R before evolution, retaining all its outputs. For every
larger patch containing its support, partial trace commutes with that
operation. Its whole-patch/reference marginals therefore continue to
match exactly across volumes. The operation is part of the declared
experiment, with its locality, gauge/reference and resource constraints.
It may change the original source state. The default numerical resource
certificate does not include an additional preparation of this kind.

An arbitrary reference extension is covered only when its compatible
whole-patch/reference marginals are explicitly supplied. In particular,
a pure canonical purification of a patch cannot generally be attached
while preserving the original correlated global source density: purity
would factor the patch/reference pair from the exterior. The finite-box
simulator may prepare its own canonical purification to implement its
mixed input. Its cutoff error can retain that reference, which is traced
when comparing with the periodic family. These are distinct statements.

## 3. Complete local polynomial records at any fixed duration

Fix a>0, m, kappa, L and any finite T>0. The actual periodic source is
`H_m=H_E+dGamma(h_m)-kappa sum_p Re Tr_M15(U_p)/15`. All these local
parameters are independent of volume. Choose any specified finite
sequence of even, gauge-invariant local instruments, or explicitly
compensated instruments with their local references. Their supports
and complete feedback policy are fixed; all outcomes, failure labels
and zero labels are retained. The default menu is the origin occupation
PVM with labels0,...,60. The preparation is the matching entire local
marginal, not an independently chosen vacuum.

For a radius-R open comparison box, keep every outgoing boundary link
as an electric spectator. Let r contain the observed cells, d=2m+1,
b=max(d,2), D=R-r-2b, and let J_root,L_rate be the owning polynomial
interaction constants. For k intermediate local instruments and s
observed cells, the complete local history bound is

```
C_R=2 J_root[2(R+d)+1]³+3 kappa[2(R+2)+1]³,
e_R <= min(1,(2k+1)2s C_R T 2^(-D)[exp(2 L_rate T)-1]). (1)
```

The same collar proof applies inside a periodic torus of side at least
`2(R+2b+d)+1`; all relevant supports have a nonwrapping coordinate lift.
The initial entire box marginal matches exactly by section1. For fixed
action, apparatus and T, C_R grows polynomially in R while 2^(-D)
decays exponentially. Thus e_R tends to zero. Any two sufficiently
large periodic complete record laws differ by at most2e_R through the
same open experiment. Finite-outcome probability vectors are complete
in total variation, so they have a unique limit; the limit is normalized
and retains every record. The same argument covers compatible retained
reference outputs from section2 in half trace norm. It claims no
closeness of the global source state.

This polynomial locality theorem is independent of Wilson gaps. Local
measurements can break translation symmetry and inject energy, so the
all-time roughness theorem for uninterrupted invariant evolution is
not reused between instrument steps. No regulator gap premise is needed
for (1). The high-kappa constant may make the radius enormous, but it
is fixed across the volume limit. This distinguishes the construction
from a planner whose kappa grows with the comparison box.

## 4. Exact logarithmic volume planning

The planner evaluates no giant exponential. Since e<4,

```
exp(2 L_rate T)-1 <= 2^ceil(4 L_rate T).
```

Let `c_R=ceil(log2((2k+1)2s C_R T))`, evaluated with exact integer bit
lengths of a positive rational. Then

```
e_R <= min(1,2^(c_R+ceil(4 L_rate T)-D)).               (2)
```

Starting at R=r+2b+1 and doubling R eventually makes the exponent at
most minus the requested bit accuracy. The code stores R and the
exponent as integers and exports the requested dyadic bound; it never
allocates2^R. The default example chooses T=1, m=0 and the fixed local
state planner's L and kappa, with24 bits allocated to (2). Large radii
are a mathematical resource estimate, not an executed large simulation.

## 5. Finite representation code for the same record experiment

The actual open marginal is supported in per-link degree at most L,
so `K0=3(2R+1)³L` suffices. Its density bra does not add ket degree.
The simulator uses `P_K H_open P_K`, retains the entire Fock carrier
and prepares that actual marginal with its own compensating purifier.
The [owning cutoff theorem](../particle_observer_outputs/FIELD_ROUTE.md#a-controlled-representation-cutoff-for-complete-observer-histories)
uses bounded interaction norm

```
B=(2R+1)³ b+3(2R+1)³ kappa,
x=BT,    step=max(4,2m+1).
```

Electric evolution preserves the degree filtration. The default local
occupation pointers also preserve degree and are complete unitaries
when their records are coherently retained. Interleave those fixed
pointer unitaries with the electric interaction picture before the
Dyson expansion. They have norm one and do not change the word-step
count. At total insertion order j, the time-ordered volume over all
evolution intervals is bounded by T^j/j!. Consequently the first
omitted total order is `M=floor((K-K0)/step)+1`, and the complete
history error is bounded by `2 sum_(j>=M) x^j/j!`. One does not restart
the initial-degree bound after each measurement.

For M>=4x and M>=1, the factorial tail satisfies

```
sum_(j>=M) x^j/j!
 <= (x^M/M!)/(1-x/(M+1))
 <=2(e x/M)^M <=2(3/4)^M.
```

Here `M!>=(M/e)^M` follows by integrating log over1,...,M, and e<3.
Since `(3/4)^3<1/2`, a complete-record bound is

```
e_K <=4(3/4)^M <=2^(2-floor(M/3)).                    (3)
```

Choose `M=max(1,ceil(4x),3(bits+2))` and
`K=K0+step(M-1)`. This proves e_K<=2^(-bits) without evaluating a giant
factorial. It keeps all occupation and pointer outcomes. A different
apparatus with nonzero coordinate-degree cost needs that cost included;
the certificate does not silently apply (3) to such an apparatus.

The finite code therefore approximates the limiting complete local
polynomial law by at most e_R+e_K. With24 bits allocated to each,
the default total is2^(-23). Preparation, numerical evolution and
compiled detector errors are additional complete-channel terms.
Source primitives can specify the exact finite operators. This
resource certificate does not execute them at its enormous dimensions.

## 6. A nonzero source-local dynamical diagnostic

The mixed Fock marginal does not make every local process static.
At m=0,a=1, the actual hopping block to one neighboring cell is
`B=-(4/17) beta P_mu^± tensor R(U)`, with
`Tr(B†B)=480/289`. For the full tracial Fock state, the CAR gives

```
J=i(c_x† B c_y-c_y† B† c_x),
<J>=0,    <J²>=Tr(B†B)/2=240/289.
```

This is pointwise in the gauge coordinates because R(U) is unitary.
The six distinct neighbors on a nonwrapping cubic patch have zero
cross-current expectation under the complete Fock trace. Electric,
onsite and plaquette terms commute with origin occupation n_x. Hence

```
<dot(n_x)²>=1440/289.                                  (4)
```

The code checks the native spin/internal singular values and the exact
two-mode CAR trace factor. Initially `rho` commutes with n_x, so a
complete initial occupation read does not alter its unconditioned state.
For the complete two-time occupation law, with all61² label pairs,

```
d²/dt² sum_(n,n') (n'-n)² Prob(n at0,n' at t)|_(t=0)
 =2<dot(n_x)²>=2880/289.                               (5)
```

On a finite code K>=K0+1, the first hopping image is contained and
the same quadratic coefficient is exact. This derivative statement
uses distinct neighbors, as in the large boxes of this package; it is
not asserted for a one-cell torus. It supplies a nonzero local temporal
response of the actual state and action, with no finite-time probability
enclosure, scattering or mass interpretation. The initial complete
number law is `binomial(60,n)/2^60`.

## 7. The mathematical boundary

The package establishes consistent finite marginals, a controlled
fixed-spacing thermodynamic limit of complete local polynomial records,
and finite-code approximation at fixed nonzero duration. The charged
purifiers remain physical components of each experiment but are not
identified across volumes. The maximally mixed Fock preparation is not
a Lorentz vacuum. No empirical calibration is removed.

Passing to exact-sign dynamics requires a local operator-response bound
on this same quantum state. A bound on local low-Wilson-mode weight is
useful but does not by itself control its electric dependence or the
many-body Duhamel commutator. Passing to a relativistic continuum also
requires uniform estimates as spacing and polynomial order change,
surviving interactions, chiral content and a matching detector dictionary.
