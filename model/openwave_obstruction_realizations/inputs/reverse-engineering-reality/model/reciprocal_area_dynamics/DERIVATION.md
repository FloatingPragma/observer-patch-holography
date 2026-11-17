# Exact thermal geometry response on the reciprocal source

The complete carrier is `C^8 direct_sum C^8`. Write `M_g=beta H_g`, with
the rational matrices derived independently from the source action in
`finite_area_response.response.matrices`. They are

```
M_g=a(r_g)D+c(r_g)(XXI+IXX),
D=diag(2 popcount(mask)+3),
a(r)=r/4+1/(7r)+3r^2/28,       c(r)=-1/(7r),
r_0=9/10, r_1=11/10,          beta=1/sqrt(14).
```

All sixteen energy labels are retained. Put `w_i=exp(-beta E_i)` and
`pi_i=w_i/Z`. The unchanged transition on populations is

```
T_(j,i)=w_j/[16(w_i+w_j)]       for opposite geometries,
T_(j,i)=0                     for distinct same-geometry labels,
T_(i,i)=1-sum_(j!=i) T_(j,i).
```

The reservoir, energy readout and four fair native proposal coins are the
specified services of `native_equilibration.reciprocal_vacuum`. A complete
collision controller and its repeated thermal preparation are not executed
here. Ideal finite transition statements do not certify those source costs.

## One metric, energy and constrained entropy family

The original action-density dictionary gives
`g(r)=diag(-r^2,r^2,1/r,1/r)`. Both branches have proper volume one;
`A_0=20/9`, `A_1=20/11`, and `Delta A=A_1-A_0=-40/99`.
The joint reference is `rho=exp(-beta H)/Z`, with conditional states
`tau_g=exp(-beta H_g)/Z_g`, geometry probability `p*=Z_1/Z`, and
conditional energies `e_g=Tr(tau_g H_g)`.

For every quantum state sigma, including coherences between geometries,
let `p=Tr(Pi_1 sigma)` and
`bar sigma_p=(1-p)tau_0 direct_sum p tau_1`. For `0<p<1`, the difference
`log(bar sigma_p)-log(rho)` is constant on each geometry block. Taking its
expectation proves the exact relative-entropy decomposition

```
D(sigma||rho)=D(sigma||bar sigma_p)+D_Bernoulli(p||p*).       (1)
```

The endpoint version uses the supported geometry block. Since
`beta Tr(sigma H)-S(sigma)=-log Z+D(sigma||rho)`, (1) proves that
`bar sigma_p` is the unique minimum-free-energy state at fixed expected
area. No stochastic projection or extra reset is added to the dynamics.

Define the source-derived number and the conditional excess energy

```
eta=beta(e_1-e_0)/(A_0-A_1),
bar E(p)=(1-p)e_0+p e_1,       E_ex(sigma)=Tr(sigma H)-bar E(p).
```

The exact enclosures prove `0.578<eta<0.579`. Since
`beta bar E(p)+eta<A>_p=beta e_0+eta A_0=:C`, the full identity is

```
S(sigma)+eta<A>_sigma-beta E_ex(sigma)
 =log Z+C-D(sigma||rho).                                  (2)
```

On the minimum-free-energy family this becomes

```
S(bar sigma_p)+eta<A>_p=log Z+C-D_Bernoulli(p||p*).          (3)
```

The first variation vanishes uniquely at p*. The second derivative in p
is `-1/[p(1-p)]`; the second derivative in expected area is
`-1/[(Delta A)^2 p(1-p)]`. This is a fixed-proper-volume variation because
the volume operator is the identity on the complete carrier.
It varies probabilities between two fixed action metrics; it is not a
local metric-field variation or a small-ball shape derivative.
The full-state functional without `beta E_ex` is different. In particular,
the same-geometry energy tangent in the one-step area package continues to
exclude its stationarity on all fixed-volume states. Equation (3) obtains
its variation space from the explicit variational minimization (1), not
from an assertion that the actual dynamics preserves it.

The number eta depends on this temperature, two selected metric values and
finite field cutoff. Two distinct area labels make their conditional mean
energies affine in area. This conversion supplies no local edge multiplicity,
regional ultraviolet subtraction, universal Newton coefficient, or
Einstein small-ball stationarity statement.

## Exact elimination of the matter populations

Use the real Hilbert space of observables with inner product
`<u,v>_pi=sum_i pi_i u_i v_i`. Let `P=T^T`. Pairwise detailed balance
`pi_i P_(i,j)=pi_j P_(j,i)` makes P self-adjoint. Set

```
f_i=1_(i in geometry1)-p*,       v=<f,f>_pi=p*(1-p*),
Qz=z-<1,z>_pi 1-f<f,z>_pi/v,
a=<f,Pf>_pi/v,      b=QPf,      C=QPQ on range(Q).          (4)
```

The range of Q has dimension fourteen. If
`c=sum_(i in geometry0,j in geometry1) pi_i P_(i,j)`, the Dirichlet-form
identity gives `a=1-c/v`. Certified values are
`c=0.08817371496...` and `a=0.63916680624...`.

For any initial population vector q, put `u_0=q/pi-1`. Detailed balance
gives `u_n=P^n u_0`. Decompose

```
x_n=<f,u_n>_pi=p_n-p*,       u_n=(x_n/v)f+z_n,
x_(n+1)=a x_n+<b,z_n>_pi,
z_(n+1)=C z_n+(b/v)x_n.                                  (5)
```

Solving the second equation and substituting into the first proves

```
x_(n+1)=a x_n+sum_(j=0)^(n-1) k_j x_(n-1-j)+F_n,
k_j=<b,C^j b>_pi/v,       F_n=<b,C^n z_0>_pi.              (6)
```

This is an identity at every finite update count. The measured area is
`<A>_n=A_0+Delta A(p*+x_n)`.
For the source matter replacement `q_i=p_g*/8`, `x_0=0`; its forcing
gives the one-step area response in the parent package. Mixture amplitude
epsilon multiplies every x, F and energy deviation exactly.

Energy belongs to the same eliminated system. For
`h_i=E_i-<E>_pi`, let `alpha=<h,f>_pi/v=e_1-e_0`. Then

```
z_n=C^n z_0+(1/v)sum_(j=0)^(n-1) C^j b x_(n-1-j),
<E>_n=<E>_pi+alpha x_n+<h,z_n>_pi.                       (7)
```

Thus conditional excess energy and every heat increment are reconstructed
without replacing the matter state by a conditional Gibbs state at each
step. In particular `k_0=0.02170170003...>0` proves that the conditional
thermal family is not invariant. A geometry-only Markov update would omit
a nonzero term of the actual source process.

## Positivity, memory dimension and convergence

Every row of P has exit sum strictly below 1/2. Gershgorin's bound and
self-adjointness therefore put every eigenvalue strictly above zero.
The exact spectral enclosures give `1<beta E_i<5`. Since `exp(4)<63`,
every cross-geometry transition exceeds `1/1024`, and every hold exceeds
1/2. For two distinct states in the same geometry, the eight paths through
the other geometry give
`(T^2)_(j,i)>8/1024^2=1/131072`. Holds cover diagonal states, and a
cross-geometry transition followed or preceded by a hold covers the other
pairs with a stronger bound. Hence

```
T^2 = (1/8192) U_uniform + (1-1/8192) R,
```

where R is a stochastic matrix. Its action contracts the oscillation of
every observable by at most `1-1/8192`. A nonconstant eigenvector of P
therefore has eigenvalue

```
0<lambda<=sqrt(1-1/8192)<q:=1-1/16384.                    (8)
```

The restriction of P to mean-zero observables has norm below q.
Compression shows `0<C<q I` on the range of Q. Consequently

```
k_j>=0,       k_j<=k_0 q^j,
sum_(j>=N) k_j<=16384 k_0 q^N,
|F_n|<=||b||_pi ||z_0||_pi q^n,
|p_n-p*|<=sqrt(v) q^n ||q/pi-1||_pi.                    (9)
```

These bounds use the same finite law and its collision count. No asymptotic
temperature or fitted relaxation time enters them.

The code encloses `b,Cb,...,C^13 b` by outward rational arithmetic and
performs exact-algorithm Gram orthogonalization in the pi inner product.
Every squared norm has a strictly positive rational lower bound. The last
is approximately `4.51237*10^-37`, with enclosure width below `1.1*10^-77`.
All fourteen directions are therefore independent. Since range(Q) has
dimension fourteen, b is cyclic for the self-adjoint operator C. Its
spectrum consists of fourteen distinct positive poles, each with positive
weight in `k_j=sum_l weight_l pole_l^j`. This proves minimal scalar
memory dimension fourteen; these are equilibration poles, not particle poles.

The centered geometry vector is cyclic for P on the fifteen-dimensional
mean-zero space: `Pf=a f+b`, and the space generated by powers of P on f
contains each `C^j b`, by induction using (4). Thus all fifteen nonconstant
relaxation modes are visible in the geometry response. No eigenvalue fit
or floating rank threshold establishes these counts.

## Complete bath records and coherent entropy production

For either native proposal law the sixteen equally weighted proposal words
are in bijection with the two geometry and eight field labels. Each energy
input has 32 proposal/bath outcomes. Keep all 512 rows, including 384
rejections or no-ops. An accepted row from i to j is reversed by selecting
i from j and swapping the recorded incoming/outgoing bath bits. Every
rejected row is its own reverse. The proposal probability remains 1/16;
the thermal qubit probabilities give

```
kappa_r/kappa_reverse(r)=w_j/w_i=exp[-beta(E_j-E_i)].      (10)
```

The code verifies the involution on every complete source label and proves
`Delta E_system+Delta E_bath=0` as an integer identity in the energy labels.
The paired auxiliary proposal word is a valid word of the same law.
The assertion is a stochastic record reversal, not time reversal of native
implementation pulses or a physical undoing of committed memory.

For a path of n updates, reverse the order and each row, and start the
reverse experiment with `q_n=T^n q_0`. At faithful q_0 the path log ratio is

```
Sigma=log(q_0(i_0)/q_n(i_n))-beta sum_(k=1)^n Q_k,
E_forward[exp(-Sigma)]=1,
E_forward[Sigma]=H(q_n)-H(q_0)-beta(<E>_n-<E>_0)
                =D(q_0||pi)-D(q_n||pi)>=0.              (11)
```

The first equality follows by multiplying (10), the integral identity by
summing the normalized reversed path probabilities, and nonnegativity
from relative entropy of the two complete path laws. Nonfaithful inputs
have the same mean identity with the usual support convention; the integral
identity requires the stated absolute-continuity condition.

For an arbitrary quantum input sigma let `Delta_E` be the complete rank-one
energy/geometry readout of this source, and q its populations. The operator
`log(rho)` is diagonal in that basis, and so is `log(Delta_E sigma)`.
The pinching trace identity gives

```
D(sigma||rho)=D(sigma||Delta_E sigma)+D(q||pi).
```

Readout leaves the source energy and expected area unchanged. After this
readout and n thermal collisions, with final state diagonal at q_n,

```
Delta S-beta Q
 =D(sigma||Delta_E sigma)+D(q||pi)-D(q_n||pi)>=0.          (12)
```

The first term is the entropy increase in the unread energy measurement,
not a derived measurement heat cost. Combining (2) and (12) proves

```
Delta S+eta Delta<A>-beta Delta E_ex = Delta S-beta Q.    (13)
```

Every record remains in the joint ledger. Formula (12) is a marginal
entropy balance when subsequent source updates ignore the retained tape;
it does not assert a decrease of complete-record distinguishability by
erasing those records. The quantity `-beta Q` is thermal log-likelihood
flow. An outgoing bath can be correlated and nonthermal, so that flow is
not asserted to equal its von Neumann entropy change. Fresh baths, the
readout, intervention apparatus and pulse-controller work remain the source's
declared services. The field-randomization preparation energy is accounted
in its parent receipt and excluded from subsequent bath heat Q.

## Arithmetic and physical domain

Characteristic-polynomial root isolation uses rational intervals of width
below `10^-100`; independent LDL congruences certify all matrix eigenvalue
counts. A 200-term exponential series on `|x|<=3` has tail bounded by
`|x|^201/[201!(1-|x|/202)]`. Every later interval operation rounds outward
at denominator `10^120`. The proof tests every Gram lower endpoint, rather
than inferring positivity from a printed decimal approximation.

The metric, energy, entropy, records and collision count concern this one
finite scalar family. The metric lapse r gives `d tau=r dt`, proper-time
generator `H/r` and proper inverse temperature `r beta` in a fixed branch.
A common model-coordinate generator does not calibrate a proper duration
per collision. The family has no spatial coarse-graining map, common local
stress tensor with the Weyl/contact theory, edge-area density or the
regional variations required by the Einstein theorem. Equations (6) and
(13) establish relaxational finite geometry and thermal balance on the
specified services, without filling those missing premises.
