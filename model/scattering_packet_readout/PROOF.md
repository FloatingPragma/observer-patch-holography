# A finite experiment for the moving contact scattering phase

Use the complete infinite cubic two-fermion source of
[moving contact scattering](../moving_contact_scattering/PROOF.md), with
`gamma=1`, `c=sqrt(3)/a^2`, `1/U=a^2 W0/2-a^3/pi`, and `0<a<=1/65536`.
The comparison clock is `J=6/a^2-H`, with free comparison `J0=6/a^2-H0`.
All eight upper-edge valleys and all four spin components are retained.
This result concerns a prepared two-particle scattering experiment. It
supplies neither a many-body ground vacuum nor a relativistic field theory.

## An actual free source packet

Let `K=aP`, `|P_j|<=1`, and write `s=B_K-E`, where `B_K` is the exact
free upper edge. Let `rho_K(E)` be the contact-vector spectral density of
the original free source. In each fiber, the normalized contact vector
chi is the constant spin singlet. Define

```
q(s)=sqrt(3003/512) [1-16(s-1)^2]^3,  3/4<=s<=5/4,
q(s)=0 otherwise,
r(p)=sqrt(3003/2048) (1-p^2)^3,      |p|<=1,
r(p)=0 otherwise,
phi(P)=r(P1)r(P2)r(P3).
```

Both `integral q(s)^2 ds` and `integral |phi(P)|^2 dP` are one, since
`integral_-1^1 (1-x^2)^6 dx=2048/3003`. Use Lebesgue measure dK for the
total-momentum direct integral and normalized torus measure for relative
momentum. The actual source vector is

```
psi_K=a^(-3/2) phi(K/a) F_K(H0(K)) chi,
F_K(E)=q(B_K-E)/sqrt(rho_K(E))
```

on the stated energy interval and zero elsewhere. The spectral theorem
proves normalization. This is a functional-calculus vector of the full
source, with its exchange antisymmetry, spin amplitudes and every valley.
There is no freely supplied phase gate or selected valley. Its preparation
depends only on the free spectral density, rather than on the interaction
scattering phase.

## Analytic source bounds on this packet window

The source coarea formula uses `u=r omega`, `|omega|=1`,

```
t=R-a^2 s/(2sqrt(3)), lambda=R^2-t^2,
R^2=sum_j cos(K_j/2)^2,
F(r omega)=r^2 f(r^2,omega),
f=sum_j cos(K_j)omega_j^2
  +[sum_j cos(K_j/2)sin(K_j/2)omega_j sqrt(1-r^2 omega_j^2)]^2/t^2,
J(u)=product_j(1-u_j^2)^(-1/2),
N(u)=1-sum_j sin(K_j/2)^2 u_j^2/t^2.
```

Extend s to `|s-1|<=1/2` and each K component to a complex disk of radius
a about a real point in `[-a,a]`. Elementary sine/cosine power-series
bounds give `3/2<|t|<2`, `|f-1|<=32a^2` for `|r|<=1/2`,
`Re(lambda)>=a^2/3`, `|Im(lambda)|<=a^2`, and `|lambda|<=2a^2`.
For example `|K_j|<=3a`, `|cos K_j-1|<=5a^2`,
`|cos(K_j/2)sin(K_j/2)|<=2a`; the last sum in f has modulus at most
`4a`, giving `|f-1|<=5a^2+16a^2/(9/4)<32a^2`.

Split the u cube into its radius-1/2 ball and exterior. On the inner ball
deform the radial contour to `r=exp(i pi/4)x`, `0<=x<=1/2`, followed
by the radius-1/2 circular arc back to the real endpoint. The analytic
square roots have their branches continued from r=0. There `|JN|<=4`.
On the ray,

```
|i x^2-lambda| >= (x^2+a^2)/9,
|F-lambda| >= (x^2+a^2)/10.
```

For the first bound use the real part when `x^2<=2a^2`, and the imaginary
part when `x^2>=2a^2`; the second subtracts `32a^2 x^2`.
On the arc `|F-lambda|>=1/5`. The same bound holds on the exterior:
for real u, `Re(F)>= (1-13a^2)|u|^2`. Its J-weighted normalized cube
measure is at most one. No endpoint of the original cube is discarded.

The ray contributes less than 9 to W: its radial absolute integral is
at most 20 and the angular factor is `4/pi^2<4/9`. The arc contributes
less than 2; the exterior contributes at most 10. Thus `|W|<21`, and
the prefactor `|t/(2c)|<a^2` gives the convenient analytic bound

```
|g_K(B_K-s+i0)| <= 32a^2.                              (1)
```

For real s the deformed contour is the boundary value approached with
positive imaginary energy: the corresponding radial poles lie below the
real contour. The deformed integral provides its holomorphic continuation
through this particular energy window. The bounds above also exclude
zeros on the deformed path throughout the complex parameter disks.
This establishes (1) as an analytic statement, rather than a bound on a
finite broadening or a differentiated principal-value integrand.

Cauchy's formula with s-radius 1/8 and K-radius a gives on the real window

```
|g_s|<=256a^2, |g_ss|<=4096a^2,
|rho_s|<=86a^2, |rho_ss|<=1366a^2, |rho_(K_j)|<=11a.
```

Here the analytic continuation of rho is
`(conj(g(conj(K),conj(s)))-g(K,s))/(2pi i)`; on real parameters it is
the positive free density. Its lower bound follows independently from
the real coarea identity

```
rho=a^2 beta sqrt(lambda) A(lambda)/pi^2,
beta>=3/4, sqrt(lambda)>=a/2, A(lambda)>=3/8.
```

The existing Jacobian bounds also give `A(lambda)<=2+48lambda<3`.
Consequently

```
a^3/100 <= rho <= a^3/2.                               (2)
```

## A finite collision time, with both wave-operator tails

The source trace-class theorem supplies complete wave operators on each
contact cyclic absolutely continuous fiber. The other reducing channels
are retained. With the convention
`Omega_+=lim_(t->+infinity) exp(itJ)exp(-itJ0)` and the corresponding
negative-time limit, `S=Omega_+^* Omega_-`. Their finite approximants give

```
S_T=exp(itJ0)exp(-2itJ)exp(itJ0)
   =exp(-itH0)exp(2itH)exp(-itH0),  t=T>0.              (3)
```

The scalar shifts cancel exactly. For any free contact packet whose
energy overlap A has two integrable derivatives and zero endpoint values
and first derivatives, twice integrating by parts gives
`|<chi,exp(-itJ0)psi>|<=||A''||_1/t^2` in each fiber. The same argument
holds for its scattered packet, with overlap `A_out=sqrt(rho) S q`.
Thus Cook's integral and fiberwise norm integration give

```
||(Omega_-(T)-Omega_-)psi|| <= U C_in/T,
||(Omega_+(T)-Omega_+)Spsi|| <= U C_out/T.
```

Completeness gives `Omega_-psi=Omega_+Spsi`. Inserting this identity
between the finite wave operators proves

```
||(S_T-S)psi|| <= U(C_in+C_out)/T.                     (4)
```

Both tails are required. Replacing the outgoing packet by the incoming
one in the second tail is not assumed. The common K-dependent edge phase
has modulus one and does not affect the fiber bound.
The finite echo evolves the full two-particle carrier. Any finite-time
amplitude in the interacting bound subspace is included in (4); no bound
state, dark channel or unsuccessful history is projected away.

For explicit constants set `u=sqrt(rho)`. Equations (1)--(2) imply
`|u|<=a^(3/2)`, `|u'|<=500sqrt(a)`,
`|u''|<=2000000/sqrt(a)`. The profile has
`||q||_1<=2`, `||q'||_1<=48`, `||q''||_1<=192`.
The existing denominator estimate `|D|>=7a^3/44`, where `D=1/U-g`,
and its reflected ratio S give

```
|S'|<=3220/a, |S''|<=11000000/a^2.
```

Indeed the second derivative is bounded by
`2|g''|/|D|+4|g'|^2/|D|^2`, with rational coefficient
`510033920/49<11000000` after using `a<=1`. Product differentiation
therefore yields
`C_in+C_out<=40000000/sqrt(a)`. Since `U<=4/a^2` and `a<=1`,

```
||(S_T-S)psi|| <= 160000000/(a^3 T).                   (5)
```

The normalization of phi is one. Consequently (5) controls its actual
moving packet, including the continuum of total momenta.

## Spatial truncation and a finite lattice

The finite-dimensional divided-difference formula gives
`||d f(A)||<=4||f'||_infinity ||dA||` for a Hermitian four-by-four matrix.
Explicitly, in an eigenbasis of A its derivative maps an entry E_ij to
`f[lambda_i,lambda_j] E_ij`, using f' on repeated eigenvalues. The mean
value formula bounds every divided difference by `||f'||_infinity`.
Thus its operator norm is at most its Frobenius norm, which is at most
`||f'||_infinity ||E||_F<=2||f'||_infinity ||E||`. The displayed factor
four is a conservative upper bound; no spectral separation is required.
The density derivatives above imply `|F|<=40a^(-3/2)`,
`|F_s|<=173000a^(-5/2)`, and
`|partial_(K_j)F_K(E)|<=368000a^(-7/2)` at fixed E.
Use `|partial_(K_j) B_K|<=2/a` and
`||partial H0||<=4/a^2`. Also `||partial_j phi||_2<20`.
The relative-momentum derivative of psi is at most `3000000a^(-9/2)`;
its total-momentum derivative is at most `5000000a^(-9/2)`.
Since `partial_(p1)=partial_K+(1/2)partial_q` and
`partial_(p2)=partial_K-(1/2)partial_q`, each of the six original
one-particle momentum derivatives has L2 norm at most `7000000/a^5`.

The original momenta p1,p2 are torus coordinates. Take K to be the
representative of `p1+p2` in the principal torus cell and set
`q=p1-K/2` modulo `2pi`. On this chart `p2=K/2-q` modulo `2pi`.
The relative-coordinate function is periodic in q. The K profile is
supported in `[-a,a]^3` and vanishes with its first two derivatives at
its edges, strictly inside the principal cell. Therefore crossing the
principal K-cell cut contributes no boundary distribution. Changes of
the q representative also give no seam, by the exact torus periodicity
of H0 and the contact vector. The derivative identities above consequently
hold as weak derivatives on the original two-particle torus, and Parseval
applies to the physical lattice Fourier coefficients.

Parseval and the coordinate moment bound show that truncation to both
particles in `[-R,R]^3` loses squared norm at most
`300000000000000/(a^10 R^2)`. When `20000000/(a^5 R)<1`, the truncated
vector has nonzero norm. Normalizing it costs at most

```
||psi_R-psi|| <= 40000000/(a^5 R).                    (6)
```

The cutoff retains both spin modes at every included site. It respects
exchange antisymmetry and is a projection of the actual source packet.

Enlarge the box by L sites in every direction for the finite experiment.
On the full N=2 occupation basis, at most twelve hopping entries leave
each configuration, with magnitude `c/2`; the absolute hopping matrix has
norm at most `6c<12/a^2`. In the interaction picture of the diagonal
contact and scalar terms, these absolute matrix elements are unchanged.
The three segments of (3) have total absolute duration 4T. A history that
sees the boundary requires at least L hops. The finite/infinite difference
on psi_R is bounded by

```
2 sum_(n>=L) (48T/a^2)^n/n! <= 2^(2-L),
provided L>=288T/a^2.                                (7)
```

Here `n!>=(n/e)^n` and `e<3` reduce the tail to a geometric series.
This bound is independent of U, whose effect is included in the diagonal
interaction picture. It does not invoke a finite-dimensional S matrix.

## Preparation, native evolution and every output

The compact vector psi_R is computable without an eigenvector oracle.
Its Fourier coefficients use only finite four-spin matrix functional
calculus and the explicit free coarea density. Positive radial derivatives
bound the coarea root away from a multiple root, so interval bisection
terminates. The angular integral uses six smooth cube-face charts of the
sphere. Subdivision with rational trigonometric and square-root enclosures
converges uniformly on these compact domains. The positive lower bound
(2) controls every reciprocal. Fourier coefficient integration and finite
normalization then give rational complex approximants at any requested
vector error. The interaction constant W0 is also computable: split out
a ball about each of its eight nodes, bound its integrable inverse-square
tail, and subdivide the nonsingular remainder.

The coefficient-approximation budget controls the normalized finite vector
defined by these rational approximants against psi_R. The separate
amplitude-preparation budget controls the complete circuit isometry against
that specified rational target. Normalization is certified before its
angles are computed, using (6) to retain a strictly positive vector norm.
These two budgets concern different maps; neither assumes an exact
continuous-amplitude state preparation.

A finite complex amplitude tree prepares those rational approximants from
blank field bits. A binary amplitude tree prepares their magnitudes; a
complete Walsh-diagonal circuit supplies their phases. Rational zero
coefficients are decidable. All retained nonzero phases and rotation
angles have computable enclosures. This is an algorithm for the specified
source vector, rather than an assumed successful state replacement.
The existing exact twelve-port source identities implement its Pauli
rotations. Every first proper-code export failure is an output event.

For V sites and M=2V field bits, an explicit loose preparation pulse
bound is `11M 2^M`. The original hopping/contact Hamiltonian has at most
16V Pauli terms and sum of coefficient norms at most `Lambda=10V/a^2`.
Controlled exponentials use the exact two-factor Pauli construction,
including scalar terms. First-order splitting with n steps on each of
the three intervals in (3) has full-Fock operator error at most
`6T^2 Lambda^2/n`. A loose pulse count is

```
N_pulse <= 11M 2^M + 768n V(M+1) + 10.                (8)
```

The dyadic port interface has per-pulse complete-isometry error at most
`29*2^(-p)`. Enclosing the computable angles separately and choosing
`29 N_pulse 2^(-p)` below its budget controls the complete native dilation.
The precision schedule can be evaluated using an upper bound on
`log2 N_pulse`; it does not allocate the enormous tape merely to report
its cost. Exact rational integrated ports refine the readable source law;
fixed binary64 port commands are not extrapolated to arbitrary precision.

The executable interfaces are the signed integrated-port exponentials in
[`source_scalar.compiler`](../observer/source_scalar/compiler.py), the exact
controlled Pauli construction in
[`native_spectral_detector.circuits`](../native_spectral_detector/circuits.py),
and the refining dyadic commands in
[`native_spin_observer.precision`](../native_spin_observer/precision.py).
Negating every port coefficient and reversing a finite word implements its
inverse proper-code unitary; the complete exports of that inverse word
are retained independently. Consequently `exp(+2iHT)` in (3) is constructed
from admitted negative-angle controls with its own cost, rather than from
passive time reversal. The underlying source simulator executes these
port exponentials and exports sequentially. It adds no undeclared idle
Hamiltonian between commands. H is the same compiled field action used by
the inherited active-control experiments; no simultaneous unrefocused
background-H evolution or passive repair is inserted in this protocol.
Hardware with such an additional drift requires a separate control proof.

Prepare a blank interference bit in |+>, apply controlled (3), and read
Pauli Y by a native basis rotation followed by the complete pointer
instrument. Its plus probability is

```
p_plus=(1+Im<psi,S_T psi>)/2.                           (9)
```

All six pointer labels, preparation stops and first-export failures are
retained. Coarsening to the event "valid Y-plus" retains every other
record in its complement. No source phase formula is compiled as the
controlled operation in (9).

## A separated complete observation law

At `a=1/65536` the inherited uniform scattering error is
`2101/32768<1/10` on this packet window. Its limiting phase satisfies

```
Im[(1-i sqrt(s))/(1+i sqrt(s))]=-2sqrt(s)/(1+s)<=-49/50.
```

The squared magnitude is at least `48/49>(49/50)^2` throughout the window.
Allocate vector/operator error at most 1/100 to each of seven operations:
Cook truncation, spatial packet truncation, coefficient approximation,
amplitude preparation, finite volume, product formula and native port
execution. Equations (5)--(8) give a finite schedule for every allocation.
Their complete-dilation errors add to at most 7/100. Therefore

```
p_plus < 1/8.
```

Deleting the interaction in the same finite echo gives the identity
exactly; with the same error ceiling its complete plus probability is at
least `1/2-7/100`. The two complete recorded laws are separated on this
one event by more than 3/10. This is a finite native instrument theorem,
with a very large resource schedule. The fine-spacing experiment is not
claimed to have been numerically or physically executed.

For sampling, 2000 independent complete histories for each law estimate
both event probabilities within 1/20 simultaneously with failure below
1/1000: Hoeffding gives `4 exp(-10)`, and the positive exponential series
proves `exp(10)>4000`. Fresh prepared copies and apparatus memory are
included services; rejected histories are counted, not recycled as fresh
samples. The executed small control tests the source echo, orientation,
pointer, untouched references and complete failure bookkeeping separately
from the analytic scattering and finite-cutoff theorem.
