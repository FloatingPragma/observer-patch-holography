# Quantitative native readout of the ordered-spin packets

The ordered-spin spectral window has an explicit finite implementation with
certified success, rejection and conditional-state errors. A compact energy
filter is constructed here, including its computable Fourier coefficients,
time truncation and native control budget. Its application uses Heisenberg
conjugations, so no absolute ground-energy calibration is required. The
result strengthens the general finite LCU observation in
[FA102](../native_joint_spin_spectroscopy/README.md): it supplies the missing
filter, a power-law success bound for the actual FA92 packets, spatial/time
cutoffs, and a fully costed coherent SELECT circuit.

The infinite-state statement assumes usable marginals of the specified
ordered state. The finite-state statement instead uses FA98's source-computable
positive-field ground. Their distinction is retained. No effective convergence
rate of those finite grounds to a specified infinite marginal, sharp mass shell,
or identification with the continuum Stone-field vectors is asserted.

## 1. Source and a completely explicit packet

The source and selected state are those of [FA92](../spin_goldstone_sector/README.md):

```
H_J = J sum_<xy> (Sx.Sy-1/4), J>0, x in Z^3,
omega = infinitesimal staggered-field-selected zero-field ground,
K >= 0 its ground GNS generator, m >= m0=1/1000.
```

Each unoriented bond occurs once. All momenta here are the recentered internal
spin momenta of that proof. Write Fourier integrals with measure
`d^3k/(2pi)^3`. Choose a computable envelope F of unit L2 norm with Fourier
support in `1/2<=|k|<=1`, and certified constants

```
Lambda >= (1/2) sum_(alpha in {0,1}^3) ||partial^alpha F||_1,
Lambda4 >= (1/2) sum_alpha ||(1+|X|_infinity)^4 partial^alpha F||_1.
```

For `0<p<=1`, put `f_p(x)=p^(3/2) F(px)` and
`A_p=sum_x f_p(x) Sx^2`. Poisson summation has no aliases: the support lies
strictly inside the Brillouin zone. Thus `sum_x |f_p(x)|^2=1` exactly, and its
momentum support is `p/2<=|k|<=p`. Iterating the one-dimensional estimate
`h|u(x)|<=integral_cell |u|+h integral_cell |u'|` gives

```
lambda_p=(1/2) sum_x |f_p(x)| <= Lambda p^(-3/2).           (1)
```

There is no unspecified envelope oracle in the concrete family. We supply

```
F(X)=C exp(3iX1/4) product_j sinc(Xj/64)^8,
C=(1/64)^(3/2)/(pi*q)^(3/2),
q=(1/15!) sum_(k=0)^7 (-1)^k binom(16,k) (8-k)^15.
```

The box-convolution formula for the Fourier transform of sinc gives
`integral sinc(y)^16 dy=pi*q`. It also puts the support in
`5/8<=|k|<=sqrt(51)/8`, strictly inside the required annulus. The normalization
is therefore exact and computable. On `|y|<=1/2`, sinc is at least `23/24`;
`(23/24)^16>1/2` proves `C<1/128`. With `delta=1/64`, elementary pointwise
bounds on `u(x)=sinc(delta*x)^8` give

```
||u||_1 <= L=16/(7delta),   ||u'||_1 <= D=208/21,
||(1+|x|)^4 u||_1 <= L4=128/(3delta^5),
||(1+|x|)^4 u'||_1 <= D4=256/delta^4.
```

For example, use `|u|<=min(1,|delta*x|^-8)`, and bound the derivative by
`8delta/3` for `|delta*x|<=1` and `16delta |delta*x|^-8` outside. These give
the explicit rational choices
`Lambda=((7L/4+D)(L+D)^2)/256` and
`Lambda4=((7L4/4+D4)(L4+D4)^2)/256` in `envelope.py`.
This envelope has finite polynomial tails rather than being Schwartz.
Smooth Fourier mollification inside the strict annulus converges in the
required L2 and spatial L1 norms, so the FA92 packet inequality extends to it
by operator-norm convergence of A_p. No new infrared estimate is imported.

For the finite cutoff `|px|_infinity<=R`, with `R>=2`, denote the truncated
operator by `A_(p,R)`. The same cell estimate, now using the weight and
`|X|_infinity>=R-1` on omitted cells, gives

```
||A_p-A_(p,R)|| <= Lambda4 p^(-3/2)/(R-1)^4.               (2)
```

## 2. An explicit compact energy filter

Let `0<d<a<=b`, let k_d be the probability density obtained by convolving two
uniform densities on `[-d/4,d/4]`, and set

```
G = 1_[a-d/2,b+d/2] * k_d.
```

This piecewise quadratic function satisfies `0<=G<=1`, is one on `[a,b]`,
and vanishes outside `[a-d,b+d]`. `filters.py` evaluates it with rational
arithmetic. Set `L=b-a+d`, `c=(a+b)/2`. Its inverse Fourier representation is

```
G(E)=integral g(s) exp(i s E) ds,
g(s)=exp(-i c s) sin(Ls/2)/(pi s) * sinc(ds/4)^2.          (3)
```

The removable value at zero is `L/(2pi)`. The pointwise bounds

```
|g(s)| <= min(L/(2pi), 1/(pi|s|), 16/(pi d^2 |s|^3))
```

imply, by splitting at `1/L` and `4/d`,

```
||g||_1 <= (2+2 log(4L/d))/pi <= Cg=(2+2 ceil(log2(4L/d)))/3,
integral |s g(s)| ds <= 16/(3d),
integral_(|s|>T) |g(s)| ds <= 16/(3d^2 T^2).              (4)
```

The second bound follows by splitting only at `4/d`. We use `pi>3` and
`log 2<1` solely for rational upper bounds. Also
`|g'|<= (b+d)L/6`, directly from the compact Fourier integral.

For the actual FA92 band choose

```
a=m0/8, b=32/m0^2, d=m0/16.
```

Then `Cg=28`; the support is positive and contained in
`(m0/16,64/m0^2)` up to harmless endpoints. Since G is one on the already
certified band, the actual ground vector satisfies

```
Psi_p=G(K/(Jp)) A_p Omega,
||Psi_p||^2 >= m0^2 p/96.                                (5)
```

The use of a nonsmooth-at-second-order spline changes no premise: (5) follows
from the proved positive spectral mass of `[aJp,bJp]`. Its L1 Fourier kernel
defines a norm-convergent quasi-local operator directly.

## 3. The instrument and its success cost

Use the site normalization `Lambda p^(-3/2)` from (1), and the time
normalization `Cg+1=29`. The target successful Kraus operator is

```
T_p = [Lambda p^(-3/2) (Cg+1)]^(-1)
      integral g(s) alpha_(s/(Jp))(A_p) ds.              (6)
```

In a ground representation, `alpha_t(A)Omega=exp(itK)A Omega`. Thus its
successful vector is exactly Psi_p divided by the displayed normalization.
In particular

```
q_p=||T_p Omega||^2 >= m0^2 p^4/[96 Lambda^2 (Cg+1)^2].   (7)
```

The fourth power includes the normalization cost of spatial LCU preparation;
the spectral weight proportional to p is not itself an instrument success
probability. The bound is constructive, not asserted optimal. Independent
attempts with fresh input copies have mean count at most the reciprocal
of (7). Recycling a rejected state is not assumed to provide a fresh copy.

For finite site/time lists, PREP has amplitudes proportional to the square
roots of the absolute complex coefficients. SELECT applies their phases and
the unitaries

```
exp(i t H) sigma_x^2 exp(-i t H).
```

Undo PREP and measure the whole label register. Its blank outcome has the
required weighted sum as Kraus operator; all other outcomes complete the
unitary instrument and are retained. A normalization upper bound is made
exact by padding with two opposite-phase identity terms. These cancel in the
blank Kraus operator but remain in the full unitary dilation.

Adding a scalar `z I` to H cancels in every conjugation before summation.
Consequently this filter needs no absolute ground-energy estimate or assumed
laboratory clock. For a finite FA98 positive-field ground, use H_B in both
conjugations. For the infinite ordered state, use H_J. A finite H_B ground is
not substituted as an exact ground of H_J.

`controls.py` independently checks a four-label LCU with complex phases,
all Kraus effects, and the untouched half of a maximally entangled input.
Removing a reject, reversing the filter phase, or replacing conjugation by
one-sided evolution fails these controls. These are exact logical controls,
not a new six-level native execution.

## 4. Finite cutoffs, coefficients and conditional errors

Truncate (3) at `|s|<=T` and use midpoint cells of width h. Let
`beta=B/J` for the finite selection source, and beta=0 for H_J. On the
normalized site block, `||[H,A]||/(J Lambda p^(-3/2))<=12+beta`:
each site touches at most six norm-J bonds, and the selected on-site field
contributes at most B to the Pauli commutator. Comparing the cell-integrated
g coefficients with midpoint coefficients gives the following unnormalized
time-filter errors:

```
tail <=16/(3d^2 T^2),
coefficient quadrature <= ((b+d)L/6) T h/2,
Heisenberg midpoint <= (6+beta/2) Cg h/p.                (8)
```

Give the rational complex coefficients a further total modulus error e_c.
`coefficients.py` implements a terminating algorithm: Machin pi intervals,
range-reduced rational sine/cosine series, and interval multiplication enclose
every midpoint coefficient; precision increases until its allotted error is
met. Thus the coefficient algorithm is explicit. `quadrature_plan` chooses
T, h and e_c with total error below the requested epsilon<1. The coefficient
L1 sum is at most Cg plus the coefficient quadrature and rounding errors,
strictly below Cg+1, so the padding in (6) is valid.

The spatial coefficients also require a finite precision budget. For rational
p the explicit envelope is enclosed by the same trigonometric intervals and
integer square roots; `envelope.coefficient_rectangle` returns intervals for
`b_x=f_p(x)/(2 Lambda p^(-3/2))`. These satisfy `sum_x |b_x|<=1`.
Choose rational complex midpoints with total modulus uncertainty alpha, and
a rational upper bound S on their L1 norm with excess at most nu. Divide all
midpoints by `max(1,S)`. The resulting exactly specified coefficients have
L1 norm at most one and error at most `2alpha+nu`: any rescaling costs at most
`S-1<=alpha+nu`. `contractive_coefficients` chooses alpha,nu each at most
one quarter of a requested e_spatial. This both enforces the padding premise
and adds at most e_spatial to the successful operator error. No decision
whether a general computable amplitude equals zero is needed; the finite
rational midpoints have decidable zeros. Their nonzero magnitudes, phases
and PREP angles have computable square-root/arctangent enclosures. The
constructive finite algorithm uses rational p; the spectral bounds hold
throughout the displayed real interval.

The success-operator error after all normalizations is at most

```
Lambda4/[Lambda (R-1)^4]
  + (tail+coefficient quadrature+Heisenberg midpoint+e_c)/(Cg+1)
  + e_spatial + e_LR + e_product + e_native.             (9)
```

Here e_native is the error of the complete finite native dilation. The other
terms compare its successful Kraus operator to (6); they are not a fictitious
one-to-one identification of all finite rejection labels with a continuous
time register. For fixed R, the inherited completely bounded FA93 cutoff
applied to each controlled local Pauli is

```
B_l(t) <= 6J (2l+1)^2 t 2^(-(l-r)) (exp(48Jt)-1),
t=T/(Jp), r=ceil(R/p).
```

The time/site LCU normalization does not increase this bound. The logarithmic
integer test in `locality_radius` certifies it without allocating an integer
of size exp(T/p). All interaction/ancilla couplings are those of the explicitly
controlled local interventions. There is no extra inaccessible reference
operation. This supplies a finite spatial cutoff for every p and error.

For an ideal input with success q>=q0, apply the success-operator error to
a purification of that input. If eta<sqrt(q0), actual success is at least
`(sqrt(q0)-eta)^2`. Normalization and partial trace give unhalved conditional
trace error at most `4eta/sqrt(q0)`. If the prepared input density differs
in trace norm by delta from the ideal input on the accessible system,
require `delta<(sqrt(q0)-eta)^2`; success remains positive and the bound becomes

```
4eta/sqrt(q0) + 2delta/(sqrt(q0)-eta)^2.                 (10)
```

Contractivity on the successful block gives exactly the delta term for the
accessible conditional output. This includes a supplied regional marginal;
it does not bound correlations with an unknown purifying system. A guarantee
for the joint output with an arbitrary retained reference requires delta to
bound the joint input trace distance, or a preparation diamond-norm bound.
Purifying the ideal input in the operator-error proof does not supply that
stronger preparation guarantee. The exact same-marginal/orthogonal-Bell-state
countercontrol in `controls.py` records this distinction. No unobserved outcome
is conditioned away in claiming the complete native error. In particular eta of
order epsilon*p^2 and delta of order epsilon*p^4 suffice under (7). An
arbitrary trace error independent of p does not certify normalized packets.

For the finite source, the successful instrument composes with FA102's complete positive Fejer
energy-clock POVM `{E_k}`. Its joint successful probabilities are
`omega(T_p^* E_k T_p)`, with total q_p; a coarse rejected label has probability
`1-q_p`. Thus every energy bin and all rejection probability remain present.
The implemented finite instrument retains its finer individual rejected
labels as well. For comparing it to this coarse target, a success-operator
error eta gives successful-block trace error at most `2eta` (purify the
input and use both Kraus norms at most one), and rejected-weight error at
most `2eta`. Therefore the unhalved L1 error of the joint energy/reject law
is at most `4eta+delta+e_detector`, where e_detector is a separately certified
unhalved complete-output error of the subsequent detector. Conditional
energy-bin error instead follows (10). No comparison discards a bin. The
subsequent spectral clock retains FA102's separately certified energy-origin
phase: the scalar cancellation proved above applies to the filter's
Heisenberg conjugations. The Fejer clock's finite bandwidth, resolution and
possible aliasing stay as specified by FA102. A uniformly calibrated infinite
spectrometer is not inferred from this finite composition.

## 5. Explicit native gate budget

For a finite graph with n sites, E edges and A overlapping unordered edge
pairs, put `C=A J^2+E J B`. First-order splitting with N steps approximates
one Heisenberg conjugation for `|t|<=t_max` with error at most
`2 C t_max^2/N`. This is the maximum over SELECT label blocks, not its sum.
Each exchange contributes four encoded Pauli exponentials including its
scalar, and each field term one. Both conjugations and the local Pauli/phase
therefore use at most `r0=2N(4E+n_field)+2` exponentials per label.

For L labels pad to `D=2^q>=L`. Expanding a q-bit control projector in Z
characters uses exactly D commuting controlled-Pauli factors. The existing
native parity circuit uses at most `8w+1` pulses for a weight-w Pauli.
Magnitude-tree rotations plus basis-state phases give PREP and its inverse
in at most `4D^2(8q+1)` pulses. Hence a deliberately conservative explicit bound is

```
P <= 4D^2(8q+1) + L D r0 [8(2n+q)+1].                 (11)
```

There is no hidden coherent-choice overhead. Tree coefficients, phases and
their precision come from the computable finite lists. Approximate every
ideal pulse angle within `2^-b`; this adds at most `P2^-b`. The FA93
source-read dyadic port commands add at most `29P2^-b` to the complete
isometry. `native_plan` chooses b for their combined `30P2^-b` budget.
It does not extrapolate fixed binary64 controls to zero error.

Ground preparation, subsequent spectroscopy, and optional complete C/E
repairs are added separately. FA98 gives a terminating source-only finite
ground algorithm; its cost is not hidden in (11). Repairs may be inserted at
completed protected-spin gate boundaries, with all coins and complement
exports retained and their full pulse cost added. Partial pulses are not
claimed repair invariant.

## 6. A certified finite example

Use the complete unchanged 256-vector FA98 witness for the open eight-site
cube, J=1, B=1/4, and the spline `(a,b,d)=(1,2,1/4)`, with normalization 5.
The exact polar certificate gives `||H-QDQ*||<=r`, `||Q-V||<=d0`, and a
ground-vector error e. Writing delta_E=r+ground_origin_error, the Fourier
first moment M1 in (4) gives an operator error `M1 delta_E` for G. Thus
the filtered-weight error is at most
`2M1 delta_E+d0(2+d0)+2e`. The energy-numerator error is at most
`(1+56M1)delta_E+28d0(2+d0)+56e`, using the exact excitation norm bound 28.
Every overlap is an integer/rational contraction; no degenerate eigenvector
is promoted to a distinguished state. The results include:

| Probe | Successful LCU probability | Conditional excitation energy |
| --- | --- | --- |
| local Z0 | [0.01636954798, 0.01636954802] | [1.1650825184, 1.1650825819] |
| uniform sum Z/8 | [0.00075962729, 0.00075962734] | [1.0922839962, 1.0922853592] |
| staggered sum Z/8 | [0.01450080916, 0.01450080921] | [1.0922846419, 1.0922847134] |

The receipt includes complementary rejection probabilities and a finite
quadrature/control plan giving conditional trace error below 1/100, provided
the finite ground input has trace error at most 10^-8. FA98's terminating
preparation supplies that finite input accuracy. This is a certified target
calculation and resource construction; the new large native tape has not
been executed. These are selected finite-source excitation energies, not
relativistic rest masses or measured particle assignments.

## Verification and remaining boundary

Run `python -W error -m pytest -q model/spin_packet_readout`, then
`python -W error -m model.spin_packet_readout --write` and `--check`.
The current runtime import closure, original integer witness and named proof
inputs are bound. Complete historical receipt chains are not falsely claimed
replayed. Both Boolean/numeric type substitutions and scope changes are rejected.

The new p^4 success and error budgets address observer access to the proved
packets. They do not supply a computable volume/field-selection rate for the
particular infinite ordered state, suppress its multiparticle continuum,
identify FA92's variational packet with a FA100 Stone vector, or derive a
sharp mass shell, Lorentz symmetry, physical species, or Einstein dynamics.
