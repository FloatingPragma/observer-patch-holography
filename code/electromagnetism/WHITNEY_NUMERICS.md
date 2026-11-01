# Whitney numerical methods

These methods evaluate the supplied finite Whitney action and neutral quantum
packet. They preserve positive kinetic reduction, projected phases and
normalization without changing the physical action or preparation. The
analytic identities below are distinct from the numerical resolution policies;
neither establishes a physical clock, a source-selected action or a quantum
propagation certificate.

The kinetic implementation is in [whitney_interacting_quantum.py](whitney_interacting_quantum.py);
the state and packet interfaces are in
[whitney_quantum_state.py](whitney_quantum_state.py) and
[whitney_quantum_packet.py](whitney_quantum_packet.py).

## Kinetic positivity controls

Set all edge coefficients to zero, all 13 complex matter coefficients to the real value
`T=1e9`, charge to `1/4`, and the scalar imaginary velocity to
`x=(12,-1,...,-1)`. The remaining velocity components vanish. The original
source action reduces independently to the one-variable minimization

\[
 \min_z \{Az^2+B(1+eTz)^2\}
 =\frac{AB}{A+B e^2T^2},\qquad
 A=169(30-10\sqrt5),\quad B=\frac{114}{5}(10+10\sqrt5/3).
\]

The exact positive answer is approximately `2.0656721888405686e-14`.
Unprotected dense Schur subtraction can return `-1.6542323066914832e-13`; a positive
determinant sign does not establish a positive metric or a valid
momentum/velocity contraction. Even at `T=1e6`, where density evaluation is
practical, the subtractive Gaussian log amplitude differs by about `4.6e-6`.
Repeating the same unstable subtraction in a second implementation would
not provide independent validation.

[The original-input controls](test_whitney_kinetic_precision.py)
derive the scalar modes and absolute determinant from the original cone,
simplex moments and icosahedral adjacency spectrum. They do not use producer
quadrature or a rounded reduced matrix as the expected answer.

## Kinetic reduction and domain

Write the positive quadrature action as `||H v + V eta||²`, where `H` contains
the 56 slice columns and `V` the 12 mean-zero gauge columns. Normalize each
column of `[V,H]`, apply Householder QR, and undo the column scales. The
trailing triangular block is a factor `F` of the minimized action; the leading
block gives its minimizing gauge velocity. This replaces subtraction of
nearly equal Gram matrices. Gaussian density uses `2 sum(log(abs(diag(F))))`
directly. A density-only calculation omits the potential, so an unused
quartic overflow does not reject a resolved kinetic density.

The scalar edge derivative pairs its endpoint values before multiplying
by barycentric weights. Uniform matter at zero edge field therefore has
exactly zero edge derivative, even at large amplitude. The vertical scalar
factor uses the Ward identity `i e Psi lambda` directly.

The spatial gradient requires the same constant-mode protection. At
`a=0`, uniform neutral matter `T=1e80` and zero potential couplings,
unprotected differentiation can give energy `9.640784334323916e128` instead
of zero. Subtracting the transported constant mode before differentiation
preserves that zero.
A nearby field `psi=(1e15+1,1e15,...,1e15)` independently checks the nonzero
energy `30-10 sqrt(5)`; the unprotected calculation gives `9.088123195321897`.
Zero mass/quartic terms are omitted before evaluating powers, and genuinely
required energies outside binary64 range are explicitly refused.

The implementation checks both the scale of the horizontal factor before
projection and the conditioning of the normalized vertical columns. Its
binary64 resolution indicator is

\[
 68\epsilon\,\frac{\|H\|_2}{s_{\min}(F)}
 \kappa_2(V_{\rm normalized}).
\]

Returning a dense metric additionally requires
`56 epsilon ||abs(F).T abs(F)||₂ / s_min(F)² <= 1e-7` and successful Cholesky
factorization. Both indicators must be at most `1e-7`. These are numerical
resolution policies, **not certified forward-error or quadrature bounds**.
There is no eigenvalue floor or projection onto positive matrices.

The distinction matters in actual controls: the uniform `T=1e6` factor and
Gaussian density are evaluated, while the unresolved dense `T=1e9` metric and
cotangent are refused. With matter `(1e12,1,...,1)`, testing only `cond(F)`
would accept a projection already wrong by about `2e-5`; the combined test
refuses it. The nearby `(1e6,1,...,1)` configuration is retained and checked
against original degree-four moments.

Binary64 entry points validate each original scalar before array coercion.
Booleans, strings, missing values and nonrepresentable exact inputs are
refused. Exactly representable integers, fractions, decimals and subnormal
configuration values remain valid inputs. This is a numerical evaluator's
contract, not a restriction on the analytic configuration space.

## Independent kinetic replay

The two receipt verifiers share one original-input monomial assembly.
They lift supplied scalars and vertices into a private mpmath context before
products, the gauge rechart or constrained reduction. Their Cholesky-checked
moment calculation is independent of the producer's quadrature/QR algorithm.
Only final reported values are converted to binary64; unrepresentable reports
are refused. The context includes 80 guard digits beyond a scale-based working
budget. This remains a finite numerical replay, not an interval certificate.

The controls exercise each action component, all 13 scalar modes, inertia,
minimizer, absolute density and public cotangent. They include real/imaginary
uniform fields, nonuniform complex fields, signed and zero charge, changes of
coordinate frame, concentrated fields, and large neutral or scaled-charge
fields whose unused potential overflows. The independent replay also rejects
indefinite reduced metrics with positive determinant and checks nonuniform
dressing monomials separately. Representation tests challenge mixed numeric
containers before conversion. The Windows/Linux workflow runs these controls
with warnings treated as errors.

A separate control challenges a tiny physical tangent hidden behind a
large pure-gauge tangent. For uniform matter one and gauge vector
`xi=(0,1,-1,0,...,0)`, supply edge velocity `D xi`, imaginary scalar velocity
`xi/4`, and an additional center velocity `1e-100`. The exact recharted
velocity is just that last entry. Unprotected binary64 gauge solving leaves
`O(1e-16)` noise; a working-precision budget based only on the largest input
also loses it. Public binary64 packet evaluation refuses that unresolved
cancellation, while direct tiny velocities through `1e-200` remain valid.
This additional normwise gate includes the gauge solve's condition number;
it is not a componentwise error certificate. Independent replay budgets from
the full original exponent span and numerator/denominator bit lengths, before
any products, and retains the `1e-100` signal. It also preserves near-one
Fraction/Decimal differences. Exact-zero comparisons use an absolute
roundoff bound rather than declaring numerical solve residuals exactly zero.

## Scale-sensitive kinetic controls

| Guarded failure mode | Required computation and control |
| --- | --- |
| For `q[30]=2e154`, width one, squaring too early gives Gaussian log amplitude `-inf` instead of about `-1e308`. | Divide coordinates by two before squaring; refuse genuinely unreportable logs. Preserve tiny widths and large balanced coordinate/width units. |
| Norm overflow can let the nonzero gradient `a=2^600 D xi` pass the Coulomb condition. | Compare the slice residual and norm in common, bounded units before squaring. The same non-Coulomb direction must refuse at ordinary and large scales. |
| Unused covariant products can overflow for a resolved neutral density. | Omit unused spatial evaluation and identically zero charge products. Large opposite scalar components and nonzero radiative coordinates remain valid for neutral density evaluation. |
| Small couplings times large field powers can overflow, and large couplings times small powers can underflow, despite a representable potential. | Accumulate positive quadrature products using separate mantissas and binary exponents. Test mass and quartic terms, real and imaginary matter, and actual reporting-range failures. |
| A damaged stiffness quadratic assigns magnetic energy `1.2158875978002469e13` to a pure gradient of amplitude `1e15`. Adding one to a boundary edge gives `1.99658875978003e14` instead of `(5-sqrt(5))/6`. | Compute original face circulations with exact rational addition, then evaluate their positive curl factors. No floor is applied to a damaged stiffness quadratic. |
| A resolution comparison can overflow for a representable large cotangent; an unreportable cotangent can become infinity/NaN. | Compare normalized scales, use real-coordinate norms, and validate every required returned quantity. Preserve large finite complex components and resolved subnormal tangents. |
| Treating a boundary edge `1e-12` as zero omits magnetic energy `4.606553370833684e-25`. | Require exact original edge-cycle closure before using the zero-field monomial model. The verifier accepts small genuine gradients and arbitrary supported tangents. |
| A numerical gauge-solve residual can underflow on a large exact pure-gradient configuration. | Construct its mean-zero gauge potential from the proved original gradient identity and derive its zero Coulomb field algebraically. This is not numerical clipping. |

Independent cross-checks use original complex simplex moments and original
curl-pair integrals, not the producer's potential accumulator. The replay
tests separately challenge domain refusal, fixed working precision,
determinant-only positivity and missing reporting checks. These remain finite
diagnostics; the exact analytic and physical scope above does not change.

## Neutral projection controls

All examples below have width and hbar one. Indices are zero based: the
first 30 coordinates are radiative; the final 26 are 13 scalar real parts
followed by 13 imaginary parts.

1. Set the center to zero, point coordinates `(0,4,8)` to `(1,1,1)` and
   momentum coordinates `(0,4,8)` to `(2**53,1,-2**53)`. Other entries vanish.
   Every input is exactly representable in binary64. An unprotected dot
   product loses the exact phase one and gives
   `3.1605377469223292e-12` instead of
   `1.7076458324454292e-12 + 2.659500810425262e-12j`, changing interference.
2. Set the center to zero, scalar point coordinate 30 to one and scalar
   momentum coordinate 30 to 256. A 256-angle quadrature gives
   `+1.3759422606592422e-11`; the original integral gives
   `-4.838248674915097e-12`. Rotating the point by 0.373 radians also changes the
   quadrature result. The nearby momentum-10 case is a resolved control.

[The projection controls](test_neutral_packet_projection.py) distinguish
phase loss and quadrature aliasing from correctly resolved ordinary inputs.
Numerical quadrature alone supplies no certified error bound.

## Exact reduction from the supplied seed

Write a point as `(x_r,y)`, its center as `(c_r,c)` and momentum as `(p_r,p)`.
Here `r` denotes the 30 radiative coordinates, and `y,c,p` have 26 scalar
coordinates. Let `J(u,v)=(-v,u)` and `R_theta=cos(theta) I+sin(theta) J`.
For positive width `s` and hbar `h`, the seed Lebesgue half-density is

\[
 g(x)=(2\pi s^2)^{-14}
 \exp\left[-\frac{|x-c_{\rm full}|^2}{4s^2}
       +\frac{i}{h}p_{\rm full}\mathbin\cdot(x-c_{\rm full})\right].
\]

Rotating scalar center and momentum together gives
`g_theta(x)=(2 pi s^2)^(-14) exp(C+U cos(theta)+V sin(theta))`, with

\[
\begin{aligned}
C&=-\frac{|x_r-c_r|^2+|y|^2+|c|^2}{4s^2}
   +\frac{i}{h}\bigl[p_r\cdot(x_r-c_r)-p\cdot c\bigr],\\
U&=\frac{y\cdot c}{2s^2}+\frac{i}{h}p\cdot y,\\
V&=\frac{y\cdot Jc}{2s^2}+\frac{i}{h}(Jp)\cdot y.
\end{aligned}
\]

In particular the constant scalar phase is `-p dot c/h`. The imaginary
coefficient in `V` is `(Jp) dot y/h`, with the displayed sign.
Expand

\[
 U\cos\theta+V\sin\theta
 =\tfrac12(U-iV)e^{i\theta}+\tfrac12(U+iV)e^{-i\theta}.
\]

The exponential series converges absolutely and uniformly on the compact
circle. Only equal powers of the two Fourier factors survive its average.
Consequently, for all such original real data, including nonzero `B` below
and complex coefficients `U,V`,

\[
 P_0g(x)=(2\pi s^2)^{-14}e^C
   \sum_{k=0}^\infty\frac{[(U^2+V^2)/4]^k}{(k!)^2}
 =(2\pi s^2)^{-14}e^C I_0(\sqrt{U^2+V^2}).
\]

This is an entire function of `U^2+V^2`; there is no square-root branch
choice in the resulting value. The last equality is the defining Bessel
series ([DLMF 10.25](https://dlmf.nist.gov/10.25)); the usual single-coefficient
integral is recorded in [DLMF 10.32](https://dlmf.nist.gov/10.32).
This derivation does not assume successful quadrature or a restricted
momentum orientation.

The [existing norm theorem](../../paper/tex_fragments/WHITNEY_QUANTUM_PACKET.tex)
gives

\[
 A=\frac{|c|^2}{4s^2}+\frac{s^2|p|^2}{h^2},\quad
 B=\frac{p\cdot Jc}{h},\quad Z=\sqrt{A^2-B^2},\quad
 \|P_0g\|^2=e^{-A}I_0(Z)>0.
\]

Cauchy--Schwarz and the arithmetic--geometric mean inequality give
`A >= |B|`. The normalized neutral half-density therefore is

\[
 f_0(x)=(2\pi s^2)^{-14}
        e^{C+A/2}\frac{I_0(\sqrt{U^2+V^2})}{\sqrt{I_0(Z)}}.
\]

This evaluates exactly the state defined in the paper. Neutrality follows
also directly: a scalar rotation of `y` rotates the pair `(U,V)`, preserving
`U^2+V^2` and `|y|^2`. The physical configuration wavefunction still requires
the existing metric factor `rho^(-1/2)`; this function returns the Lebesgue
half-density, not that wavefunction or a propagated state.

## A normalization consequence

For `p=h Jc/(2s^2)`, one has `A=B`, `Z=0` and `U^2+V^2=0`.
Also `p dot c=0`, and the scalar part of `C+A/2` is exactly
`-|y|^2/(4s^2)`. The normalized scalar projection is thus the centered
Gaussian for every center size. This proves why rejecting every state whose
separately reported projection norm underflows would be incorrect. For
`c[0]=1e100`, `p[13]=.5e100`, `s=h=1`, the amplitude at the origin is the
ordinary value `(2 pi)^(-14)`, although the squared projection norm cannot
be returned reliably as binary64. Its standalone reporting API must still
refuse that norm. The normalized-state API evaluates the combined expression.

## Numerical contract and independent controls

Pointwise entry points retain the original finite, exactly binary64-
representable real-scalar input contract. They reject malformed vectors,
Booleans, masked values, nonfinite scalars and inexact coercions before
array conversion. Width and hbar must be positive. Displacements, dot
products, `U^2+V^2`, `A^2-B^2` and `Re(C)+A/2` are computed from the original
inputs as exact fractions. No numerical zero floor is used.

High-precision arithmetic evaluates the final normalized expression
before conversion. Its precision budget covers large absolute phases and
cancelling logarithms, with successive higher-precision evaluations required
to agree relatively to `1e-30`. The final complex binary64 result must agree
with that calculation relatively to `1e-12`. This is a normwise numerical
stability policy, **not an interval certificate or proof of an exact node**.
Unresolved zeros and unreportable outputs raise an error. The public seed
logarithm additionally requires its returned unwrapped phase to have at most
`1e-12` radians absolute conversion error; a small relative error in a huge
phase alone would not suffice. The projection evaluates its phase internally
and can remain evaluable when that unwrapped binary64 report is unavailable.

The legacy `nodes` argument remains validated as an integer at least 16, but
has no effect on the closed-form evaluation. Raising the angle count is not an accuracy intervention.

Controls compare against integration of the original rotated seed and its
Gaussian overlap, and against independently simplified cases. The general
overlap integrand is complex when `B` is nonzero; its average is the positive
squared projection norm. The controls include
complex scalar centers and momenta, nonunit width and hbar, interference,
near Bessel zeros, rotation, coordinate permutations, normalization
compensation, caller precision contexts and valid/reporting-range pairs.
The controls run with warnings treated as errors on Linux.
Their numerical integrations are independent diagnostics, not certified
interval enclosures.

## Scalar norm and radius evaluator

The scalar evaluator preserves both the projection probability and the scalar
nodal squared radius. With zero center, scalar momentum `1e10`, and width
and hbar one, the radius is approximately 25; subtractive evaluation can
instead erase it. With scalar center `(1,0)` and momentum `(0,1e8)` in one
complex plane, the norm retains the finite `exp(-1/2)` suppression.

For scalar parts `x,p`, write `X=|x|²/(4 sigma²)`,
`Y=sigma² |p|²/hbar²`, `A=X+Y`, `B=p·Jx/hbar`, `D=X-Y`,
and `K=4XY-B² >= 0`. The producer forms these invariants over exact
rationals from the supplied real scalars, before array coercion. Then
`z=sqrt(D²+K)`, `A-z=B²/(A+z)`, and, for negative `D`,
`z+D=K/(z-D)`. The zero case is handled separately. The equivalent radius
is `2 sigma² (13+(z+D)-z(1-I1(z)/I0(z)))`. These identities remove the
large subtractions; a private mpmath context budgets extra precision for
the remaining Bessel ratio. The two public observables share this evaluator.
Masked, boolean, complex and nonfinite inputs are rejected. Integer,
Fraction and NumPy real scalars retain their original values.
NumPy integers and integer components inside a Fraction are converted
to unbounded Python integers before rational arithmetic; wrapping a NumPy
integer in Fraction alone can retain fixed-width overflow.

Every returned nonzero scalar must retain relative `1e-12` accuracy when
converted to binary64, or that readout explicitly refuses the range. This
is a reporting criterion, not a certified error enclosure. An unreportable
probability does not prevent a separately representable normalized radius.
The pointwise projection evaluates the closed form of the circle integral,
a modified Bessel function of the complex cosine and sine coefficients,
derived [above](#exact-reduction-from-the-supplied-seed).
No angle count controls its accuracy: the legacy `nodes` argument is
validated and has no effect. The normalized amplitude is returned even when
the separately reported projection norm is unrepresentable.
Its binary64 seed, rotation and phase-space callers validate the original
entries too. They reject masked/boolean data and scalars whose conversion
would change their value, including an integer displacement erased above
`2^53`. This explicit caller limitation does not restrict the scalar
norm/radius evaluator's exact-rational input arithmetic.
The seed and the projection form displacements, dot products and the
combined exponent as exact fractions of those entries. A private mpmath
context evaluates the normalized expression; successive precisions must
agree to relative `1e-30`, the binary64 result to relative `1e-12`, and a
returned seed phase to `1e-12` radians. Unresolved zeros and unreportable
amplitudes raise an error instead of returning NaN, zero or an underflowed
Gaussian sample. This is a numerical stability policy, not an interval
certificate.

The independent receipt replay imports neither the packet producer nor the
interacting coefficient evaluator. It uses positive one-dimensional integrals
instead of Bessel functions, with a rescaled
Gaussian integration coordinate at large `z`. Its scalar comparisons have
no absolute tolerance that could accept erasing a small positive value.
The phase-space replay first authenticates the reported numerical coordinates;
the scalar integrals then use those coordinates, so roundoff in a second
gauge solve cannot hide a wrong small charge. The integral representation
and large-argument behavior are consistent with [DLMF 10.32](https://dlmf.nist.gov/10.32)
and [10.40](https://dlmf.nist.gov/10.40).

## Additional phase and range controls

Rounding a displacement before the phase dot product is also unsafe. Set two
radiative point coordinates to one, the first center coordinate to `2^-54`,
and their momenta to `(2^54,-2^54)`. The exact phase is `-1`, while subtracting
the first center in binary64 loses that phase. The control computes
the original displacement and expected Gaussian independently. Both seed and
projected evaluations must preserve that phase.

Additional retained controls use independently summed Fourier coefficients
for a nearly chiral scalar seed with `c=(C,0)`, `p=(+/-1/C,C/2)` and
`C=2^500` or `2^1023`. The surviving complex circle coefficients remain of
ordinary size while the separate center and momentum terms are enormous.
A further phase control combines the largest finite binary64 momentum with
the smallest positive binary64 hbar. The normalized projection remains
representable even though its unwrapped phase cannot be returned as a
binary64 value. These are independent numerical controls, not additions to
the analytic hypotheses or physical claims.

## Reproduction and evidence scope

From the repository root, with the pinned requirements:

```bash
python -m pytest -q -W error \
  code/electromagnetism/test_whitney_kinetic_precision.py \
  code/electromagnetism/test_whitney_kinetic_inputs.py \
  code/electromagnetism/test_whitney_kinetic_replay.py \
  code/electromagnetism/test_whitney_kinetic_scale.py \
  code/electromagnetism/test_whitney_public_scale.py \
  code/electromagnetism/test_whitney_quantum_packet.py \
  code/electromagnetism/test_neutral_packet_observables.py \
  code/electromagnetism/test_neutral_packet_projection.py
```

The live packet receipt is produced by
`python code/electromagnetism/whitney_quantum_packet.py` and independently
checked by `python code/electromagnetism/verify_whitney_quantum_packet.py`.
The state, trial-history, packet and real-continuum receipts have their own
canonical producers and verifiers; history pins the state receipt, so any
source-binding refresh must follow that dependency. Frozen evidence is separate.

The packet receipt stores phase-space, projection-norm and radius diagnostics;
its builder does not call the pointwise seed or projected-amplitude evaluator.
Those pointwise functions therefore require their own direct tests.
Independent controls use original simplex moments, exact rational invariants,
Fourier coefficients, scalar integrals and known special cases. They challenge
signs, phase, normalization, precision budgets and refusal boundaries rather
than treating agreement between two rounded calculations as a proof.

The dedicated workflows run the numerical controls on Linux.
Their scope is the declared finite computation. A passing numerical replay
supplies neither an interval enclosure nor a physical apparatus error budget.
