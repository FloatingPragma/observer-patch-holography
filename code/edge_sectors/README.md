# Finite gauge-sector readouts

This package evaluates the **supplied** finite Hamiltonians in the Standard
Model gauge paper. It selects neither a physical Hamiltonian nor an OPH
coupling. It distinguishes the positive sector probabilities, a fitted
heat-kernel parameter, and the prediction for a sector excluded from that fit.

```sh
python code/edge_sectors/heat_kernel_holdout_validation.py
python -W error -m pytest -q code/edge_sectors
```

## Reproduced defects

At main `b8714c98`, the original-input controls produced six scientific
failures and six passes before the implementation changed:

| Input | Old result | Independent result |
| --- | --- | --- |
| S3, `h=1e12`, sign sector | `1.204775757793917e-53` | `1.2056327160495836e-53` |
| S3, `h=1e15`, sign sector | zero; the printed comparison then crashed | `1.205632716049383e-65` |
| Z2, `h=1e6`, nontrivial sector | `3.125277814319816e-14` | `3.125000000000879e-14` |
| Z2, `h=1e7` and `1e8` | roundoff-scale Fourier differences | approximately `3.125e-16` and `3.125e-18` |

The S3 eigensolver controlled absolute vector error, which did not resolve
the tiny sign component. The Z2 Fourier transform subtracted overlaps near
one and clipped negative results. Neither operation preserved positive rare
probabilities. The live S3 paper table also disagreed with its declared
Hamiltonian at ordinary couplings: at `h=.5` it printed
`(.909,.0013,.089)` instead of approximately
`(.962927,.00024760,.036825)`. The corrected table is checked against the
original three-by-three matrix, independently of the producer.

## Exact abelian reduction

For `n>=2`, `h>0`, `K=1`, `Gamma=5`, the eight-link configuration Hamiltonian
has nonpositive off-diagonal entries and connected single-link hopping.
Its ground vector is unique and strictly positive. Each Gauss-star and
noncontractible electric-cut symmetry is a permutation commuting with the
Hamiltonian, so it fixes that vector: a positive vector cannot acquire a
nontrivial unit-modulus phase. In the electric Fourier basis the ground state
therefore has zero divergence and zero winding.

Such electric fields are precisely plaquette boundaries over Z_n, including
composite n. Integrating the divergence-free field on the dual torus gives
plaquette potentials; zero winding makes the integration path independent.
The potentials are unique modulo a common constant. Fixing `c00=0` and writing
`(a,b,c)=(c01,c10,c11)` leaves `n^3` states, compared with `n^8` link states.
In link order `(00x,00y,01x,01y,10x,10y,11x,11y)` their fields are

```
(-a, b, a, c-a, b-c, -b, c-b, a-c) mod n.
```

The specified restricted star measures `q=b-c mod n`. Each magnetic plaquette
shifts one potential, or all three when the fixed potential is shifted and
the representative is restored. Opposite shifts coincide at n=2 and **add**.
Removing the common `-20-8h` energy leaves diagonal

```
4h [sin²(pi a/n) + sin²(pi b/n) + sin²(pi(b-c)/n) + sin²(pi(c-a)/n)]
```

and hopping `-1/2` for each of the eight translations. The readout is a sum
of positive squared amplitudes in each charge sector, not a difference of
nearly equal expectation values. Z5 uses 125 states instead of 390,625.

In the large-h branch, write the shifted matrix as
`[[0,-b.T],[-b,A]]`. When A is strictly diagonally dominant it is a positive
M-matrix. Taking the origin amplitude to be one gives
`z=(A+sI)^-1 b` and `s=b.T z` for ground energy `-s`. The right side decreases
in s; the root is bracketed by zero and `b.T A^-1 b`. Solving for the root
relative to this upper bound avoids an absolute root-finder floor. These
positive tail solves resolve the retained Z2 examples through `h=1e154`.
The other branch uses a sparse Hermitian eigensolver.

This algebraic reduction is exact; its floating-point implementation is not
an interval certificate. The coupling must be exactly representable in
binary64. Both branches require positive amplitudes and a componentwise
eigen-equation residual at most `1e-8`, measured against positive incoming
amplitudes. The original returned amplitudes are squared and summed by charge
over exact rationals, normalized exactly, then converted once to binary64.
This preserves a representable sector even when individual squared terms
would underflow. Conversion must preserve every positive sector to relative
`1e-12`; otherwise the readout raises. Unresolved tails and nonfinite data
also raise instead of yielding zero sectors. These diagnostics are not a
proved relative population-error bound. `h=0` is refused because the full abelian
model then has degenerate topological ground sectors; the code does not
choose a ground state silently. Normalized printed comparisons additionally
refuse a log-gap of magnitude at most `1e-8`, where the nearly uniform
numerical weights do not justify that normalization. Even above that cutoff,
the held-out numerator can be a roundoff-dominated subtraction: at
`h=5.708163265306122e-9`, the former report printed a positive residual while
an independent 65-digit original-incidence calculation gave approximately
`-1.09016217468756855e-9`. Reporting now also requires the difference of the
fitted and measured log gaps to exceed
`64 eps (1 + abs(fitted_gap) + abs(measured_gap))` in magnitude, where eps is
binary64 machine epsilon. This numerical cancellation guard refuses both
unresolved signs and apparent zero mismatches. It is not an interval bound
for the upstream solve; the probability API remains available independently.
Resolved small residuals use scientific notation instead of printing zero
percent at the ordinary display precision.

The generic fitted-time and held-out prediction helpers also evaluate the
original binary64 inputs in private 100-digit arithmetic and require final
conversion to relative `1e-12`. A nonzero subnormal is not automatically a
resolved result: the audit reproduced a 6.2% rounded Z5 held-out prediction
at `h=1e61`, and an 11% rounded fitted time for nearly equal weights divided
by `1e308`. These comparisons now refuse insufficient output precision;
the Z5 sector probabilities remain separately available. Resolved subnormal
predictions, including the Z5 report at `h=1e59`, still work. This conversion
check does not certify upstream eigensolver accuracy or logarithms by intervals.

## Complete S3 comparison within the declared model

In character order `(triv,sign,std)`, with four electric links,

```
H(h) = [[0,0,-1], [0,24h,-1], [-1,-1,12h-1]], h>=0.
```

The positive ground vector is proportional to `(1,x/(24h+x),x)`, where its
energy is `-x` and

```
F(x,h) = x+12h-1 - 1/x - 1/(24h+x) = 0.
```

F is strictly increasing in x, giving a unique positive root. Exact rational
bisection starts in a bracket with endpoint ratio at most six, so tiny roots
need no exponent-dependent search. Supplied integers, Fractions, Decimals
and binary real scalars retain their original values. The returned bracket
has rational endpoint signs; probabilities and logarithms use a private
100-digit context after 256 bisections. This does not make the logarithms
certified intervals. Every nonzero binary64 output must survive conversion
with relative error at most `1e-12`, or the evaluator refuses that readout.
The probability-only API remains available when normalizing a fitted-time
diagnostic is unresolved.

Let `z=x+12h` and `y=x(x+24h)`. The secular equation gives
`(z-1)y=2z`. Under the stated `d_R exp(-t lambda_R)` convention, fitting the
standard sector predicts `p_sign_pred=p_std²/(4p_triv)`, hence

```
p_sign / p_sign_pred = ((z-1)/z)² < 1,
log(p_sign / p_sign_pred) = -2 log1p(1/(z-1)).
```

The latter identity preserves small discrepancies even when rounded
probabilities coincide. From `z²-2z/(z-1)=144h²`, z increases from two to
infinity. Thus the measured/predicted ratio increases strictly from 1/4
toward one: **no finite coupling in this model gives the exact fitted law**.
The standard-sector diffusion time `log(2/x²)/3` is positive precisely above
`h*=(1+sqrt(3)-sqrt(2))/24`. It vanishes at h*, where the normalized log
residual and log-ratio are undefined; below h* the fit is outside positive
diffusion. On the positive-time branch, the absolute raw discrepancy
decreases while the fitted time increases, proving the decay of the
normalized residual and log-ratio excess. These are complete statements
about this fixed model, not refinement or continuous-group results.

## Independent validation and impact

The controls include full original-link Z2/Z3 Hamiltonians; separately
enumerated electric-link Z5 matrices; composite Z4 and tiny/large coupling
checks; and 350-digit diagonalization of the original S3 matrix. They check
each probability and comparison component separately, the entire printed
S3 table, real command-line output, malformed data, genuine zero/negative
fitted-time boundaries and explicit precision refusals. Isolated bad
implementations must fail scientific assertions. The finite-gauge workflow
runs the entire package on Linux; its tests reject omitted,
collection-only, filtered or ignored-failure dispatch.

Ordinary Z2/Z3/Z5 values retain their declared model and agree with independent
replay. The corrected S3 numbers and fixed-model result replace the erroneous
live table and finite-trend wording. The corresponding live claim and
novelty/falsification rows are synchronized; their assumptions, dependencies,
gates and `declared_structure` class are unchanged. The gauge paper PDF and
release manifest are rebuilt. This package has no frozen scientific receipt;
no frozen target, custody record, SU3 table, book claim, physical input or
absolute-coupling prediction is created or refreshed. Work belongs to the
standing [evidence audit #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
