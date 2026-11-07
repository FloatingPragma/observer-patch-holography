# Born representation from a countable effect domain

The [native record construction](../native_born_bridge/PROOF.md) identifies a
specific gap: protocols with finite discrete descriptions form a countable menu, whereas
the finite effect-valuation theorem quantifies over every complex matrix
effect. The theorem below reduces that mathematical quantifier to an exact
countable domain. It does not establish that the native source realizes the
domain or assigns equal probabilities to equal effects across implementations.

## Exact theorem

Fix a positive finite dimension `d`. Let `E_Q` be the Hermitian matrices over
the Gaussian rationals `Q(i)` satisfying `0 <= A <= I`, with the usual complex
positive-semidefinite order. Suppose a function `v:E_Q -> R` satisfies:

1. `v(A) >= 0` for every `A` in `E_Q`.
2. `v(I) = 1`.
3. `v(A+B) = v(A)+v(B)` whenever `A`, `B` and `A+B` belong to `E_Q`.

Then `v` is 1-Lipschitz in operator norm. It has a unique continuous
extension to all effects in `M_d(C)`. This extension is an additive effect
valuation, and there is a unique density matrix `rho` such that

```math
v(A)=\operatorname{Tr}(\rho A)\quad(A\in E_{\mathbb Q}),\qquad
\rho\succeq0,\quad\operatorname{Tr}\rho=1.
```

Conversely, every density matrix defines such a valuation on `E_Q`. Neither
continuity nor a trace probability rule is assumed for `v`. The complex
matrix representation, its positive order and the three all-domain conditions
above are hypotheses. The proof is mathematical prose; the referenced full
finite-effect theorem is machine checked, while this dense-domain extension
has no Lean implementation in this package.

## Positivity supplies the continuity bound

Additivity at `0+0` gives `v(0)=0`. Applying it to `A+(I-A)` gives
`0 <= v(A) <= 1`. If `A <= B`, then `B-A` is a rational effect and
`v(B)-v(A)=v(B-A)>=0`, so `v` is monotone.

For every positive integer `n`, all partial sums of `A/n` are effects;
therefore `v(A/n)=v(A)/n`. If `r=m/n` is nonnegative rational and `rA`
is an effect, its partial sums are also effects, giving `v(rA)=r v(A)`.
This includes `r>1` whenever `rA<=I`.

Take `A,B` in `E_Q` and rational `delta>=||A-B||`. The rational effects

```math
\frac{A}{1+\delta}\preceq\frac{B+\delta I}{1+\delta}\preceq I
```

give, by monotonicity, compatible additivity and rational homogeneity,

```math
\frac{v(A)}{1+\delta}
\le\frac{v(B)+\delta}{1+\delta}.
```

Interchanging `A` and `B` yields `|v(A)-v(B)|<=delta`. Taking rational
`delta` down to `||A-B||` proves the claimed 1-Lipschitz bound. This argument
uses order comparisons, not a commutativity assumption about `A` and `B`.
Its constant is independent of `d`; the number and cost of programs needed
to realize an approximation are separate questions.

## Density and compatible approximation

Gaussian-rational Hermitian matrices are dense among complex Hermitian
matrices in operator norm: approximate the finitely many real diagonal and
real/imaginary off-diagonal coordinates by rationals while preserving
Hermitian symmetry. Entrywise convergence implies operator-norm convergence
in fixed finite dimension.

Let `A` be any effect. Choose a Gaussian-rational Hermitian `B` with
`||B-A||<=delta`, for positive rational `delta`. Then
`-delta I <= B <= (1+delta)I`, and

```math
A_\delta=\frac{B+\delta I}{1+2\delta}\in E_{\mathbb Q},\qquad
\|A_\delta-A\|
\le\frac{2\delta}{1+2\delta}.
```

The error bound follows from `||I-2A||<=1`. Thus `E_Q` is dense even at
the singular boundary of the effect cone. Its countability follows from
the finite number of rational coordinates.

Lipschitzness makes `v(A_n)` Cauchy for every sequence of rational effects
approaching `A`, and makes the limit independent of the approximating
sequence. Define that limit to be `v_bar(A)`. The extension is 1-Lipschitz,
nonnegative, normalized, and unique among continuous extensions.

It remains to establish compatible additivity; approximating each effect
independently without controlling their sum would be insufficient. Let
`A,B,A+B` be effects. Choose Gaussian-rational Hermitian `R_delta,S_delta`
within `delta` of `A,B`, respectively. Set

```math
\widehat A_\delta=\frac{R_\delta+\delta I}{1+4\delta},\qquad
\widehat B_\delta=\frac{S_\delta+\delta I}{1+4\delta}.
```

Both numerators are positive. Their sum is bounded above by
`(1+4delta)I`, because `A+B<=I`. Hence both matrices and their sum belong
to `E_Q`. They converge respectively to `A,B,A+B`. Exact additivity on
`E_Q` and continuity therefore give
`v_bar(A+B)=v_bar(A)+v_bar(B)`.

The extension satisfies precisely the hypotheses of
[`EventAlgebra.finite_busch_gleason`](../../Lean/EventAlgebra/FiniteBuschGleason.lean).
That theorem supplies the unique density matrix. Conversely, positivity
and linearity of the trace pairing give all three conditions for every
density matrix. If two density matrices agreed only on `E_Q`, their
continuous trace pairings would agree on all effects by density, and the
same finite theorem gives equality. This completes both directions.

## What the executable checks

[`dense_effects.py`](dense_effects.py) checks exact rational instances of the
order sandwiches, single-effect and compatible-pair regularizations, and
complex tomographic coordinates. Its positive-semidefinite test uses exact
Schur elimination, retaining null directions. The finite fixtures support
the proof's algebra; they do not verify a quantifier over every effect or
every source program.

A retained counterexample makes the finite-data limitation explicit. On
qubit effects define

```math
w(E)=f(\operatorname{Tr}(E)/2),\qquad
f(t)=3t-6t^2+4t^3=\frac{1+(2t-1)^3}{2}.
```

Since `f'(t)=3(2t-1)^2`, this assignment lies in `[0,1]` for every effect,
is normalized, and obeys `w(I-E)=1-w(E)`. It agrees with `rho=I/2` on
every rank-one projector, hence on any projective tomographic menu. But
for `E=F=I/8`,

```math
w(E+F)=7/16,\qquad w(E)+w(F)=37/64,
```

with defect `-9/64`. A density reconstructed from finitely many projective
probes therefore cannot certify additivity of the source's assignment on
the countable domain. The repository's stronger, machine-checked
[committed-fixture countermodel](../../Lean/EventAlgebra/OperationalAdditivityBoundary.lean)
has the same logical role and remains applicable.

## Remaining source obligations

To apply this theorem to native observer records, an implementation must
supply the complex effect representation and an exact realization of the
specified countable menu, or a separately controlled approximation theorem.
It must prove that two programs representing the same effect have the same
record probability, and that every compatible sum is represented by a
disjoint union in an admitted complete instrument. All failure, timeout and
complement outcomes belong to that instrument.

The proof derives a continuous mathematical extension from exact
all-domain conditions. It does not infer those conditions from a finite
number of successful measurements or assume that a gate set's density
alone establishes implementation equivalence. Approximate gate synthesis
requires explicit error bounds and a compatible limiting record law.
Finite observer capacity also requires its own uniform error/cost analysis.
Fixed bounds on program description, memory and execution admit only a finite
part of the menu. Applying the countable theorem through increasing resource
ceilings requires compatible valuations across those stages; a single finite
calibration is insufficient.

This result reduces the effect-coverage requirement to a countable domain.
It selects no dimension, source state, physical vacuum, spacetime clock,
particle spectrum or relativistic field theory. No canonical physical
premise is discharged by the finite fixtures or by this conditional theorem.
