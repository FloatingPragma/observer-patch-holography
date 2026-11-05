# Complete native action of the corrected regional modular generator

## One source, state and regional algebra

Use the specified open `2 x 2 x 1` Weyl/contact action with eight
canonical anticommutation modes, `c=1/8`, `U=1` and `beta=1/32`. Its
complete field dimension is 256. The exact source preparation executes
all 256 classical column records, with every first export abort retained.
Let `Sigma` be its unnormalized success field marginal and
`sigma=Tr_out Sigma` the signed restriction to modes `(0,1,2,3)`.
This regional algebra is the full `M16`; its exterior has dimension 16.
The joint field-and-record state need not be faithful even when this
field marginal is faithful.

The [certified source control theorem](../certified_reference/README.md)
and its [finite entropy attachment](../prepared_regional_entropy/PROOF.md)
give an exact rational floor and operator bound

```
sigma >= m I,     m>1/20,
||K-P8-a I|| <= E < 2e-9,       K=-log sigma.            (1)
```

Here `P8` is the complete eighth-order regional logarithm polynomial of
the same action, and `a` is the Gibbs normalization scalar. The checker
recomputes the preparation and logarithm error bounds and verifies the
literal source certificate's full source manifest. It does not repeat
the 65,535-node preparation-angle proof.

The actual reduced logarithm has noncommuting boundary corrections.
`P8` is not replaced by `beta` times the isolated two-site Hamiltonian.
The polynomial's exact Pauli decomposition has 39 nonzero terms; all
are retained, including its scalar identity term. Its matrix equality
is independently reconstructed with exact rational complex arithmetic.
The decomposition provider lists highest-bit-first labels, whereas the
native compiler takes little-endian labels. Every label is reversed and
the full matrix reconstruction checks this convention.

`P8` commutes with regional number and parity. Every individual Pauli
term has an even number of `X,Y` factors and is even in fermionic
parity. Individual split factors need not commute with regional number.
The complete exact target has that number symmetry; no exact number
conservation is claimed for a finite product approximation or for
individual native processor pulses.

## Modular sign and scalar gauge

Use the observable convention

```
alpha_s(A)=sigma^(is) A sigma^(-is)=exp(-isK) A exp(+isK).
```

Its state unitary is `W_s=exp(+isK)`. The implemented polynomial target
is `U_s=exp(+isP8)`, with its scalar polynomial phase retained in the
source word. Anti-Hermitian Duhamel and (1) imply

```
||U_s-exp(-isa) W_s|| <= |s| E.                        (2)
```

The remaining scalar in (2) cancels from every channel and observable
conjugation. This is a gauge comparison of unitary lifts, not an equality
of the two absolute logarithms. No normalization by `Tr sigma` is needed:
multiplying sigma by a positive scalar changes only this irrelevant phase.
The source success probability and abort complement are retained.

The supplied modular parameter is dimensionless. A physical second,
an electronvolt, or a geometric diamond charge is not specified by (2).

## Native symmetric product and its finite error

Write `P8=sum_j b_j P_j` using the exact real rational Pauli coefficients.
For a step `h=s/n`, execute each `exp(+ih b_j P_j/2)` in listed order,
then execute the same factors in reverse order. Repeat this palindrome
`n` times. Each Pauli factor is the existing source compiler's literal
basis-change, parity-CNOT, `Z` rotation and uncompute word. There is no
matrix primitive or free basis rotation. The identity term incurs its
actual source pulse. Zero modular requests retain the full schedule.

Let `B=sum_j |b_j|`. Both the exact one-step exponential and the
palindromic product have identical derivatives through order two at zero.
One may see the second-order identity directly from the ordered product,
or use the palindrome's `S(-h)=S(h)^(-1)` and its first derivative
`S'(0)=iP8`. Differentiating this relation twice gives
`S''(0)=S'(0)^2`.

Every factor is unitary for real `h`. Expand its third product derivative
by the multinomial rule. The factor generators have total norm at most
`B`, so the norm of the product's third derivative is at most `B^3`
at every real parameter. The exponential has the same bound. Taylor's
integral remainder therefore gives

```
||S(h)-exp(+ihP8)|| <= B^3 |h|^3/3,
||S(s/n)^n-U_s|| <= B^3 |s|^3/(3n^2).                 (3)
```

The second line is a unitary telescoping estimate. It holds at finite
parameter, with no asymptotic truncation or eigensolver premise. Its
bound is dimension independent once the source-specific `B` is supplied.
The retained request uses `s=1/2`, `n=32`, `B<0.051`, giving an exact
rational splitting allowance below `5.4e-9`.

## Literal port controls and complete output

The source is the existing exact six-level response over `Q(sqrt(5))`.
Each emitted binary64 angle and each of its twelve integrated port
coefficients is interpreted as its exact rational value. An independent
Pauli template checks every kind, target, fixed basis angle, variable
Pauli angle and CNOT routing. The inherited outward interval kernel
bounds exact multiples of pi and every full-source generator difference.

If `G_e` is the emitted anti-Hermitian source generator and `T_e` its
ideal proper-code generator, then

```
||exp(G_e)-exp(T_e)|| <= ||G_e-T_e||.
```

The full six-dimensional row-norm enclosure includes coupling out of
the intended two- or four-dimensional code. Scalar angle rounding and
port rounding are counted separately; their sums give `d_lit`. This is
an all-input bound, unlike a state-preparation estimate on one blank.
The code checks `d_lit<1e-8` on the admitted literal-word domain and
records the tighter actual rational value.

For each pulse, coherent load, source action and complete export form
a Stinespring isometry. The export contains every complement row and
satisfies the identity resolution on all six processor inputs. If
`V_e` denotes its successful code block, the word's successful map is
`V_L ... V_1`. At each first failure, all corresponding complement rows
are applied to the same successful prefix. The failed word halts and
retains the remaining field and reference. The identity

```
W*W+sum_f F_f*F_f=I                                   (4)
```

follows by telescoping the per-pulse resolutions. The exact ideal code
has zero failure amplitudes; these zero slots remain part of the common
comparison channel. Isometric telescoping bounds the entire literal
word, including every failure and arbitrary untouched reference, by
`d_lit`. In particular its mathematical total abort probability is at
most `d_lit^2`, without conditioning on successful runs.

Combining (2)--(4) gives a complete Stinespring comparison, up to the
scalar gauge, and hence a half-diamond channel bound

```
d_total=d_lit+B^3 |s|^3/(3n^2)+|s|E < 1e-8.            (5)
```

For the retained request the rational total is approximately
`6.385e-9`. An observable of norm at most one has the corresponding
complete Heisenberg error at most `2 d_total`, for any chosen norm-one
extension to the retained failure outputs. No division by success
probability occurs. The exact source and coherent load/export services
are the admitted simulator specification; their agreement with physical
apparatus is not proved by (5).

The raw host calculation executes the full six-level native pulse for
all sixteen regional input columns. Its frozen source tape and prefix
maps expose all first-failure Kraus maps on request. Independent matrix
exponentials check the polynomial target and modular sign. Those floating
residuals are diagnostics with separate admissions, not a proof that the
host exponential is correctly rounded. No polar correction or source
success normalization is performed.

## Continuation of the actual prepared source and observer records

Write `v_j` for every unnormalized successful field column of the
original native preparation, including its recorded probability. With
the region on the lowest four occupation bits, the continued success
columns are exactly

```
v'_j=(I_exterior tensor W) v_j,       j=0,...,255.       (6)
```

The signed CAR restriction is the same factorization used in (1).
Every actuator first-failure map similarly gives
`(I_exterior tensor F_f) v_j`. Its aggregate probability is
`Tr(sigma sum_f F_f*F_f)`. The original preparation abort outputs halt
before this operation and remain unchanged. The classical preparation
record indices are retained, together with field coherences involving
untouched quantum references and all first-failure labels. There is no source reset or
replacement by an ideal Gibbs array between preparation and continuation.

The full-preparation checker executes the original probability word,
all 255 adaptive readout events, all 256 conditional preparation words,
and (6) with the independently executed all-input regional map. It
retains the original preparation object and its failure reconstructors,
and exposes every actuator failure on all record columns. This use of
the all-input map is exact linear channel composition; it does not
substitute an assumed logical target for the native operation.

The default cold receipt replay executes the all-input actuator and
reverifies the mathematical preparation attachment; it does not rerun
the full conditional preparation. The separate `--prepared-check` and
the independent integration test execute those full preparation words.
Their numerical comparison forms `-log sigma` from the actually executed
unnormalized regional density and retains both preparation and actuator
failure probabilities.

Since the target `W_s` preserves sigma, (5) bounds the corresponding
regional stationarity error. On the full plaquette, a regional unitary
can change correlations and need not preserve Sigma. This does not
intertwine the regional modular group with the global one, and it does
not invalidate the finite modular-inclusion obstruction at a boundary.

The exact even polynomial target commutes with the exterior CAR algebra.
The approximate complete channel inherits its exterior/reference bound
from (5); no exact locality or parity statement about individual native
pulses is substituted for that bound.

## Counted resources and scientific boundary

The retained regional request has 37,568 pulses, the same number of
proper-code loads and exports, four regional working bits and four
unchanged exterior field bits. The eight preparation record bits and
every inherited preparation failure label are retained. The existing
source processor is reset by its specified export service. The actuator
requires no new source initialization, measurement, random choice or
erasure of the preparation record. Its exact pair-code count and complete
failure-row count are reconstructed from the emitted tape. The original
preparation costs remain those of its owning certificate and are not
absorbed into the actuator count.

For several admitted requests, complete channel errors add by
contractivity while records and quantum references are retained. The
exact regional target obeys its modular group law. Finite compiled
products approximate that law within the summed errors; they are not
declared an exact group of finite words.

This supplies an operational, corrected regional modular action on the
same finite interacting source and prepared reference. It supplies no
geometric cap, conserved Lorentz stress identification, regional area
normalization, positive Newton constant, fixed-volume generalized-entropy
stationarity, spatial continuum limit or physical time calibration.
