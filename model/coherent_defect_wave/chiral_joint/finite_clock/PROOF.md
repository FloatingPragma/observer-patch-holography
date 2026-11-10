# Exact finite clocks with controlled interacting records

## Same Hamiltonian and declared finite service

Let H be a finite self-adjoint Hamiltonian on its complete physical Hilbert
space, with `||H||<=B`. A nonorthogonal rational basis uses its positive Gram
matrix G, so self-adjointness means `H^T G=G H`. The scalar time t is expressed
in the declared source clock. Neither an energy shift nor a replacement
generator is made.

The architecture permits exact algebraic finite matrix operations. This is
the same type of declared functional-calculus service as the finite Cayley
clock in [source-state response](../../../particle_field_state_response/PROOF.md).
The calculation below supplies its explicit rational operation and error.
It does not assert an efficient local circuit on the full cubic carrier.
Implementing this service within a narrower hardware instruction set requires
that implementation and its resource bound.

For `X=tH`, define

```
a = I-X^2/10,
b = X/2-X^3/120,
R33(X) = (a+i b)^(-1)(a-i b).
```

All factors are polynomials or inverses of polynomials in the unchanged H.
They commute, and

`D=a^2+b^2=I+X^2/20+X^4/600+X^6/14400 >= I`.

Consequently the inverse exists on every finite self-adjoint source, and
`R33^dagger R33=I`. Every exact symmetry or Gauss constraint commuting with H
commutes with R33. The source energy of H is conserved. The nonlinear
eigenphase of R33 is not identified with a replacement physical energy.

## Global operator error without a small-time assumption

Set `p(z)=1+z/2+z^2/10+z^3/120`. Repeated integration by parts gives the exact
identity

`exp(z)p(-z)-p(z) = -z^7/720 integral_0^1 exp(uz) u^3(1-u)^3 du`.

The code verifies this identity symbolically using the recurrence for each
integral of `exp(uz)u^k`. The apparent divisions by z in that recurrence
extend continuously to z=0. The weight integral is exactly `1/140`.
For imaginary z, `|exp(uz)|=1`, while the positive denominator identity gives
`|p(-z)|>=1`. Thus

`|p(z)/p(-z)-exp(z)| <= |z|^7/100800`.

The finite spectral theorem, with `z=-itE`, proves

`||R33(tH)-exp(-itH)|| <= |tB|^7/100800`.

For n equal substeps, telescope the difference of products of unitary
operators. The complete bound is `|tB|^7/(100800 n^6)`. The bound holds on
the full space, hence on every prepared state, including entanglement with
a retained reference. There is no projection-dependent omitted remainder.

For any normalized initial and final states, the amplitude error is at most
this operator bound. A complete projective measurement retains its
orthogonal complement. Its outcome total variation is also bounded by the
pure-state trace distance, at most the operator error. Rare-event positivity
can require the stronger amplitude comparison rather than subtracting this
error from an extremely small probability.

## Neutral propagation on the actual Peter–Weyl source

Consume the [complete Peter–Weyl experiment](../peter_weyl_links/PROOF.md),
including all five matter representations, both spin preparations and the
physical fermionic reference. The unchanged Hamiltonian has bound B=1896.
At `t=151/4641746545968`, its exact exponential has an opposite-spin neutral
arrival amplitude with leading magnitude

`L=151 t^5/3840`

and full Taylor remainder at most L/2. The one-step rational clock has error
`delta=(1896 t)^7/100800`, and exact rational comparison gives

`0 < delta/L < 1/4` (approximately 0.02351544).

The finite-clock arrival amplitude therefore lies between L/4 and 7L/4 in
absolute value. Its complete exclusive arrival probability is enclosed by

`[L^2/16,49L^2/16]`.

The complementary probability lies in `[1-49L^2/16,1-L^2/16]`. The source's
local packet-and-star detector contains this exclusive target, so its
complete local arrival probability has lower bound L²/16. Its complement
is retained. Same-spin arrival remains exactly zero because the rational
clock preserves the same rotation and seam symmetries. The source and its
reference are not replaced by a normalized successful subensemble.

This is a finite service certificate on the huge cubic source. The code
does not materialize or invert that entire Hilbert-space matrix. The exact
functional-calculus theorem applies to its declared finite H, whose full
norm and source moments are established by the consumed proof. All consumed
inputs and their transitive pins are checked before the receipt is accepted.

## All five matter representations with the same finite clock

The [compensated family transport](../family_transport/PROOF.md) experiment
uses one source, one static dual-native reference and time `t=1/1508` for
all five native representation blocks. Its prepared local arrival amplitude
is at least `5t/24`; the preparation-off local amplitude is at most `t/8`.
These estimates follow from complete first and second source powers, with
all interaction and reference states retained.

Two R33 substeps of the same H have full state error
`epsilon=(1896t)^7/(100800*2^6)<t/96`. Therefore the prepared local amplitude
is at least `19t/96`, while the preparation-off local amplitude is at most
`13t/96`. The two probabilities have the strict same-clock contrast

`Pr_prepared(F_r)-Pr_unprepared(F_r) >= t^2/48 > 0`.

Each setting retains the complete no-arrival outcome, and the static
reference remains part of the state. The preparation-off input is a control
input; it is not identified with the stationary interacting vacuum. The
same two-substep clock and bound apply to every native representation.

## Complete rational execution on the periodic source

A separate concrete execution uses the complete 72-state physical block of
[matter-generated ring dynamics](../matter_induced_ring/PROOF.md). Its
actual Hamiltonian, positive Gram matrix, initial compensated e^c state and
all 67 eliminated matter alternatives are retained. No energy-dependent
Schur kernel is used for time evolution. This source has six-state links;
its spectrum is not transferred to the Peter–Weyl119 source.

At t=1/64, solve the two rational systems

```
D Re(psi_out) = (a^2-b^2) psi_in,
D Im(psi_out) = -2ab psi_in.
```

The code checks both full residuals and exactly verifies the final G-norm
equals one. The two normalized loop projectors have ranks one and one;
their orthogonal complement has rank seventy. Their three probabilities
are exact nonnegative rationals summing to one, with both loop probabilities
strictly positive. The receipt includes each loop's complex amplitude as
two rationals and the probability of every other state. The projectors
extend to the other conserved blocks with zero loop success and identity
complement. Those blocks have zero initial weight.

The finite source/reference preparation is the declared preparation in the
ring proof. A three-state pointer records the projective instrument via
the controlled cyclic shift. Every outcome and its unnormalized conditional
state is retained. The record probabilities are the squared norms of these
branches, with no success conditioning. The general bound above compares
this complete finite clock with the unchanged source exponential; a numerical
generalized-eigenvector calculation is an independent test, never the
rational probability certificate.

## Domain

These are exactly unitary finite source clocks and complete observer records.
Clock approximation error is separated from spectral, continuum and unit
identification questions. The result supplies no many-cell stationary
vacuum, relativistic mass shell, Standard Model identification or GeV scale.
No empirical mass, coupling or time standard enters the calculation.
