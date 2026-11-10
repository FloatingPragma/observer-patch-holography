# Complete position and field records after a finite energy clock

The source, reference, preparation and autonomous Hamiltonian are exactly
those of the [finite charged clock](../affine_gauss_charged_clock/PROOF.md).
The detector accepts a specified set of its twenty-qubit phase labels.
It implements no exact ground projector. Conditional on this accepted
record, the actual source density operator rho has ground fidelity
`<g|rho|g> >= f`, where the parent proves `f>0.99999999439512`.
Every rejected phase label and the pump failure remain records.

## From a ground certificate to an actual detector probability

Let P be any orthogonal source projector, `p=<g|P|g>`, and
`q=tr(rho P)`. Purify rho and choose an ancillary unit vector whose
product with g has overlap `sqrt(<g|rho|g>)`. Splitting this overlap with
P and I-P, then applying Cauchy-Schwarz to the two terms, gives

```
sqrt(f) <= sqrt(p q)+sqrt((1-p)(1-q)).
```

Writing p and q as squared sines solves this elementary inequality. The
sharp lower bound is zero when `p<=1-f`; otherwise it is
`[sqrt(f p)-sqrt((1-f)(1-p))]^2`. The sharp upper bound is one when
`p>=f`; otherwise it is `[sqrt(f p)+sqrt((1-f)(1-p))]^2`.
Both envelopes are monotone in p on their applicable intervals. Thus
an interval `[p_lo,p_hi]` for the actual ground gives a rigorous interval
for q by evaluating the lower envelope at p_lo and the upper at p_hi.
The implementation encloses every square root using integer square-root
inequalities on rational endpoints. No floating-point square root decides
a probability bound.

This retains the actual detector error. Reporting the certified ground's
probabilities as exact probabilities of the accepted finite-clock state
would omit that error. A global trace-distance estimate is also valid but
can be unnecessarily loose for a rare event; the binary fidelity envelope
uses the size of that event explicitly.

## Twenty-six physical joint outcomes

In the complete one-charge carrier, let Pi_(x,0) select source vertex x
and transverse vacuum, and Pi_(x,1) select that vertex and every nonvacuum
transverse occupation through the declared total cutoff. There are thirteen
vertices. The twenty-six orthogonal projectors sum to the identity on the
physical one-charge sector. They preserve the longitudinal Gauss embedding
and the retained compensator. Their position and transverse-number factors
commute with each other, but need not commute with the Hamiltonian.

The [full residual verifier](../affine_gauss_full_spectrum/PROOF.md) admits
the complete 70928-component binary trial v and bounds its sine angle to
the actual unique ground by s. For each of the twenty-six projectors,

```
|<g|Pi|g>-||Pi v||^2/||v||^2| <= s.
```

All trial norms are summed as exact rational squares of the binary64
entries. No eigenvector component or transverse occupation is dropped.
The fidelity envelope above then bounds the actual post-clock probabilities.
A separate application to each of the thirteen position projectors and
to the total nonvacuum projector gives tighter marginal bounds. Individual
interval endpoints are not renormalized; they enclose one complete law.

## Ordered instrument, reference and autonomous time

Let D_z create the normalized compensated source wave packet. The actual
pump prepares it with probability 1/4. Let K_k be the implemented rational
finite-clock Kraus operators. They retain their common autonomous source
evolution. After the phase record, perform the joint projective readout
only for an accepted label. Its complete original-attempt probabilities are

```
pump failure:                      3/4;
rejected label k:                  ||K_k D_z|vac>||^2/4;
accepted k, vertex x, flag b:       ||Pi_(x,b) K_k D_z|vac>||^2/4.
```

All phase labels, all twenty-six joint outcomes and every zero-probability
number branch remain. The effects sum to one by completeness of the clock
and the joint readout. A twenty-six-state pointer records the joint outcome;
on the complementary number domain its controlled interaction is the
identity and the existing number flag identifies that branch. Its logical
capacity adds ln(26) nats. The finite controlled readout is
`U_read=sum_j Pi_j tensor S_26^j`, with identity on the complementary
number sector. Here S_26 shifts the initially blank pointer and j labels
the twenty-six outcomes. Let L be a Hermitian logarithm with
`exp(-iL)=U_read` and `||L||<=pi`. During one source-time unit, declare
the interaction `V(t)=exp(-itH_raw)L exp(+itH_raw)`. Its complete propagator
is exactly `exp(-iH_raw)U_read`: the autonomous source is never paused.
The pointer has zero Hamiltonian, so this final common source evolution
does not change its outcome probabilities. The recorded position refers
to the start of this readout interval, not a simultaneous final-time
position eigenstate. The source energy change is bounded by twice the
raw source norm; the controller supplies that work.

This interaction-picture finite-matrix control is a declared apparatus
service, not a compiled passive repair circuit. It preserves the full
Gauss and number constraints. Synthesis, controller and environment
workspace are separate from the logical pointer capacity.

The receipt's intervals aggregate the accepted k values; multiplying their
conditional probability bounds by the parent acceptance interval gives
original-attempt bounds. The phase labels themselves are not erased.
Ground fidelity is invariant under any additional autonomous waiting time,
so these bounds hold uniformly before the joint readout. They do not imply
that the accepted state is stationary, that position was definite before
measurement, or that its post-measurement evolution is unchanged.

This result supplies actual finite-clock charge/field observations with
complete error accounting. The finite boundary problem has no assigned
Standard Model identity, relativistic mass shell or GeV mass.
