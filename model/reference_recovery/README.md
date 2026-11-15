# Source implementation of finite reference recovery

This package compiles the adjoint of the existing
[reference observation](../reference_observation/README.md) into native source
pulses. It implements a bounded finite service for a **supplied faithful
reference and CAR restriction**. It does not recover an arbitrary lost global
state, select the reference, or establish a physical clock.

For restriction Kraus operators \(K_a\), global reference \(\rho\), and
\(\rho_A=\sum_a K_a\rho K_a^*\), the recovery Kraus operators are

\[
 R_a=\rho^{1/2}K_a^*\rho_A^{-1/2},\qquad
 V=\sum_a |a\rangle\otimes R_a.
\]

Faithfulness and trace preservation give \(V^*V=I\). Tracing the retained
environment gives the Petz channel, and applying it to \(\rho_A\) gives
\(\rho\). Its adjoint is the same reference observation map. A global
occupation-pointer readout after recovery therefore implements its regional
POVM. Both pointer outcomes remain present; conditioning on a pointer result
does not reproduce the original global reference.

## Compilation and complete outputs

`compile_isometry` eliminates the two or four admitted input columns by
Gray-row QR. The initial Gray permutation, diagonal phases and inverse
two-row rotations are all compiled with the existing native source compiler.
This constructs a unitary completion on at most six work bits. It does not
introduce an arbitrary matrix gate. All other input columns see that same
completion; the unused input bits are declared initialized blanks. Zero-angle
requests count as pulses. The bounded admission limit is one million pulses.

`execute_complete` executes the full six-dimensional source primitive for
every pulse and retains each complement Kraus branch at every first export
failure. Separate complement branches add incoherently. The literal word,
initial columns and full primitive matrices are immutable snapshots, so each
failure vector can be reconstructed by `run.failure(event)`. Streaming their
Gram matrices reduces host storage; it neither discards outputs nor
normalizes successful runs. The last two entries in the four-column weight
array are padding for a two-code (`J`) pulse, whose complement has dimension
two. Source service counts include a load and export for every pulse.

For retained success map \(A\), target isometry \(V\) and failure operators
\(F_j\), let \(G=\sum_j F_j^*F_j\). The reported quantity
\(\sqrt{\|A-V\|^2+\|G\|}\) bounds the complete dilation error, including an
untouched reference. The failure Gram is a column Gram, not a reduced density
matrix; its trace on identity input columns is a sum over the complete input
basis, not a normalized failure probability. On a normalized prepared input
its trace is the actual failure probability. Floating residual admission is
a numerical control, not an outward-certified hardware error bound.

`nested_word` composes recovery services only when their intermediate
dimensions and faithful references agree. It retains both environments. The
same-reference tower identity is tested against direct recovery on every
matrix unit. This is not a regional retraction or modular covariance theorem.

## Interacting-field instance

The default instance uses all four CAR modes of the two-site Weyl/contact
field, coordinate speed `1/8`, contact `U=1`, and the faithful Gibbs reference
`beta=1, mu=0`. Recovery from the two-mode region uses six work bits and
348,697 native pulses. `reference_preparation` also compiles the regional
reference purification, recovery and a global-mode pointer, retaining the
regional purifier and recovery environment on nine initialized bits. Thus
this route does not assume an already prepared correlated global input.
The preparation, design parameters and initialization are still supplied
resources. This faithful reference is distinct from the pure half-filled
ground used for the conditional particle candidates.

The tests cover arbitrary complex isometries, full input/reference columns,
an independent matrix-power construction of the Petz channel, nested
recovery, complex references, both pointer outcomes, injected source leakage,
separate nonadjacent export branches and malformed primitives. Run:

```sh
python -m pytest -q -W error model/reference_recovery
python -m model.reference_recovery.recovery
```

The result supplies CP-02's bounded source implementation of the reference
map. Common native full-field reference selection, scalable compilation,
physical initialization and clock calibration are open. No canonical
physical premise is discharged by this package.
