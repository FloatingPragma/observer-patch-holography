# Modal Maxwell factorization boundary

This note records a laptop-scale theorem package for the V3
electromagnetism lane, issue #733. It does not promote the OPH light-signal
surface to physical electromagnetism.

## Result

The existing free-vector construction provides, for chart scale `a` and
momentum `k`, a nonnegative exact FZ-12 frequency `omega(a,k)`, a scalar
spatial action `omega(a,k)^2`, and a two-dimensional transverse fibre when
`k` is nonzero.

On that fibre define the modal operator

```text
C_(a,k)(v) = [omega(a,k) / sqrt(k.k)] (k x v).
```

Its complex Fourier-curl form is `D_(a,k) = i C_(a,k)`. The first-order
Maxwell-shaped candidate used below is

```text
G(E,B) = (D_(a,k) B, -D_(a,k) E).
```

`Lean/Screen/ModalMaxwellFactorizationBoundary.lean` proves exactly that:

- `k.C_(a,k)(v) = 0`, the modal `div curl = 0` identity;
- `C_(a,k)^2(v) = -omega(a,k)^2 v` for transverse `v` and nonzero `k`;
- `D_(a,k)^2(v) = omega(a,k)^2 v` on the complex transverse fibre;
- the complex spatial action is exactly the scalar extension of the existing
  real `photonSpatialAction`, not a separately inserted kinetic symbol;
- the opposite-sign paired generator `G` squares to minus the complexified
  FZ-12 spatial action on both amplitudes, and therefore obeys the committed
  second-order oscillator equation;
- `G` lands in a complex transverse fibre for each of its two amplitudes;
- a same-sign mutation `(D B, D E)` instead squares to the positive spatial
  action and fails the oscillator equation on a tested amplitude whenever
  the spatial action on that amplitude is nonzero.

The proof is a composition of the existing exact FZ-12 frequency identity
with the standard three-dimensional vector triple-product identity. It adds
no project axiom and runs no simulation.

This construction removes no registered premise. PR-20, PR-21, and PR-22
remain declared inputs to the composed OL-F1 surface; defining and factoring
a candidate generator does not source-produce the selected action or realize
its oscillator and physical-frequency readings. PR-53 and PR-54 still own the
open physical photon/field attachment and the source gauge-field, current, and
action attachment.

## Adversarial control and what is still open

The theorem `sameSignCurlMutation_fails_wave` is the adversarial control. It
keeps the same momentum, FZ-12 frequency, complex Fourier-curl block, and two
transverse fibres, but changes only the relative off-diagonal sign. The
resulting second-order sign is wrong. Thus the opposite-sign pairing is an
exact factorization of the existing modal oscillator, not merely a suggestive
notation.

More importantly, the coefficient `omega(a,k)/|k|` is momentum dependent.
The construction is an exact modal/pseudodifferential factorization, not a
proof of a local position-space curl operator. The remaining gap includes:

- a source-produced dynamics and identification of the two amplitudes as
  electric and magnetic fields;
- assembly of modal fibres into one real field with the required reality
  pairing between `k` and `-k`;
- existence and locality of the assembled position-space operator;
- a U(1) gauge potential, gauge quotient, and Maxwell action;
- a typed bridge from the rational seam-incidence Gauss receipts to the
  modal divergence used here;
- a conserved physical current and its inhomogeneous coupling;
- continuum/refinement control, Lorentz covariance, and laboratory readout.

Thus this result supplies only the exact pointwise/modal algebraic
factorization that was previously absent. Local and source-produced Maxwell
dynamics remain open. It does not discharge any registered premise, promote
an observation row, satisfy the positive closure contract, or close issue
#733.

## Independent exact replay

`code/electromagnetism/modal_maxwell_factorization.py` replays the complex
vector identities over exact Gaussian-rational arithmetic at an off-axis
transverse mode. Its tests verify the opposite-sign factorization, the
same-sign mutation, a wrong-normalization mutation, and the zero-frequency
boundary. The replay deliberately records that it does not independently
compute the FZ-12 symbol; Lean supplies the universal composition with that
symbol. This is a finite algebra check, not an empirical campaign.

## Novelty and flagship assessment

The vector-calculus core and Fourier Maxwell sign pattern are standard. The
useful contribution is their exact composition with the OPH-specific FZ-12
symbol plus a machine-checked sign mutation. This is modest formal progress
for detailed electromagnetism/issue custody, not a novel physical Maxwell
derivation, and does not merit scarce flagship space.
