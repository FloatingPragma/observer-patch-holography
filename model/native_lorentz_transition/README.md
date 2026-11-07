# Native null-direction transition test

Twelve exact source directions suffice to check whether an invertible
linear chart transition is conformal Lorentz. The package reuses the
native null-tomography kernel and decoder, classifies time orientation
and parity separately, and propagates the conformal factor through
composition. The source vectors are derived from the registered response
Gram. The proposed linear map and its physical null interpretation are
inputs.

[`PROOF.md`](PROOF.md) states the real-matrix theorem, the finite-error
bound and the source-attachment requirements. `transition.py` implements
exact `Q(√5)` controls and rejects floating-point inputs. Its
`conformal_form(B)` helper returns the factor for any symmetric form
passing the twelve zero tests; positivity requires the additional
invertible-congruence hypothesis. `classify(L)` checks that hypothesis,
and `proper_future(L)` checks both orientation signs.

Run from the repository root:

```sh
python -m pytest -q model/native_lorentz_transition/test_transition.py
python -m model.native_lorentz_transition.certificate --check
```

The exact certificate includes a native port-permuting rotation, an
algebraic boost, a rational rotation and a dilation. Spatial parity,
time reversal and singular null collapse are explicit countercontrols.
These controls verify a conditional finite implication. They supply
neither observer-selected chart transitions nor a physical spacetime
metric or dynamics.
