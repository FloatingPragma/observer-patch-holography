# Quadratic twists of the actual source clock

The ordered-spin source has an operator-norm quadratic expansion under
each of its existing canonical Q and P controls. Its remainder is O(p²)
on the same long-wavelength clock. The quadratic coefficient is an actual
bond observable with a complete finite Pauli-score instrument.

The [proof](PROOF.md) derives an explicit local commutator bound that
transports macroscopic coefficient phases through the full source time
histories. They can be retained as central coefficients in an enlarged
history representation. Scalarity in the chosen selected state and a
closed equation on the original canonical fields are separate statements.
The aggregate bond-energy coefficient is scalar without ergodicity:
the source ground inequality and complete finite-region reset give an
explicit surface-to-volume fluctuation bound and finite-box energy bracket.
Exact source force checks and a positive two-oscillator control distinguish
quadratic control response from a single-frequency dynamical law.

The source action, state and clock are unchanged. No measured masses,
dispersion fit or independent Gaussian vacuum is supplied.

```sh
python -W error -m pytest -q model/spin_twist_closure
python -W error -m model.spin_twist_closure --write
python -W error -m model.spin_twist_closure --check
```
