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
The [selected finite-ground equation](STATE_CONTROL.md) sharpens the
regional energy-density squared error to `36J^2/L^2` and bounds its
complete excitation tail. Opposite ordered source phases have identical
energy statistics and distinct local spin readouts, so this energy
estimate does not select a full regional state.
Exact source force checks and a positive two-oscillator control distinguish
quadratic control response from a single-frequency dynamical law.

The [scalar polarization limit](POLARIZATION_LIMIT.md) removes the actual
extra directional current in the Q-field commutator. Its weighted test
is a discrete derivative with zero integral, which kills the entire
mean-ergodic current component. The source energy bound transports this
vanishing through every fixed finite time history. Separate boost and
unbounded-generator limits are not inferred.

The source action, state and clock are unchanged. No measured masses,
dispersion fit or independent Gaussian vacuum is supplied.

```sh
python -W error -m pytest -q model/spin_twist_closure
python -W error -m model.spin_twist_closure --write
python -W error -m model.spin_twist_closure --check
```
