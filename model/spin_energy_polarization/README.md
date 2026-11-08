# Source energy polarization with growing cutoffs

The actual ordered cubic spin source admits a measured energy-polarization
change under its canonical field preparations. On one fixed source, state
and clock, this change has a controlled spatial limit and the exact affine
translation law expected of an energy first moment. The result uses the
existing source energy and current identities. It supplies a concrete part
of a boost construction, while retaining the field-action and state-covariance
requirements needed for a physical boost.

The [proof](PROOF.md) establishes:

- exact cutoff independence on compact preparation supports;
- explicit mesh and cosine errors for the actual weighted energy change;
- the source spatial-shift law, including its common energy clock;
- an exact bond-current remainder in the candidate scalar-field action;
- Schwartz weighted-energy tails and complete charged-preparation truncation.

The field-limit readout is a source energy-polarization functional. Its
identification with a continuum boost operator still requires control of
the actual propagated source current correlations and energy transport into
the limiting representation. The existing finite-rotation/one-boost theorem
consumes operational covariance; this construction does not assert it.

```sh
python -W error -m pytest -q model/spin_energy_polarization
python -W error -m model.spin_energy_polarization --check model/spin_energy_polarization/results.json
```

The rational mesh controls use a translated compact C3 test. Results are
reported per unit of the actual source coefficient kappa. Its value is
not assigned by those controls. No large native tape, empirical species,
mass shell or Lorentz speed is inferred.
