# Nuclear spin gradients and measured energies

A positive square of the difference between native spin currents excludes
two unwanted odd nuclear channels below their two-nucleon thresholds while
preserving five disclosed mass, axial and binding calibrations. The full
2,520-state singlet spectrum retains a separate spin-three bound channel
relative to its two-Delta threshold. This is a finite construction comparison,
with no complete physical nuclear identification.

The [proof](PROOF.md) derives the native operator, complete sectors and
fragment bounds. The [receipt](results.json) records calibrated GeV energies,
both excluded channels, all other levels and preparation failures.

Run `python3 -m model.openwave_spin_gradient_nuclei --check`.

The source's exact calibrated actions are declared in
[calibration_realizations.json](calibration_realizations.json), with the same
five empirical inputs and all seven retrospective gradient controls. The
producer authenticates the input and calibration-recipe hashes and reconstructs
the complete native action from those rational coefficients. `fit_parameters`
retains the numerical fitting recipe; evaluating a declared source does not
refit it. These remain calibrated agreements, not independent predictions.
