# Neutrino masses and oscillations

A neutral family tensor calibrated to six measured oscillation parameters gives Dirac masses of 0.001,0.00871206 and 0.05013980 eV for the stated lightest-mass choice. The actual native weak current produces complete neutrino and antineutrino records. Their free propagation agrees with the calibrated vacuum oscillation law within a total-variation bound below 4.1×10⁻⁹ over the specified phase domain. A single common scale multiplying the charged-lepton hierarchy cannot fit both measured mass splittings.

The masses and mixing are calibration inputs; physical beam preparation, matter effects and decay rates require their own source and apparatus comparison. See the [proof](PROOF.md) and [complete records](results.json).

Run `python3 -m model.openwave_neutrino_calibration --write` or `--check`. The calibration API and `neutral_family_tensor()` in `calibration.py` expose the native charged-family-frame tensor for explicit installation in other carrier constructions.
