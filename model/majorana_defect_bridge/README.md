# Neutrino masses in the same field action

The optional Majorana interaction uses the same native currents, Higgs
doublet and forty-five-component matter representation as the charged
action. Its full matrix and momentum spectrum retain every field. The
calibrated example has one zero neutrino mass and two masses set by supplied
oscillation splittings; it is a conditional calculation, not a mass prediction.

There is also a concrete [native fixed-repair readout](../common_source/native_repair_mass_attachment.py).
Its specified law fixes identical up, down and charged-lepton Yukawa matrices,
CKM=I and a zero Majorana coefficient. The receipt reconstructs that exact
matrix and checks its excluded charged-lepton ratio windows. Reading this
source specification fixes its model; it does not make that spectrum fit
nature or select the separate empirical Weinberg coefficient.

The exact mixing graph also constrains defect assignments. A massive
Majorana field can carry only an order-two additive charge. With full
mixing, excluding the zero-mass neutrino forces that charge to vanish for
all three neutrinos. These constraints leave composite and nonadditive
mechanisms separate from integer winding.

The [proof](PROOF.md) states the current, mass, mixing and topology domains.
The [receipt](results.json) retains the calibrated inputs and exact enclosures.
Physical identification requires a defect-to-field map and a controlled
interacting theory in the same reference and clock.

```sh
python -m model.majorana_defect_bridge --write
python -m model.majorana_defect_bridge --check
python -m pytest -q model/majorana_defect_bridge
```
