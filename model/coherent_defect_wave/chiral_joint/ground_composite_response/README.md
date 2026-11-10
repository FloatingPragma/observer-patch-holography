# Composite response of the quantum ground

Applying the native neutral composite to the actual interacting ground
doublet gives nine positive energy lines. The lowest is the complete
local physical gap, approximately 0.0895553. Exact probabilities cover
all twenty-six states in each accessible charge sector, including the
energy-one states. The [proof](PROOF.md) and [receipt](results.json)
retain every probe setting, failure and complementary outcome.

The input is the selected normalized ground projector, with an explicit
fermionic reference in the probe apparatus. Spatial transport and a
many-cell stationary vacuum require their own source calculation.

Run from the repository root with the pinned `model/requirements.txt`:

```sh
python -m model.coherent_defect_wave.chiral_joint.ground_composite_response --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/ground_composite_response
```

`--write` regenerates the receipt. The certificate uses exact matrix
identities and rational root isolation. The test suite also compares
the weights with a separate numerical generalized eigensystem.
