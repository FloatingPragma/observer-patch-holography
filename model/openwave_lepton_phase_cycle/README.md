# Coherent lepton mass readout

Four ancillary phase rotations isolate the native analyzer's charged
pairing on the full fermion space. A symmetric finite product then reads
the same three rest coefficients with a complete-output error below
10^-12 using 15229 coherent cycles, compared with the conservative
108 trillion reset-round budget. Every final ancillary outcome is retained.

The [proof](PROOF.md) preserves the supplied Higgs background, phase,
electron calibration and physical-identification boundary. This is a
resource reduction for a declared rest-mass apparatus, not a correction
to the muon ratio or a physical lepton identification.

```sh
python3 -m model.openwave_lepton_phase_cycle --write
python3 -m model.openwave_lepton_phase_cycle --check
python3 -m pytest -q model/openwave_lepton_phase_cycle
```
