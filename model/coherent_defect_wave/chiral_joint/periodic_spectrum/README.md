# Charged energies on the periodic quantum source

Two cells exchange charged excitations through two quantum links while both
Higgs memories interact with the leptons. The complete neutral gauge sector
has 3424 states and a unique ground. A compensated charge-exchange detector
sees equal addition and removal thresholds, `0.97009612853 +/- 3e-8` in
the shared clock's units, with a certified positive weight near `0.22249697`.
All preparation failures and detector alternatives are retained. Joint
source/reference translation phases are also recorded; they do not determine
an intrinsic particle dispersion or a GeV mass.

The [proof](PROOF.md) fixes the active source, state, clock and apparatus.
The [receipt](results.json) records the complete physical spectra, readout
intervals and source hashes. This calculation uses the declared reduced
lepton/link/scalar system; the interacting full-family ground is a separate
calculation.

From the repository root:

```sh
python -m model.coherent_defect_wave.chiral_joint.periodic_spectrum --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/periodic_spectrum
```

Replay checks the stored rational eigenbasis proposals with exact arithmetic.
Use `--generate-witnesses` explicitly to replace those proposals, followed
by `--write` to regenerate the receipt.
