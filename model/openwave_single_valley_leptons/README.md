# Three charged-lepton minima

A declared nearest-neighbour Higgs-Wilson term gives the complete native
charged source exactly one momentum minimum per family. The three
colourless charged branches have electron-anchored rest energies
0.00051099895069, 0.1056594144959 and 1.7769849708134 GeV.
All quarks, neutral modes and high-energy corner states are retained.
The phase and electron calibration are explicit; interacting physical
identification and the neutral continuum are separate requirements.

The [proof](PROOF.md) derives the unique minimum and complete-record error
bound. The [receipt](results.json) contains every side-four spectrum and
complete arrival law. Reproduce with:

```
python -m model.openwave_single_valley_leptons --write
python -m model.openwave_single_valley_leptons --check
python -m pytest -q model.openwave_single_valley_leptons
```
