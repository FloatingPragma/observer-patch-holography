# Reproducible OpenWave experiments

Three headless tests execute the same complete native 180-coordinate field
inside a pinned OpenWave checkout: stable free propagation, the three charged
rest coefficients, and the finite lattice mass shell. Each compares every
source-coordinate amplitude with an independent matrix-symbol solution,
retains preparation failure and rejects a changed frequency. The analyzer
and electron calibration retain their target exposure; these tests provide
no unexposed physical prediction or official OpenWave status.

Run from RER with its scientific dependencies installed:

```
python3 -m model.openwave_admission_packet --write --openwave-root /path/to/openwave --python /path/to/openwave/python
python3 -m model.openwave_admission_packet --check --openwave-root /path/to/openwave --python /path/to/openwave/python
```

The OpenWave checkout must be clean at
`ef55175d04737e73c92f33f9f9cf3d17c593da54`. Its untracked
`openwave/xperiments/oph_admission` overlay contains the portable original
Taichi engine, exported source and three-test runner. Unrelated existing
overlay files are rejected. No tracked upstream file, cell or submission is
changed. [The mathematical domain](PROOF.md) and
[input/prediction audit](prediction_audit.json) delimit admission claims.
