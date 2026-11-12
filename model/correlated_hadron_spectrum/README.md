# Correlated hadron spectrum

The retained quenched Wilson correlators are fitted with their complete
time covariance and a fixed conditioning rule. Every one of the 32 paired
stream deletions is refitted with its own covariance. All 48 channel, sink,
window and kappa fits remain visible; the original diagonal-time fit is a
separate replayed control.

The [receipt](results.json) carries the paired chiral and GeV scale chain,
joint covariance, fit-choice sensitivity and original unresolved nucleon
label. The rho-anchored pion remains a calibration identity. The separate
string anchor exposes scale dependence. No sea-quark, continuum, volume or
physical source-clock error is bounded.

```sh
VECLIB_MAXIMUM_THREADS=1 python3 -m model.correlated_hadron_spectrum --write
VECLIB_MAXIMUM_THREADS=1 python3 -m model.correlated_hadron_spectrum --check
python3 -m unittest model.correlated_hadron_spectrum.test_fit -v
```

The [proof](PROOF.md) states the covariance algebra, all fixed numerical
rules and the distinction between first-stage correlated fitting and the
retained approximate chiral regression.
