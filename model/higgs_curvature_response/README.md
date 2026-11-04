# Higgs curvature, stress and entropy variation

A curvature coupling of the same Higgs doublet changes both its stress
tensor and the Wald area coefficient. The [proof](PROOF.md) derives their
common normalization and a corrected generalized-volume variation on small
balls, with an explicit sixth-order remainder. Exact controls retain
nonzero gauge and Yukawa exchange, the full Higgs mass kernel and every
polynomial moment. The [receipt](results.json) binds these calculations to
the retained matter action. A physical entropy state, stationarity at the
specified generalized volume and a common curved-source limit are separate
inputs; the result implies a conditional linearized scalar-tensor response.

```sh
python -m model.higgs_curvature_response --write
python -m model.higgs_curvature_response --check
python -m pytest -q -W error model/higgs_curvature_response
```
