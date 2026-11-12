# A coupling measured at different scales

The declared colour gauge ensemble gives a coupling that increases with loop
size at one unchanged bare coupling. Four independent cold/hot productions on two
volumes supply the raw records. An exact finite free-field normalization
removes the geometric scale trend. Finite-volume effects are explicit;
the result supplies neither physical QCD running nor a particle-core scale.
See [the derivation](PROOF.md) and [results](results.json).

Analysis: `python3 -m model.openwave_scale_coupling --write` or `--check`.
Tests: `python3 -m pytest -q model/openwave_scale_coupling`.

The production environment uses Python 3.12.6 and [pinned dependencies](requirements.txt).
Generate the raw records from the repository root:

```sh
python -m model.openwave_scale_coupling.production --side 6 --beta 6 --seed 3101 --start cold --samples 400 --warmup 200
python -m model.openwave_scale_coupling.production --side 6 --beta 6 --seed 3102 --start hot --samples 400 --warmup 300
python -m model.openwave_scale_coupling.production --side 8 --beta 6 --seed 3103 --start cold --samples 400 --warmup 300
python -m model.openwave_scale_coupling.production --side 8 --beta 6 --seed 3104 --start hot --samples 400 --warmup 400
```

Every command writes its measured records, final configuration and metadata
under `data/`. Analysis and production are separate commands.
