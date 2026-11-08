# Measurements with a shared outcome menu

The existing nine-bit source controls implement up to eight successive weak
measurements on the same 66-dimensional source. Their outcomes form one
complete menu, including a quantum residual and every failure. The effects
can be noncommuting; the successful state changes are derived explicitly.
Each stage uses two admitted source words. Passing stored-primitive checks
gives a complete menu error below 1/50000. The
[proof](PROOF.md) distinguishes this bounded construction from universal
measurement access, arbitrary precision and a growing physical carrier.

Run the exact controls and reproducible receipt from the repository root:

```sh
python -m pytest -q model/native_coexistent_effects
python -m model.native_coexistent_effects.packet --write
python -m model.native_coexistent_effects.packet --check
```
