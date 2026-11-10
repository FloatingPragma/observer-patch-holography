# Weak-current records

The actual finite native source has a positive Higgs-dressed weak-current
expectation and a nonzero inelastic current signal. Its energy-selected
charged state also has a positive current-density insertion. The
calculation retains the original physical preparation errors, all
detector failures and the charged reference energy.

[The proof](PROOF.md) fixes the source, current and ordered instruments.
[The receipt](results.json) includes a nonzero-phase probability above
0.03793 per original preparation attempt and a positive charged-bin
current mean. These records do not identify a W boson or a weak
scattering cross section.

From the repository root:

```sh
python -m model.native_weak_current_response --check
python -m pytest -q -W error model/native_weak_current_response/test_response.py
```
