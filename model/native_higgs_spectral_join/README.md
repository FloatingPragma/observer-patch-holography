# Elementary and Higgs probes at one native energy

The actual eight-cell charged source has a positive overlap between its
elementary e^c probe and the quantum-Higgs screened L probe in the same
energy bin. The full operator equation and the native ground bounds give
this result without a Higgs expectation value or a transferred two-cell
spectrum. A complete interference experiment has preparation-inclusive
joint-record contrast above 0.07247.

The [proof](PROOF.md) defines a source-only common-energy selector and
distinguishes it from the maximum elementary-residue selector. The
[receipt](results.json) checks every native current, raw lepton number,
parity and full-domain operator identity. Weak interactions, a relativistic
interacting pole and measured mass identification are outside this result.

```sh
python -B -m model.native_higgs_spectral_join --write
python -B -m model.native_higgs_spectral_join --check
python -B -Werror -m pytest -q model/native_higgs_spectral_join
```
