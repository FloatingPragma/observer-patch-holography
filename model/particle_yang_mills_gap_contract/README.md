# From Euclidean correlations to a spectral gap

Exponential decay of a dense set of local correlations excludes low-energy
spectrum of the same reconstructed Hamiltonian. A surviving finite-energy
excitation also keeps its threshold finite. The [proof](PROOF.md) gives
this conditional bridge without assuming equality of finite and continuum
spectra, and checks the exact clock scaling of the selected Wilson source.
Exact finite transfer matrices test the criterion and its failure modes.
The required continuum decay and Yang--Mills reconstruction are unproved
for that selected source.

Run `python -m pytest -q model/particle_yang_mills_gap_contract`, then
`python -m model.particle_yang_mills_gap_contract --write` and `--check`.
