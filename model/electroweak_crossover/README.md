# Electroweak crossover of the declared coupling chain

The static long-distance modes of the hot electroweak theory are described by a three-dimensional SU(2) gauge
theory with a Higgs doublet. Its two parameters `x` and `y` follow from the four-dimensional couplings and the
temperature by dimensional reduction, and its lattice version has an exact continuum relation. This package
simulates that lattice theory, validates the code against published results, and places the couplings of the
declared chain in its phase diagram.

| Quantity | Value | Status |
| --- | --- | --- |
| `x` of the declared chain at the transition | 0.2961 | map imported, inputs declared |
| kind of transition | crossover | endpoint `x_c = 0.0983(15)` imported |
| pseudo-critical temperature | 161.13 +- 0.13 +- 1.6 GeV | measured; second error from the reduction |
| Standard Model lattice result (not blind) | 159.5 +- 1.5 GeV | imported |
| Higgs mass below which the chain gives a first-order transition | 72.5 GeV | imported endpoint, this map |
| cosmic time at the transition | 9.02e-12 s | radiation era, `g* = 106.75` |
| one-doublet row, assumption A | 145.11 GeV | Higgs and top masses assumed |
| one-doublet row, assumption B | 149.33 GeV | Higgs and top masses assumed |

The chain's couplings are declared forms. Their derivation is open. The package takes them as inputs.

From the repository root:

```sh
python -m model.electroweak_crossover --write
python -m model.electroweak_crossover --check
python -m pytest -q -W error model/electroweak_crossover
python -m model.electroweak_crossover.simulate --n 24 --beta-G 8 --x 0.296 --y -0.09 --sweeps 20000 --out run.npz
```

- [PROOF.md](PROOF.md) states each result with its label and citation.
- [lattice.py](lattice.py): the lattice theory in numpy, with the relations between lattice and continuum parameters.
- [kernels.py](kernels.py): compiled kernels of the production runs. They need numba. Without it `simulate` uses `lattice.py`.
- [reduction.py](reduction.py): the map from four-dimensional parameters and temperature to `(g3^2, x, y)`.
- [reweight.py](reweight.py): multi-histogram reweighting with block jackknife.
- [production.json](production.json): the lattice measurements of the production runs. The raw series are not stored.
- [certificate.py](certificate.py) and [results.json](results.json): the derived numbers and 44 named checks.
- `test_lattice.py`, `test_reduction.py`, `test_results.py`: gates. They run in under a minute.

The production used 53 runs at the chain's point and the validation runs listed in `production.json`.
