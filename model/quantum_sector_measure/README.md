# Native quantum electric-sector measure

The actual 60-dimensional quantum link has an invariant ten-dimensional
Gauss-singlet carrier. Its complete Cayley/position channel fails closure
on the two weak SU(2) electric labels, including a coherent-state witness
and a diagonal within-sector witness. Retaining family and endpoint labels
closes the population law. The executable also solves the full stationary
quantum family, including coherences between two equal-dimensional blocks.

The resulting centre balance is `1-w_2-w_6`. It is 3/5 for a declared
tracial preparation, and zero in the limit of the native d=6 input.
This fixes conditional inputs to the seam scale formula; it supplies no
physical electroweak or gauge-coupling prediction.

Run from the repository root:

```
python3 -m model.quantum_sector_measure --write
python3 -m model.quantum_sector_measure --check
python3 -m pytest -q model/quantum_sector_measure
```

[PROOF.md](PROOF.md) states the source, invariant restriction, complete
records, exact stationary family and comparison boundary. `results.json`
binds the source imports and includes all original three-round records.
