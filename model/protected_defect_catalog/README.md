# Protected objects before particle names

The [catalogue](CATALOG.md) distinguishes source memory charges, conserved
spin sectors and prescribed topological textures. It records how each is
prepared, protected and measured, together with its exact classification
domain. The producer replays complete bare-source and ordered-cover charge
classifications and exhausts all 65,536 Pauli-product charges of the selected
spin cube. The same source also has an exact integer spin charge; a local
creation has two nonzero charge-resolved outcomes with certified weights.
The ideal complete charge instrument retains all 256 records. No Standard
Model particle or physical mass is identified.

The [proof](PROOF.md) gives the graph classification and the remaining
source-to-particle requirements. A complete catalogue over unspecified source
architectures is not claimed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.protected_defect_catalog --check
python -m pytest -q -W error model/protected_defect_catalog
```
