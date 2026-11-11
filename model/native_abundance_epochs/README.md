# Native defect populations and persistence

The native integer-load repair and six-class register generate an exact population history on a closed cell
with three ports and two faces. Every attempt is retained. The class shares remain normalized, total load and
face charge are conserved, and the combined load/register cost decreases. Under independent uniform initial
loads from zero to five and zero initial flux, the persistent occupied-face share is exactly `1007/13824`,
about `7.28443%`.

The [proof](PROOF.md) derives the law and its finite-epoch bounds. The [receipt](results.json) retains all 129
population cuts from zero through 128 attempts, the exact terminal probability and the source hashes. Attempt
counts have no physical temperature or time calibration; the cell has no physical particle dictionary.

```sh
python3 -m model.native_abundance_epochs --write
python3 -m model.native_abundance_epochs --check
python3 -m pytest -q -W error model/native_abundance_epochs
```
