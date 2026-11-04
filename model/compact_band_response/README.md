# Work and entropy on the compact interacting defect

Changing the source hopping moves its isolated defect band. The interaction
acts on that band, so its force includes the derivative of the localization
projector. [The proof](PROOF.md) derives this term, a volume-independent
bound, the same localized work readout on every memory record, and thermal
response using the compact particle operator. The complete source calculation
retains every complementary mode. Relativistic stress, a physical metric and
regional Einstein entropy balance require further construction.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.compact_band_response --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/compact_band_response
```
