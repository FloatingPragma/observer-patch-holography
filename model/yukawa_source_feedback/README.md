# Couplings from ordered feedback

The exact memory theorem supplies coherent ordering and noise for the
correlated twelve-port feedback. Its next order gives an explicit algebraic
correction to the family coefficient. The [proof](PROOF.md) uses the owning
all-lag theorem and derives this coefficient without fitting a mass; the
[receipt](results.json) checks them against the full finite channel.

These are consequences of the declared logical feedback pulses and clock.
Identifying their coefficients with a fermion-Higgs vertex requires a
source-to-field attachment.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/yukawa_source_feedback
OPENBLAS_NUM_THREADS=1 python -m model.yukawa_source_feedback.response --check
```
