# Preparing a measurable correlation

Fair native repairs prepare a complete four-read correlation experiment from
an arbitrary quantum state, with an explicit error bound. The calculation
controls only the low-degree observable needed by the experiment, avoiding
the exponential dimension in a full-state preparation bound. All raw records
and resource stops are retained. The result uses a declared graph, occupation
and attempt clock; identifying a physical vacuum and particle spectrum
requires their common field interpretation. See the [proof](PROOF.md),
[implementation](preparation.py) and [receipt](results.json).

```sh
python -m model.native_correlation_preparation.packet --write
python -m model.native_correlation_preparation.packet --check
python -m pytest -q -W error model/native_correlation_preparation
```
