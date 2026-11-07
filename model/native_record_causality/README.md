# Independent repairs

Two fixed repairs on disjoint native endpoints commute, including their
complete quantum output and every exported or retained record, when the
command record slots are compared in the same order. The exact proof covers
both full native carriers and arbitrary entangled references. Adaptive
feedback and early memory timeouts provide explicit counterexamples to a
broader claim. This finite observation equivalence supplies no spatial
dimension or relativistic field limit. See the [proof](PROOF.md) and the
[replay certificate](results.json).

```sh
python -m model.native_record_causality.packet --check
python -m pytest -q model/native_record_causality
```
