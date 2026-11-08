# Local experiments at growing capacity

The complete local repair record, including its visible coin and quantum
disturbance, has a controlled comparison across source sizes. A declared
common clock readout and matching local geometry are required. Exact
finite-size allowances combine the native dependency cone, marked event
counts and canonical state restriction, without enumerating the global
state space. Physical capacity, the selected density and a unique spatial
limit require their own derivation.

See the [proof](PROOF.md), [exact evaluator](limit.py) and
[source-bound receipt](results.json).

```sh
python -m model.native_local_history_limit.packet --check
python -m pytest -q model/native_local_history_limit
```
