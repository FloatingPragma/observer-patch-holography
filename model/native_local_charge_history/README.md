# Reading local charge transport

Complete native E records determine each decoded local occupation history
and its exact boundary flux. All snapshots use one native attempt counter;
the same complete-output error controls their joint law. Zero outcomes,
coverage erasures and the empty/full ambiguity remain explicit. The
[proof](PROOF.md) derives the operator continuity identity with arbitrary
reference correlations. A source occupation charge becomes a physical
particle charge only after its field and interaction maps are established.

```sh
python -m pytest -q -W error model/native_local_charge_history
```
