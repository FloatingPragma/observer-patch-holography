# Local observations at large capacity

The native fair federation's stationary fixed-number state has an exact local
restriction that can be computed without constructing the whole carrier.
Every local number sector is retained. The calculation reuses the
[canonical-sector readout](../observer_time_fields/canonical_sector.py),
extends it to arbitrary finite windows and integer populations, and checks
nested restrictions and finite-size error bounds. It does not select a
physical capacity law or identify port count with cosmological capacity.

The [source-family interface](PROOF.md) states that domain. A separate
[dynamic product restriction](PRODUCT.md) preserves the weighted native
repair/hold observer, complete records and arbitrary correlated references
when unused cells are appended. Its coherent proper-code services and a
connected spatial field are outside this embedding.

Reproduce the exact receipt with:

```sh
python -m model.native_capacity_restriction.packet --check
python -m pytest -q model/native_capacity_restriction
```
