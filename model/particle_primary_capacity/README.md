# Primary capacity in nats

The [proof](PROOF.md) counts the actual Wilson carrier and bounds its
Gauss-constrained physical dimension. The smallest ambient member has
approximately447.41nats; physical capacity is smaller and not counted
exactly. Capacity, prepared-state entropy, correctable record capacity and
classical simulation memory are different quantities.

`capacity.py` authenticates the source and uses its complete degree-one
119-state link space. All higher-cutoff values remain bounds. No mass
capacity law or cosmological capacity is inferred.

```sh
python -m model.particle_primary_capacity --write
python -m model.particle_primary_capacity --check
python -m pytest model/particle_primary_capacity
```
