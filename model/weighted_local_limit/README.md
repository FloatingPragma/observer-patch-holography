# Local observations at large capacity

For the declared weighted repair family, a fixed set of cells approaches one
product reference as the total number of cells grows, for declared mean
occupations of two or four per cell. The [proof](PROOF.md)
gives an explicit bound without simulating the whole carrier, retaining every
local occupation sector. Its [certificate](results.json) includes capacities
larger than 10^122 binary sites. The density and weighted extension are supplied
choices; physical capacity, preparation and connected observation dynamics
require separate derivations.

From the repository root:

```sh
python -m model.weighted_local_limit.packet --check
python -m pytest -q -W error model/weighted_local_limit
```
