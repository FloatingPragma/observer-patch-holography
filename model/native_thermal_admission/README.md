# Native equilibrium and its obstruction

The unchanged 1,296-state population law has 46 recurrent classes and a
45-dimensional simplex of stationary measures. Every full-support Gibbs
law at finite positive temperature is excluded in the fourteen nontrivial
load sectors by exact negative cost drift. The two endpoint sectors have
identity dynamics. Recurrent-only thermal weights are optional inputs,
not selected native abundances.

The [proof](PROOF.md) gives the general drift argument and complete class
classification. The [receipt](results.json) retains all 7,776 original
attempt records, sector witnesses and source bindings. Physical thermal
abundances require a separate bath, temperature and clock attachment.

```sh
python3 -m pytest -q model/native_thermal_admission
python3 -m model.native_thermal_admission --write
python3 -m model.native_thermal_admission --check
```
