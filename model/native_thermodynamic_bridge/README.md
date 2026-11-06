# Native repair and finite thermodynamics

On the supplied binary total-two sector, the actual fair-coin endpoint repair
is exactly uniform conditional resampling on each seam fibre. All thirty
native kernels preserve one uniform reference. Their connected configuration
graph makes that common stationary law unique for positive random scans and
for repeated covering scans. The same reference maximizes state entropy, and
each native row maximizes conditional transition entropy on its allowed fibre.

The [proof](PROOF.md) derives the exact entropy increase for fixed schedules,
the additional mixing-information term for random scans, and the supported
fluctuation identity. Exact rational logarithmic certificates check these
identities without floating entropy calculations. The complete endpoint and
coin/current tape retain the initial classical information.

The binary-sector restriction matters: one source repair of the uniform
78-state startup law decreases entropy by `log(3)/26`. Physical energy,
temperature, time, erasure costs and cooling remain separate. The uniform
reference has constant modular energy.

```sh
python3 -m model.native_thermodynamic_bridge.certificate --check
python3 -m pytest -q -W error model/native_thermodynamic_bridge
```
