# Native charge and Gauss flux

The complete interacting eight-cell source transports a negative unit
native charge with a retained positive reference at the origin. On each
of its six actual arrival outputs, the full Gauss constraint forces a
unique one-dimensional link-flux sector. The arrival matter-plus-scalar
charge and its enclosing outward flux both equal minus one. All 256
cuts, both periodic routes to each neighboring site, and every failure
output are retained.

The calculation uses the actual 119-state endpoint matrices and source
Hamiltonian. Central hypercharge telescopes across cuts; a mutation
control shows why the weak contribution to electric charge cannot be
cancelled on arbitrary nonarrival link states. On the attached singlet
arrival sector, electric charge and hypercharge agree exactly.

The [proof](PROOF.md) gives the operator identities and complete joint
readout. The [receipt](results.json) contains exact spectral multiplicities,
source-bound route probabilities and source hashes. The calculation
supplies no electromagnetic continuum normalization or free-electron
identification.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_native_charge_flux --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_native_charge_flux --check
python3 -m pytest -q model/openwave_native_charge_flux
```
