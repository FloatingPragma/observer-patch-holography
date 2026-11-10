# Native charge arrival and energy

The same finite native charged source supports a position record followed
by a finite energy readout. Arrival at time 1/4 and acceptance near the
actual admitted charged ground has unconditional probability above
4.50e-6 at the apex and 1.13e-6 at an adjacent boundary site, with every
preparation and detector failure retained. A coherent spatial preparation
followed by the energy clock has accepted ground infidelity below 2.14e-8.

The [proof](PROOF.md) keeps the reference cost, measurement backaction and
actual Hamiltonian explicit. The ground belongs to a selected native
spectral bin; the full charged-sector ground and a relativistic particle
interpretation are outside this result.

From the repository root:

```
python3 -m model.native_abelian_energy_records --check
python3 -m pytest -q -W error model/native_abelian_energy_records
```

The [provenance](PROVENANCE.md) specifies the fresh source replay and
immutable research comparison used by the [receipt](results.json).
