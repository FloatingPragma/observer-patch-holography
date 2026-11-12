# Energy-filtered charged arrivals

The actual e^c energy-bin preparation has a strictly positive native
charge-arrival probability at time zero after its energy filter.
An additional local atomic-energy tag makes the first projector tail
computable from the same native pole residues. Complete source bounds
retain preparation errors, all reference systems and all seventeen
observer outcomes.

The time-zero signal is virtual dressing resolved by measurement. It
does not demonstrate motion. The same source clock supplies explicit
bounds on subsequent probability changes, with the initial tail retained
as the baseline. The Hamiltonian, hopping coefficient, charge reference
and ground preparation are unchanged.

See the [proof](PROOF.md), [receipt](results.json), and
[standalone paper result](PAPER_RESULT.tex).

```sh
python -m model.openwave_energy_resolved_transport --write
python -m model.openwave_energy_resolved_transport --check
python -m pytest -q model/openwave_energy_resolved_transport
```
