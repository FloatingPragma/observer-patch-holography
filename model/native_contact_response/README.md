# Two-body contact correlations

The native repair source has a finite two-body mode below the threshold of
its independent-pair comparison. Complete repair records observe a strictly
positive share of that mode. Its contact term and spectral weight are
derived from the original source law. Connecting the cell to the larger
source produces a proved leakage, so the isolated-cell mode requires a
dressed spatial construction before it can represent a persistent particle.
The [proof](PROOF.md) gives the response, its readout and the persistence
criteria without assigning a physical mass.

```sh
python -m model.native_contact_response.packet --check
python -m pytest -q -W error model/native_contact_response
```
