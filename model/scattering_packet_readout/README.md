# Reading an interacting scattering phase

A finite source experiment can distinguish the moving two-particle contact
interaction from free evolution. An explicit incoming packet, collision time
and lattice cutoff connect the scattering theorem to a complete interference
measurement. Every preparation failure, detector outcome and export failure
is retained. The [proof](PROOF.md) gives a probability separation greater
than 3/10 and a finite resource schedule, although that schedule is enormous.
The fine-spacing experiment is not claimed executed.

A separate five-qubit control executes the original source pulses on all
32 quantum rows and checks both time orientation and the free control.
Exact port certificates bound its literal compilation errors; complete
numerical arrays are regenerated during replay. The [receipt](results.json)
separates this execution from the analytic scattering and cutoff theorem.
The result concerns a nonrelativistic two-particle family; physical species,
a common many-body vacuum and Einstein dynamics require additional results.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.scattering_packet_readout --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/scattering_packet_readout
```
