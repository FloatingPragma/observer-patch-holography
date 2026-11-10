# A local detector for the native quantum current

A five-state pointer reads the original gauge-invariant current through an
explicit coupling at one source cell. The source Hamiltonian keeps evolving
during this short pulse. An exact 400-state construction and a full-source
error bound preserve the nonzero negative response of both prepared
excitations, with every preparation failure and pointer outcome retained.

The construction declares a local polynomial-current/pointer coupling. Its
short duration requires a large, explicitly bounded control strength and
switching work. The preceding preparation and energy-bin services retain their
matrix-control assumptions. No autonomous repair compiler or particle mass
identification follows. The [proof](PROOF.md) gives the source clock, complete
instrument and resource bounds; [results](results.json) bind the executed code.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -W error -m model.native_local_current_readout --check
.venv/bin/python -B -m pytest -q -Werror model/native_local_current_readout
```
