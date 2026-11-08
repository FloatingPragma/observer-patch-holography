# Local energy and the counted transfer clock

The fair binary source supplies a positive transfer, an interacting exchange
generator and exact local energy currents. On its specified faithful
reference, the counted transfer Hamiltonian converges to a common local
exchange clock for fixed local observables and complete fixed finite detector
menus, with explicit finite-size bounds.

Unread repair has a separately measured positive energy input. The
[proof](PROOF.md) retains this distinction, preparation assumptions and the
additional requirements for relativistic stress and physical time.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_local_energy_clock.packet --check
python -m pytest -q -W error model/native_local_energy_clock
```
