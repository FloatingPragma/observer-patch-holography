# Charged spin excitations on the selected-field source

The declared three-dimensional spin Hamiltonian has an isolated lattice
quasiparticle band for sufficiently small positive `J/B`, by an external
weak-coupling theorem with its source hypotheses checked. Opposite
sublattices give conjugate internal-charge branches. Exact source matrices
fix the first nonzero hopping at second order. A complete local instrument
creates charge while retaining failure. A finite residual bound certifies
a short native preparation of its interacting ground. [The proof](PROOF.md) distinguishes these
statements from a physical particle identification.

| Source creator | Conserved internal charge | Spectral readout | Physical species |
| --- | --- | --- | --- |
| Local spin flip on the even sublattice | -1 | The isolated band `E(k)` in the proof, with nonzero original overlap at sufficiently small `J/B` | Unassigned |
| Local spin flip on the odd sublattice | +1 | The same folded band by exact source symmetry | Unassigned |

Both creators are even in the original CAR grading. These are local spin
excitations; a protected topological-defect interpretation, physical
electric charge and Lorentz spin have not been derived. The conjugate
branches have equal source energies, without identifying a measured mass.

The same ground-relative source has gap at least `B/8` and propagation
speed at most `72J`. At physical spacing `a`, bounded rest energy and bounded
`J/B` force the physical propagation speed to vanish as `a` tends to zero.
This excludes the specified local massive relativistic refinement of this
regime. It does not exclude the separate ordered-spin limit, a finite-cutoff
comparison, or a different justified source dictionary.

The replay uses the literal `native_ordered_ground.source.matrix`, all
virtual states needed for the exact second-order identities and the full
onsite Fock completion of the abstract instrument. The literal native pointer
agrees on the protected spin code; its out-of-code extension is different.
The replay uses no measured masses.
Neither an explicit infinite-volume smallness threshold nor a finite-coupling
dispersion remainder is certified. On the eight-site cube at `J=1,B=400`,
the product-to-ground half-trace error is at most `1/197`. The replay
executes all 167 preparation, pointer and feedback pulses on 131072 work
rows, retaining 608 native first-failure slots. The analytic compiler bound
keeps the complete-output error from the exact-ground instrument below
`1/100`. Floating execution diagnostics check the implementation; they do
not replace the analytic error bound. This finite certificate does not
evaluate the infinite-volume smallness threshold. The longer adiabatic
plan at `J=1,B=100` is retained for comparison and is not executed.
No physical species is
assigned and no absolute or relative physical mass is predicted.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/high_field_particle
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.high_field_particle --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.high_field_particle --check
```
