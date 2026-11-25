# A magnetic spin doublet on one charged orbital

A vacancy on the twelve ports of a guarded Interval Echo cell has a two-dimensional lowest band when a scalar phase service supplies one quantum of magnetic flux through the sphere. On this band, perpendicular spatial half-turns that preserve the applied field anticommute. The detector sends every accepted launch to its dark exit when it compares the two orders. The source has one vacancy and no two-component spin register. The half-turns combine port permutations with the phase corrections independently determined by the link field. Guarding the cell and supplying the flux are essential experimental conditions. [Protocol](PROTOCOL.md) · [Complete records](results.json) · [Proof and compiler](PROOF.md)

| Flux quanta | Lowest band dimension | Gap in contact-generator units | Pairs giving negative order sign |
| --- | ---: | ---: | ---: |
| 0 | 1 | 1.381966 | 0 of 15 |
| +1 | 2 | 1.763356 | 15 of 15 |
| −1 | 2 | 1.763356 | 15 of 15 |
| +2 | 3 | 2.000000 | 0 of 15 |

The applied phase is equal on all twenty oriented faces. Local gauge corrections follow from the requirement that each actual half-turn return every applied link coefficient to itself. Changing the phase convention or the face used to set up the calculation gives the same order sign. Bare spatial turns give the positive sign at all four settings. Because the compensated half-turn order has the same sign on the entire twelve-port space, the independent spectral measurement is needed to establish that the source occupies a doublet. A controlled order measures a relative phase; a single half-turn squared also depends on the chosen overall phase of that operation.

Two sharp input addresses reach the doublet using twelve coherent accept-or-refuse clocks, with a reference frequency set to the calculated lowest energy. Each clock uses the finite native contact sequence. The total admissions are 16.642% and 16.649%; their conditional overlaps with the exact doublet exceed 99.995%. The score includes the first refusal at every appointment and both phase exits at every terminal port. Each acquisition has 262,144 simulated starts. On the supplied exchange-strength-three control the gap is 5.290067 while the one-quantum spatial sign is preserved. [Complete detector and phase service](measure.py)

## Which particles could inhabit the doublet?

The one-quantum twelve-port rotation action contains three inequivalent bands of dimensions two, four and six, each exactly once. Only one copy of the lowest spin representation fits in this sector. In a common root-face-zero phase convention the heavy, light and middle intensity modes of the independent Koide detector have respective ground-doublet weights 13.599%, 14.854% and 21.546%. A coherent combination of the three modes is entirely dark to this doublet. The three distinct magnetic energies therefore do not constitute copies of one charged-lepton spin representation. A different charge register or additional material degrees of freedom could give other multiplicities. [Charged three-line detector](../../../interval_echo_lepton_register/README.md) · [Integer band certificate](bands.py)

Disable the background service and open the guards to the twenty-cell federation. At one native time unit the extremal return probabilities of normalized doublet inputs are 48.889% and 77.081%. An independent 84-contact word exchanges two matching guarded cell templates, preserving every transit address. Aligned doublets give a positive permutation sign in 262,144 of 262,144 starts. The antisymmetric orbital control gives a negative sign at the same contacts. This experiment separates the magnetic spatial spin read from fermionic exchange of the composite.

For the ordinary single-defect free law the material current across a written edge equals a difference of its private coin. A phase multiplied onto its hopping amplitude can be absorbed by a rephasing of that coin. The full spectrum of a fully written single-defect board is therefore invariant under an independently assigned twist on every edge. On a guarded written triangle the two signed preparations give unchanged port probabilities under a face twist of $2\pi/20$. Blank contacts give an absolute change of 0.031973 in a port probability under the same diagnostic external twist. Coherent rearming of a collar, a second mobile defect, ring exchanges and direct phase-reference probes use different operations. [Exact local columns and receiving records](passive.py) · [Rotation proof in Lean](../../../../../Lean/SourceLaw/MagneticVacancy.lean)

The Python function `charged_source(owner=0)` returns the two columns, the projector and the compensated half-turn pair. The two columns are the Gram–Schmidt frame of the spectral projections of sharp ports 0 and 1. All addresses are local. An experiment in a different cell must keep the link and phase chart and the physical local phase routes consistent. The complete registers, preparation and readers can be replayed from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.flux_doublet.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/flux_doublet
```
