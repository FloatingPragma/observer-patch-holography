# An observer can reopen a written contact

One native echo on the empty written contact (152,166), before the declared pair acquisition, erases its neutral memory and restores its native free interaction. The same echo on an initially blank, permanently severed comparator writes a neutral memory while its external switch keeps that contact off. Both use the original two loads, acquisition and complete 31-outcome detector. The predictions below are ideal finite-instrument probabilities, including source exhaustion. No detector trials are sampled. [Design and comparator](DESIGN.md) · [Calculation](experiment.py) · [All outcome probabilities](results.json)

| Strength, exposure 37 | Native one-echo home probability | Permanent-cut home probability | Total variation of all 31 outcomes |
| --- | ---: | ---: | ---: |
| 0 | 0.518592986107 | 0.995196057737 | 0.565097070335 |
| 1/2 | 0.851854289714 | 0.995945068415 | 0.436551931551 |
| 8, the registered strength | 0.995714986100 | 0.996109885986 | 0.065691080528 |

Two native echoes restore the neutral memory, giving the no-echo distribution. Those equal arms follow from the identical reduced operators and are propagated once. At strength 8 the full address distribution distinguishes the interventions much more strongly than home occupancy alone. The primary design gate, total variation above 0.01 at strength zero and exposure 37, passes. A raw-board pair graph and a separate Chebyshev propagator cross-check that appointment; their source acquisition vector comes from the declared protocol. Exposure is the supplied appointment coordinate, with no independent proper-time calibration. The calculation is a retrospective supplement to the row-04 package.

There is also an exact distinction at every real strength. Put $F=364/729$ and $G=\sum_{j=0}^{7}F^j$. For the unchanged, unconditioned accepted source, erasure increases the coefficient of squared exposure in quiet probability by $49G/729=0.133727901844\ldots$. Source exhaustion is the same constant contribution in both arms. Every possible undetected kinetic destination is retained in this coefficient; diagonal contrast does not contribute. This establishes a local dynamical response to the memory operation without an all-time loss claim.

## The eleven marks cannot translate freely

The seam-memory projector onto all zero-contrast labels has dimension 18, retaining both coins and all nine endpoint words. The exact native local generator annihilates that subspace for every strength. Hermiticity gives the reverse zero product. Every other term acts as the identity on that seam's memory register, so the full free generator commutes with its projector. A neutral mark at a particular seam therefore stays there at every free exposure, including when the loads interact with other seams or carry quantum references. This is an operator statement on the declared free Hamiltonian; addressed echoes are interventions outside that flight.

The eleven-memory material cannot translate all eleven marks to another location under this free law. Loads can move relative to the fixed marks. The intervention distinguishes immutable external cuts from physical memory-controlled contacts, while a programmable material trap with the same operations predicts the same records. The particle-versus-trap interpretation and OpenWave candidate decision are unresolved. Persistence is the literal row-04 test; this obstruction imposes no additional criterion on it.

## Reproduction

The first command recomputes the exact memory and short-exposure certificates and the primary appointment, comparing it with the archive. The full replay covers all nine appointments and takes longer. Verification preserves the archive; `--record PATH` explicitly creates a first record and refuses an existing path.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_04.memory_intervention.experiment --headline
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_04.memory_intervention.experiment
OPENBLAS_NUM_THREADS=1 python -B -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_04/memory_intervention
```
