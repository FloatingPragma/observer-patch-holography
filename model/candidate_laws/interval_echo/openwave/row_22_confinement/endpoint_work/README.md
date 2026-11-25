# How much work can two local endpoints force?

An observer prepares the same extensive texture with an endpoint in two different owners. The spectrum of each local term is measured on fresh copies. The change in mean energy of the entire board agrees with the change measured at the exchanged owners' boundary in all 18 comparisons. Separate minimum-energy inequalities give an upper bound independent of separation when the admissible states differ only by relocating an endpoint. The [face-parity experiment](../parity_sources/) tests what happens when the surrounding faces are also constrained.

## Two inequalities for relaxed endpoints

For a common family of normalized states let E(H) be its attained minimum energy. Add any two bounded source terms A and B. For each state the expectation of A lies between its least and greatest eigenvalues. Evaluating both minima on one another's minimizers puts E(H+A)−E(H) in this same interval with or without B present. Subtraction and interchange of the two sources give

    |E(H+A+B)−E(H+A)−E(H+B)+E(H)| ≤ min(width(A), width(B)).

This covers arbitrarily many other loads and correlated memories allowed to relax. It places a bound of 33 exchange units on the connected ground work when either additive source switches a single native seam at contrast eight. It does not impose a bound on changes in trial-state energy between different preparations.

For hard constraints on an endpoint let a transposition exchange the twelve port spaces and all thirty internal seam memories of the old and new owners. Other owners and external seam memories retain their addresses. Suppose this permutation takes the *entire* family of allowed states with the old endpoint bijectively onto that with the new endpoint. The other source occupies a third owner. In the sum of local Hamiltonians the exchanged internal terms cancel exactly. Every external term incident on either owner can change its expectation by at most its spectral width. In the written record basis its exchange part is a parity projection commuting with the squared contrast. Blank records have zero energy there. The full native spectrum is contained in {0, 8, 9, 32, 33}. At most six external seams touch the two owners. For a local port label present as a gateway on just one owner, a private reference cost can change by at most 32. There are at most six such labels. Hence every trial state's energy changes by at most

    33 × incident external seam count + 32 × differing gateway label count ≤ 390.

Applying the trial comparison in both directions bounds the absolute difference of the *two constrained minima* by the same number. A common gateway mask gives an upper bound of 198. Subtracting a second independently varying one-source minimum costs up to a second such bound. Physical transposition during a flight is unnecessary. [Lean proofs of the variational steps](../../../../../../Lean/Screen/EndpointGroundWork.lean)

The bound is uniform in board size and occupancy. It allows finite intervals of linear growth and permits a confining sector whose flux constraints change throughout the interior. It provides neither a colour assignment for the endpoints nor a physical conversion to GeV/fm.

## What the observers read

A 64-position phase detector measures a uniformly chosen native term or a known private port cost on each independent preparation. Multiplication by the number of eligible terms gives the respective mean total energy. Each nonempty cohort comprises 2^42 simulated categorical draws. The two independent work estimates use all terminal bins. The enrolled joint 99% Hoeffding bounds cover both acquisitions and the zero-contrast controls. The term selector, controlled timing, sharp load preparation and Born probabilities are supplied. Vertex distances below count owner steps from the fixed source. [Protocol](PROTOCOL.md) · [Measurement program](measure.py)

| Owners | Distance | Local labels | Whole-board work at contrast eight | Boundary expectation |
| --- | ---: | --- | ---: | ---: |
| 20 | 4 | Matched | 64.995 ± 0.059 | 65 |
| 20 | 5 | Differ | −32.480 ± 0.059 | −32.5 |
| 80 | 10 | Matched | −32.517 ± 0.235 | −32.5 |
| 80 | 11 | Differ | 65.004 ± 0.235 | 65 |
| 320 | 22 | Matched | 0.005 ± 0.939 | 0 |
| 320 | 23 | Differ | 32.431 ± 0.939 | 32.5 |

The homogeneous endpoints above have zero net private-reference work. A separate positive control fills precisely the unmatched gateway ports before relocation. At 23 steps the whole-board work was **159.557 ± 0.939**. The incident seams contributed **−32.500** and the unmatched private terms contributed **192.000**. Omitting these private readings gives the wrong sign. Setting the contrast to zero removes their contribution on all three boards. These measured energies belong to sharp, textured preparations. The minimum-energy statements above follow from the variational proof. [Original counts](records.json) · [Private-reference control counts](private.json) · [Independent sharp-state tests](test_work.py)

Run from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_22_confinement.endpoint_work.measure
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_22_confinement.endpoint_work.measure --private-control
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_22_confinement/endpoint_work
(cd Lean && lake build EndpointGroundWork)
```
