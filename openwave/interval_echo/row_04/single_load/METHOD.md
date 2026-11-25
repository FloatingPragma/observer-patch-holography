# M10.5 · One load released on an open federation

On the blank one-load sector the free Hamiltonian in exchange units is

$$
H=\tfrac12 L+6\kappa I,\quad
B=Q\exp(-iL_{\mathrm{cell}}/2),\quad
\rho=\sum_{m=0}^{47}r_m B^{1024}|j_m\rangle\langle j_m|(B^{1024})^\dagger .
$$

Here the row reads whether a prepared particle state persists in free evolution. The candidate is one material excitation in a bounded acquisition. The home cell has twelve ports and three outward gateways. During an attempt the operator Q removes the three clicked gateway components. A click places the load at a known port and the operator routes it by swaps to the next sharp launch. The incoming weight of attempt m is the product of previous failure probabilities. The unnormalized density above therefore includes every accepted attempt at its actual weight. Its trace deficit is the exhausted-source bin. During a pulse the contrast is 5κ times the identity. During flight it is 6κ times the identity. Both phases cancel from observed probabilities.

For home port i at exposure t the observable probability is the diagonal element of the freely propagated density at i. Two compatible unsigned seam reads at i resolve that port. A geometric cell can report twelve home ports, accepted exterior, or exhausted source. The native echo, copied seam records and inverse echo implement these probes on the one-load sector. One terminal read uses a fresh ensemble. The external exposure uses exchange units.

## Why the excitation can persist

Let A be internal adjacency and G read the gateway coordinates. The kernel of the stacked matrices G, GA, GA² and GA³ is invariant under A because its minimal polynomial divides (A−5I)(A+I)(A²−5I). Each vector in the kernel has zero gateway amplitude. Extending it by zeros gives identical internal flight on the fully connected federation. Its projection survives every clear query. The protected part of a sharp launch j has weight P(j,j), where P is the rational orthogonal projector onto that kernel. Summing these contributions with lower bounds on the incoming attempt weights gives the unconditional floor. The lower bound remains valid for every exposure. A port measurement after preparation and random swaps during flight each change the subsequent home records.

| Evidence | Result |
| --- | ---: |
| Certified floor at every free exposure, all three boards, any complete port sweep order | 0.999940596796 |
| Lowest coherent home probability, 11,520 finite 20-cell settings | 0.999940596806 |
| Highest random-swap home probability on that grid | 0.348844205608 |
| Smallest coherent advantage over a release-port measurement on that grid | 0.342859559000 |
| Historical cell 13, exposure 37, coherent / dephased / uniform sharp / random swaps | 0.999941 / 0.570144 / 0.388503 / 0.151119 |

The grid covers every home, every cyclic starting address in each direction and all 24 integer exposures from 15 through 38. An independent raw-edge implementation checked one prespecified setting per cell, including exterior population and exhaustion. The largest error among its 1,120 comparisons was below 1.3e-12. All 645,120 development bins are retained with the source. The exact coherent bound covers all 420 cells of the three registered boards. The control bounds of 0.916448 after dephasing and 0.559219 for random swaps apply on the 20-cell confirmation window. Each arm uses the same home detector.

## Reproduction and scope

The platform script builds the cell and full-board Laplacians from the original federation file. It uses a dense matrix exponential, a matrix power for the 1,024 clear rounds and explicit retry weights to check the 56 historical pilot bins. The frozen author source also checks native echo columns, resource limits, integer enclosures and all prospective settings. These are simulation-generated record probabilities. The one-load source filter and the localization target were known during design. The scoped claim is free persistence after the apparatus is removed. A species, electric charge, spin, inertial mass and physical time scale require additional observables. The first official M10.5 task commit still determines the prospective replay. Its outcome belongs in a subsequent receipt.

| Formula or test | Implementation |
| --- | --- |
| Native gauge of the one-load generator, signed queries and observed home ports | [Gateway preparation at the public parent](https://github.com/FloatingPragma/observer-patch-holography/blob/09f062f77f7f540f2783f128ae6314c5d8b818ea/model/candidate_laws/interval_echo/openwave/row_04_compact/gateway_acquisition.py) |
| Invariant projector and exact incoming-weight bound | [Integer certificate](https://github.com/FloatingPragma/observer-patch-holography/blob/09f062f77f7f540f2783f128ae6314c5d8b818ea/model/candidate_laws/interval_echo/openwave/row_04_compact/gateway_bound.py) |
| Finite development grid and historical pilot | [Grid summary](survey.json) · [Pilot receipt](pilot.json) |
| Raw board to the four pilot arms | [Platform script](openwave.patch) · [Public board](https://github.com/FloatingPragma/observer-patch-holography/blob/09f062f77f7f540f2783f128ae6314c5d8b818ea/model/executable_model/carriers/FEDERATION_LEVEL0.npz) |
| Held selection, gates and time order | [Versioned task](../../../../model/candidate_laws/interval_echo/openwave/row_04_compact/m10_5/TASK.md) |

From the OpenWave checkout, give the script the public board file.

```sh
python openwave/xperiments/m10_oph/research/scripts/m10_5_one_load.py --board /path/to/FEDERATION_LEVEL0.npz
```
