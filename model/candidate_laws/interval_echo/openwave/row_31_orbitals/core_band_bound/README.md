# Which orbital levels can an empty written core support?

For one amount-one load on a blank host, **every pattern of stationary empty writes excludes the energy interval (−3,0)** in native exchange units. For an amount-two load the excluded interval is (−27,0). This holds on each registered federation, irrespective of core shape, number of writes or refinement. Energy is measured from the uniform free load on the same board. For a localized core this is the blank-host dissociation reference. An extensive patterned medium uses the same reference as a comparison. [Exact local identities and finite-board checks](measure.py) · [Enclosures and multiplicities](results.json)

Echoing an empty seam in vacuum prepares a definite written memory and leaves its endpoints empty. Vacuum allows the writes to share endpoints. Inject one load of amount a in {1,2}. On a blank seam, conjugating the registered operator on all 99 states gives half the exchange Laplacian and the native mismatch cost 8a² at an occupied endpoint. On a written seam the entire local generator kills each relevant material input (empty, left occupied, right occupied), including all amplitudes into other memory states. The verifier checks these identities entry by entry over the exact field ℚ(√2). Empty writes on other seams leave the local calculation valid. The registered degree term compensates for a port's missing neighbors in the blank reference.

Let s(x) count written seams at port x. Let L be the graph Laplacian of the remaining blank seams and d(x) the original degree. After subtracting the uniform blank energy 48a², the one-load Hamiltonian is

$$H=L/2-8a^2\operatorname{diag}(s).$$

The inequalities 0 ≤ L/2 ≤ diag(d−s) follow by summing squared differences and squared sums of neighboring amplitudes. Thus H is between the two diagonal matrices −8a² diag(s) and diag(d−s−8a²s). Since d≤6, the latter is bounded above by diag(6−(8a²+1)s). Both outer diagonals sort by decreasing s. The eigenvalue min–max principle bounds each ordered eigenvalue. **Exactly as many levels belong to each band as there are ports of that written degree.**

| Written seams at a port | Amount-one energy band | Amount-two energy band |
| ---: | ---: | ---: |
| 0 | [0, 6] | [0, 6] |
| 1 | [−8, −3] | [−32, −27] |
| 2 | [−16, −12] | [−64, −60] |
| 3 | [−24, −21] | [−96, −93] |
| 4 | [−32, −30] | [−128, −126] |
| 5 | [−40, −39] | [−160, −159] |
| 6 | {−48} | {−192} |

The number of negative levels equals the number of touched ports. Adding more writes can add bound states, while none enter the excluded interval. As the distinct memory patterns form orthogonal invariant blocks, a coherent superposition can mix their energy measures but cannot fill their common gap during free evolution. Multiplicities in the table are conditional on a definite written pattern.

For a hydrogenic sequence anchored to the measured ground binding at the blank threshold, Eₙ=E₁/n². The lower bounds E₁≥−48 and E₁≥−192 imply that every amount-one prediction for n≥5, and every amount-two prediction for n≥3, lies in its respective gap if there is a bound ground. In particular increasing the number of empty written seams with board size cannot yield a sequence accumulating at the blank threshold. Time-dependent writing, populated memory states and additional mobile matter fall outside this argument. The OpenWave criterion asks whether the candidate produces hydrogen-like levels; this result settles this stationary-core construction.

The executable checks 70 full-cell spectra and verifies the combinatorial bounds on representative blank, single-seam, intersecting, intercell and completely written patterns on 20, 80 and 320 cells. Run `OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_31_orbitals.core_band_bound.measure` from the RER root.
