# A public material read of core erasure

A native contact transfers the core memory’s blank-or-written value to the position of the orbital load. On the 20-cell board, after one attempted erasure and 64 units of native flight, the observer read a blank core on **74,603 of 4,028,408** admitted launches. On 80 cells the corresponding count was **74,590 of 4,028,686**. Each setting began with 4,194,304 simulated launches. [Native detector and routing](measure.py) · [Compressed joint public records](results.json.gz)

After the terminal flight, measure the position of the single load. The measured position selects a path of blank seams to the first core endpoint. Every routed π contact is a complete native swap on the blank one-load sector at κ=8. The paths avoid the core seam and cover all 240 or 960 possible addresses. If the load starts on the second endpoint, the detour reaches the first without using the core contact. The receiver uses the same registered read apparatus after either the native or the exchange-only flight.

The material input to the core is thus (1,0), independent of the position record. For blank memory the native generator on (1,0),(0,1) is 8I+(I−SWAP)/2. Evolving for π moves the load to the second endpoint. For the vacuum-written memory the generator kills both inputs and the load stays at the first endpoint. The measured effect for a flip is exactly the blank projector on this two-memory subspace; the other public load position reads written. The program checks the complete 99-state map, all nine public outcomes, and the off-diagonal elements of the two-memory effects. These facts apply also to a coherent superposition correlated with the orbital position. The initial sharp position read can disturb the orbital and preserves the marginal memory probabilities. This detector covers the blank and vacuum-written subspaces of the enrolled source. The clock, address read, controller, isolated routing and terminal occupations are supplied.

| Operation on the core before native flight | Flips / admitted on 20 cells | Flips / admitted on 80 cells |
| --- | ---: | ---: |
| Leave the write in place | 0 / 4,028,039 | 0 / 4,028,629 |
| One native echo | 74,603 / 4,028,408 | 74,590 / 4,028,686 |
| Two successive native echoes | 0 / 4,027,993 | 0 / 4,028,828 |
| Supplied complete erase | 4,028,677 / 4,028,677 | 4,027,977 / 4,027,977 |

For one echo the simultaneous 99% intervals on the admitted blank fraction are 1.8259–1.8782% and 1.8255–1.8777%, respectively. Omitting the terminal contact gives zero flips even on a blank core; a half-time contact registers half of the blank population in the ideal model. The receipt includes both flight laws, delays 0 and 64, all four interventions and all three receiver settings. All position/terminal-word counts are retained with explicit sparse-zero semantics and a source-refusal bin. No observer has direct access to the private collar. Every bank is drawn from integer weights on the 2^-40 grid of the [release sampler](../measure.py), with source refusal as the exact complement; the largest recorded L1 distance of a sampled law from its computed probabilities is 4.7×10⁻¹⁰.

From the RER root run `OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_31_orbitals.native_core_release.public_readback.measure`.
