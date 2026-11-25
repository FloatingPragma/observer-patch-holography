# OW-22: separating colour sources on the intact federation

**Status: HOLD.** The simulated receivers resolve an attractive colour-dependent contact. When the sources are separated by two to 22 owner steps, their connected energy is exactly zero in the tested fixed-source space. A fourth, freely hopping load mediates measurable attraction at two and three steps. Its additional binding approaches zero at long separation. A different preparation gives a linear energy cost of **17 native energy units per changed cell** in an intact polarized domain. The spectators that select the bulk phase account for the cost. OW-22 asks for a linear inter-charge potential near 1 GeV/fm. Neither the identification of that potential nor an independent conversion to those physical units is established here. [Frozen protocol](PROTOCOL.md) · [Complete records](records.json)

## What two separated receivers measure

A holds one binary load in the lower colour triplet of a native cell. B and C are adjacent owners, each holding one such load, with their joint colour in the antisymmetric antitriplet. These are the measured one-cell triplet modes used in the three-load colour experiment. The source here prepares their joint state explicitly. Its two settings are the normalized epsilon tensor and an orthogonal octet state formed by cycling A's three colours. They have the same single-source marginals and the same BC marginal. Addressed mode projectors identify the three copies, total load and local populations. An exchange of the two colour labels tests BC antisymmetry. The joint determinant projector distinguishes the two preparations. At a contact the laboratory frames identify the two contacting port projections with the same first colour. These choices define the colour comparison; physical gauge transport between distant frames is a separate question. [Three-load source and reader](../row_13_baryon/README.md) · [Native current construction](../../currents.py)

The energy measurement uses separate copies. Each of the three external seams leaving A can report swap-odd, both occupied, or the remaining outcome. A controlled swap with a reference arm is a parity analyser; the commuting double-occupation flag splits the even result. The registered local unitary at a π pulse realizes the swap on the blank seam for all three contrasts used here. For each seam we compare preparations ABC, A alone, BC alone and blank. The local density operators on the other seams cancel in this combination, including on the bridge inside BC. Each setting keeps the actual federation and its outside ports. There are 134,217,728 simulated shots per seam setting, including quiet trials. Counts at different settings come from independent copies. Source control and calibrated reference arms are apparatus. [Native binary-sector generator](../../../contrast_memory_spherical_v1/sectors.py) · [Simulation](measure.py)

| Measurement on 320 owners | Exact connected energy | Simulated readout, one standard error |
| --- | ---: | ---: |
| Singlet at one owner step | 0 | −0.000015 ± 0.000049 |
| Octet preparation at one step | −17/96 | −0.176961 ± 0.000149 |
| Both, two through 22 steps | 0 | Largest absolute reading 0.000088 |

For example, a straight line through the first two octet readings predicts +0.177009 at distance three. The held detector gives +0.000004. The two smaller federations show the same contact and plateau through four and ten steps. On a shared bridge the local port marginals are each 1/12. Double occupancy of the two aligned ports is forbidden by the singlet and has probability 1/96 in the specified octet. The connected expectation is (-(1+2\kappa)\langle n_i n_j\rangle). Exchange alone gives an octet contact of −1/96. Reversing contrast gives +15/96. The term counts permit all three reductions. Independently, signed-quadrature pulses on the full 99-state seam distinguish the laws. On the symmetric one-load input the respective bright-Y probabilities for contrasts +8, 0 and −8 are 0.146447, 0.5 and 0.853553. [Shot records, including pulses and held forecasts](records.json)

## Does relaxing the colours change the potential?

For the fixed-source calculation, allow every state in the nine-dimensional space of A's triplet and BC's antitriplet. Compress the full generator to this space with one load in each named owner. The single-cell energies include their three connections to the rest of the federation. A hop between two owners leaves the fixed-number subspace. The compressed odd projector on an external bridge consequently equals ((n_i+n_j-2 n_i n_j)/2). Subtracting the disconnected operator for the same source positions leaves only (-(1+2\kappa)n_i n_j) on each A-to-BC bridge. At distance two or more there are no such bridges. This proves equality of the joint and disconnected operators, hence also of all their eigenvalues, even after minimizing over arbitrary admitted colours. The paired BC self-energy stays in both operators. [Full compression](spectrum.py) · [Independent comparison with explicit swaps on three-load bit words](test_confinement.py)

At contact the interaction has operator norm 17/32 on all three federations. On the 320-owner board the minimum falls by 0.271235 when all nine colours can vary, by 0.302338 when the octet flag is held and by zero for the one-dimensional singlet. The first table gives first moments of prepared states. These subsequent numbers are constrained minima. Enforcing the source triplet modes and their locations would require external clamps or repeated projections. This calculation leaves relaxation into other spatial orbitals and sectors with extra material open.

There is a size-independent bound even when the three binary loads move freely in blank memory. In unshifted units their native Hamiltonian is (144 I + K -16 A), where (A) counts occupied seams. Three distinct loads have at most three occupied seams. The exchange operator (K) is half a graph Laplacian on three-load configurations with maximum degree at most 18, and therefore (0\le K\le 18 I). All three-load energies lie in **[96, 162]** on every carrier of maximum port degree six with the registered reference completion. This excludes asymptotically unbounded positive energy growth in this fixed three-load sector. A finite interval of a dressed potential, or a string carried by a growing number of background degrees of freedom, is compatible with the bound. [Native sector definitions](../../../contrast_memory_spherical_v1/sectors.py)

## Can a mobile load mediate attraction across empty cells?

Yes, in a finite dressed sector. One additional binary load hops freely among all 924 ports outside the three source cells on the 80-owner federation. Local term energies and hopping amplitudes come from the registered law. No string or sequence of displaced markers is prepared. The apparatus fixes the source locations and either clamps their joint singlet or admits arbitrary joint colours. In the scalar comparison the missing source colours can be held on inert external references, preserving the same marginal at each remaining source. All four reference runs use the same three blocked owner regions. Native exchange terms at their boundaries retain the diagonal cost of a hop toward a blocked cell. [Dressing protocol](DRESSING_PROTOCOL.md) · [Native compression and phase reader](dressing.py) · [Dressing receipt](dressing.json)

In each arm, subtract its source-only minimum and the mobile load's constant reference energy. Let the resulting ground energy be E(X) for sources X. We simulated a two-quadrature phase measurement on externally prepared lowest-energy states at unit clock time. The reference branch cancels the stated offsets. Each quadrature uses 16,777,216 simulated shots. The additional interaction available from sharing a single mobile load is E(ABC) − min(E(A), E(BC)). The independent calculations include an eigenpair residual bound of 10⁻¹⁰. These are constrained equilibrium energies. The experiment supplies energy selection and does not measure a spontaneous equilibration time.

| Separation on 80 owners | Additional ground energy, singlet held | Independent phase read, one standard error |
| --- | ---: | ---: |
| 2 | −0.109367 | −0.109385 ± 0.000338 |
| 3 | −0.018137 | −0.018420 ± 0.000343 |
| 4 | −0.003046 | −0.003715 ± 0.000343 |
| 5 | −0.000269 | −0.000883 ± 0.000343 |
| 10 | −0.00000000248 | −0.000559 ± 0.000343 |

At two steps, 15.23% of the computed mobile probability lies on the accessible ports across A's boundary, and 15.23% across the BC boundary. At ten steps these weights are 0.00013% and 18.38%. The mobile load can reduce its energy at both separated sources at short range. At long range essentially all of that binding gain is already available from the best isolated source. The ordinary alternating sum E(ABC) − E(A) − E(BC) + E(empty) is still +0.079217 at ten steps, because each isolated reference can bind its own copy of the same single mediator. Extrapolation from one and two steps predicts −0.215805 at the first held distance of three; the measurement gives −0.018420. A descriptive fit to the last five measured placements, selected after acquisition, gives a slope of −0.000001 ± 0.000109 native units per step.

On 20 owners we also let the nine admitted source colours mix and entangle with all 204 mobile positions. The extra ground-energy shift is −0.007888 at two steps, −0.000000484 at three and below 10⁻¹⁰ in magnitude at four. With the singlet held on this same board the two-step shift is −0.109251 at registered contrast, −0.000191 for exchange alone and +0.001767 at reversed contrast. A direct four-load bit-word test checks the compressed Hamiltonian on a coherent source and mobile state. These finite results identify a short-range effect of the mobile load. A persistent flux tube would require evidence beyond the one-mediator sector.

## Where a linear energy cost does appear

In a separate sharp preparation, put one marker in each owner along a geodesic and place a spectator across its side bridge. Compare two degree-six positions of each interior marker. In the aligned phase it touches its spectator. In the other phase it points forward along the route. Flipping a connected interval leaves the number of markers and spectators fixed. Occupancy at independently addressed ports locates both phase boundaries. The energy detector reads the full set of affected exchange and occupancy terms using 1,048,576 simulated shots per seam. A sharp binary word has mean shifted energy (sum_{p\;\mathrm{occupied}} d_p/2-(1+2\kappa) A). Each flipped site loses precisely one occupied seam while keeping its degree. Its energy cost is 17 at registered contrast, 1 for exchange alone and −15 when contrast is reversed.

| 320-owner preparation | Energy result in native units |
| --- | ---: |
| Connected interval with spectators, fit to lengths 1 and 2 | 17.00529 ± 0.00341 per cell |
| Predictions at held lengths 3 through 12 | Largest deviation 1.87 standard errors |
| Same connected interval, spectators removed | −0.00101 ± 0.00371 per cell |
| Same two outer wall records 12 cells apart, interior permitted to realign | 34.00452 ± 0.00763 total |

The last preparation has two additional interior walls. The aligned gap permits each outer wall to belong to a separate one-cell domain. Watching only the outer two records cannot tell whether the entire interval has remained in the excited phase. In the connected preparation the spectators explicitly lift the phase degeneracy by a bulk energy density. These sharp initial-state measurements establish an imposed domain cost. Equilibration, a colour charge on the walls and a source-induced colour flux tube each require different evidence. The earlier written-junction and belt-wall experiments use different preparations. [Belt-wall program](../quark_winding/belt_walls/) · [Written-junction program](../quark_winding/written_junction/)

## Replay

From the repository root with NumPy, SciPy, NetworkX and pytest installed:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_22_confinement.measure
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_22_confinement
```

The default command reconstructs the three native federations, regenerates the detector counts, constrained spectra, one-mobile-load comparisons, controls and held forecasts, and compares them with both receipts. Use `--write` to regenerate it. Scores, physical-unit calibrations and upstream submissions are separate from this research result.

[Relocating a local endpoint on an occupied board](endpoint_work/) compares independent whole-board and boundary energy measurements and bounds the relaxed separation work for local soft sources and transposed hard constraints.

[Two face parities and their static energy](parity_sources/) shows when a string is forced, how additional odd faces screen the endpoints, and how much work two such local pins can require on the fully active half-filled board.
