# A charged source emits two energy-resolved neutral pairs into filled cells

Prepare three vacancies and two excesses in owner 0 of the twenty-cell Interval Echo federation. The other nineteen owners start with one unit on every material port. After free flight, owners 6 and 4 can each receive an excess and a vacancy in a locally identified 50-unit orbital, while owner 0 retains a negative charge in its lower orbital. At exposure 128 the three observers jointly record this channel in **3,503,157 of 4,194,304** simulated launches. Their independently calibrated energies are **8.47062, 50 and 50** in exchange units relative to the filled material. The source and the outgoing preparation both have mean energy **108.47068931**, with full energy variance **1.00085617**. Selecting an individual formation line, measuring a weak current, or assigning a lepton lifetime would require additional observations. [Simulated records](results.json) · [Emitter and reader](measure.py) · [Physical routes and pulses](apparatus.py)

## What arrives at the two receivers?

One fourteen-port bench has two routes, both starting in the same owner. On the eight-port route a vacancy occupies a two-port orbital at the left end. A vacancy and an excess occupy the fourth and fifth sites, next to the bridge. These three defects share one connected Hamiltonian. The six-port route holds a second vacancy and excess on sites two and three. The rest of both routes starts filled. Each receiving orbital is on the first two ports across its bridge. Its two possible vacancy–excess orders can have even or odd relative amplitude. The regional reader resolves the even and odd states, and records a third outcome for a miss. The charged reader resolves the two vacancy energy orbitals and its own miss. All 27 joint outcomes enter the protocol. The two routes have separate free generators and the two prepared neutral pairs are part of the initial inventory. [Full 99-column inventory conservation](../neutral_ejecta/README.md)

| Source preparation and free contact | Matching three-way arrivals at 128 | Starts |
| --- | ---: | ---: |
| Lower charged orbital, two odd pairs | 3,503,157 | 4,194,304 |
| Lower charged orbital, two even pairs | 3,199,116 | 4,194,304 |
| Upper charged orbital, two odd pairs | 22,006 | 4,194,304 |
| Reverse the odd-pair experiment | 3,504,171 | 4,194,304 |
| Close the first bridge | 0 | 4,194,304 |
| Halve contrast during flight with fixed source and readers | 30,276 | 4,194,304 |

A charge on the first route changes the emission probability. With the lower charge prepared, 3,800,782 starts at exposure 128 give the same lower charge and odd first pair, regardless of the second owner. Filling the surplus vacancy gives 1,714,922 odd first pairs in an equally large bank. The first pair therefore interacts with its charged surroundings. The second receiver can operate independently. With the first bridge closed it still records 3,866,653 odd pairs. A source with a single vacancy and no seeded pairs yields zero two-neutral arrivals for every free exposure, by conservation of the three individual load counts.

At exposure 256 only 24,324 starts give the outgoing channel that accepted 3,503,157 at 128. Preparing those outgoing orbitals and waiting 128 gives the near-equal reverse count in the table. The fixed real Hamiltonian gives exactly reciprocal transition probabilities between these real source and daughter vectors. These are finite, reversible emissions of initially present material pairs. Their return and the measured energetics distinguish this mechanism from spontaneous heavy-lepton relaxation. [Resonance into preloaded matter](../resonant_delivery/README.md) · [Source whose upper band exceeds the free daughter energies](../daughter_shell/README.md)

## Local energy measurement

Each region calibrates its charged vacancy or its neutral pair on fresh copies of the route with the other tested matter removed. In the vacancy orders on the end bond, the compression of the native generator is
\[
H_C=\begin{pmatrix}8.5&-0.5\\-0.5&17\end{pmatrix}.
\]
On either interior pair orbital it is
\[
H_N=\begin{pmatrix}49.5&-0.5\\-0.5&49.5\end{pmatrix}.
\]
The charged observer estimates the mean energy from the vacancy probability on the second port and the antisymmetric probability on its bond, as \(8+8.5p_2+p_-\). The neutral observers measure \(49+p_-\), where the 49 comes from three sharp squared load differences and two boundary-exchange terms. Their calibrated preparation counts give charged levels 8.47061718 and 17.02933502, and neutral levels 49 and 50. The two neutral levels distinguish even and odd exchange of the vacancy and excess. The difference of 1 is an orbital level spacing.

In both the initial and the resolved outgoing subspaces the complete Hamiltonian compresses to the sum of the three regional operators. The proof checks every integer half-unit matrix entry in the two local product bases, including the bridge bonds. The eigenvectors of this compression are finite packets of the full generator. For the primary preparation the sharp full-energy variance is 1.00085617 in squared exchange units. Thus the calibrated common sum is the full mean for the prepared orbitals, and reading all three orbitals resolves the specified outgoing channel. A single full-Hamiltonian energy for an individual emission would require another instrument. [Native word and compression checks](bound.py)

## Can the neutral pair take the charged orbital’s excitation energy?

The local gap of the charged reader is \(\sqrt{293}/2=8.55862138\). Start with the charge and first neutral pair in their upper levels. For **every flight time**, the probability to find the charge in its lower orbital and the first pair in either calibrated receiving orbital is **at most 0.033790**. This bound includes every possible second-wing record. At exposure 128 this lower-charge channel had zero counts among 4,194,304 starts.

Let \(V\) embed the two daughter modes with the lower charge, \(P=V^*\), \(K=V^*HV\), and let \(E\) be the mean energy of the raised-charge source \(s\). The lowest singular value of \(E-K\) equals the charged gap \(g\). For any free propagator \(U\), expand \((E-K)PUs\) into \(PU(E-H)s+(PH-KP)Us\). The two terms have norms bounded by the initial energy spread \(\sigma\) and the operator norm \(\epsilon=\|HV-VK\|\). Consequently \(\|PUs\|^2\le(\sigma+\epsilon)^2/g^2\). In the four sharp vacancy-order and pair-order words, leakage out of the chamber has Gram matrix \(\operatorname{diag}(1/2,1/2,3/4,3/4)\) at both source and destination. The numbers follow from the two paths out of the pair and the additional path out of the second charged site, each with hopping magnitude 1/2. Rotating the charged basis gives \(\sigma^2=0.74914383\) and \(\epsilon^2=0.50085617\). Other destinations, other pair energies and more extended source dressings have different bounds. [Lean energy-width theorem](EmissionBound.lean) · [Exact compression and leakage test](bound.py)

The bound also covers a less selective detector. Keep the lower charged orbital and accept **every nontrivial neutral word** anywhere in the three-port receiving owner’s active region. This gives a stronger all-time upper bound of **2.556284%**. Direct simulated evolution at 128 gives probability \(2.57\times10^{-8}\), with zero successes in 4,194,304 starts. Extending the neutral chamber treats paths that stay inside it as ordinary free propagation and reduces its leakage. In all six neutral order words the exact twice-Hamiltonian has weighted absolute row ceilings at most \(2(4001/80)\), using weights 200 on the near adjacent pair, 5 on the separated pair and 1 on the far adjacent pair. The only neutral escape hop crosses the receiver entrance with amplitude magnitude 1/2. The independent charged escape contributes \((1-17/\sqrt{293})/8\) to the squared leakage. [Exact matrix, six receiving words and all exposures](results.json)

For any length of homogeneous receiving path terminating at one open end, with the same charged window, entrance and initially filled receiver, the all-time ceiling is **3.26785%**. On adjacent neutral-pair words the potential is at most \(6\kappa=48\), with at most three exchange bonds. On separated words it is at most \(4\kappa=32\), with at most four exchange bonds. The diagonal exchange and the absolute off-diagonal row sum together contribute at most one per such bond. Thus the neutral spectral ceiling is 51 at any path length. The single entrance gives the same upper escape norm. This extension concerns the declared route and charge tag; additional exits, different host interactions or other receiving matter change the energy budget.

## Preparing the bench and reproducing the experiment

The bench uses the actual 1–82 and 3–48 intercell bonds. Before loading, a vacancy visits the nine private degree contacts to write zero contrast. Equal unit loads write 44 exterior physical contacts silent. Four temporary openings beside filled buffer ports provide access to fresh degree-five phase references. The executable verifies the full registered echo undoes a unit-pair write and checks the zero-contrast kernel for all endpoint loads. It also reconstructs the unit-filled preparations, the addressed routes for five initially stored defects and the local pulse frames. The two registered Hamiltonians have dimensions 168 and 30. Their sum on the product space has dimension 5,040. The complete free evolution keeps every physical contact in its prepared state or its invariant zero-contrast sector.

A routed full native exchange lends the addressed load to a blank private reference for a phase dwell and returns it. Private quarter phases around a partial dimer exchange rotate either type of two-port orbital into its energy basis. The source seals the openings while their buffers carry unit load. At the terminal read each observer first checks its buffer. A filled buffer permits the inverse unit-pair write and reference access, and every other buffer observation is retained as a miss. Signed port measurements then distinguish the two energy labels. Starting inventory, bench geometry, detector pulses and exposure are supplied laboratory settings. The route lengths and appointments were explored for resolvable transfer before collecting this bank; the counts test transport at those settings. No lepton mass input enters the Hamiltonian or the channel decoder.

From the repository root, run with a single numerical thread.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.emission_channels.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/emission_channels
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_11_spectrum/emission_channels/EmissionBound.lean)
```
