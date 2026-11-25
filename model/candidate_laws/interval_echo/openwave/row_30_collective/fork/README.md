# What happens when a strand wall reaches a three-way junction?

Two observers in different outgoing cells can detect walls together after an incident wall is prepared at a junction. Their joint records also distinguish how the incoming arm was prepared when every other cell's reading is discarded. On the 20-cell screen with two dimers per arm, one orientation gave **10,470,504** joint detections in **268,435,456** launches. The same apparatus with the zero-wall input gave **71,370**. The blank reference hand placed the exposure at approximately two exchange radians. The result is a material interaction across the fork. Extending the strand's measured single-particle Klein–Gordon law directly to all neighboring wall sites predicts far too many exclusive turns. [Design and follow-up](PROTOCOL.md) · [Complete simulated counts](counts.npz) · [Count summaries](results.json)

## Three neighboring owners measure the same fork

The center is the native seam between ports 6 and 7 of cell zero. A gateway dimer spans a different cell boundary on each arm. A written connector joins that dimer to the center. The following connector on each arm sits wholly in the neighboring cell. All other seams outside the prepared dimer and connector network have idle memories. The preparations use the blank dimer, coherent connector and idle collar of the [strand construction](../README.md). The supplied fork addresses, collars, one load per dimer, initial orientations and terminal analyzers are apparatus. During flight every contact has the registered law. The 20-, 80- and 320-cell boards realize the same local fork with different addresses. Comparing them checks the embeddings and counting reproducibility. It does not increase the number of inequivalent junctions tested.

Two neighboring dimer endpoints have a wall when their loads disagree. All dimer orientations can be read simultaneously. Equal endpoint loads on every connector define two zero-wall words related by flipping every dimer. Starting from each word, reverse all dimers on arm A. This produces one wall on the connector nearest the junction in A. The source is a sharp wall at the fork; it has no prepared direction of momentum. B and C have identical local preparations in the incident and zero-wall runs. Tests use both incident orientations and keep all wall counts, including many-wall events.

With three dimers per arm at four exchange radians, the complementary incident input had **98,846,623** B-and-C coincidences, against **50,089,922** on its zero-wall counterpart. Both counts include every other wall configuration, giving an increase of about **18.16 percentage points**. At two radians the other orientation on these longer arms gave **26,305** versus **3,904**. Timing comes from a separate exchange hand over one quarter of each flight. Full simulated orientation counts at every appointment, the two phase-analyzer settings, and the independent clock bins are archived. With contrast removed during flight, initially disjoint dimer orientations evolve independently. The remote-only joint probability is then exactly the same for either source, even when differences survive after selection on the entire fork.

The exclusive counts measure a different event. In the two-dimer complementary run at two radians, **10,451,992** launches ended with exactly the two outgoing walls one link beyond the junction and no other walls, against **43,594** zero-wall launches. Opposite quarter-turn measurements on the central dimer distinguish a coherent transition from the input orientation to the orientation with both outgoing walls at the junction. On this source their all-start probability difference reaches about **−0.647** at one exchange radian. Coincidences alone would not establish an uninterrupted history of either wall, or entanglement between the outgoing owners.

## An exact obstruction to a direct single-wall turn

In the prepared orientation space, each active blank dimer flips its own bit with matrix element −h, where h=1/2. The coherent connector contributes diagonal energy 4κ(1−ab) on its endpoint loads a,b. The idle collars contribute zero. Degree references are constant on each degree-matched dimer. Up to a scalar, the entire Hamiltonian is therefore a sum of independent bit-flip operators −hX and a diagonal interaction. This restriction follows from the full native seam columns, including their components outside the selected words. [Generator and record simulator](measure.py) · [Independent contact projection and finite matrix checks](test_fork.py)

On a tree every full wall word has two preimages, a dimer word and its complement. Let L≥2 be the number of dimers in each arm. To change the incident wall on the first A connector into a lone wall on the first B connector requires reversing either **all 2L dimers of A and B** or **the center and all L dimers of C**. A step along A requires one dimer flip. Turning one incident wall into the two near-junction walls in B and C also requires only one flip, at the center. Moving both walls one link farther needs three flips, one at the center and one on each of the two outgoing gateway dimers.

If two specified orientations differ in d places, the first nonzero off-diagonal term of the free propagator is (iht)ᵈ. No lower power can reverse those bits. At power d there are d! possible flip orders and no diagonal insertions. Thus the amplitude to turn into a lone outgoing wall starts at order t³ on the short fork and t⁴ on the longer fork. A direct step along A and the near-junction split both start at order t. This comparison holds for every diagonal interaction on these dimer words. The contrast influences subsequent times and the detector probabilities.

A finite-time upper bound follows without a series truncation. For a set S of required reversals let Pₛ project onto states in which these bits are reversed from their initial orientations. Evolve inside Pₛ by its compressed Hermitian Hamiltonian. The only terms feeding it from outside are flips of one of its required bits. The driven Schrödinger equation gives

$$
\|P_S\psi(t)\| \le h\sum_{i\in S}\int_0^t \|P_{S\setminus\{i\}}\psi(u)\|\,du
\le (ht)^{|S|}.
$$

The last bound follows by induction. Summing the squared bounds for the two disjoint terminal orientations bounds the probability of an exclusive turn by \((ht)^{2L+2}+(ht)^{4L}\), and by one if this expression exceeds one. For fixed ht<1 increasing the arm length suppresses this one-wall turn arbitrarily strongly on such a family of forks. The executed boards realize lengths two and three.

## Test against the scalar dispersion measured on a strand

Take one wall site per connector. Join two sites when their connectors touch the same dimer, and let Lg denote that graph's ordinary Laplacian. A direct single-particle continuation of the straight strand has positive-frequency generator √(m²I+bLg). On the infinite line it gives the measured lattice dispersion with squared spatial wave number 4sin²(k/2). Independent strand meters supply m=15.00051 and b=16.00156 in exchange and dimer units. The comparison identifies the sharp incident material wall with a sharp site of this scalar theory. Its forecast uses the actual graph and the independent timer. Both parameters range over ±1% and each timer over its simultaneous count interval. The calculation bounds the change of the complex amplitude over this full box before squaring. [Forecast, calibration and fresh acquisitions](scalar.py) · [All follow-up counts](scalar.json)

At the earliest nonzero appointment the fresh banks gave **zero or one exclusive turn per 268,435,456 starts** across both outgoing directions, both input orientations and both arm lengths. Every scalar forecast exceeds **0.000615 per launch**. The simultaneous counting upper limit on any of the eight turn readings is **7.13 × 10⁻⁸**. All four straight-arm readings overlap their scalar forecast intervals. This excludes this particular graph continuation. A different junction condition, a different initial dressing, or a field that allows interactions and extended interfaces requires its own observer test.

Intervals invert the binomial lower and upper tails directly. A Bonferroni allocation across at most 50,000 reported marginal probabilities gives at least 99% simultaneous count coverage. Parameter variation in the scalar test is an explicit robustness allowance. The public [OW-30 criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L136) names emergent massive dispersion for matter waves and specifies no spatial dimension. The fork tests extrapolation beyond a strand. It supplies neither an accepted OpenWave state and decoder nor a public held acquisition on the proposed packet. [Scope of the row decision](../../submission/ROW_30_DECISION.md)

From the repository root, replay both acquisitions and the independent checks with

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
P=model.candidate_laws.interval_echo.openwave.row_30_collective.fork
.venv/bin/python -B -W error -m "$P.measure"
.venv/bin/python -B -W error -m "$P.scalar"
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_collective/fork
~~~
