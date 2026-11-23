# Row 4 · Free stability of a localized Interval Echo pair

## Two loads on a connected board, independent of coupling

On the 20-cell federation, a cell-13 observer starting with one load on 156 and one on 157 can prepare a pair whose probability of detection on an original home seam is at least **99.4771841% at every free exposure, for every real contrast strength**. This includes failed preparations. The terminal measurement is the M10.3 instrument with all 30 seam addresses and quiet. The observer prepares eleven neutral seam records before release. These records persist during flight. The full board of 240 ports stays connected and the positive-interference combination of the same two pair configurations can reach another cell. [Protocol](PROTOCOL.md#connected-release-from-the-historical-sharp-pair) · [Source operators and full-board certificate](open_pair.py)

The launch seam fixes a degree-six endpoint $a=157$ and a degree-five endpoint $b=156$. Among the neighbors of $b$, only $c=166$ has the degree of $a$ and is distinct from and nonadjacent to $a$. Their second common neighbor is 163. A neutral record on an empty contact makes its exchange and contrast columns vanish for every subsequent load word. Four native sharp swaps allow the same two loads to write eleven of the thirteen contacts leaving the three-port path. The two contacts to 163 stay blank. The final detector also operates on written contacts. Its echo, flag and inverse echo measure the same binary occupancy on both blank and neutral records. [Native 99-state checks and preparation route](corridor_pair.py) · [Open contacts](open_pair.py)

Write $X=ab$, $Y=ac$ and $Z=bc$. During preparation the observer addresses the two path contacts and the degree reference at $b$. The contrast on these three words is $2I$. Twice the exchange is the path Laplacian with diagonal $(1,2,1)$ and off-diagonals $-1$. Its eigenvectors include $D=(X-Z)/\sqrt2$ with eigenvalue 1 for twice the exchange. For the bright combination $(X+Z)/\sqrt2$ the survival amplitude under a pulse of length $2\pi/3$ followed by an occupied read at $b$ is $1/3$. The corresponding dark amplitude is $e^{-i\pi/3}$. The common strength-dependent phase cancels from every record. Three occupied reads from $X$ leave the unnormalized vector $(-13X+14Z)/27$ up to phase. A first vacancy leaves the sharp separated pair $Y$ and a single addressed swap restores $X$. The first-vacancy probabilities in an attempt are $4/9$, $4/81$ and $4/729$. With failure probability $f=364/729$ and eight attempts, the exhausted fraction $f^8=0.0038635883$ contributes to quiet. The acquisition uses at most 24 pulses and local occupancy queries and seven retry swaps in addition to the initial writes and swaps. Flight has no intermediate query. [Acquisition and complete instrument](open_pair_release.py)

During release the complete native generator has $H D=(1+3\kappa)D$. The two paths into each of $ac$ and $b\,163$ cancel. The detector contains both components of $D$, so the joint probability is bounded below by its initial, unconditioned dark weight $(1-f^8)/(1+9^{-3})$. More generally $R$ occupied reads and $M$ attempts give loss at most $9^{-R}+2^{-M}$. The orthogonal bright state has probability $t^8/18432+O(t^{10})$ of reaching a different cell at short exposure. A competing Poisson swap law also loses home detections. At exposure one the single jump into $Y$ with no other jumps has conditional probability at least $(e^{-1}-e^{-2})/2$ in a prepared ensemble. [Exact flight and path counts](open_pair.py) · [All 1,116 outcome probabilities, in trillionths](connected_records.json)

| Strength, exposure 37 | Prepared coherent pair | Dephased pair | Random swaps | Sharp unfiltered pair |
| --- | ---: | ---: | ---: | ---: |
| 0 | 0.995196 | 0.652907 | 0.123196 | 0.655440 |
| 1/2 | 0.995945 | 0.926296 | 0.123196 | 0.929889 |
| 8 | 0.996110 | 0.986455 | 0.123196 | 0.990281 |

With all twelve boundary contacts of a different three-port path written, a closed-corridor control gives an all-strength, all-time joint floor of $2975/3072$. The closed path provides a check of the instrument and a comparison for geometric confinement. In the connected experiment the eleven written contacts also constitute deliberately prepared material structure. The protected state is a product: port 156 holds one load with certainty and the other load occupies the antisymmetric mode of ports 157 and 166. The contrast interaction contributes a common phase. Two of the written contacts, (152,166) and (157,170), are external links of cell 13, so the cell exchanges load with the other cells through (159,226) alone. The result concerns two loads localized relative to that structure. It supplies neither its charge and mass nor evidence of a freely translating particle on an initially uniform board. The historical sharp, blank-memory data below retain their coupling-selection history. [Closed source and analytic time law](corridor_release.py) · [OpenWave review](https://github.com/openwave-labs/openwave/pull/745#pullrequestreview-6005365373)

The independent reader reconstructs the 28,680 states from the hash-checked board and propagates with Chebyshev polynomials and Poisson uniformization. Its maximum absolute difference from the archived bins is below $1.8\times10^{-11}$. The receipt rounds individual probabilities to a trillionth. The numerical comparison allows $2\times10^{-10}$ for propagation. Both are separate from the exact bound on the protected weight. [Independent reader](open_pair_check.py)

Reproduce the source alone with the first command. The second emits all nine settings and four arms. The third recomputes the same probabilities independently. The last command checks the obstruction for the blank-memory preparation.

```sh
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_04.open_pair_release
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_04.open_pair_release --full
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_04.open_pair_check --full
python -B -O -W error -m model.experiments.interval_echo_consensus.pair_window
```

## Release on blank memories

An initially adjacent pair of unit loads was released on the complete spherical federation with every seam active. After 150 native time units on 20 cells, a one-cell observer could still detect the pair on one of its 30 internal seams with probability **98.3881%**. Only **2.8306%** lay on its starting seam. The material distribution moved within the patch while the constituents overwhelmingly stayed together. Removing the contrast interaction lowered the same local detection probability to **14.9559%**. [Frozen test](PROTOCOL.md) · [Measurements, detector counts and addresses](results.json)

| Sphere and last exposure | Preparation or interaction | Adjacent anywhere, lowest over tested times | On a seam in the original cell, at last time | Both in original cell, at last time |
| --- | --- | ---: | ---: | ---: |
| 20 cells, 150 | Sharp adjacent release | 99.0688% | 98.3881% | 98.7525% |
| 20 cells, 150 | Native partial exchange before release | 95.0280% | 94.4304% | 95.8774% |
| 20 cells, 150 | Exchange-only control | 19.1675% | 14.9559% | 25.3091% |
| 80 cells, 32 | Sharp adjacent release | 99.0545% | 98.3465% | 98.6924% |
| 80 cells, 32 | Native partial exchange before release | 95.0696% | 94.5248% | 96.4974% |
| 80 cells, 32 | Exchange-only control | 24.2755% | 21.9452% | 34.4417% |

The table gives finite-unitary probabilities over every preparation, including those recorded as quiet. A separate sample had **16,093/16,384** native patch detections at time 150, against **2,474/16,384** with exchange alone. Each detector uses a fresh ensemble. It supplies no stabilizing intermediate measurements. At that exposure, the patch radius enclosing 90% of preparations was zero cell steps for the native sharp and perturbed pairs and four steps for exchange alone. These radii use the simulated position distribution. The one-cell detector uses the registered native echo on the 30 internal seams and can directly register either a seam address or quiet. [All appointments and detector outcomes](results.json) · [Source and terminal instrument](../../../../experiments/interval_echo_encounter/README.md)

## What stabilizes the pair

The registered native interaction on a blank collar contributes an exchange and a penalty proportional to the squared difference in endpoint load. Interval Echo includes the common degree-reference terms. With two distinct units, their total contrast is 10 on an occupied seam and 12 when separated. After subtraction of 12κ times the identity, the full binary-sector Hamiltonian is H = K − 2κP, where P detects an adjacent pair. The measurement uses the registered κ = 8. All internal and external exchanges operate throughout flight. The only supplied interventions are initialization, the optional initial partial exchange and the terminal read. Evolving with κ = 0 preserves the same load count and admits the large losses shown above. [Native law](../../../contrast_memory_spherical_v1/law.py) · [Exact sector restriction](../../../contrast_memory_spherical_v1/sectors.py)

The isolated low-energy band is invariant under free evolution. The registered all-size bound gives Pr(adjacent at time t) ≥ 1 − 40/(2κ−6)² = **0.60** for every real t and every normalized release supported on occupied seams. Its proof uses bounded degree, the native cell incidence and the external matching. This is a statement about constituent separation. Spatial containment of the pair in its original cell has the finite exposure and board limits in the table. The two screen sizes have different gateway placements and observation horizons. Repeating the local result on them gives a size check at their common appointments. [Binding proof](../../../contrast_memory_spherical_v1/PROOF.md)

Before release the sharp source places one unit at port 0 and one at port 1. A separate check used the 36 uniformly weighted, initially separated two-unit preparations on cell zero of the 20-cell sphere. After one unit of unassisted evolution, a sharp read at the specified seam prepares precisely the released state with probability **0.01891418%**. The remaining **99.98108582%** of those starts are quiet at that seam. The high persistence probabilities in the table apply after a successful preparation. The perturbed arm moves 4% of the preparation probability into a nonadjacent configuration by a single native partial exchange. [Independent encounter experiment](../../../../experiments/interval_echo_encounter/README.md) · [Present source and branch calculation](./__main__.py)

## Scope for an OpenWave submission

**The dynamic free-persistence test passes for the adjacent pair released on blank memories at the registered κ = 8.** The submission concerns a charge-unassigned composite on a discrete quantum federation. Its electromagnetic coupling, spin and identification with the electron have yet to be established. Derrick's continuum scaling theorem concerns a different class of field functionals. Here a spectral gap explains the measured relative binding. The simulation runs in the RER executable. Upstream validation and placement are separate decisions. [Public criterion](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md) · [Candidate specification](../../SPEC.md)

## Reproduction

Run from the RER repository root with its model dependencies installed. The first command rebuilds both boards, evolves all three arms and verifies the archived receipt. The second compares the sparse evolution against a direct dense cell exponential and validates the native preparation and readout columns.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04
OPENBLAS_NUM_THREADS=1 python -B -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_04
```

The binary basis contains 28,680 configurations on 20 cells and 460,320 on 80 cells. A zero of energy has been removed. Uniform Chebyshev propagation uses the enclosure [−2κ,12] and an analytic Bessel-series tail below 10⁻¹² at each appointment. For any s > 0, the contour representation gives |Jₖ(z)| ≤ exp(z sinh s − ks). Summing this bound geometrically controls the omitted tail without access to a desired trajectory. The exact single-cell exponential and a direct small-board sparse exponential provide separate numerical comparisons. The receipt records model probabilities, simulated finite counts and the supplied native exposure for each setting. [Implementation](./__main__.py) · [Independent checks](test_release.py)
