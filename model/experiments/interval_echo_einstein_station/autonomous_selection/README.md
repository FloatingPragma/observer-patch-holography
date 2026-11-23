# Does free flight select the measured boundary-capacity equilibrium?

The registered free two-load dynamics preserve distinct interior entropy gradients at identical probabilities for every total-energy reading. This holds with every exterior interaction open. A full-energy-matched stationary preparation gives an interior sharp-record entropy slope of **−log(13/11) = −0.167054 nats per unit probability transferred** on the twelve five- and six-cell shapes. Every reduced exterior record and the boundary completion capacity have zero slope along this preparation change. The unbiased member of the same family has zero interior slope. The experiment uses local endpoint observers and includes every release. [Preparation and readout](PROTOCOL.md) · [Regional flight records](results.json) · [All-seams-open design](OPEN_HOST_DESIGN.md) · [Inclusive open-board counts](open_host.json)

## A preparation with the same full energy distribution

On a fixed board let D denote the number of allowed sharp two-double-load words. In each of three enrolled cells, two independent orbitals in the quadratic cell harmonic vanish at all three gateways. Their kinetic energy is three, and each single double load carries a constant offset of 192 in the full federation. Put one such orbital in the first cell and one in the second cell, giving a four-dimensional product space A. Replacing the second cell by the third gives a disjoint product space B. Both particles occupy different cells, so possible cross-cell collision terms multiply zero gateway amplitudes. Acting with every registered hopping and diagonal term gives energy 390 for every product. The full-board stencil checks this cancellation including every possible hop out of the region. Its largest column residual on the three boards is below 5.7 × 10⁻¹⁴.

The density at the start of each separate release is

$$
\rho_\alpha=\frac{I+\alpha(P_A-P_B)}{D},\qquad |\alpha|<1.
$$

A concrete preparation draws from a complete orthonormal word menu. When the two loads start in different cells, each cell uses its own complete local basis, which includes the two dark orbitals and ten orthogonal completions. Native amount-two splitters and degree phases prepare all twelve local basis inputs. When both start in the same cell, use the sharp pair basis. Uniform probabilities on this menu give the identity density divided by D. Multiply the probabilities of the four A words by 1 + α and the four B words by 1 − α. This is a preparation choice. The free Hamiltonian subsequently selects neither the sign nor the magnitude of α.

All these states commute with the free Hamiltonian. Because A and B have the same energy and the same rank, every energy-projector probability is identical to the uniform state's value. Their reduced exterior states are identical as well. Each biased product has empty exterior and gateway rails, and tracing the interior difference yields zero. This holds at every generator parameter and after averaging with any common distribution of appointments. An active probe which changes the preparation or traverses the region defines another experiment.

## What the endpoint observers see

On the three cells selected in each of the twelve shapes the rank-two orbital projector has diagonal 1/3 on six rails and zero elsewhere. In the cell quadrupole the six antipodal port pairs form a regular simplex. Three distinct forbidden gateway rows leave coefficients on the other three antipodal pairs summing to zero. Their two-dimensional projector has diagonal 2/3 on each pair, split equally between the two ports. The A and B product projectors therefore each have 36 sharp words with diagonal 1/9. These are completely interior words and can be read using terminal routed occupation flags.

At α=3/4 a word in the positive class has probability (1+α/9)/D and a word in the negative class has probability (1−α/9)/D. Transferring one unit of probability uniformly from the negative class to the positive class gives the Shannon derivative −log(13/11). Within each shape, both signs and the unbiased preparation give the same reduced exterior law, including correlations between boundary ports. The exterior conditional completion bits consequently have zero derivative. Fixing the full energy law, total load and the static region does not force the previously measured entropy balance for this autonomous family. These states are faithful finite pair ensembles. Their preparation does not select a physical vacuum.

With interactions parked solely at the regional cut, fresh observers measured slopes between −0.16779 and −0.16537 nats, with standard errors 0.00061–0.00074. With the entire 20-, 80- or 320-cell federation open, the law is uniform plus the same bias on its much larger pair space. A local observer records an A-class event, B-class event or the inclusive complement. A pair partially or wholly outside the region always contributes to the complement. At 2³⁴ independent simulated trials per setting the measured slopes were −0.17405 to −0.16484 nats, at least 33.65 standard errors from zero. The transferred fraction for this single pair of subspaces is 4|α|/D. [Compiled preparations and regional source](measure.py) · [Open-board source and readers](open_host.py)

## Fixed ancestral region

Use all descendants of the same five original cells, giving 5, 20 and 80 source cells in successive federations. Allocate equal numbers of distinct cell pairs to A and B and include all four equal-energy orbital products per pair. Under the registered generator these mutually orthogonal products remain stationary. Repeated antipodal gateway addresses in some descendants change the individual projector diagonals. Group all sharp observations by their actual signed projector weight with their measured multiplicities. [Regional protocol and calibration correction](OPEN_HOST_DESIGN.md) · [Inclusive counts](ancestral.json) · [Source](ancestral.py)

| Descendant cells | Total probability transferred at α=3/4 | Exact-law entropy slope, nats | Measured slope, nats |
| ---: | ---: | ---: | ---: |
| 5 | 0.0005230 | −0.16705 | −0.16720 |
| 20 | 0.0006191 | −0.16154 | −0.16141 |
| 80 | 0.0006431 | −0.15952 | −0.15950 |

The perturbation occupies a finite fraction of the source law on this family. If the host has N cells and the ancestral patch N/4 cells, each sign has floor(binomial(N/4,2)/2) cell pairs. Thus the transferred probability converges to 1/1536 at α=3/4. The enrolled shape has fixed ancestral angular fraction. These address counts do not identify a shrinking proper ball.

## Independent transient and low-energy tests

Two different sharp-word mixtures on each of the twelve regional shapes were launched with the same mean native energy as the uniform pair state. They have different spectral populations. Fresh copies received a single terminal read at exposure 1, 4, 16 or 64. Pooling seventeen uniformly selected appointments between 64 and 128 gave interior sharp-record total variation **0.424–0.712** from the monitored equilibrium. Calculating the time average from resolved energy projectors at a 10⁻¹⁰ level-grouping tolerance gave **0.355–0.678**. Widening the grouping tolerance to 10⁻⁸ changed complete-word TV by at most 0.0171. Very small splittings can affect a numerical infinite-time extrapolation. The exactly stationary dark preparations avoid that extrapolation entirely.

The lowest pair energies of the same regions range from 259.96 to 290.94. The resolved lowest clusters contain one to five levels at the 10⁻⁸ grouping tolerance. They are separated from the next level by at least 0.00578, with maximum column residual below 7.4 × 10⁻¹³. For the cluster's equally weighted mixture, select the interior word set on which its probability exceeds the uniform pair law. The smallest eigenvalue of this event restricted to the cluster bounds its probability for *every* density supported there. That event witnesses at least 0.928 interior TV from the monitored equilibrium on every shape. A fixed-charge sector ground and an externally monitored infinite-temperature word law answer different equilibrium questions. Neither measurement supplies a geometric vacuum, Newton's constant or the area identity.

Run all measurements locally from the repository root into fresh output paths.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_einstein_station.autonomous_selection.measure --output /tmp/regional-auto.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_einstein_station.autonomous_selection.open_host --output /tmp/open-auto.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_einstein_station.autonomous_selection.ancestral --output /tmp/ancestral-auto.json
python -m pytest -q model/experiments/interval_echo_einstein_station/autonomous_selection
~~~

## Half filling with the full energy distribution fixed

The registered 12-port cell has 924 sharp states with six mobile doubles. Six actual graph automorphisms permute its three cut ports as the symmetric group on three objects. Let S and T be two conjugate reflections in this group. Each preserves every internal seam and the full-host degree correction, hence commutes with the regional Hamiltonian. If R conjugates S into T, it also commutes with H and with every energy projector. Consequently every operator

$$
 X=p(2H-cI)(S-T)
$$

commutes with H and has zero trace separately on each eigenspace. Here c centres the integer pair generator. Requiring the partial trace of X onto all three cut ports to vanish gives linear equations in the polynomial coefficients. Using the full 924-word matrices, the first degree with a nonzero interior word bias is four. The coefficient vector and its maximum absolute row sum L are deposited as integers. The symmetry identities, partial trace constraints and nonzero interior diagonal are checked in exact arithmetic. The finite operator comparison states (I ± 3X/(4L))/924 are positive, stationary, and have identical probabilities for every energy outcome and identical full cut density matrices. [Half-filling design](DENSE_PROTOCOL.md) · [Finite calculation](dense.py) · [Complete held counts](dense.json)

At the plus preparation, the interior sharp-word entropy derivative per unit of probability transferred between observable interior words is **−0.215675 nats**. Separate simulated endpoint counts gave **−0.215402 ± 0.000209 nats**. The uniform reference has zero slope. The polynomial density comparison specifies the coherent initial state. A native preparation circuit for that density is outside this calculation.

An independently executable sharp input provides a persistent record on the same physical cell. Let O be the mean diagonal of X conditional on the observed interior word, rescaled so that its largest absolute value is one. Two ordinary initial mixtures draw sharp words with probabilities (1 ± 3O/4)/924, and then undergo undisturbed free evolution. They also have identical full energy laws. Taking a diagonal and pooling by interior word both commute with relabeling by R, so O itself is a difference of conjugate operators whose energy-projector traces cancel. Initially their cut distributions agree. The terminal detector records each interior word and reports its O score. Fresh measured score differences were +0.000724, −0.000115 and +0.000433 at exposures 1, 4 and 16. These are instantaneous records.

For the infinite-time mean, let D_H be dephasing between distinct energy levels, retaining coherences within each degeneracy. The two terminal score means differ by 1.5 Tr[D_H(O)²]/924. The conserved X and Hilbert–Schmidt Cauchy–Schwarz give

$$
 \frac{1.5}{924}\operatorname{Tr}[D_H(O)^2]
 \ \geq\ \frac{1.5\,\operatorname{Tr}(OX)^2}{924\,\operatorname{Tr}(X^2)}
 \ >\ 0.0001980 .
$$

The score comes entirely from the interior sharp record. Thus the time-average record of these physically prepared inputs cannot converge to one common interior law even when their complete energy distributions match. The exact retained cut-density comparison pertains to the stationary polynomial family. The sharp-launch experiment imposes no matching condition on its terminal cut densities.

Replay this supplemental test with `python -m model.experiments.interval_echo_einstein_station.autonomous_selection.dense --output /tmp/dense-autonomy.json`.
