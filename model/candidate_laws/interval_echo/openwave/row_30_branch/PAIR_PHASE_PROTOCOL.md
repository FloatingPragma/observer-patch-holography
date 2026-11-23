# Reading a pair's relative phase in a different cell

Use the blank-memory 20-cell Interval Echo board, all 630 seams, the original degree references and two unit loads. The historical choice of strength 8 used persistence data. The earlier frequency and adjacent-arrival results, including the cell 0 to cell 4 measurement at exposure 256, are known. No phase-transfer measurements enter the design. The question is whether remote observers can distinguish two joint phases even after accounting for how often a pair reaches their analyzer.

In the common port labels of each cell, take the lexicographically first ordered four distinct ports $(u,v,w,x)$ with contacts $uv,wx,uw,vx$. Use the same labels at every station, without aligning the detector to an arrival. One load starts at each of $u,w$. On separate copies offer the two disjoint rails $uv,wx$ native exchange exposures of $(\pi/2,3\pi/2)$ for the plus word and $(\pi/2,\pi/2)$ for the minus word. Read only whether the signed contrast on $uw$ is zero by echo, copying one binary flag, and inverse echo. Retain failure as a terminal record. On success the two disjoint seam words $X=uw,Z=vx$ have amplitudes proportional to $(1,\pm1)$. Both signs have the same preparation success, joint occupancy probabilities and single-load reduced density operator. No source label is sent to the receiver.

At each appointment use a fresh preparation. Locally query unit absolute contrast on each rail, then zero contrast on $uw$. Record the first failed query. The remaining two-dimensional subspace is spanned by $X,Z$. Offer both rails an exposure of $\pi/2$ and read zero versus nonzero contrast on $uw$. These two counts resolve the real off-diagonal element between disjoint pair addresses. Count source refusal, all three analyzer refusals and both final outputs per enrolled launch. On separate copies report the full-board probability of an adjacent pair, a separated pair or a double occupation as an audit of association. No measurement acts during flight.

Use sources 0, 2 and 3. For each use its own cell and all three neighboring cells at exposures 0, 16, 64 and 256. Source 0 is the reference panel. The six directed neighbor routes from sources 2 and 3 at 256 are the held-out phase panel. Freeze their predicted outcome probabilities before generating counts. Propagate the complete two-load sector with strength 8 and with exchange alone. A third flight evolves two independent identical bosons under the one-load generator, including double occupancy and both assignments at a coincidence. Use the same entangled input, source success, apparatus and clock in each arm. Dephase the two input words before flight as a sign-blind control. Also dephase the two admitted receiver words before its pulses as a control for unequal admission rates. A parked flight has zero distant admission.

Give each sign on each held-out native route $2^{26}$ independent simulated launches. Use NumPy PCG64 with seed sequence [41029, source, receiver, sign index], plus first. Sample the six-outcome distribution from a separately assembled propagation after freezing the predictions. Reject a numerical reproduction if any predicted bin differs by more than $2\times10^{-9}$ or a retained distribution fails normalization at that tolerance. Report both the per-launch signal and the fraction of equal flags among receiver admissions. Conditional on the two admission totals, apply two-sided Fisher exact tests to the two final counts for each sign. For detection require the predicted direction and $p<0.01/6$. Report every route including zero admissions. Call the phase reproducible across sources if at least one route per held-out source passes. Treat the other times, arms and sizes of the phase contrasts as comparisons. Test binding separately by a uniform-in-time lower bound on adjacent probability above one half for both native launches, and by reporting adjacency in every competing flight. Coincidences at separated appointments describe fresh copies. They do not supply a continuously monitored trajectory.

All durations, the four ports, the controls, the statistical threshold and the sample budget are apparatus choices. The comparison uses native exposure and record addresses. It asks for neither a mass in physical units nor an OpenWave cell decision.

## What the receivers record

The four local port labels are $(u,v,w,x)=(0,1,5,9)$. Half the enrolled preparations return the source refusal. On acceptance the state is $(|u,w\rangle\pm|v,x\rangle)/\sqrt2$. The two words have no common occupied port. Tracing over either load removes their cross term exactly. Copying only the source's binary contrast flag preserves that coherence within the accepted state.

At the receiver let $a,b$ be the unnormalized amplitudes on its two admitted words, including the source success amplitude. The two terminal probabilities per enrolled launch are $|a-b|^2/2$ (equal) and $|a+b|^2/2$ (unequal). Dephasing before release gives each sign the average of the two undephased six-bin distributions. Dephasing after receiver admission preserves its rate and splits its two final outcomes equally. The replay checks the first control against incoherent sharp-word flights and the second against the actual local pulse. [Native instrument and predictions](pair_phase.py) · [Independent propagation and tests](pair_phase_check.py)

Each table entry gives equal flags / all local admissions out of 67,108,864 simulated preparations of that sign. The Fisher test compares the flags within admissions. The six-route cutoff is 0.001667. Both source cells have two passing routes. At this budget the reciprocal routes between cells 2 and 3 miss the cutoff. The receiver can distinguish ensembles of preparations. A single remote launch rarely returns either flag.

| Source → receiver | Plus | Minus | Two-sided p | Predicted direction detected |
| --- | ---: | ---: | ---: | --- |
| 2 → 1 | 122 / 590 | 1,350 / 1,987 | 2.23 × 10⁻⁹⁴ | Yes |
| 2 → 3 | 222 / 407 | 301 / 658 | 0.00554 | No |
| 2 → 9 | 3,025 / 3,315 | 2,740 / 3,721 | 4.89 × 10⁻⁸⁶ | Yes |
| 3 → 2 | 360 / 579 | 284 / 530 | 0.00419 | No |
| 3 → 4 | 674 / 1,006 | 11,771 / 12,101 | 5.59 × 10⁻²⁰⁶ | Yes |
| 3 → 8 | 280 / 364 | 215 / 345 | 2.83 × 10⁻⁵ | Yes |

Predicted remote admissions range from 5.42 to 181.1 per million enrolled launches across these routes and signs. On 2 → 1 the equal flag probabilities per million launches are about 1.47 and 19.30. Among admissions their fractions are 0.164 and 0.675. Independent bosons also carry phase on this route (fractions 0.476 and 0.860), and exchange alone gives 0.532 and 0.961. The distinction supplied by the native interaction is the high probability that the two loads remain adjacent in a terminal position measurement. [All laws, source cells, exposures and outcomes](pair_phase_predictions.json) · [Held-out counts](pair_phase_counts.json)

Within the declared coherent dynamics, the conditional change also constrains the two-word transfer matrix. Let $A$ be the free transfer between the source and receiver code spaces multiplied by the common source success amplitude. On the orthogonal signed inputs its columns have squared norms $q_+,q_-$. Let $f_+,f_-$ denote their equal-flag fractions on admission. For the two normalized output vectors the binary statistical distance $|f_+-f_-|$ is bounded by their pure-state trace distance. Taking the determinant of their Gram matrix gives

$$
\sigma_{\min}(A)^2 \ \ge\ \frac{q_+ q_-}{q_++q_-}(f_+-f_-)^2 .
$$

Here the product of the squared singular values equals the Gram determinant and their sum is $q_++q_-$. For 2 → 1 the predicted probabilities give a lower bound of $1.77\times10^{-6}$. A channel that always emits one fixed receiver state and varies only its intensity cannot give these fractions. The replay also evaluates the smaller singular value directly. This bound uses the predicted instrument, one coherent projected propagation per admission and a supplied phase convention. The sampled tests above concern two real input phases.

## Binding during free propagation

Write the shifted generator as $H=K-16P$, with $K$ the positive exchange Laplacian and $P$ the projector onto adjacent pairs. The accepted source has $\langle P\rangle=1$ and $\langle K\rangle=17/4$ for both signs in all three cells. Energy conservation therefore gives $\langle P\rangle_t\ge 47/64$ at every exposure. Including the source refusal gives at least $47/128$ per enrollment. At exposure 256 the actual adjacent probabilities for the six source/sign settings after acceptance range from 98.99% to 99.36%. Exchange alone gives 15.45% to 26.63%, and independent bosons give 4.70% to 5.92%. These probabilities refer to the full released ensemble. They place no corresponding lower bound on an individual detected pair's history.

There is also an exact obstruction to rigid motion. In the adjacent subspace a single swap connects two occupied seams only when they belong to a triangle. Every triangle is internal to a cell. The adjacent-only generator thus decomposes into twenty connected blocks of thirty internal seams and thirty individual bridge states. Its free evolution gives zero transfer between cells. In the full generator both loads can pass through a bridge. Starting on an internal seam incident on the departing bridge port and ending on an internal seam incident on its far port takes four moves. The first load must cross first and the second must cross last. The two intermediate moves can occur in either order. The integer move matrices through order four verify this count for all 1,500 ordered shortest crossings and exclude other intercell internal-seam endpoints at these orders. Each shortest path contributes $(-1/2)^4$ to $H^4$. For those individual seam words the leading amplitude is $t^4/192$, independent of the diagonal contrast. Thus a strong terminal adjacency signal is compatible with transfer through separated configurations; a projection enforcing adjacency throughout would remove the transfer. [Exact crossing calculation](pair_crossing.py)

## Replay

From the repository root, set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` for the extended calculations. The forecast uses the registered Chebyshev expansion. The second implementation constructs the exchange graph from the hashed board, advances it by sparse matrix exponentials, and recomputes the receiver flags separately. It checks 1,728 receiver probabilities and 288 association probabilities before drawing the six multinomial count tables. The data are simulated counts and the appointments use native exposure units.

```sh
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_30_branch.pair_phase
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_30_branch.pair_phase_check
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_30_branch.pair_crossing
```
