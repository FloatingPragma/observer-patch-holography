# Boundary entropy and independently measured transmission at finite filling

In an equilibrium of mobile double loads at half filling, adding transmitting boundary contacts by changing the shape of a region produces an interior entropy change whose magnitude approaches **1.38629 nats per available contact**. The same experiment at other fillings gives different coefficients. At fixed shape and half filling, fluctuations in how many boundary contacts are occupied change the interior entropy while leaving the mean transmitting aperture fixed. [Acquisition protocol](PROTOCOL.md) · [Native contact calibration and counts](measure.py) · [Full results](results.json)

The five-cell cycle and open path have equal numbers of ports and different numbers of boundary seams. Adding the same sixth cell gives a second equal-volume pair. Each microscopic subdivision quadruples the number of ports and doubles the number of cross-boundary seams. All seam addresses were checked on the eight deposited meshes from 20 to 327,680 cells. One survey follows the descendants of the same base patches. Another takes the parent patches at levels 0 through 3 and their descendants on host levels 1, 3, 5 and 7. In the spherical subdivision chart these parent patches shrink, while the number of microscopic boundary contacts grows. Their spatial proper radius requires an independent identification.

| Finest shrinking patch family | Two loads | Quarter filling | Half filling | Three-quarter filling |
| --- | ---: | ---: | ---: | ---: |
| Five-cell parents | 0.00073 | 0.74307 | 1.37413 | 2.22964 |
| Six-cell parents | 0.00062 | 0.74417 | 1.37705 | 2.22285 |

Entries are nats of conditional boundary capacity gained per additional available contact. The half-filled prediction differs from the two measured transfers by 0.88% and 0.67%. Launch volume was measured by uniform address membership. A separate probe sampled a uniform regional address on a fresh preparation. A vacant rail admits a distinguishable tag. If the tagged address faces the exterior, a calibrated half-contact transfers it to a blank exterior rail with probability one half. Ineligible launches and absent contacts count as quiet. Twice the observed click frequency times the measured launch volume estimates the transmission aperture. This gated probe measures available contacts. Free flight through arbitrary occupied matter is a different experiment.

At fixed half filling on a *single* host, compare the aperture cubed divided by the squared launch volume for a smaller parent with that for the base parent of the same shape type. The measured ratios for five-cell cycles are **2.001, 4.003 and 8.009** for parent levels 1, 2 and 3. Separate positive area and volume calibration constants common to the board cancel from each comparison. The exact ratios are 2, 4 and 8 because area proxy and volume proxy scale as one and two powers of the parent length. In a similarity-preserving three-dimensional interpretation both proxies would give the same cubed-area to squared-volume ratio for the two shapes.

## Which entropy the apparatus reads

Inside a region of \(V\) rails choose \(N\) occupied sharp addresses uniformly without replacement, leaving exterior contacts parked during free evolution. The blank native binary sector preserves \(N\), so the resulting identity density on that sector is freely stationary. This preparation fixes uniformity. For \(m\) cut rails and \(K\) occupied cut rails, the conditional number of distinct exterior words completing any given interior word is \(\binom mK\). Under this particular uniform law, their conditional entropy is the logarithm of that number. The information observer reads the occupancy class and scores it using these known multiplicities. The largest alphabets require the exact equiprobability result to interpret the samples. The transmission and volume observers use independent preparations and counts. Holding \(V,N\) fixed gives \(H_{\rm int}=\log\binom VN-C_m\), where \(C_m=\mathbb E\log\binom mK\).

Let \(f=N/V\) and \(h(x)=-x\log x-(1-x)\log(1-x)\). There is a uniform finite error bound

\[
0\le m h(f)-C_m\le \log(m+1)+(V-m)/(V-1).
\]

For \(x=k/m\), the \(k\)-term is a mode of the binomial law with parameter \(x\). Its probability is between \(1/(m+1)\) and one, giving the corresponding two bounds on \(\log\binom mk-mh(x)\). Averaging, the concavity gap \(h(f)-\mathbb E h(K/m)\) is the mean binary relative entropy against \(f\). The inequality \(\log y\le y-1\) bounds this by \(\operatorname{Var}(K/m)/[f(1-f)]\). Sampling without replacement gives \(\operatorname{Var}(K/m)=f(1-f)(V-m)/(m(V-1))\). Thus for equal-volume shapes with a fixed ratio of distinct cut sizes tending to infinity, their capacity difference per additional available contact approaches \(h(f)/(1-f)\). In the fixed-two-load case the log capacity per contact vanishes. The limiting half-filling slope is an information-to-transmission conversion for this apparatus.

## What a second tag adds

At half filling, move probability from the central \(K\) class in equal amounts to two classes at distance \(w=\max(1,\lfloor\sqrt m/2\rfloor)\). This keeps total load and mean transmitting aperture fixed while changing the cut-occupancy variance by \(w^2\) per unit of probability transferred. Two disjoint simultaneous tag pulses measure the second factorial moment of the number of available sites. For an unbiased ordered pair of distinct sites their joint click probability conditional on \(K\) is \(p^2(m-K)(m-K-1)/[m(m-1)]\), where \(p=1/2\) is measured on the isolated seam. At fixed mean \(K\) the resulting second-moment slope is exactly the variance slope. For even \(m\) the interior entropy change is \(\log[\binom m{m/2}/\binom m{m/2+w}]\). Taking logarithms of the \(w\) successive binomial ratios gives \(2w^2/m+O(m^{-1})\) for \(w=O(\sqrt m)\). This predicts an entropy change from the coincidence counts independently of the interior counts. [Coincidence protocol](COINCIDENCE_PROTOCOL.md) · [Complete two-tag counts](coincidence.json)

On the finest shrinking five-cell patch (\(m=80\)), two-tag counts predicted an interior change of \(0.3752\pm0.0138\) nats per unit transfer. The independently calculated entropy slope was 0.39569. For a transfer among five classes at equal spacing with weights proportional to \((1,-4,6,-4,1)\), every joint two-tag flag has exactly zero slope. Each flag probability is a quadratic in \(K\); the signed transfer annihilates even cubic polynomials. The interior entropy slope on the \(m=20\) shrinking cycle is 0.007420, measured as \(0.007361\pm0.000403\) nats. At \(m=80\) this control gives a smaller predicted 0.001965 nats and its held sample has a standard error of 0.000837. Higher moments of the same exterior completion distribution carry information omitted by the joint transmission read.

The five-class transfer also fixes mean native preparation energy. At fixed \(K\) the sharp-word expectation of the parked-region Hamiltonian sums diagonal seam costs and a degree correction. A pair of equal-class ports contributes the sampling-without-replacement probability of unlike occupations; ports in different classes contribute the corresponding product probabilities. Each term is a polynomial of degree at most two in \(K\), so its five-class slope vanishes exactly. The three-class variance transfer generally changes that energy. At the largest fixed-footprint resolutions, finite coincidence and five-class samples become too noisy for these small derivatives. The exact outcome polynomials and signed log multiplicities give the population tests at those resolutions.

These measurements concern sharp-word entropy, an addressed aperture and a stationary preparation. A geometric interpretation requires independently reconstructed proper area and volume, a physically selected vacuum and control of the absolute first-variation error at the Einstein small-ball scale. [Autonomous free evolution](../autonomous_selection/README.md) · [Material-clock geometry](../proper_ball_ruler/README.md)

Use the eight files accepted by the repository's pinned refinement reader. Run both commands from the repository root into fresh output paths.

~~~sh
python -m model.experiments.interval_echo_einstein_station.density_area.measure --carriers /path/to/carriers --output /tmp/density-area.json
python -m model.experiments.interval_echo_einstein_station.density_area.coincidence --output /tmp/two-tags.json
~~~
