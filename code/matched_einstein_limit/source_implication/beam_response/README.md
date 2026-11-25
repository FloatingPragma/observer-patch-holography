# Imaging the same finite source through a native gateway

A collar correlation changes the angular image of a material pulse passing between two cells. For the first arrivals, the surviving angular response is exactly proportional to a connected correlation of the two stored heads. The other angular coordinate disappears. A finite exposure resolves both coordinates. Calibrating the change in signed image area on ordinary matter then gives 53.3–57.0% of the full modular prediction for a vacancy and the opposite sign for an excess. The source, fixed-volume adjustments and ordinary modular calibration are those of the same prepared 80-cell equilibrium. [Reference, realizability and volume constraint](../full_algebra/PROOF.md) · [Population measurements](results.json) · [Exact-source enclosures](certificate.json)

## The two screens

The intercell gateway is the edge between ports 1 and 38. The neighbors within each cell form a five-cycle. Starting at the marked contacts, the launch order is \((0,5,9,8,7)\) and the detector order is \((39,40,47,46,42)\). This is one direction of transmission. Repeat the experiment for an excess or a vacancy without combining the species or reversing the setup. On a fresh source launch at one of the five ports, release the native free flight for the specified exposure, and measure which receiver pixel is occupied. The sixth outcome is a miss. An ensemble of sharp launches measures all 25 unconditional hit probabilities \(P_{ji}\). The native flight, including both head parities, determines every probability. [Measurement](measure.py)

Assign the five bearings \(b_i\) to the vertices in this cyclic order. The readout uses the first sine and cosine harmonics rounded to the rational weights in the receipt, with \(\sum_i b_i=0\). The observer can introduce two independent angular biases by using launch probabilities \((1+2u\cdot b_i)/5\) for small real \(u\). The detected image position is the mean receiver bearing on hits. The total hit fraction, reported independently, fixes its normalization. At \(u=0\) its derivative is the two-by-two matrix
\[
 J=2\left(\frac{B^T P B}{W}-\bar b_{\rm out}\bar b_{\rm in}^{T}\right),
 \qquad W=\mathbf{1}^TP\mathbf{1}.
\]
Here \(B\) has the bearings as rows. The bars are means with the receiver and sender marginal weights of the same joint hit table, respectively. The signed image-area gain is \(a=\det J\). Reversing a screen orientation changes the sign of \(a\) and both intervention slopes together. Their transfer ratio is independent of that convention. The reference gains for both species at exposure four are positive, so the results there also hold for the unsigned image area in a neighborhood of the reference.

## The first image is a one-coordinate correlation measurement

Each sender–receiver pixel pair has a unique three-hop route at the first nonvanishing order. Its encounters are sender to the local gateway, the blank bridge, and gateway to receiver. After removing the common third-order amplitude normalization, the blank hops act as the identity on the stored heads. Only sender 0 sees the left marked ingress and only receiver 39 sees the right marked egress. In each parity, let \(L\) and \(R\) be their positive intensity operators. They act on distinct heads. The right operator is independent of the left parity. The four possible intensity entries are then \(1\), \(l=\langle L\rangle\), \(r=\langle R\rangle\), and \(c=\langle LR\rangle\); only the marked-to-marked entry has the last value. The code checks these identities against the exact integer third powers of the native generator for both probe signs. At exposure \(t\), multiply the entries by \(t^6/2304\) to obtain the leading probabilities. The omitted terms are higher order on the fixed finite board.

Substitute these entries into \(J\), with marked bearing \(b_0=(1,0)\) on either end. With \(W=16+4l+4r+c\) in the reduced units,
\[
 J=\frac{50}{W^2}\big(\langle LR\rangle-\langle L\rangle\langle R\rangle\big) b_0 b_0^T .
\]
The proof uses \(\sum_{i\ne0} b_i=-b_0\) on both screens. Thus the one visible angular gain is an independent optical measurement of a head correlation. At the nominal equilibrium it is 0.007365 for vacancy and 0.007109 for excess. Replacing the preparation by the product of the same marginals makes it exactly zero. The leading two-coordinate image area and all its source derivatives are zero for every preparation. At the amplitude level any coherent input supported on the four unmarked pixels with coefficient sum zero is also annihilated by the leading transmission matrix. There are three independent complex inputs of this kind.

The path condition is local and has a check across the supplied carrier sizes. At every one of the 30, 120 and 480 intercell gateways on the 20-, 80- and 320-cell federations, each of the 25 pixel pairs has exactly one walk of length three and none shorter. The same rank bound holds whenever that local condition and the two-head marking are retained. A fixed sequence of such initial-arrival relays using sharp intensities, and any state-independent stochastic reshuffling of their pixel labels between relays, also has at most one transmitted angular coordinate. The first link has only two distinct input columns. Normalizing any mixture of the two transmitted output vectors places the result on an affine line. Additional paths, heads or measurements during flight can change this conclusion.

## Finite-time transfer at the ordinary calibration

At a finite appointment use the full propagator on 960 ports and both stored heads. The numerical scan covers six exposures from one quarter to eight native exchange units. Exposure four was selected for the interval and counting tests after this scan. The ordinary change comes from an unread left sign operation, and unpairing replaces the heads by the product of their marginals. The separate survey's right-hand rotation keeps its volume fixed in both cases. Both changes increase the full head entropy. For unpairing the entire entropy change comes from the collar. Measure \(a'_o,a'_p\) from the optical image and freeze \(a'_o/S'_o\) as the ordinary coefficient. The held ratio is \(a'_p S'_o/(a'_o S'_p)\). A single scalar response of this image area to the full modular charge requires a ratio of one.

| Probe at exposure 4 | Uniform-launch hit fraction | Signed area gain at the nominal reference | Ordinary area slope, nominal | Unpairing area slope, nominal | Held ratio, certified |
|---|---:|---:|---:|---:|---:|
| Vacancy | 0.00128044 | 0.000145012 | +0.000828624 | +0.000257897 | [0.53303, 0.56987] |
| Excess | 0.00133439 | 0.0000000234497 | −0.0000000434896 | +0.0000000116872 | [−0.49952, −0.45406] |

The signed reference areas are strictly positive across the enclosed equilibrium for both species. In the changed-law control the squared-contrast strength is set to zero while keeping the same graph, launch, receivers and head preparations. Both species retain a two-coordinate image, and the computed unpairing response is below \(10^{-15}\) at the center using the original survey's matched direction. The control does not assume that direction would hold a recomputed control-survey volume fixed. [Complete scan and changed law](results.json)

The interval replay uses the bounded integer Chebyshev propagator. Let \(M=\widetilde B^TP\widetilde B\) with \(\widetilde B\) the three columns consisting of the constant and the two bearing coordinates. Then \(a=4\det M/M_{00}^{3}\) by the Schur complement. Separate fixed translations of the two screens preserve this formula and improve interval conditioning. The calculation bounds all nine Hermitian moment operators, the true source and volume-adjuster directions, the derivative of the determinant, and the independently measured entropy slopes. Each amplitude block has operator error below \(10^{-19}\); positive Gram formation and the angular weights add less than \(10^{-16}\) to every reported moment or raw direction expectation. This is separate from statistical uncertainty. [Interval calculation](certificate.py)

The synthetic count experiment uses \(2^{72}\) shots per sender, species and preparation, drawn in multinomial batches. Each row of six integer counts includes misses. Five fresh banks measure the reference, the two intervention endpoints and the two opposite quarter rotations needed for the volume derivative by parameter shift. With simultaneous 99% sampling coverage, the excess ordinary slope is in \([-4.393,-4.304]\times10^{-8}\) and its unpairing slope is in \([1.151,1.187]\times10^{-8}\). The vacancy slopes are in \([0.00082262,0.00083471]\) and \([0.00025521,0.00026036]\). These intervals use the numerical reference as the sampling population. They include the independently certified range of matching angles. The replay verifies the frozen count file by SHA-256, checks both the frozen and newly generated counts against every live detector probability at simultaneous 99% coverage, and reconstructs the signed intervals directly from each count table. [Raw counts and confidence intervals](counts.json) · [Sampling procedure](counts.py)

The measured area is an angular image-area gain between material five-pixel screens. For an interpretation as spacetime focusing one needs the relation between those pixels and neighboring null rays, their cross-sectional separation in proper units, an affine calibration, and control of diffraction and shear. A normalized centroid from a low-transmission angular channel carries different physical information. These records establish both an optical signature of collar correlations and a quantitative failure of its direct scalar identification with universally calibrated light focusing.

From the repository root run

```sh
export OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:.
python -m matched_einstein_limit.source_implication.beam_response.measure
python -m matched_einstein_limit.source_implication.beam_response.certificate
python -m matched_einstein_limit.source_implication.beam_response.counts
python -m pytest -q code/matched_einstein_limit/source_implication/beam_response
```
