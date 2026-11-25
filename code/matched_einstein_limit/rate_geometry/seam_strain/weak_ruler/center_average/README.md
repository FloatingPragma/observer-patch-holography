# Does averaging native radar balls recover entropy balance?

Each observer in this calculation fixes the volume of their ball to the same value. On the five-station radar filling, averaging their areas gives the fourth power of radius required for a finite curvature reading. The coefficient misses the prediction calibrated on an ordinary rail change by a factor greater than 59 on a change of head pairing. There is also a physically accessible family that keeps the averaged area and enclosed count exactly fixed while the measured regional entropy increases. Both statements use one preparation protocol. The continuum of observer positions and the identification of the twelve-port aperture with their source region require separate justification. [Native effects and interval results](native.json) · [Five-station source](../hinge.json)

## What the observers measure

On identically prepared copies of the 80-cell reciprocal-strain prototype, evolve the vacancy pump for four ticks. Read its spent mobile, both heads and strain on the twelve-port region enrolled by the hardware aperture. The number meter projects onto the enrolled mobile modes. To measure the energy of this region, give each bond joining an inside and an outside port half weight. In exchange units the resulting regional operator is \(H_R=(PH+HP)/2\), where \(P\) projects onto the inside and \(H\) is the generator that moves the mobile. Regional energy adds over disjoint partitions. A shift of the energy origin by \(c\) adds \(cN_R\). The signed regional entropy response comes independently from the negative logarithm of the reference regional state on its occupied support. [Regional algebra and entropy](../regional.py) · [Energy and number calculation](native.py)

On other copies of the same preparation, retire the mobile and send weak monitor probes from the five radar stations. Pool both charges and reciprocal orientations at 16 inspections per tick, including missed returns at the deadline. Fill the four measured tetrahedra as in the original experiment. Three tetrahedra meet at each interior radial edge. Their dihedral angles sum to the full cone angle α. [Radar protocol and filling](../hinge.py)

Enroll observer centers uniformly by proper volume inside narrow cylinders about the middle third of each of the four interior edges. Choose one small sampling volume \(W>0\) for all nearby preparations. A common transverse cutoff \(d\) is fixed by \(W=d^2\sum_i\alpha_i L_i/6\), where \(L_i\) is the full edge length. Choose \(W\) small enough that the cylinders are disjoint and lie inside the star. Their ends specify which centers to sample. An observer measures their entire ball, including the part outside the enrollment cylinder. These centers are a specified continuum extension of the five actual stations.

Let \(\ell\) label the Euclidean volume \(4\pi\ell^3/3\) and let \(C\) be the coefficient of the mean area loss multiplied by \(W\). At a center farther than \(\ell\) from an edge, an individually volume-matched ball has Euclidean area. Centers within \(\ell\) occupy a fraction of order \(\ell^2\), and an affected ball loses area of order \(\ell^2\). Integrating over the center and enforcing its proper volume gives exactly \(W(4\pi\ell^2-\langle A\rangle)=C\ell^4\) for all sufficiently small radii on these cylinders. Fixing only the ensemble's mean volume gives a different coefficient. [Derivation and independent integrations](GEOMETRY.md)

## Prediction and controlled comparisons

Use the ordinary rail switch alone to fix the conversion from regional entropy response to area loss. Freeze it and independently vary the fraction of pairs prepared with both original marginals but independent heads. For each unit increase in this fraction, the prediction for the individual-volume area coefficient is between 0.00002898 and 0.00003033. The radar reconstruction gives 0.0018187 to 0.0019239. Calibrating on regional energy instead predicts less than 0.00000000480. The aperture, enrollment rule, clocks and normalization are held throughout. [Ordinary calibration and held predictions](native.json)

A second experiment removes dependence on the choice of coupling. Use the ordinary rail mixture and the angle of a right-head marginal rotation as two independent controls. The determinants of the two control responses exclude zero. The finite-dimensional implicit function theorem supplies smooth curves on which both specified readings are exactly constant for sufficiently small positive re-pairing. The reference rail weights lie inside the permitted interval and the reference head rotation admits both signs. The following are outward enclosures of first derivatives along these physical curves, per unit re-pairing at the reference. Energy is in native exchange units and entropy is in nats. [Source-level calculation](native.py) · [Independent checks](test_native.py)

| Two readings held fixed | Independent response | Enclosed derivative |
|---|---|---:|
| Individually matched mean area and regional count | Regional energy | −0.2690 to −0.2393 |
| Same pair | Regional entropy | +0.01354 to +0.04251 |
| Common-radius mean area coefficient and regional count | Regional energy | −0.253670 to −0.253662 |
| Same pair | Regional entropy | +0.02406 to +0.03317 |
| Regional energy and regional count | Individually matched mean area coefficient | +0.0017486 to +0.0018577 |

Fixing the individual area coefficient fixes the entire averaged area at each sufficiently small prescribed volume. Fixing regional count makes the energy contrasts invariant under every additive choice of energy origin. Any proposed conversion of this regional entropy to a shrinking-ball source with a common nonzero gain fails the matched-area comparison. Identifying a different regional algebra with those proper balls would require independently enrolling its modes and measuring its response.

A finite replay uses a re-pairing probability of one in ten million. Solving the two geometric and number conditions before opening the energy and entropy records gives a right-head turn of 0.00024524212 radians. Directly assembling the positive input state gives a regional energy change of −2.53725 × 10⁻⁸ exchange units and an entropy change of +2.83630 × 10⁻⁹ nats. The secant slopes, −0.253725 and +0.0283630, fall inside the two independent derivative enclosures. Doubling the quadrature used for the area solve gives a coefficient residual below 10⁻²⁵ in the numerical difference calculation. The implicit-function proof supplies the exact matching claim. [Finite input and readouts](finite.json) · [Independent nonlinear replay](finite.py)

## Would a different radial weight fix the source?

Sort the twelve aperture ports by graph distance from port 0. The four classes contain one, three, three and five ports. Read number and half-bond energy independently on each. Two positive lotteries of native marginal rotations and re-pairings give **exactly the same eight readings** when admixed in equal proportions to the paired reference. The derivative of the difference between their individually volume-matched area-loss coefficients lies between 1.9601 × 10⁻¹¹ and 2.7589 × 10⁻¹¹ per unit common admixture. The input direction was constructed to match the eight readings, so its small numerical size bears the cost of realizing it as a positive lottery. Its strictly positive lower bound gives distinct geometric coefficients at sufficiently small positive admixtures. [Interval measurements and mixture weights](shells.json) · [Construction](shells.py)

Any function of those eight identical readings gives the same result for the two sources, including arbitrary signed radial weights and every common energy origin. This argument is specific to this twelve-port partition by graph distance. A ball region chosen by an independently established proper metric could have other modes, boundary couplings or coherences. The separate finite-angle area calculation describes exactly what such a readout would have to reproduce.

## A necessary bound on isolated hinges

There is a stronger pointwise test. An observer on the interior of a hinge has fixed-volume area \(A=4\pi\ell^2(\alpha/2\pi)^{1/3}\). The first hinge on the measured star gives \(-0.003763158 < (\partial A/\partial t)/\ell^2 < -0.003763156\) when \(t\) increases the re-paired fraction. An entropy slope bounded by a constant times \(\ell^4\) cannot cancel this response at any fixed positive area coupling for all small radii. The finite fourth-power ensemble average and this second-power individual response occur on the very same star.

More generally, suppose the same ball remains in a product hinge neighborhood under refinement. Write \(a,s\) for its area and modular variations and assume \(|s|\le S\ell^4\) and \(|s+\eta a|\le\epsilon\ell^4\) with \(\eta>0\). The derivative of the exact area in α has magnitude at least \(2\ell^2/3\) on \(\pi<\alpha<2\pi\). The mean value theorem gives the same lower bound times the angle contrast for finite changes in that range. The triangle inequality then gives

$$2\eta\,|\Delta\alpha|\le 3(S+\epsilon)\ell^2.$$

Thus a refining hardware ensemble retaining isolated hinges must suppress their source-dependent angle changes at least this fast if the modular slope stays bounded and the coupling stays positive and bounded below. With fixed point locations for finitely many local material generators, a center on an open edge can also be chosen a positive distance from every generator. Its ball eventually has the trivial interior algebra and zero entropy response in that dictionary, even though the geometric area responds. A refinement where many hinges or material modes enter each ball calls for its own volume, entropy and uniform error calculation. A weighted second-moment bound quantifies the loss of uniformity when a small fraction of centers carries the whole mean. [Lean bounds](../../../../../../Lean/ObserverPatchHolography/EinsteinBranch/StationarityReadout.lean) · [Observer sampling requirements](GEOMETRY.md#weak-deficits-and-observer-sampling)

The field equation uses direction-dependent source and curvature readings on a common family of proper balls. The scalar comparison here rules out the specified pairing of averaged spatial geometry with the twelve-port aperture. Proper clocks, directional sectional responses, the physical shrinking-ball algebra and conservation of its total stress remain inputs to the conditional theorem. [Full hypotheses and error rates](../../../../PROOF.md#geometric-normalization-and-the-three-errors)

## Reproduce the measurements

From the RER root with the model Python environment active, run

~~~sh
export PYTHONPATH=code:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m matched_einstein_limit.rate_geometry.seam_strain.weak_ruler.center_average.native --write
python -m matched_einstein_limit.rate_geometry.seam_strain.weak_ruler.center_average.finite --write
python -m matched_einstein_limit.rate_geometry.seam_strain.weak_ruler.center_average.shells --write
python -m pytest -q -W error code/matched_einstein_limit/rate_geometry/seam_strain/weak_ruler/center_average
~~~
