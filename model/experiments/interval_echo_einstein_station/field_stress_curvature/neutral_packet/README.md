# A propagating neutral field predicts the endpoints of two material probes

A receiver made from the same field as the source reads a phase that every passive ring-energy and current history misses. On two fresh held mixtures, this reading predicts all 96 labelled terminal probabilities for the previously calibrated vacancy and excess probes. Reversing the phase labels misses 43 of 48 in each mixture. The test uses the one-cell, 458-word neutral sector of the written-field extension. Its scope is the two stated preparations and the two existing charged instruments. [Neutral pulse protocol](PROTOCOL.md) · [Held transfer protocol](TRANSFER_PROTOCOL.md) · [Earlier charged survey and energy records](../README.md)

## Interference read

The first local head phase at (4,11) prepares $c g+i\sigma s v$ at the time the charged survey starts. Here $g$ is the uniform ground, $v=e^{-iH/2}Z_{4,11}g$, $c=\cos(\pi/8)$, $s=\sin(\pi/8)$, and $\sigma=\pm1$. Each positive ring projector annihilates $g$. Every passive sequence of ring-energy/current observations therefore gives exactly the same law for the two phases, as shown in the [parent calculation](../README.md#which-source-data-agree-exactly). To interrogate the field, a separate cohort receives the phase $c+i\tau s Z_{10,11}$, $\tau=\pm1$, using the same unitary operation as the source preparation on a different head. Let $U=e^{-itH}$ be the uninterrupted free evolution after the second pulse and $h_f$ the detector's energy projection. The amplitude reaching its click port is

$$
 i\sigma cs\,h_fUv+i\tau cs\,h_f U Z_{10,11}g
 -\sigma\tau s^2\,h_f U Z_{10,11}v.
$$

The no-pulse vacuum term vanishes. In the four-setting mixed difference $D_f=\sum_{\sigma,\tau}\sigma\tau P_{\sigma\tau,f}$, the last term and both individual intensities cancel. The surviving coefficient is $8c^2s^2=1$, so

$$
 D_f(t)=\operatorname{Re}\langle Uv,h_f U Z_{10,11}g\rangle.
$$

Independent single-packet ring counts for $v$ and $Z_{10,11}g$ bound this overlap by the geometric mean of their click rates. The two probe signs on the vacuum have equal ring probabilities at every delay. The source coherence becomes readable because the second head phase can mix the ground and excited subspaces. This is active field spectroscopy, with energy supplied by the applied pulses. If the detector does not contain the probe head, $h_f$ commutes with $Z_{10,11}$ and $h_f Z_{10,11}g=0$ at zero delay. Later changes then require propagation from the probe. This conclusion concerns the finite-time field response; no distance or speed is assigned to a hop.

Each counted setting uses $2^{27}$ simulated, aggregated binary preparations. Across 640 banks, Hoeffding's inequality and a union bound give 99% simultaneous coverage with raw probability radius $0.000209306$ and four-setting radius $0.000837224$. All twenty faces and four delays were included. There are 18 faces that omit the probe head, of which 12 share no probe vertex. The numbers of resolved contrasts on the 18 faces are 0, 2, 4, and 7 at delays 0, $3/16$, $1/2$, and 1. For the 12 vertex-disjoint faces the counts are 0, 0, 0, and 2. At delay one the disjoint face (2,3,4) gives $0.00332907$, with interval $[0.00249184,0.00416630]$ rounded outward. All vacuum controls include zero and every measured contrast agrees with the independently propagated matrix element. [All ring and single-packet click counts](counts.json) · [Finite evolution and inferential checks](measure.py)

## A held test of the field-to-probe calibration

Choose face (2,3,4) after the first experiment, then freeze the held protocol before collecting the new banks. The two unknown preparations are classical mixtures of the positive and negative *source phases*. All four predicted and measured instruments use separate randomized copies of each source. A neutral receiver at delay one and a material probe at dwell $3/16$ need not interrogate the same copy. The first preparation time, head phase, initial field, angular launch ports and respective evolution laws agree with the archived pure standards. The neutral field uses the ring Hamiltonian. Introducing a vacancy or excess accesses the compensated matter-field Hamiltonian of the parent experiment. The two Hamiltonians apply to different probe operations on copies of the same starting field. Each party can return any of the twelve labelled ports, including the launch port. The simulations assume the supplied written-field extension, coherent addressing, and the calibrated exposures.

With measured neutral differences $a_0,a_1$ on the two pure sources, the true difference on a mixture with positive fraction $p$ is $(1-p)A_0+p A_1$. Each measured difference is accurate within $2e$ when both raw rates are accurate within $e$. Consequently the observed held difference lies within $4e$ of $(1-p)a_0+p a_1$. Inverting this inequality determines the interval for $p$ from neutral bits alone. Each predicted charged bin is the same convex combination of *its own* two observed calibration probabilities. We sweep the inferred range of $p$, then allow the charged calibration and held bin errors. A union bound covers all 248 marginals used across both held mixtures and the earlier $p=1/2$ control. Each neutral setting has $2^{38}$ independent simulated Bernoulli trials; each charged multinomial has $2^{32}$ trials. One charged terminal histogram has twelve correlated bins and contributes eleven independent degrees of freedom. The 96 bins count the whole set of outcomes, without treating their tests as independent.

| Source read | Positive fraction, simultaneous 99% | Held bins covered | Bins missed after swapping pure standards |
|---|---:|---:|---:|
| First mixture | [0.24469, 0.25535] | 48/48 | 43/48 |
| Second mixture | [0.74475, 0.75541] | 48/48 | 43/48 |

The hidden preparation fractions were $1/4$ and $3/4$. The decoder does not take these fractions as input. The greatest probability shortfall of the swapped forecast exceeds 0.07260 in each mixture after all allowances. The earlier independent equal-mixture control covers 48/48 bins when its known half weight is used. The source readings, pure charged calibrations and held charged readouts concern the neutral field and its sign-dependent influence on two materials. The one fraction is identifiable because the pure neutral standards have distinct measured responses. For unrestricted density matrices a single difference observable leaves many states with the same expectation; no claim of arbitrary-source reconstruction follows. The neutral read requires one more freely evolved exposure after its probe introduction. The predictions compare separately prepared cohorts and do not imply an earlier operational forecast in a single history. [Held mixture counts](transfer_counts.json) · [Register-only decoder and changed assignment](transfer.py)

**A concrete ambiguity outside the family.** The vacuum has exactly zero signed neutral response. In the pure-source calibrations the negative-phase difference is strictly below zero and the positive-phase difference strictly above. Convexity therefore supplies a genuine mixture with exactly zero difference. Inverting the two calibrated means with their errors places its positive fraction in $[0.0123641,0.0176876]$. That whole range misses the independently deposited vacuum charge histograms in **46 of 48** bins, by more than 0.06469 probability in the worst one. The retrospective challenge has a separate simultaneous 99% bound on its four neutral and 144 charged marginals. The observable that aliases vacuum and mixture is the signed difference used by the predictor. Keeping the two pulse probabilities separately tests this ambiguity.

On the same two pure sources the earlier survey measured opposite-phase changes of $+0.092388$ and $+0.011183$ in vacancy and excess transport curvature. The positive held result supplies a physical field read for forecasting what these two calibrated material instruments will register. It does not make the directional ring current a full stress tensor, identify the transport distance with a proper metric, or calibrate a positive Newton constant against ordinary matter. The independent [three-dimensional chart test](../../field_geometry_equilibrium/README.md) missed held distances by 35–96%. The experiment has one finite mesh, no sequence of fixed-proper-volume balls and no calibration relating the measured ring energies to an Einstein tensor. The committed counts are synthetic draws from the specified finite Born probabilities.

## Reading the vacuum fraction independently

Each of the three calibrated preparations gives a point in the plane of the two neutral click rates. Independent mixtures occupy the triangle with those vertices. The calibration uncertainty bounds the determinant between $2.8047\times10^{-5}$ and $2.8154\times10^{-5}$, keeping the triangle nondegenerate. Its infinity-norm condition number is at most 10.286. For each held cohort the decoder uses only its two neutral counts to enclose the three preparation fractions. Over the three-source family the second rate resolves the vacuum ambiguity of the signed difference.

In the two earlier held cohorts the upper vacuum fractions are 0.000239 and 0.000171. In the held interior cohort the inferred negative-phase fraction is 0.24720–0.25271, the positive fraction 0.12168–0.12838 and the vacuum fraction 0.62435–0.62568 (rounded outwards). They enclose the independently prepared proportions of 1/4, 1/8 and 5/8. The charged interval forecasts overlap the independent counts in all 144 bins across the three cohorts. Cyclically permuting the three charged calibrations excludes 47, 46 and 44 of the 48 bins, respectively. The joint coverage bound includes 938 marginals and exceeds 99 percent. This prediction uses the calibrated three-source family and linear mixing of its terminal probabilities. [Frozen protocol](HULL_PROTOCOL.md) · [Independent counts](hull_counts.json) · [Fraction and endpoint calculation](hull.py)

Each probability has a finite-count Bernstein interval. For $k$ successes in $n$ trials and $L=15$, with failure probability at most $2e^{-15}$ the error in the expected count is bounded by $20+\sqrt{400+30k}$. This follows by replacing $np$ with $k+|k-np|$ in the two-sided Bernstein inequality. Rounding this radius upward in integer counts covers all 938 marginals simultaneously with probability above 99%, including all 640 settings considered when the receiver was chosen. For each held cohort the decoder intersects the probability constraints with the triangle of nonnegative weights summing to one. The resulting polygon contains the true weights on that event. The minimum and maximum of each charged endpoint forecast occur at its vertices, with the appropriate upper or lower calibration probabilities. Independent copies provide separate held charged histograms.

Run from the repository root with single-threaded BLAS.

~~~sh
p=model.experiments.interval_echo_einstein_station.field_stress_curvature.neutral_packet
python -m "$p.measure"
python -m "$p.transfer"
python -m "$p.hull"
python -m pytest -q --pyargs "$p"
~~~
