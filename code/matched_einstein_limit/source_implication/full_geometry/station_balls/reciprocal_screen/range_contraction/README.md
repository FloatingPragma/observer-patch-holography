# Does an energy write contract the distances reported by clocks?

Four observers exchange signals on the selected 80-cell Interval Echo preparation. Each station records the time of return, including a preset time for a signal that times out. Independently prepared copies supply an energy reading and a product-detector entropy measurement. The ordinary energy write raises entropy and reduces the volume inferred from those clocks. Replacing a correlated right-hand record by an independent draw also raises entropy and reduces that volume, while preserving ordinary parked energy. Both interventions are specified before clocks are read. The number of hardware ports enters none of the reported lengths or volumes. [Simulated counts and intervals](results.json) · [Native clocks, positive preparation and equilibrium certificate](../../../../full_algebra/PROOF.md)

## What the observers measure

The baseline is the published positive lottery on the two parked records, with preparation weights rounded to a fixed dyadic grid. The ordinary intervention is the existing unread reversal-sign operation. It sends each of the two parity blocks to half their sum and raises the stored exchange energy by the change in odd weight. No new object is inserted into the clock survey. The correlation intervention replaces the right record with one drawn from an independent copy of the same source. Its joint state in each parity is the product of the original left marginal in that parity and the full right marginal. It preserves each marginal and the odd weight exactly.

Each of the twelve directed routes produces an original-denominator histogram of completed round trips at ticks 2 to 11, with all remaining trials in the capped bin at tick 12. The effective range is half the round-trip time. The clock tick is half an exchange-time unit. The four clock gains are estimated from ratios of the directed means. The six gain-corrected lengths determine the surface area and volume of the surveyed tetrahedron. These are estimates from records in a particular radar protocol. Independent runs measure the two local Pauli outcomes in each of nine settings together with reversal parity. Their observed frequencies determine entropy and the parked energy above baseline.

| Preparation | Parked energy above the common offset | Entropy (nats) | Area from clocks | Volume from clocks |
| --- | ---: | ---: | ---: | ---: |
| Equilibrium preparation | 0.202432–0.202441 | 1.77727–1.77819 | 14.29817–14.29951 | 2.76290–2.76405 |
| Ordinary energy write | 0.499996–0.500005 | 1.96712–1.96804 | 14.22354–14.22487 | 2.73713–2.73826 |
| Independent right record | 0.202432–0.202442 | 1.88816–1.88908 | 14.27082–14.27217 | 2.75500–2.75618 |

The table rounds the simultaneous 99% counting intervals outwards. All three panels use identical settings. Each count bank is a multinomial simulation with $2^{38}$ starts per route or Pauli setting. Integer propagator bounds control the return-time probabilities. Hoeffding bounds on all 36 timing means and 96 tomography coefficients give the joint coverage. The clock uncertainty propagates through the gain correction and tetrahedron formula with interval arithmetic. The entropy uncertainty uses the trace-distance continuity bound on the two four-dimensional blocks taken together. The reference and written sources are separate ensembles. These data exhibit contraction without a volume-matching feedback operation, an assumed pressure or a fit to a target curvature. [Acquisition and reconstruction](measure.py)

## A finite entropy–length theorem

Let $X_*$ be the exact, independently calibrated interior stationary preparation. Write $r_j^*>0$ for its twelve directed range means and $A,V$ for the area and volume of its radar tetrahedron. At this point the finite theorem gives $dS+\eta\,dA-\lambda\,dV=0$ on every trace-one head direction. All twelve clock means are affine in the source. In nats define
\[
 w_j=r_j^*\big(\eta\,\partial_jA-\lambda\,\partial_jV\big)_*,
 \qquad b_{ij}=w_{ij}+w_{ji}.
\]
The certified interval derivatives on the original equilibrium give $14.71<w_j<23.76$ and $236.8196<\sum_e b_e<237.5429$. Both bounds include the range, ordinary-coupling and volume-multiplier uncertainties. They hold for the exact reference even though the simulated lottery approximates it.

Write $a_i=\sum_j(\log r_{ij}-\log r_{ji})/4$. The clock length is defined by $2\log\ell_{ij}=\log r_{ij}+\log r_{ji}-a_i-a_j$. Thus $\log\ell=L\log r$ is a fixed linear map. Because $A,V$ depend only on the six lengths, the chain rule gives $w=L^\mathsf{T}b$. Indeed scaling both directions of a single edge equally leaves all gains fixed and changes just its log length, so its coefficient is $w_{ij}+w_{ji}$. Euler homogeneity also gives $\sum_e b_e=2\eta A_*-3\lambda V_*$.

For **any** state $X$ on this source algebra with these same clock effects, put $z_j=\log(r_j(X)/r_j^*)$ and $K_*=-\log X_*$. The entropy differential at the reference acts as $\operatorname{tr}(K_*\,dX)$. Since this action and the twelve range effects are linear on trace-zero displacements, the tangent equation gives the *finite* identity $\Delta\langle K_*\rangle=-\sum_j w_j(e^{z_j}-1)$. Combining it with the relative entropy identity gives
\[
 S(X_*)-S(X)
 =D(X\Vert X_*)+\sum_jw_j(e^{z_j}-1-z_j)
    +\sum_e b_e\log\frac{\ell_e(X)}{\ell_e^*}.
\]
Everything on the right except the last sum is nonnegative. Let $B=\sum b_e$ and $R=\prod_e (\ell_e(X)/\ell_e^*)^{b_e/B}$. Then $R\le e^{-[S(X)-S(X_*)]/B}$. **An entropy increase forces this calibrated average of the observers' six log lengths to decrease.** If all six lengths are fixed, the exact reference uniquely maximizes entropy, and the deficit is at least $D(X\Vert X_*)$. On a shape-preserving change the survey volume ratio is $R^3$. For an arbitrary change of shape the measured volume requires all six lengths. The two unconstrained interventions above give its sign experimentally. The proof uses the interior stationarity and the measured signs of the twelve coefficients. A boundary information capacity by itself supplies neither. [Lean proof of the weighted bound and its strict consequence](RangeContraction.lean)

## Which gravity test this completes

This is a finite response of the observers' survey to a prepared change in a measured exchange energy. The clock signals can scatter and the inferred tetrahedron can change shape. A universal mass response requires an independently characterized moving or freely formed source and an observer-calibrated proper ruler across probes and survey sizes. In fact the twelve range effects already have rank twelve on correlation changes which hold ordinary energy and both head marginals fixed. At the interior reference a small independent parity-energy write can therefore be compensated in all twelve raw means. By the same right inverse it can be combined with either sign of a small uniform fractional range change. Positivity of the native source lottery holds on a neighborhood. Hence the sign of the parked energy increment alone cannot determine the sign of the survey volume under *every* admissible preparation. For an explicit shape control, increase both log times of edge 01 at unit speed and decrease both of edge 02 at unit speed. All four clock gains stay fixed. The interval derivatives give $V'(0)>0.7494$ and $S'(0)=b_{02}-b_{01}>16.95$ nats. The rank-twelve correlation directions realize these range derivatives at fixed parity. Starting instead with the trace-preserving block shift $(-I_4/4,+I_4/4)$ and correcting its range derivatives in those same directions gives the same positive entropy and volume slopes with an increasing ordinary parity energy. The exact stationary point has a positive preparation neighborhood, so the signs persist for a sufficiently small positive step. This is a different preparation from the ordinary write in the experiment. [Rank and preparation proof](../../../../full_algebra/reuse/INVISIBLE.md)

To test the local Einstein coefficient, enroll two species of clocks and a transported angular frame on the same source. Record proper transverse areas in twelve signed ray sectors at matched affine distances, over nested observer-defined regions. Separately measure the ordinary and correlation contributions to the regional null stress. Freeze the coupling on ordinary changes, then forecast all held correlation changes and alternate probes without recalibration. Keep shot, ruler and angular errors small enough to resolve the area response at order $r^4$ and bound the residual across the refinement. The existing directional inversion supplies the test in every null direction once these records are available. [Protocol, error budget and controls](../../DIRECTIONAL.md)

Run from the repository root with one numerical thread.

~~~sh
PYTHONPATH=code OPENBLAS_NUM_THREADS=1 python -m matched_einstein_limit.source_implication.full_geometry.station_balls.reciprocal_screen.range_contraction.measure
PYTHONPATH=code OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error code/matched_einstein_limit/source_implication/full_geometry/station_balls/reciprocal_screen/range_contraction
~~~

## The same writes viewed from two additional stations

Enroll ports 5 and 6 against the four original ports $(0,1,38,39)$ on each of the same three prepared stocks. For each of the 30 directed routes, monitor the destination at half-tick appointments, then the return to the emitter. Give all histories unfinished at tick 12 the full timeout score of 3. The two charge signs have equal weight, and the source parity weights come from each stock. The simulator draws one multinomial of $2^{50}$ starts per stock and route from bounded integer first-return effects. The count record retains the ten successful return bins and the terminal timeout bin.

For each stock fit the original four endpoint gains from their own reciprocal clocks. Enroll each outside clock against those same four without refitting the anchors. Score the means in $[0.5,3]$. A Hoeffding half-width of $2.5\sqrt{\log(18000)/(2\cdot2^{50})}$ and a $10^{-10}$ numerical allowance give simultaneous coverage of at least 99% for all 90 means by a union bound. The native timing bound applies per directed route since every interrogation concerns one destination and one emitter with four coin columns each. In the binary64 readout each successful bin sums at most 32 parity-block matrix products. The input and contraction roundoff and the correction of a tiny negative bin to zero fit within the added allowance. Propagate all intervals through the logarithmic gain decoder and the distance Gram determinants before differencing preparations. The same confidence event covers every reported volume difference.

The ordinary write shrinks 12 tetrahedra and expands 3. Decorrelation at fixed ordinary energy shrinks 7 and expands 8. Every sign is resolved on the simultaneous event. Using station indices in the port order above, quartet 0123 has ordinary volume change $[-0.025794,-0.025756]$ and decorrelation change $[-0.007902,-0.007864]$. At quartet 0345 both changes are positive, $[0.000170,0.000222]$ and $[0.017062,0.017114]$. At quartet 1245 the ordinary write increases the volume by $[0.000471,0.000521]$ while the decorrelation change is $[-0.026872,-0.026822]$. All changes use the same survey volume units.

For each preparation the five leading principal minors of the six-station displacement Gram matrix based at station 0 are positive throughout the measurement intervals. Its Euclidean rank is five. Each of the 15 four-station subsets has its own positive tetrahedral volume. The six-station distances exclude a shared flat three-dimensional Euclidean embedding. Adding subset volumes as disjoint spatial cells, or interpreting any of the clocks as a probe-independent proper ruler, requires an additional physical construction.

[Native count acquisition and decoder](held_stations.py) · [Timing histograms, interval volumes and signed changes](held_stations.json) · [Independent Cayley–Menger volume checks](test_held_stations.py)
