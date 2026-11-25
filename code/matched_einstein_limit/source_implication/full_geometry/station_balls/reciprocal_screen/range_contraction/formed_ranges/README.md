# Do both radar probes see the same change when two formed heads are separated?

A freely propagating excess encounters the two heads on the 80-cell Interval Echo board. After four Hamiltonian time units, trace out its position and retain the even-parity two-head stock $F$. An independent control $I=F_L\otimes F_R$ pairs a left and a right head drawn from separate copies. This holds both single-head densities and reversal parity fixed. The acquisition probes these stored preparations with a fresh vacancy or excess. It measures how previously formed correlations scatter a timing signal. A concurrent experiment with the forming load in flight requires a different protocol. [Formation column and normalized head stock](../../../../../full_algebra/native_precision.py)

## Calibration that the comparison cancels

For each charge $q$, stock $s$ and directed station pair measure the mean $r^q_{ij}(s)$ of the capped return. The emitted probe visits the destination and then returns, with position checks every half time unit. Tick $n$ has score $n/4$ and all unfinished histories have score 3 at tick 12. Define $D^q_{ij}(s)=\sqrt{r^q_{ij}(s)r^q_{ji}(s)}$ and the relative shape response $U^q_{ij}=D^q_{ij}(F)/D^q_{ij}(I)$. Consider the common-range hypothesis

$$
D^q_{ij}(s)=c^q_{ij}\,g^q_i(s)g^q_j(s)\,d_{ij}(s).
$$

Here $d_{ij}(s)>0$ is the same separation for both probes. The factors $c^q_{ij}>0$ may differ for every edge and species and are held fixed across the stock change. All endpoint gains $g^q_i(s)>0$ may vary independently with stock and species. The stable edge factors cancel in $U^q_{ij}$. Under the hypothesis, the discrepancy $h_{ij}=\log(U^{+}_{ij}/U^{-}_{ij})$ has the form $a_i+a_j$. Thus on each quartet $i,j,k,l$ all three disjoint-edge matching sums must agree,

$$
h_{ij}+h_{kl}=h_{ik}+h_{jl}=h_{il}+h_{jk}.
$$

These equalities are sufficient for a positive common edge response up to the stated factors. The complete six-station graph has 15 undirected edges. Its endpoint map has rank six. The 45 differences between matching sums have rank nine and annihilate that map, so their common kernel is exactly its image. Equality would establish algebraic compatibility on these two stocks. A Euclidean embedding and physical proper lengths require additional measurements. [Exact integer rank and independent calibration test](test_formed.py)

## Native observations

The 30 directed routes use the same ports $(0,1,38,39,5,6)$ for both species. The integer transfer kernels enclose every true four-coin subblock and incorporate the pre-timeout first-hit recurrence. We evaluate all 120 directed means, one per route, stock and charge, using exact rational contractions and outward intervals. The source column has vector error below $10^{-20}$ before final fixed-point rounding. Normalization, contraction of trace distance under partial trace, and the triangle inequality give an allowance of $10^{-18}$ in trace norm for the normalized formed head stock and twice as much for the product. Every capped effect has norm at most 3. The code propagates the native effect error and the larger stock allowance through every reported quantity. The receipt contains mathematical intervals on expectations. They have no finite-trial sampling term. [Native timing error bound](../../../../../isometric_collar/certify.py) · [Directed readings and all 45 comparisons](results.json)

All 45 comparisons exclude a common shape response, and 18 have opposite response signs. They measure nine independent linear constraints. Write $C^q$ for a chosen two-edge-product ratio, and $T^q=C^q(F)/C^q(I)$. The largest absolute discrepancy uses quartet indices $(1,2,3,5)$, the physical ports $(1,38,39,6)$, with $C=D_{12}D_{35}/(D_{15}D_{23})$. The enclosing intervals are

| Quantity | Lower | Upper |
| --- | ---: | ---: |
| $\log T^{-}$, vacancy | 0.000384312747 | 0.000384312758 |
| $\log T^{+}$, excess | 0.001396186444 | 0.001396186456 |
| $T^{+}/T^{-}$ | 1.001012385805 | 1.001012385824 |

For any six endpoint corrections $a_i$, each matching difference involves four residuals $h_{ij}-a_i-a_j$. Hence the maximum absolute residual over edges is at least one quarter of the observed absolute log discrepancy, greater than $0.0002529684$. At indices $(2,3,4,5)$, or ports $(38,39,5,6)$, the log response of $C=D_{23}D_{45}/(D_{25}D_{34})$ is $+0.0001021881$ for a vacancy and $-0.0001370449$ for an excess. The two signs refer to the same two stored preparations.

The result excludes this common multiplicative range interpretation of the unpooled six-unit mean clocks on the formed-stock comparison. Histograms acquired at multiple deadlines can reconstruct the same prescribed pooled target from either charge. A correction depending on material along a route, stock and probe, or a separate physical ruler, would be a different measurement. Assigning shape change to null curvature or gravitating mass requires the corresponding source and geometry measurements. [Existing four-station shared radar estimators](../../../../../shared_radar/README.md)

From the repository root, replay the native proof with

~~~sh
PYTHONPATH=.:code OPENBLAS_NUM_THREADS=1 python -m matched_einstein_limit.source_implication.full_geometry.station_balls.reciprocal_screen.range_contraction.formed_ranges.measure
~~~
