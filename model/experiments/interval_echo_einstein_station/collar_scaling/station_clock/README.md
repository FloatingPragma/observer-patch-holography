# An all-resolution count of clocked collar admissions

A protected collar word can be visible to an observer only if its station is in the observed region. On the spherical midpoint family, clock each adjacent station transfer in seven addressed appointments. For a corner observer and the shrinking annuli below, the one-word admission probability divided by squared relative radius tends to **27/16**. For two distinct uniformly placed words the probability divided by the fourth power tends to **729/256**. The routed signals and region are independent of the word correlations. Their modular contrast has a strictly positive fourth-order coefficient, while any area reconstructed solely from these route records has zero response. [Finite surveys](results.json) · [Source laws and readout](../README.md)

## One clock on the same prepared source

Retain one read station at port zero per cell, both original written words, and the protected port four. During a crossing, send a sharp load from the station to the outgoing gateway, across the external seam, and from the other gateway to its station. The accessible internal template has diameter three. Each crossing therefore takes at most seven swaps. Fill any remaining appointments by servicing a spare blank contact disjoint from the stored tag and protected word. Its two zero loads are a stationary vector of the native contact. The separately calibrated single-swap hand advances by one duration π per appointment. All transfers start and finish at a station. Restrict the route search to these complete crossing services. A path with h cell crossings takes exactly 7h appointments; the shortest such path measures the dual graph distance. The single- and double-load probes transfer at the same hand reading. Doubling exchange on the blank two-load sector preserves the single-load hand but prevents transfer of the double-load probe. [Native transfer and changed-law check](measure.py) · [Underlying hand](../clock_balls.py)

This is an additional clock protocol on active Interval Echo. The original fastest protected route permits direct gateway-to-gateway services and variable crossing durations. Its registered four-board results stay measurements of that instrument. On the three bundled boards the program checks every directed crossing against the actual labelled contacts, counts the padding and compares the entire rooted distance histogram to a separate subdivision. At higher levels the calculation uses the stated midpoint-incidence continuation. It applies to any labelling with three distinct accessible gateways per cell. Every limit below refers to this declared family of finite instruments.

## The path count

Write m for the subdivision frequency. A face of the original icosahedron contains m² microscopic triangles, each indexing a station. There are N=20m² stations on the sphere. Select a microscopic corner triangle touching a distinguished original vertex as home. Five original faces meet there. In each of these sectors index an upward triangle by nonnegative i,j, with its vertex nearest the corner at lattice distance k=i+j along the triangular grid. The downward triangle paired with that upward triangle exists when k≤m−2. Give them respective levels 2k and 2k+1. At each available level t the number of triangles in one sector is floor(t/2)+1.

A dual step inside a sector changes the level by one. Each triangle can reach its own corner triangle in exactly its level many steps by descending. Crossing a radial side between sectors keeps the level and changes the sector by one. The five corner triangles themselves form a cycle. Consequently a station in sector i at level t has distance t+min(i,5−i) from home as long as the route cannot leave the five sectors. Any route that leaves first reaches level 2m−2 and then crosses the opposite edge. Balls of radius at most 2m−2 are therefore counted exactly by

$$
 B(R)=\sum_{j\in(0,1,1,2,2)}F(R-j),\qquad
 F(x)=\lfloor(x+2)^2/4\rfloor \quad(x\ge0),\qquad F(x)=0\quad(x<0).
$$

The sum counts sector offsets with multiplicity. In particular B(2k)=5k²+4k+1. For R divisible by four the population of the annulus (R/2,R] is exactly 15R²/16+R. This calculation fixes its coefficient before reading either word.

The relative outer radius is ℓ=R/H, where H is the greatest station distance from home. There is an independent global calibration. Give the original icosahedron vertices their graph distances 0,1,2,3 from the distinguished vertex. The five nearest vertices have height one, the next five height two and the opposite vertex height three. Multiply by m and interpolate linearly on each original face. Every microtriangle has heights (k,k+1,k+1) or (k,k,k+1). Assign these types potentials 2k and 2k−1 respectively. Across a dual edge the potentials differ by at most one, including at original seams. From home to a triangle at the opposite vertex their difference is 6m−3. Hence H≥6m−3.

Conversely, each original face has a vertex at original graph distance at most two from the origin vertex. Follow at most two original edges to that vertex and then cross the target face. Each leg within one face takes at most 2m−2 microscopic steps by the sector count. At the beginning and at each change of face, at most two steps around a corner pentagon select the required sector. Thus H≤3(2m−2)+6=6m. In particular H/m→6. The finite calculation measures H independently by breadth-first distances; the limit requires no fitted distance standard.

## Consequence for the finite theorem

On the odd subdivision levels choose m=2·4ʲ and R=2·2ʲ, so R²=2m. Then ℓ tends to zero, the number of admitted stations grows, and 1/(Nℓ²) and 1/N both tend to zero. From j≥1, put t=1/R and b=H/m. The single admission divided by ℓ² is (15/16+t)b²/20. For the two ordered distinct stations the exact normalized probability is

$$
 \frac{n(n-1)}{N(N-1)\ell^4}
 =\frac{(15/16+t)(15/16+t-t^2)b^4}{400(1-t^4/5)} .
$$

Taking t→0 and b→6 gives the two limits. The Lean file checks the finite rational identity and its limit. The sector and distance arguments above connect those algebraic parameters to the actual clock schedules. [Lean normalization](../../../../../Lean/ObserverPatchHolography/EinsteinBranch/RoutedAnnulus.lean) · [Original pair-admission theorem](../../../../../Lean/ObserverPatchHolography/EinsteinBranch/CollarAdmission.lean)

Multiply the pair limit by the independently determined connected modular contrast of the original shared and fresh tags at the reversed-tag reference. The result is about 0.00941 nats per fourth power of relative radius. It is positive even using the rigorous lower bound 101/33696 nats on the conditional contrast. A strictly positive coupling places this measured collar contribution at leading order. No area depending only on this clock survey responds to this change of tags. Removing a discrepancy of that order would require an independently measured material or geometric term at the same order. This instrument resolves two-dimensional screen admission. Its clock surface has no established identification with a proper three-dimensional ball or a universal gravitating stress.

## What follows for the original fastest-swap clock

The original, variable-duration clock also gives a rigorous result for its *annuli*. Consider any continuation with the protected midpoint carrier and any regime with R→∞, R/H→0. The original closed-ball estimates give cR²≤M(R)≤CR², uniformly in the center when R is sufficiently large and R≤m. To obtain a lower bound for an annulus, take a shortest protected path to a farthest station. At floor(3R/4) swaps choose the station in the current cell. It is at distance at most three swaps from that point. Its entire closed station ball of radius floor(R/8) then lies in (R/2,R] for R≥64. The upper bound follows by containment in the outer ball. Since H is comparable to m, this supplies upper and strictly positive lower bounds for n/(Nℓ²) on the original annuli. Compactness now gives a cofinal subsequence with a positive limiting single coefficient; the existing pair theorem yields its square there. The same argument works for any fixed finite set of centers on a common subsequence. The value and uniqueness of that coefficient for the full fastest-swap sequence remain open. [Uniform closed-ball bounds and their assumptions](../README.md#which-metric-do-these-clocks-determine-at-large-scales)

To regenerate both complete-clock controls and all graph counts, run the following from the RER root.

~~~sh
python -m model.experiments.interval_echo_einstein_station.collar_scaling.station_clock.measure --write
python -m pytest -q model/experiments/interval_echo_einstein_station/collar_scaling/station_clock
~~~
