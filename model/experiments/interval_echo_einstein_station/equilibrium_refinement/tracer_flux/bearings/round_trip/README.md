# Timing a free tag to a second observer and back

Twelve observers recover a common radial scale on the same two-resident source. A tag moves freely between readings, first reaching a second observer and then returning to its sender. In units calibrated by an independent material exchange, the three limiting one-way ranges are **11, 14 and 15**. At two finite shutter settings, independent counts distinguish all three ranges on each of three resident stocks. The three leading directions recovered from the transit times agree with the previously measured angular frame to a subspace-angle sine below **0.040**. Taken as lengths themselves, the clock readings require more than **2.19%** distortion to fit any three-dimensional Euclidean chart. [Full model counts and error bounds](results.json) · [Angular experiment](../README.md)

## Two independent observers, clock and source

Use the 20-cell host, the same parked exterior contacts and the admitted 55-word ball of the [two-resident source](../../../CLOCK_BALL.md). The tag has load one and the two residents each have load two. Attempt insertion at the chosen sender. Keep every occupied-site failure. During a dwell the entire three-mobile system evolves with the native Hamiltonian. At the tick, local material comparisons read and undo the sharp registers as in the [original tracer instrument](../../README.md). The remote reader's first tag detection begins the return leg. The first subsequent detection at the sender completes it. The tag and both residents continue from their actual conditional state. There is no reflection command, route choice or replacement word at the far end. Source preparation, monitoring and the deadline are supplied.

Repeat on fresh copies with the sender and receiver interchanged, using the uniform reference, the split stock and the spread stock. The latter two have equal one-site occupancies and equal incoming native energy distributions. The 132 ordered pairs are all measured on every stock at both original short shutters, $t=1/512,1/256$. An isolated native dimer independently measures its exchange probability $q=\sin^2(t/2)$ with $2^{44}$ starts at each shutter. The unit of the reported radial scale is $5q$ per tick, the short-dwell probability for any of five possible exchanges. Bare elapsed oscillator time is also retained. [Acquisition](measure.py)

## Uniform radial calibration

Write $P_{yx}$ for the block of $e^{-itH}$ transmitting a tag from $x$ to $y$ while allowing both residents to move. A resident state $\rho$ with the tag known at $x$ gives next-location probabilities $\operatorname{tr}(P_{yx}^*P_{yx}\rho)$. Each of five neighboring transmissions is $it/2$ times a resident permutation plus $O(t^2)$. At any other site the block is $O(t^2)$. These estimates hold uniformly on the finite resident space. Measuring the current location at every tick therefore gives the simple, symmetric five-neighbor walk in the short-dwell limit after rescaling time. The actual sampled protocol also reads the residents at each tick.

For each receiver $b$, set $h_b$ to $0,11,14,15$ on the four graph-distance classes from $b$. The integer adjacency matrix $A$ satisfies $(Ah_b-5h_b)(x)=-5$ off $b$. Thus $h_b$ is the ordinary mean number of jumps until the receiver is reached. At a finite shutter, the independently calculated free propagator gives an operator bound, for every $x\ne b$,

$$\left\|5q I+\sum_y (h_b(y)-h_b(x)) P_{yx}^*P_{yx}\right\|\le 5q\epsilon.$$

The maximum $\epsilon$ over all ordered pairs is below 0.003917 at the shorter dwell and 0.014901 at the longer one. Summing the conditional drift up to the first visit, initially with a bounded stopping time and then taking its limit, gives

$$\frac{h_b(x)}{5q(1+\epsilon)}\le E[T_b\mid x]\le\frac{h_b(x)}{5q(1-\epsilon)}.$$

The bound holds for every admitted resident density, including coherent densities. Apply it again to the actual resident state at the far end for the return. Consequently every stock has the same limiting mean $h_b(x)$ for *half* the round trip in exchange units. In oscillator time, halving the dwell instead **increases** the expected duration by a factor in $[1.9625,2.0378]$ at any of the three separations. Repeated observation slows the physical motion. A common physical speed cannot equate both bare mean times to the same fixed distance.

For the numerical enclosure the code centers the native Hamiltonian by a scalar and bounds its degree-24 Taylor propagation in column norm. Absolute row and column sums bound the sparse multiplication errors. A bound on each complex Hermitian block includes the Taylor errors in both factors and the floating dot products. A control with zero seam contrast on the same alphabet and graph retains the calibration to relative drift below 0.000080 at the longer dwell. This identifies the role of native exchange and regular coordination. [Drift and dynamics certificate](dynamics.py)

## All-start acquisition and an independent radial frame

Before each dwell the shutter can expire with probability $\lambda=t^2/16384$, implemented by 32 or 30 independent fair bits. Each admitted attempt ends with completion or expiration. A bank of $2^{34}$ independent simulated starts includes both outcomes and the blocked launches. The recorded rates use fresh forward and reverse banks. A discounted first-entrance equation gives the probabilities for the simulation. The second leg begins in the actual measured word. Native transition errors enter through the verified residuals of both linear systems. The resolver norm is at most $H_{\max}=15/[5q(1-\epsilon)]$. The largest certified error in a timeout probability is below $3.49\times10^{-9}$.

For the untruncated two-leg duration $T$, the measured expiry probability is $p=1-E(1-\lambda)^T$. At every stage of the trip the remaining expected time is at most $2H_{\max}$. Thus $E[T(T-1)/2]\le 2H_{\max}ET$, and

$$p/\lambda\le ET\le p/[\lambda(1-2\lambda H_{\max})].$$

The censoring correction is below 0.15%. Simultaneous 99% binomial intervals over the timing standards and conditional timeouts enter the full range analysis. On all six panels, all ranges in the first class fall between 10.95 and 11.05 exchange units, the second between 13.95 and 14.06 and the third between 14.95 and 15.06. These are finite data in the model.

Let $D$ be the symmetrized half-round-trip mean table and $J=I-\mathbf1\mathbf1^T/12$. The observer obtains the three leading eigenvectors of $-JDJ$ directly from the clocks, before any comparison to the angular data. In the common limit $(I-A/5)h=\mathbf1\mathbf1^T-12I$, where the columns of $h$ are the hitting times above. Therefore

$$-JhJ=12(I-A/5)^+.$$

Its three largest eigenvalues are $12/(1-1/\sqrt5)=21.7082\ldots$, followed by 10 with multiplicity five and $12/(1+1/\sqrt5)$ with multiplicity three. The leading projector is precisely one quarter of the angular Gram independently reconstructed from short one-way arrivals. Its eigenvalue is a new radial calibration from clock durations. For a measured interval table, the maximum row sum of its componentwise uncertainty bounds the operator error. The measured spectral separation bounds the sine of the angle between the population and angular subspaces by $\eta/(11.7082-\eta)$, where $\eta$ bounds the distance to the ideal operator. All six panels give less than 0.040. This result supports a common three-dimensional description of the **long-wavelength component** of these finite transit records.

## Can the whole clock table be proper distance?

Any set of three-dimensional Euclidean points with pairwise distances $d_{ab}$ has squared-distance Gram $-J(d_{ab}^2)J/2$ of rank at most three. Across the six panels the corresponding fourth eigenvalue from actual symmetrized transit times is greater than 40.52 in squared exchange units. Even the least centered eigenvalue is above 22.99. If a three-dimensional fit moves every distance by at most a fraction $e$, its squared Gram differs in operator norm by at most $(2e+e^2)\max_a\sum_b d_{ab}^2/2$. The 99% input intervals imply $e>0.0219$ on every panel. This fractional discrepancy survives every common conversion of clock time to length.

There is a local curved-space test as well. Suppose the measured lengths are minimizing distances in one geodesically convex normal three-ball of proper radius $R$ and $|\mathrm{sec}|\le\kappa$ there. In parallel frames the Jacobi equation gives $\|d\exp-I\|\le\sinh(x)/x-1$ on transverse vectors for $x=\sqrt\kappa R$. Indeed iterating the radial Jacobi integral equation bounds each difference term by $x^{2n}/(2n+1)!$; radial transport is exact. For $x\le1/3$, $\sinh(x)/x\le(1-x^2/6)^{-1}\le54/53$. Every length on this ball is then between $52/53$ and $54/53$ times its normal-coordinate length. The minimizing curve lies inside the ball, and the straight normal-coordinate segment is also admissible. Its Euclidean endpoint distance consequently differs from the proper distance by at most a fraction $1/52$. This is below the measured lower bound. Thus these particular transit readings, multiplied by any common speed, exclude a common geodesically convex ball with $\kappa R^2\le1/9$ at the quoted coverage. Curvature in spacetime, a common ballistic ruler, and the stress sourcing it still require independent measurements.

From the repository root, run

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_einstein_station.equilibrium_refinement.tracer_flux.bearings.round_trip
python -m "$pkg.measure"
python -m pytest -q --pyargs "$pkg"
~~~
