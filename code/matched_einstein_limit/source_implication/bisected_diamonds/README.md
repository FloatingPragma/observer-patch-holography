# Null Ricci from three bisected clock intervals

A geodesic clock passes through the chosen event at proper time zero. Count the markers in the causal interval between its readings at −h/2 and h/2, and in the two intervals obtained by splitting its worldline at zero. Interval membership requires both transmissions, from the first endpoint and to the second. Write P for the full count divided by the flat four-volume πh⁴/24 and S for the average of the two counts, each divided by π(h/2)⁴/24. Counts are expressed in calibrated volume units. These measurements use causal order, proper clock durations and a common scalar marker intensity. They do not require a transverse optical ruler.

## Geometric inference

Use signature (−,+,+,+) and a convex normal neighbourhood in a smooth four-dimensional metric. Let m be the positive, twice differentiable marker intensity with respect to proper four-volume. Choose three centered timelike geodesics with tangents u in one orthonormal time-space plane. Each half-diamond is centered a proper time h/4 away from the common event. With H the covariant Hessian of m and b = −mR/180 + □m/60 at this event, the leading mean counts are

```
q(u) = m Ric(u,u)/30 + b + H(u,u)/40
P(u) = m + h^2 q(u) + o(h^2)
S(u) = m + h^2 [q(u)/4 + H(u,u)/32] + o(h^2).
```

The curvature coefficient comes from the small geodesic diamond volume expansion. In normal coordinates at its midpoint, the flat diamond has second moments h²/60 along time and h²/30 along each spatial axis, from integrating balls of radius h/2 − |t|. The first moments vanish in flat space. Corrections to the domain and measure from curvature have relative order h², so coupling to a linear change in m contributes at most order h³ in the normalized count. The quadratic Taylor term gives the two Hessian contractions in q. Taylor expansion of m along the midpoint clock gives the additional term H(u,u)/32 in S. Smoothness is required across all three clock directions. [Gibbons and Solodukhin, *The geometry of small causal diamonds*](https://doi.org/10.1016/j.physletb.2007.03.068)

Consequently P − 2S/3 equals m/3 + h²[m Ric(u,u)/36 + 5b/6] + o(h²). For k=(1,n) with |n|=1 use u=(1,0),(5/3,4n/3),(5/3,−4n/3) in that order. The weights w=(−9/16,81/160,9/160) sum to zero and recover X(k,k) from X(u,u) for every symmetric form X on the plane. This gives

```
Ric(k,k) = lim 36 Σ w [P(u) − 2 S(u)/3] / (m h^2).
```

The common central m can be measured from the leading rest count. An error of at most ε in each of the six normalized reads contributes at most 67.5ε/(mh²) for exact normalization. Consistency needs ε=o(h²), or absolute volume error o(h⁶), uniformly on the chosen finite set of directions. Clock-duration errors o(h³) suffice for the leading normalization error with correctly identified geodesics and endpoints. One must also control any difference between the actual clock trajectories and those geodesics. Substituting a consistent estimate of positive m in the final division preserves the limit. The formal coefficient and directional error identities have a [Lean proof](../../../../Lean/ObserverPatchHolography/EinsteinBranch/BisectedDiamonds.lean). The [source-net error theorem](../../../../Lean/Geometry/SourceNetVolumeError.lean) controls flat diamonds for its assigned weighted cover. Applying the curvature measurement to Interval Echo requires its regional clock and count errors on shrinking curved diamonds.

Read the regional ordinary plus collar energy density ρ(u) separately in each of the same three clock frames. The one comparison per ray is

```
E(n) = Σ w [36(P(u) − 2 S(u)/3)/(m h²) − κ ρ(u)].
```

If ρ(u) approaches an independently normalized T(u,u), this converges to Ric(k,k) − κ T(k,k), regardless of whether Einstein dynamics holds. With exact central normalization, full and half count errors at most ε (including expansion remainders), and directional density errors at most δ, the error bound is 67.5ε/(mh²) + (9/8)|κ|δ. It allows correlated errors. On the three-frame design an independent moving-frame energy reading detects a local pressure term invisible at rest. The preexisting nine-ray theorem converts vanishing residuals to tensor closure up to the trace. It requires a common frame, ordinary coupling, and conserved stress. In a Poisson comparison with an expected order N counts per full diamond, normalized fluctuations are of order 1/√N. The condition Nh⁴ → ∞ makes their contribution to the curvature estimate vanish in probability. The common native sampling law and clock errors have to be controlled for an application to Interval Echo.

## Controlled four-dimensional measurement

For U=t+x, W=t−x take the metric (1+cU²)η and intensity 1+νU² near the origin. A Lorentz null rotation fixing U moves each timelike geodesic to constant transverse position. The proper-time equations then give U+cU³/3 = Bs and W=s/B, with B = u⁰+u¹. The conformal factor preserves the light cones. Integrating over W and the transverse disk at fixed U gives the differential marked volume

```
[π (ΔW)^2 / (4 ΔU)] (U−Ua)(Ub−U)(1+cU^2)^2(1+νU^2) dU.
```

The curvature at the origin, independently obtained by differentiating the conformal metric, is Ric(k,k)=−2c(1+nₓ)². Five-node Gaussian quadrature is exact for the degree-eight polynomial integrand in real arithmetic. The replay also draws Poisson counts for three disjoint pieces of each full diamond. It retains the separate counts, seed, exposure and propagated sampling standard deviation. Different direction panels use independent samples. The record is a synthetic geometric comparison.

Along +x in flat space with ν=1/2, the full-diamond-only estimate is +3 and the bisected estimate is exactly zero in expectation. For c=−1/4 with the same mark, the true curvature is 2. At h=1/32 the full-only mean is 5.00241, the bisected mean is 2.00160 and the sampled read is 2.00130 with standard deviation 0.00150. For c=+1/4 the true value is −2, the full-only mean is +1.00065 and the bisected mean is −1.99883 at the same duration. Finite-radius bias and counting noise enter separately. A marker concentration can thus mimic curvature even with arbitrarily many counts from a single scale and direction design.

## Causal order and duration on the stationary 80-cell source

The sharp launch test regards a port as reachable if it has positive detection probability at arbitrarily small positive flight time. In each charge and parity sector the integer off-site hopping matrix has 5,036 identity blocks and four marked directed blocks on two undirected edges. The factors on one edge act only on the first internal head; the other edge acts only on the second. All are nonnegative and nonzero in the native basis. The program checks these facts block by block and verifies connectivity. For ports at graph distance d, the first off-site term of the flight amplitude has order d. Its matrix is a nonzero scalar times the sum of shortest-path products. Each path uses a marked edge at most once, products on distinct factors are nonzero, and the contributions have a common phase. A faithful head preparation assigns strictly positive probability to this amplitude for sufficiently small positive times. Thus each of the 960 ports can be reached arbitrarily early, on either probe. Analyticity of the finite propagator makes any subsequent zeros isolated. A uniform port-time census for an interval defined by positive probability therefore counts 960h, up to sets of zero time measure. Its duration-doubling ratio is two. This conclusion concerns exact probability support, including tiny tails.

A second test uses the measured reciprocal order-three relay lengths, permits waiting at every port and places both tips at port zero. At clock duration h with one marker per port per clock unit the interval volume is Σₓ max(h−2d(0,x),0). The relay costs are enclosed over the certified source interval. Monotonicity bounds every shortest distance and hence each full/half count. The raw integer bounds are in the receipt. The doubling ratios at durations 2, 4, 8 and 16 are approximately 6.224, 10.458, 9.500 and 6.662. On this fixed graph the ratio is exactly two below twice the cheapest incident travel time. A smooth four-dimensional proper-volume measure with positive scalar intensity approaches sixteen at small geodesic duration. None of the six sampled scales verifies that scaling. The continuum test would use a joint limit with many resolved cells per shrinking interval. The tested source supplies neither that limit nor the three boosted, synchronized, geodesic clock reads. A native stress and curvature equality therefore cannot be inferred from this census. [Relay length certificate](../common_front/certificate.py) · [Counts and bounds](results.json)

Run `PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.bisected_diamonds.measure --write` at the repository root.
