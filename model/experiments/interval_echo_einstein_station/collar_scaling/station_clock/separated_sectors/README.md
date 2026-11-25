# The regional entropy of two correlated face sectors

A protected observer can read the entropy response from the cyclic order of four native seam energies on two different stations. In a shrinking corner ball the signed entropy slope is **0.00133857665090 nats per fourth power of the relative radius**. Every sum of regional word observables with sublinear range and bounded total strength at each station has vanishing slope at that order. The two statements concern the same input, including preparations that place both words outside the ball. [Records](results.json) · [Native station survey and path counts](../README.md)

## A preparation independent of the observation radius

Use the seven-appointment station clock and two distinct uniformly placed ordered words of the protected collar experiment. Give each original icosahedral face one fair writer tag. Use the home face A and the smallest indexed face B at distance two from A around its corner pentagon. They have indices 0 and 3 in the midpoint atlas. The other eighteen tags are independent of these and each other. In the reference the two selected tags have opposite values. With probability t replace the tag on B by an independent fair tag, before drawing the word addresses. Each word is written independently given the tag on its own face. Fix these two faces at every resolution and for every radius. Trace away the preparation tags and retain the addressed orthogonal word slots, including their blank flags. The geometry survey avoids the slots and is constant along t. [Word preparation and route](../../PROTOCOL.md) · [Enumeration of native word probabilities](../../../shared_variance.py)

Write Q for the conditional joint word table with fresh tags, P for the conditional table with reversed tags, and D=Q−P. The table at one address on each selected face is P+tD. Each row and column sum of D vanishes. A pair on the same face keeps its shared tag throughout. All other pairs have an unchanged law. For any region containing a stations from A and c from B among N stations, the probability of a placement affected and fully admitted is

$$
 p=\frac{2ac}{N(N-1)}.
$$

At fixed addresses, seeing only one of the two words leaves the local distribution unchanged. Since the occupied slots are orthogonal and the address lottery is independent of t, the derivative of the complete regional von Neumann entropy at t=0 is exactly p times the entropy derivative of P+tD.

## A second read using the measured seam energies

A single word contains two ordered seam energies from 0, 32 and 33 in exchange units. Give the ordered pair score +1 on the cycle 0→32→33→0, −1 on its reversal, and 0 on equal energies. For a placement with one word in each admitted sector take the product of their scores. Assign zero to every other trial. This observable requires a comparison of the two addressed energy records. Its definition uses the calibrated energy labels alone.

On an eligible placement the 9×9 native word tables give a score derivative of 1/36. For the entropy, Q is a product and the fixed marginals give ΣD log Q=0. In units of 1/2592 the positive entries of D contribute 18 units at each of P/Q=35/36 and 3/4. Its negative entries contribute −36 units at P/Q=13/12. Hence the conditional entropy derivative and the relation between the two independent scores are

$$
 \kappa=-\sum D\log(P/Q)=\frac{\log(169/105)}{144},\qquad
 S'(0)=\frac{\log(169/105)}4\,\frac{d}{dt}\mathbb E[\mathrm{score}].
$$

The sign and gain follow from the native writer counts. On the archived independent fresh and shared word cohorts the conditional signed score difference is **0.02777628 ± 0.00005283** at simultaneous 99% coverage, compared with 1/36. Shared minus fresh equals fresh minus reversed entry by entry in the native population. Each cohort contains 2³⁴ simulated trials. The regional prediction multiplies that conditional slope by the separate station placement probability. The comparison of the addressed energies was chosen after inspecting the word law. [Cohorts](../../../shared_variance.json) · [Readout and grouping](measure.py)

## How much of the response can a sum of local terms carry?

Let R be the clock-ball radius in completed station hops, H the eccentricity of the corner and ℓ=R/H. On the midpoint refinement use R²=2m with subdivision frequency m, N=20m² and H/m→6. Up to radius 2m−2 the counts in the two faces on the whole ball are F(R) and F(R−2), where F(s)=floor((s+2)²/4) for s≥0 and is zero otherwise. Their leading term is R²/4 each. On the annulus (R/2,R] subtract the corresponding inner counts. The leading term on each face is then 3R²/16. The resulting limits are

| Region | p/ℓ⁴ | Orientation slope/ℓ⁴ | Entropy slope/ℓ⁴ in nats |
|---|---:|---:|---:|
| Whole ball | 81/200 | 9/800 | 0.00133857665090 |
| Annulus | 729/3200 | 81/12800 | 0.00075294936613 |

There is a geometric lower bound on the separation of the two words. Give the two outer vertices of A height m and all other original vertices height zero. Interpolate linearly on every original face and assign each microtriangle its largest vertex height. Adjacent microscopic vertices have integer heights differing by at most one. Adjacent microtriangles share two vertices, so their assigned heights differ by at most one even across original edges. The height on B is identically zero. In A a triangle at sector indices i,j and k=i+j has height k+1 if upward and k+2 if downward; its distance from home is 2k or 2k+1 respectively. Hence for any x on A in the ball and any y on B, d(x,y) is at least the height at x, which strictly exceeds half of d(home,x). In the annulus every such pair has separation greater than R/4. [Independent shortest-path and height tests](test_sectors.py)

For any set of annulus stations of diameter at most R/4, its entire reduced word state is constant along t. Condition on each ordered pair of written addresses. An affected pair has at most one word in the set; either partial trace of D vanishes. This proves equality of the states for every t. Any weighted sum of observables on such sets has identically constant mean, with arbitrary signs and position-dependent weights.

The whole ball includes the vertex where the sectors meet. Let a regional observable be a sum of possibly overlapping operators O_X on word sets X of diameter at most b, with incidence strength J=supₓ Σ_{X containing x} ‖O_X‖. If both affected addresses x in A and y in B occur in one X, then d(home,x)≤2b and d(home,y)≤3b. Write B(h) for the number of stations within h hops of home. There are at most 2 B(2b) B(3b) ordered placements of this kind. On each such placement the derivative of a term's mean in absolute value is at most ‖D‖₁‖O_X‖=‖O_X‖/36. The sum of the norms of terms containing both addresses is at most J. Thus

$$
 \left|\frac{d}{dt}\left\langle\sum_X O_X\right\rangle\right|
 \le \frac{2J B(2b)B(3b)}{36N(N-1)} .
$$

The corner count gives B(h)=O((1+h)²), and N(N−1)ℓ⁴ is asymptotic to a positive constant times R⁴. For b=o(R) and uniformly bounded J the ratio of this slope to ℓ⁴ tends to zero. The estimate also covers refinement-dependent normalizations when J(1+b)⁴/R⁴→0. It applies to any calibration of the local terms satisfying that bound, including a fixed ordinary-matter conversion factor.

This result quantifies the spatial range of an operator reproducing the collar modular variation on the protected word algebra. A finite-range assignment of energy or stress on this algebra misses a leading term that the two-station comparison measures. The construction uses an externally selected, anticorrelated reference and a clock ball on the two-dimensional station graph. Identifying a boost stress on proper three-dimensional balls requires a reference, metric and coupling supplied by a physical source. Other material degrees of freedom, or a transformation that represents widely separated station addresses as nearby in the observer's proper metric, fall outside the range estimate.
