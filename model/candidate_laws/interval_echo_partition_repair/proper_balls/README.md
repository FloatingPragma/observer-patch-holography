# Compatible words on ten volume-matched smooth balls

A single equilibrium distribution of endpoint records balances the first-order entropy and proper-area changes for ten supplied observers on each of three scales. The two geometries give exactly the same fifteen calibrated, native six-station ranges. Their true small-ball areas differ for moving observers. The construction supplies the two smooth metrics, a common rest frame for the material records and the number of endpoint completions as a function of the measured area. Its scope is a finite comparison for the partition-repair proposal. [Two compatible station geometries](../../../../code/matched_einstein_limit/source_implication/full_geometry/station_balls/README.md) · [Proper-area integrals](../../../../code/matched_einstein_limit/source_implication/full_geometry/station_balls/timelike_areas/README.md) · [Tomographic frames](../../../../code/matched_einstein_limit/frame.py)

## Which records are in each proper ball?

The metrics are ultrastatic. On the round station chamber their spatial tensors differ by the sign of the traceless perturbation in $B_\pm=\operatorname{diag}(2\pm1,2\mp1,2)$. For an observer moving along unit vector $a$ with Lorentz factor $\gamma$, write $L=I+(\gamma-1)aa^T$. On the observer's geodesic rest slice the chamber coordinates are $x=Ly$, $dt^2=(\gamma^2-1)(a\cdot dy)^2$. The induced metric is $H(y)=L^Th(Ly)L-(\gamma^2-1)aa^T$. It obeys $H(y)y=y$, so radial distance from the common origin is $|y|$. The chamber formula applies throughout the sampled balls. For each sign and frame we solve
\[
V_g(r_g)=\int_{|y|<r_g}\sqrt{\det H_g(y)}\,d^3y=4\pi\ell^3/3
\]
and integrate the tangential area density on $|y|=r_g$ to obtain $A_g$. The area calculation uses geometry alone.

A static material worldline at spatial coordinate $x$ pierces this ball precisely when $|L^{-1}x|<r_g$. Its stored letter persists in time, so observers intersecting a given worldline at different times read the same record. Rest, the six signed axial frames and the three diagonal frames have seven distinct projected shapes. The opposite axial velocities have the same spatial cut in these static sources. The three volume radii are $1/256,1/1024,1/4096$, in range units. All 30 selections read restrictions of **one population of 42 worldlines**, with two endpoints per distinct boundary.

For each boundary and each of 384 fixed directions $v$, compute the radial interval between the two metric crossings $r_g/|L^{-1}v|$. Take the direction maximizing the smaller gap from this interval to the origin and to any other boundary interval, divided by the lower crossing radius. Place the endpoints one quarter of that gap inward and outward. Across all 21 boundaries the smallest relative gap exceeds 0.1116. Direct tests against both metrics confirm that every cut contains exactly one endpoint of its designated pair. Every other pair is wholly included or wholly excluded. All endpoint memberships are identical across the metric switch. The archived coordinates suffice to repeat this check independently of the direction search. [Worldline coordinates, radii and whole-pair indices](results.json) · [Geometry and nested restriction checks](test_balls.py)

## One finite word rule and its entropy derivative

Let $e$ be the designated pair for one cut at scale $\ell$. Fix $\eta=1$ nat per squared range unit for the entire survey. Let $M=\ell^{-5}$ and supply the integer
\[
d_{eg}=\lfloor M\exp(\eta A_g)\rfloor .
\]
Choose $N_e=\operatorname{lcm}(d_{e+},d_{e-})\max(d_{e+},d_{e-})^2$. Conditional on metric tag $g$, a uniform pair index in $\{0,\dots,N_e-1\}$ produces an endpoint $i$ and a shift $k$ by division and remainder with divisor $d_{eg}$. The endpoints store $i$ and $j=(i+k)\bmod(N_e/d_{eg})$. Because $N_e/d_{eg}\ge d_{eg}$ this enumeration has exactly $N_e$ distinct legal pairs, with exactly $d_{eg}$ compatible pairs per observed value at either endpoint.

All worldlines carry the same fair, readable tag $g$. Conditional on $g$ the pair indices are independent and uniform. The local words consist of this tag and the letters on the included worldlines. If $I$ is the set of wholly enclosed pairs, the conditional marginal is uniform on
\[
\Omega_g=(N_e/d_{eg})\prod_{i\in I}N_i
\]
records. Changing the tag preserves each wholly enclosed pair's contribution even when its legal letters change. For **every differentiable normalized perturbation** of the common source distribution on this legal support, the entropy derivative at the fair reference is $\delta S=-\sum_g\delta p_g\log d_{eg}$. This follows by differentiating the full local entropy at local probabilities $1/(2\Omega_g)$. Variations within each tag's words drop out at first order, including correlations between pairs. The area derivative for the same source is $\delta A=\sum_g\delta p_g A_g$, and $\delta V=0$ for each enrolled ball.

For the exact proper areas, $0\le\log M+\eta A_g-\log d_{eg}\le-\log(1-1/M)\le2/M$ when $M\ge2$. Consequently
\[
|\delta S+\eta\delta A|\le (2/M)\,\|\delta p_g\|_1
\le 2\ell^5\,\|\delta p_{\rm global}\|_1 .
\]
All ten errors divided by $\ell^4$ vanish for bounded probability tangents as the radius decreases. The equilibrium can be selected by uniform redraws of the tag and pair indices, or by symmetric mixing of their finite product, carrying the indices through a tag switch. This construction specifies no local means of distributing a changed geometry tag to all the material stations.

These estimates apply to every *finite* truncation of a shrinking survey. With successive radii separated by a factor of four, all smaller crossings eventually lie strictly inside the next scale because the largest boost is three and $r_g/\ell\to1$. At a fixed scale distinct projected shapes have different quadratic forms $L^{-2}$. A generic direction avoids their finitely many pairwise coincidence sets, providing relative separation which persists for sufficiently small radii. Thus further boundaries may be enrolled with the same formula on a sequence of finite boards. A common unrenormalized entropy for an infinite product of internal records is outside this argument. Unenrolled centers, observers and time-dependent fields require additional geometric or dynamical input.

## Held metric transfer and controls

Increase the $+$ branch probability by $t$ from one half and decrease the $-$ probability by $t$. Keep the pair indices uniform, without changing the word rule. The following values divide the signed slopes by $\ell^4$. The same slopes are also the finite symmetric secants through $t=\pm1/8$, because the tag entropy is even in $t$ and the conditional entropy and surveyed area are affine. These are exact distributional formulas evaluated with numerical areas. No observer trial bank is sampled here.

| $1/\ell$ | Frame | $\delta A/\ell^4$ | $(\delta S+\delta A)/\ell^4$ | With the two branch counts exchanged |
|---:|:----|----:|----:|----:|
| 256 | $x$ | −8.974476 | +0.000368 | −17.949320 |
| 1024 | $x$ | −8.938478 | −0.000108 | −17.876847 |
| 4096 | $x$ | −8.936235 | +0.000037 | −17.872508 |
| 256 | $xz$ | −20.344748 | −0.000978 | −40.688518 |
| 1024 | $xz$ | −20.120957 | −0.000105 | −40.241810 |
| 4096 | $xz$ | −20.107115 | +0.000047 | −40.214278 |

Exchanging $x$ and $y$ exchanges the two metrics. The $y$ and $yz$ area slopes reverse the displayed signs. Rest, $z$ and $xy$ are invariant under that exchange, so their true area transfer vanishes. All seven geometries enter the same board. Biasing the reference tag odds from $1:1$ to $2:1$ adds $-\log2$ nats to the unnormalized residual for this unit transfer. The distinct responses to the label exchange and the bias test both the geometric attachment and the equilibrium assumption. [All signed and control slopes](results.json)

For a general finite change from the fair reference law $P_R$ to $Q_R$, the discrepancy contains the local relative entropy, as follows from the reference log probabilities.
\[
\Delta S+\eta\Delta A=\sum_g\Delta p_g\,[\log M+\eta A_g-\log d_{eg}]-D(Q_R\Vert P_R).
\]
For the tag-only transfer from $1/2$ to $5/8$, every region reads the tag, so $D=(5/8)\log(5/4)+(3/8)\log(3/4)\simeq0.031584$ nats at every scale. The same cost appears for $3/8$. It cancels from the centered secant in the table and persists in each one-sided change from the reference. This example rules out a finite-amplitude, vanishing fourth-order deficit on the surveyed family.

## A geometry-blind straight-line count has the wrong sign

Consider a fixed homogeneous, isotropic population of complete straight lines in the static material coordinates $x$. Count each line hitting the interior of the projected convex ball $r_g L(B_1)$. Equivalently, count its two boundary intersections away from tangencies. For each direction the penetrating lines are indexed by offsets in the orthogonal projection of that body. Dilating by $r_g$ multiplies every projected area by $r_g^2$. Even granting an independent flat-space calibration for each moving frame, the resulting area estimate is therefore $4\pi r_g^2$. Exact knowledge of the sector-specific proper-volume radius has been allowed in this construction.

Write $R_g$ for the scalar curvature of the induced rest slice at the center. The normal-coordinate volume density is $1-\operatorname{Ric}_{ij}y^iy^j/6+O(|y|^3)$. Integration on the ball gives $V_g(r)=(4\pi r^3/3)[1-R_gr^2/30+O(r^4)]$ for the smooth chamber (odd terms integrate to zero). Since $dV_g/dr$ is the area, at the same prescribed volume the two quantities have expansions
\[
r_g=\ell[1+R_g\ell^2/90+O(\ell^4)],\qquad
A_g=4\pi\ell^2[1-R_g\ell^2/30+O(\ell^4)] .
\]
Consequently the line-count response divided by the true signed area response tends to $-2/3$ whenever $R_+\ne R_-$. The direct $xz$ survey at $1/4096$ gives $+13.404798$ for the line estimate and $-20.107115$ for the proper area after division by $\ell^4$. At all three scales the $+$ projected $xz$ ball strictly contains the $-$ ball. Thus even an inhomogeneous fixed **nonnegative** measure on complete straight lines cannot decrease its hit count where the true area decreases. In each sector a proper-geodesic line measure adapted to the metric can change the distribution of paths and their offsets. That requires a physical rule giving the material probe access to the relevant geometry. This example addresses coordinate-line counting on the specified metric pair. [Measured line and area slopes](results.json) · [Tests on all four responsive frame shapes](test_balls.py)

The calculation integrates the excess over flat volume and area separately, then combines each answer with the analytic flat term at 55-digit decimal precision. This preserves the fourth-order area difference and the integer floors at the finest radius. Angular and radial orders 9 and 14 give identical multiplicities for all 42 metric-dependent areas. The maximum change in $M A$ is below $8.8\times10^{-9}$, while the smallest numerical distance of $M\exp(A)$ from an integer is above 0.0202. The largest relative proper-volume root residual in the correction integrator is below $4.4\times10^{-19}$. A separate full-matrix quadrature verifies both the target volume and the area at every archived proper radius. These checks measure numerical agreement. The analytic $2/M$ bound uses the ideal exact survey areas and needs no quadrature assumption. A device supplying approximations $\widehat A_g$ can add up to $\eta\max_g|\widehat A_g-A_g|\,\|\delta p_g\|_1$ to that bound. The observed order comparison supplies no rigorous enclosure for this device error. [Stable integrals and compatible-word indexing](measure.py) · [Independent metric, entropy and floor tests](test_balls.py)

The station ranges admit both geometries. The per-cut area-to-alphabet assignment, the shared geometry label and the large calibrated registers are explicit costs of this construction. A native Interval Echo source would have to reconstruct the balls, realize those changing completion counts and supply independently calibrated regional field and collar stress on the same variations. The gravitational coupling and the OpenWave gravity and dark-matter criteria remain open on this evidence.

Reproduce the surveys and independent tests from the repository root.

```sh
python -m model.candidate_laws.interval_echo_partition_repair.proper_balls.measure
python -m pytest -q model/candidate_laws/interval_echo_partition_repair/proper_balls
```
