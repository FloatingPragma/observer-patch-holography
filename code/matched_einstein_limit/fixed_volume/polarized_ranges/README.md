# Reading curvature from distances between paired stations

Two passing radial checks have a precise consequence and a limit. On the timed, two-resident Interval Echo source, ten sets of three orthogonal axes pass an independent material midpoint test at two different radii. Every set fails a further distance check. Both ends of an axis one step from the center are three steps from both ends of an axis two steps from it. These four lengths would force a smooth minimizing geodesic through the center to branch. An experiment on smooth curved products shows how the same comparison recovers the three sectional curvatures of a chosen orthogonal triad when the endpoints really are paired geodesic displacements. [Native and smooth readings](results.json)

## Four cross ranges

Work inside a convex normal neighborhood of a smooth three-dimensional spatial slice of bounded geometry. At event $o$ prepare both ends of each independently transported axis $v_i$, at positions $\exp_o(\pm r_i v_i)$. The $v_i$ are unit vectors. Unequal radii are allowed, with $c h\le r_i\le h$ for a fixed $c>0$ under refinement. Let $d_{\sigma\tau}$ be the proper distance between end $\sigma$ on axis $i$ and end $\tau$ on axis $j$. Ultrastatic comoving observers can read this spatial distance as half the light-clock round trip. Direction-dependent propagation or nontrivial lapse requires a separate physical calibration in a different spacetime.

Write $\mathcal R_{ij}=R(v_i,v_j,v_j,v_i)$ with positive value on a positively curved section. Normal coordinates at $o$ give
\[
 d(\exp_o u,\exp_o v)^2
 = |u-v|^2-\tfrac13 R(u,v,v,u)+P_5(u,v)+O(h^6).
\]
Here $P_5$ is homogeneous of total degree five. The coefficient follows by substituting the quadratic normal-coordinate metric correction into the energy of the straight segment between $u,v$. Variation of this segment contributes only from sixth order, by stationarity of its Euclidean energy and the cubic path correction. Equivalently use the Taylor expansion of squared geodesic distance. Simultaneously reversing $u,v$ cancels $P_5$. Averaging all four sign choices cancels the tangent inner product as well. Therefore
\[
 X_{ij}=r_i^2+r_j^2-\tfrac14\sum_{\sigma,\tau}d_{\sigma\tau}^2,
 \qquad \frac{3 X_{ij}}{r_i^2r_j^2}=\mathcal R_{ij}+O(h^2).
\]
This reading needs no independently measured angle between the two axes. For oblique axes it reads sectional curvature times the squared sine of that angle. The trace requires a calibrated angular spanning set. If $\sum_i w_i v_i v_i^\top=I/3$, bilinearity gives $\mathrm{Scal}=9\sum_{ij} w_iw_j\,\mathcal R_{ij}$. For an orthonormal triad this is twice the sum of its three distinct pair readings. The identity follows by contracting the two pairs of curvature indices with the two copies of $I/3$. Subsequently $A|_{V=4\pi h^3/3}=4\pi h^2-(2\pi/15)\mathrm{Scal}\,h^4+O(h^6)$. On a slice of zero extrinsic curvature at $o$, $\mathrm{Scal}=2G(u,u)$ in its rest frame. This route can supply the geometric side of a held stress or entropy test at the independently measured ordinary coupling. [Distance expansion, Pennec (2019)](https://arxiv.org/abs/1906.07418) · [Fixed-volume geometry](../README.md)

If each length has error at most $\epsilon$, the error in $X$ is at most $8h\epsilon+3\epsilon^2$, since radial lengths are at most $h$ and cross lengths at most $2h$. With the measured radii also in the denominator, the additional relative denominator error is $O(\epsilon/h)$ for $\epsilon\ll ch$. Thus $\epsilon=o(h^3)$ suffices for curvature error tending to zero, and for fixed-volume area error $o(h^4)$. With two independent source preparations at signed amplitudes $\pm a$ and bounded third source derivatives, an area derivative on that tangent needs control after division by $a$. For instance $h^2+\epsilon/h^3=o(a)$ with $a\to0$ suffices.

The four signs matter. Take Euclidean axes $e_1,e_2,e_3$ and radii $(\ell,2\ell,2\ell)$. Tilt the first two endpoint pairs by replacing their directions $\pm e_i$ by $\pm e_i\cos\ell+e_3\sin\ell$. All center ranges remain exactly the same. The centroids of the two pairs now lie on the same side of the center. Their cross mean gives a spurious scalar $6(\sin\ell/\ell)^2\to6$ in flat space. Measuring each pair's diameter independently detects a shortfall $2r_i(1-\cos\ell)$. Accurate center ranges alone leave this effect unconstrained.

## Independent curved clock comparison

The smooth simulator places comoving mirrors on a product of a two-surface of sectional curvature $K=-1,0,+1$ and a straight line, with unit lapse. An orthonormal triad is rotated so that its three sections have different curvature. The radial lengths are $\ell,2\ell,2\ell$ for $\ell=1/4,\ldots,1/32$. Cross ranges use the constant-curvature law of cosines for the surface factor and exact addition of squared distance in the line factor. Each mirror clock reports the nearest integer round-trip tick of size $\ell^5/4096$. The decoder receives only the radial and cross tick counts. The simulator's curvature and rotation are withheld. The attributions of sectional curvature follow separately from the product's curvature tensor.

For the positive product at $\ell=1/32$ the three readings are $0.012347,0.790368,0.197593$, against exact limits $1/81,64/81,16/81$. The scalar reads $2.000615$; the negative and flat controls give $-1.999385$ and $-0.00000007$. The largest scalar error across the three metrics falls from $0.03986$ to $0.000616$ over four refinements. The finite-radius bias is separate from the interval for rounding, whose total width is below $0.000005$ on this panel. The displaced flat control instead gives $+5.998047$. [Clock generator and interval decoder](smooth.py)

## Native positions and independently attained transfer lengths

The monitored one-cell source has 66 two-resident sharp words. The enrolled two-hop ball admits 55 of them. Keep all 66 as intermediate states. An earlier one-tag acquisition independently resolves the 30 physical port adjacencies even on coherent resident stocks. Its recovered angular Gram matrix in the canonical unit-sphere representation has diagonal $1$, antipodal entries $-1$ and remaining entries $\pm1/\sqrt5$. A sharp word $(a,b)$ has measured load position $2(u_a+u_b)$ in that representation. Five distinct preparations map to the origin. Each of two shells supplies ten opposite pairs and five unpaired sites. These are angularly reconstructed positions. [Sharp-word placements and anonymous tag counts](../../../../model/experiments/interval_echo_einstein_station/equilibrium_refinement/tracer_flux/bearings/README.md)

Measure a second length between the position preparations. Drive one allowed resident-vacancy exchange per appointment. A displacement from word $x$ to word $y$ needs at least the shortest-path distance in the two-resident word graph because a single appointment changes at most one occupied site across one adjacent contact. Scheduled swaps attain this distance. It is also the first nonzero amplitude order for unmonitored propagation under the native resident Hamiltonian, whose permitted off-diagonals are all $-1/2$. The latter identity concerns the order at short time, independently of the time of first detection. In fresh transfer trials the receiving observer verifies both residents by reversible local comparisons. The return undoes the route and restores the initial word. The same isolated material exchange hand calibrates the appointment unit. All 990 different endpoint pairs among the 45 nominated words were tested in this manner, using 1,980 terminal comparisons. These are clocked preparations and transfers; they supply no prior proper-length interpretation. [Native implementation](native.py) · [Word propagation theorem](../../finite_entropy_region/monitored_order/README.md)

For **each** of the five origin words, four inner and seven outer pairs have equal distances to the origin and their mutual distance equals twice that radius. Test all triples of these eleven pairs for orthogonality using the exact integer coefficients of $1$ and $1/\sqrt5$ in their angular Gram products. Exactly two triples pass per origin. Every triple has independently measured radii $(1,2,2)$ in exchange steps. For example at origin $(1,2)$ one such triad is
\[
 (1,10)\leftrightarrow(2,9),\qquad (1,9)\leftrightarrow(2,10),
 \qquad (5,11)\leftrightarrow(6,8).
\]

For this triad the four cross ranges for the first two axes are $(3,3,3,3)$. The other two combinations give $(2,2,2,2)$ and $(3,3,3,3)$. Direct substitution into the smooth formula returns the three contractions $-3,\ 3/4,\ -3/16$, hence the same spurious scalar $-39/8$ as in the nine other surviving triads. In a smooth metric with attained minimizing distances, choose one minimizing segment from the first endpoint to the origin. Concatenation with a minimizer to either endpoint of axis two attains length $1+2=3$, hence gives a minimizing geodesic. Its continuation at the origin with that unit tangent is unique, so the two distinct length-two endpoints would coincide. Their measured separation is four. This exact obstruction is independent of a sectional-curvature approximation or any assumed angle in a physical metric.

## How much must a proper ruler change?

Suppose the observed ratios persist on a shrinking sequence in bounded curvature. In normal coordinates centered on the reported origin, lengths normalized by the shorter native radius have Euclidean distances up to $o(1)$. Distinguish two ways of calibrating the endpoints.

**Proper antipodes.** If the corrected pairs are genuine opposite geodesic displacements from the same center, the RMS of their four cross distances is $\sqrt{r^2+s^2}+o(1)$, regardless of their angle. For a fractional correction of at most $e$ on each individual measurement this requires $3(1-e)\le\sqrt5(1+e)+o(1)$. Consequently $\liminf e\ge(3-\sqrt5)/(3+\sqrt5)=0.145898\ldots$. Orthogonal Euclidean axes with radii 1 and 2 multiplied by $6/(3+\sqrt5)$ attain the bound. The corrected radii and diameters are scaled up by $1+e$, and the four cross ranges are reduced by $1-e$.

**Freely placed endpoints.** Let $a,b,c,d$ denote all four candidate endpoints in Euclidean coordinates centered at the reported origin. The nominal center ranges are $1,1,2,2$, the diameters are $2,4$, and the four $a,b$ to $c,d$ cross ranges are $3$. The identity
\[
12(\|a\|^2+\|b\|^2)+6(\|c\|^2+\|d\|^2)
-4\|a-b\|^2-\|c-d\|^2
-2\!\sum_{x\in\{a,b\},y\in\{c,d\}}\!\|x-y\|^2
=\|2a+2b+c+d\|^2
\]
places a constraint on all ten readings even if the two corrected diameters miss the center. The right side is nonnegative. Substituting the individual error bounds gives $9(1+e)^2\ge13(1-e)^2$ for $0\le e\le1$. The sharp minimum over arbitrary Euclidean placements is thus $e_*=(\sqrt{13}-3)/(\sqrt{13}+3)=0.091673\ldots$. To attain it, take $(\pm1,0,2/3)$ and $(0,\pm2,-4/3)$ and multiply every coordinate by $6/(3+\sqrt{13})$. The corrected center ranges have ratio $1+e_*$ to the measured lengths. Every other distance has ratio $1-e_*$. Multiplying the unscaled example by $0.9$ instead gives a rational-coordinate control within ten percent of every registered distance. These bounds apply asymptotically to shrinking configurations because the Euclidean normal-coordinate error tends to zero after normalization. The one-board experiment does not determine which correction a physical refinement implements. [Kernel-checked identities and range bounds](PolarizedRanges.lean)

From the repository root run with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -m matched_einstein_limit.fixed_volume.polarized_ranges.measure
python -m pytest -q code/matched_einstein_limit/fixed_volume/polarized_ranges
~~~
