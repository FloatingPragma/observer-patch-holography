# Observer centers near an interior radar hinge

In the five-station filling, three Euclidean tetrahedra meet around each of four interior edges. In the interior of such an edge the metric is locally the product of a line and a planar cone of total angle `α`. The finite radar survey gives the lengths used in this filling. The calculation below concerns a specified ensemble of observer centers on its interior. [Five-station survey and native data](../hinge.json)

Write `δ = 2π − α` for the angle deficit and assume `π < α < 2π`. A ball centered exactly on the edge has the fraction `α/(2π)` of the Euclidean ball's volume at every sufficiently small radius. A center off the edge eventually sees an ordinary Euclidean ball. The ensemble can have a finite small-ball curvature coefficient because the centers within distance `r` of a unit length of edge occupy volume `α r²/2`.

## One common radius for the ensemble

Take a unit length of product hinge, put `β = α/2`, and place a center a transverse distance `s` from the edge. On unfolding its cone sector around the edge, the metric ball is a Euclidean ball centered at `(s,0,0)` intersected with `|φ| ≤ β`, up to faces of measure zero. The shortest angular difference is at most `β < π`. A straight chord in that unfolded sector attains the distance. The axial coordinate is independent. Because `β > π/2`, the missing part can meet the ball only for `s < r`.

Let `F(α) r⁵` denote the integral of the volume deficit over centers per unit length. Integrate the missing part over its two angular sectors, its transverse radius `ρ`, its axial chord, and the center volume element `α s ds`. At `r=1`, substituting `ρ=st` and integrating out `s` gives

\[
 F(\alpha)=\frac{8\alpha}{15}\int_\beta^\pi\int_0^\infty
   \frac{t\,dt\,d\phi}{(1-2t\cos\phi+t^2)^2}
 =\frac{2\alpha}{15}\bigl[\cot\beta+(\pi-\beta)\csc^2\beta\bigr].
\]

The inner integral is the derivative with respect to `φ` of `−[cot φ + (π−φ)csc²φ]/4`. The expression in brackets tends to zero at `π`. Restoring the radius by dilation proves the fifth power exactly on the product. The integrated area deficit at this radius is `5F r⁴`, by differentiation or the radial coarea formula.

For a probability density `w` of observer centers, put `C = Σₑ F(αₑ) ∫ₑ w ds`. Assume that `w` is continuous and bounded, with compact support away from the exterior boundary. For a fixed finite, nondegenerate star the mean volume is `4πr³/3 − C r⁵ + o(r⁵)` and the mean area is `4πr² − 5C r⁴ + o(r⁴)`. On each edge tube, rescaling the transverse coordinates gives these limits by dominated convergence. The regions within a constant times `r` of a vertex have center volume `O(r³)`; all ball volumes and areas there are `O(r³)` and `O(r²)` by dilation of the finitely many incident cones. They contribute `O(r⁶)` and `O(r⁵)` respectively. The fixed star has uniform tube and dilation constants.

The *one common radius* at which the mean volume equals `4πr³/3` is `r + C r³/(4π) + o(r³)`. At that radius the mean area is `4πr² − 3C r⁴ + o(r⁴)`. Thus the per-length coefficient before weighting by observer density is `3F`.

## Each observer fixes their own volume

For an individual ball of actual radius `R`, put `u = s/R`, and let `b(u)` and `p(u)` be its volume and area divided by `4πR³/3` and `4πR²`. Let `η = δ/2` be the half-angle of the missing wedge and `n` a uniformly distributed Euclidean direction. From the observer a ray enters that wedge at distance `s/m`, where `m = max(0, −nₓ − |nᵧ| cot η)`; for `m=0` it stays outside. Integrating the ray up to its first wedge crossing gives

\[
 p(u)=\Pr(m\le u),\qquad b(u)=\mathbb{E}\min(1,(u/m)^3),
 \qquad p=b-u b'/3.
\]

The integrand at `m=0` is one, also at `u=0`. Thus `1/2 < α/(2π) ≤ p ≤ b ≤ 1`. Both fractions increase with the offset and with the full cone angle. The identity for `p` follows by differentiating the volume in `R` with the center fixed.

For a target volume `4πr³/3`, the observer takes radius `R = r b( u )⁻¹ᐟ³`. Then `x=s/r=u b⁻¹ᐟ³` and `dx/du=p b⁻⁴ᐟ³>0`. This determines the adjusted ball uniquely. The individual area relative to `4πr²` is `p b⁻²ᐟ³`. Integrating the loss against `α s ds` gives the per-length coefficient

\[
 I(\alpha)=4\pi\alpha\int_0^1 u p b^{-5/3}(1-p b^{-2/3})\,du .
\]

The averaged area at the individual target volume is `4πr² − r⁴ Σₑ I(αₑ)∫ₑ w ds + o(r⁴)`. Radii obey `R≤2¹ᐟ³r`, so the same localization and vertex estimates apply. The estimators differ because adjustments near an edge can be a finite fraction of `r`. For example at `δ=2.6225`, a center located at `s=r/2` takes `R/r∈(1.04375,1.04377)` and reads `A/r²∈(11.0137,11.0138)`. Outside the tube, `s≥r` gives the Euclidean radius and area.

### Interval evaluation

Integrate the vertical coordinate of the uniform sphere direction at fixed azimuth. For `0≤θ≤η` put `a(θ)=sqrt(max(0,1−(u sinη/sinθ)²))` with its endpoint interpreted by continuity. The area and volume deficits relative to their Euclidean values are `π⁻¹ ∫a dθ` and `π⁻¹ ∫a³ dθ`, respectively. The latter follows from integrating `1−(u/m)³` over the directions with `m>u`. With `y=u sinη`, `h=sqrt(1−u²)`, `λ=atan2(sinη h,cosη)`, and `ψ=atan2(h,u cosη)`, the integrals are

\[
 \pi(1-p)=\lambda-y\psi,\qquad
 \pi(1-b)=\lambda-\tfrac12 y(3-y^2)\psi+\tfrac12 u y\cos\eta\, h.
\]

For the first identity change variables to `cosθ=sqrt(1−y²)cos v`. For the second subtract `y²∫a csc²θ dθ` from the first integral and set `t=cotθ`. The remaining integral is a circular segment. The unit tests compare both expressions to a separate angular integration.

The certificate evaluates these formulas by outward interval arithmetic. In each of 512 equal radial subintervals the left endpoint at the smallest angle bounds `p,b` below, and the right endpoint at the largest angle bounds them above. Separate bounds on all positive factors of the integral for `I`, with exact integrated weight `∫u du`, give the following enclosures. The angle bands contain the full interval enclosures for the four reference deficits of the native weak monitor at 16 inspections per tick. [Coefficient implementation](cone.py) · [Numerical and interval checks](test_cone.py) · [Native angle intervals](native.json)

| Deficit `δ` in radians | Common-volume area coefficient `3F` | Individual-volume area coefficient `I` | Linear Regge coefficient `4πδ/15` |
|---|---:|---:|---:|
| 2.575 to 2.580 | 1.6345 to 1.6478 | 1.7648 to 1.8152 | 2.157 to 2.162 |
| 2.620 to 2.625 | 1.6601 to 1.6731 | 1.7956 to 1.8469 | 2.194 to 2.200 |

### Derivative for a perturbed hinge

For a change of the measured full angle, write a subscript `α` for differentiation at fixed `u`. With `v = u cosη` and the same `h,ψ`, differentiating the angular integrals gives

\[
 2\pi p_\alpha=h-v\psi,\qquad
 4\pi b_\alpha=2h^3+3v^2h-3v(h^2+v^2)\psi .
\]

The first right-hand side is `∫₀ʰ t²/(v²+t²)dt`. The derivative of the second right-hand side with respect to `h` is `6h(h−vψ)`, and both expressions vanish at `h=0`. Therefore each grows with `h` and decreases with `v`. As either `u` or `α` increases, the two angle partials decrease. This supplies independent endpoint bounds on every integration bin. The formulas extend continuously to `u=0,1`, including the value zero at the outside of the tube.

Since `x dx=u p b⁻⁵ᐟ³ du` integrates to one half, put `G=p² b⁻⁷ᐟ³`. Then `I=2πα−4πα∫₀¹uG du` and

\[
 I'(\alpha)=2\pi-4\pi\int_0^1 u G
 \left[1+\alpha\left(2p_\alpha/p-7b_\alpha/(3b)\right)\right]du .
\]

Differentiation under the integral is justified by the positive lower bound on `p,b` and continuous bounded partials on each compact angle panel. The function `individual_slope_interval` encloses this derivative on the entire requested angle interval. For example, 4,096 bins at `α=2π−2.6216±10⁻⁸` give `−0.70154 < I′ < −0.69347`. The sign refers to opening the cone at a fixed observer-density weight per unit length of hinge. Changes of hinge length or observer weights contribute separately to a survey of the full star.

## Weak deficits and observer sampling

As `δ→0`, the bracket in `F` is `δ/3 + O(δ³)`. Hence `3F = 4πδ/15 + O(δ²)`. The individual coefficient has the same first-order term. Indeed `1−b` and `1−p` are uniformly bounded by `η/π`. The integrand of `I/(4πα)` equals `u[(1−p)−(2/3)(1−b)+O(δ²)]` uniformly. The integrated leading terms are the common-radius area deficit `5F` minus twice the volume deficit `F`. The resulting first-order coefficient is the smooth fixed-volume area coefficient `2π/15` times the integrated Regge scalar per unit length `2δ`. At the finite native angles the interval table distinguishes all three coefficients.

A continuous probability density of centers assigns a positive shrinking-tube population to the edge. A fixed finite menu of stations does not realize that sampling law as `r→0`. If all centers keep positive distance from the hinges and vertices, their balls eventually have zero deficit. A center on an isolated hinge with a fixed positive sampling weight contributes an area deficit of order `r²`. Centers near a physical exterior boundary may produce a term of order `r³` in the mean area, which is why the support assumption matters. A refining observer ensemble must establish the sampling and error bounds through the shrinking tube. The spatial scalar average requires additional directional clock, stress and regional entropy data to test physical null balance.

For independent draws of centers with the declared density and exact individual ball readings, an affected center occurs with probability `O(r²)`, and the absolute area deficit is bounded by a constant times `r²`. The variance is therefore `O(r⁶)`. For `N` centers, the root mean square sampling error after division by `r⁴` is `O(1/(r√N))`. Thus `Nr²→∞` suffices for convergence in probability. Where the ensemble has nonzero hinge weight, a sample with bounded `Nr²` retains a positive probability of missing every affected center. Finite accuracy of the area and volume readouts adds a separate source of error. A boundary that merely chooses which centers to enroll permits their balls to cross it inside the physical star.
