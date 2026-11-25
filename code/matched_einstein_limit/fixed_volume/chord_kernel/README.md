# The small-ball weight from two exit clocks

From a point inside a ball, send two unit-speed probes in opposite directions along the same geodesic and clock both first exits. The product of the two times divided by the diameter gives the radial weight of the Einstein small-ball energy integral exactly in a flat ball, in every direction. Bounded curvature changes its integrated value by $O(\ell^6)$ for bounded energy density in three spatial dimensions. The calculation supplies a way to observe the **shape of the weight** once the observer has a proper ball, geodesics and clocks in common units. [Simulation and clock programs](measure.py) · [Integer readings and comparisons](results.json)

## What the clocks measure

Let $\ell$ be the proper radius of a strongly convex normal ball, $x$ an interior point, $r=d(o,x)$, and $a,b$ the forward and backward first-exit distances. Take both rays from $x$ with equal and opposite unit tangents. In Euclidean space their signed intersection parameters are the two roots of $t^2+2\langle x,v\rangle t+r^2-\ell^2=0$. Thus

$$w(x,v)=\frac{ab}{2\ell}=\frac{\ell^2-r^2}{2\ell}.$$

Every direction gives the same answer. A measurement with a prescribed direction needs no angular population model. A product of two indicator records also measures the weight. Draw independent uniform appointments $U,V$ on $[0,2\ell]$. The probability that both appointments precede the respective exits is $ab/(4\ell^2)$. Multiply by $2\ell$ of exposure to obtain $w$. Reusing one appointment would give a minimum of the two distances. On a seven-slot diameter, the executable counts all 49 ordered pairs. The conversion of a clock time to a spatial distance assumes independently measured unit speed on the rest slice.

There is a uniform estimate at grazing incidence. Parametrize the geodesic by $t\in[-b,a]$ and put $F(t)=d(o,\gamma(t))^2/2$. For sectional curvature $|K_{\rm sec}|\le C$ and $C\ell^2$ sufficiently small, $1-cC\ell^2\le F''\le1+cC\ell^2$. To see the bound, on a radial geodesic of length $r$ the transverse Hessian of $d^2/2$ on an endpoint vector $v$ is $r$ times the minimum index form among fields going from zero to $v$. For any competing field $W$, parallel transport and Cauchy–Schwarz give $\int |W|^2\le(r^2/2)\int |W'|^2$ and $\int |W'|^2\ge |v|^2/r$. The lower index bound is $(1-Cr^2/2)|v|^2/r$. The linearly increasing parallel trial field gives the upper bound $(1+Cr^2/3)|v|^2/r$. The radial Hessian is one with zero mixed term.

Interpolate the two equal endpoint values $F(-b)=F(a)=\ell^2/2$ to $t=0$. Subtracting quadratic functions with the lower and upper second-derivative bounds gives

$$(1-cC\ell^2)ab\le\ell^2-r^2\le(1+cC\ell^2)ab.$$

The measured weight has relative error $O(C\ell^2)$ uniformly. Its flat leading bound is $O(\ell)$ and the volume is $O(\ell^3)$. Thus for bounded stress the integrated difference from the proper radial kernel is $O(C\ell^6\|T\|_\infty)$. In flat three-space $\int_Bw\,dV=4\pi\ell^4/15$. The Einstein modular coefficient additionally multiplies this by $2\pi$ to give $8\pi^2\ell^4/15$. That normalization and the identification with the entropy derivative of the same state require independent evidence.

Clocking the first boundary crossing to within $h$ per leg changes the weight by at most $2h+h^2/(2\ell)$ because each true leg is at most $2\ell$. The integrated clock error is $O(h\ell^3)$ for $h\le\ell$. Thus $h=o(\ell)$ suffices at fourth order. This leg error bound applies to sampling a true crossing in time. A coordinate displacement of a general boundary can produce a larger individual time shift near tangency. A volume quadrature also requires its own cell measure and stress regularity. For the cubic Riemann sums the extended flat weight is Lipschitz and zero on the boundary. The cell-size error for bounded Lipschitz stress on the rescaled ball is therefore $O(h\ell^3)$.

## Actual first exits on a scheduled cubic register

A logical $9^3$ register fits on the data patch of the 80-cell federation. The chart assigns 729 unmarked native ports to cubic stations. The router transposes two addressed ports using positive native seam contacts and restores the intermediate data ports. A separate blank dimer hand advances on every microperiod, including padding, and an echoed sharp pointer verifies every arrival of the hand. Each enrolled logical step uses 65 microperiods. A sharp reading at the destination verifies the tag at each station. Inspection samples the end of each programmed step. Paths through intermediate circuit addresses do not define logical boundary hits. A separate quadrature check reads one isolated transfer time as $\pi$. The input tag, rail chart, slot padding, rim membership, interrogation switches and clock conversion are supplied. [Native contact, route and pointer construction](../../../../model/experiments/interval_echo_einstein_station/proper_ball_ruler/README.md)

The first-outside aperture has radius 3.5 mesh steps. From the centre, both legs along each coordinate axis took four steps. From $(1,-1,2)$, the opposite leg counts were $(2,4),(4,2),(2,6)$ along the three axes. The differences between the three products at this resolution are retained. Turning the launch by a right angle after its first step changes the central legs to $(5,5)$. Changing only aperture membership to $x^2+y^2+4z^2<3.5^2$ gives $(2,2)$ along the third axis with the same route and clock services. Removing exchange from the data generator during the same positive contact schedule preserves every blank load word. The working hand registers one complete logical step and the data arrival check fails. The circuit demonstrates a guided ruler whose interpretation depends on the station chart and launch rule.

For finer cubic comparison grids the calculation sums exact integer first-exit counts over **every** interior launch point along one axis. The transverse lattice coordinates $j,k$ put its first exterior stations at $\pm\lceil\sqrt{n^2-j^2-k^2}\rceil$ for $n$ points per radius. The mass errors relative to $4\pi\ell^4/15$ are 21.63%, 11.30%, 5.73% and 2.94% at $n=8,16,32,64$. Holding the decoder fixed, the two test profiles $x_{\parallel}^2/\ell^2$ and $x_{\perp}^2/\ell^2$ give 0.12516 and 0.12382 at $n=64$ in units of $\ell^4$. Both continuum values are $4\pi/105\simeq0.11968$. These are geometric counts on finer grids. The addressed native pulse replay uses the $9^3$ chart.

## Is the signed curvature visible in the clocks?

The second comparison uses the five supplied ultrastatic space forms with constant spatial sectional curvature $K=-1,-1/2,0,1/2,1$. The clock emitter advances the geodesic of given initial bearing and finds the first clock tick at which its radial position has reached the rim. It uses the space-form cosine law for propagation. The perpendicular decoder receives the radial and exit tapes and the supplied right angle. With $C_K(r)=\cos(\sqrt K r)$ in the spherical case and $\cosh(\sqrt{-K}r)$ in the hyperbolic case, the perpendicular half-chord $d$ obeys $C_K(\ell)=C_K(r)C_K(d)$. Expanding through fourth order gives

$$d^2=\ell^2-r^2+\frac{K}{3}r^2(\ell^2-r^2)+O(K^2\ell^6),\qquad
\widehat K=\frac{3[d^2-(\ell^2-r^2)]}{r^2(\ell^2-r^2)}.$$

This statistic is more demanding than the integral weight. With $r/\ell$ bounded away from zero and one, clock quantization contributes $O(h/\ell^3)$ to $\widehat K$. The simulation uses $h=2^{-18}\ell^5$. The three launch radii are $\ell/4,\ell/2,3\ell/4$. Oblique bearings with radial cosine $1/3$, $0.8$ and $1$ check the product weight independently of the perpendicular estimator. A separate count of the surface area at matched proper volume supplies $G_{uu}$; on these ultrastatic slices $G_{uu}=3K$. The area and chord meters share the supplied space form and scale, and use different data. [Area experiment and count protocol](../README.md)

| Radius | Largest relative weight error | Largest sectional error | Largest contraction difference between chord and area |
|---:|---:|---:|---:|
| 1/2 | 0.050060 | 0.067941 | 0.141103 |
| 1/4 | 0.011908 | 0.016148 | 0.034551 |
| 1/8 | 0.002941 | 0.003988 | 0.008470 |
| 1/16 | 0.000733 | 0.000994 | 0.002041 |

## A curvature test with three unequal sectional readings

A supplied product of a curved two-surface of curvature $K=\pm1$ and a straight line has different sectional curvatures in different planes. Shoot from radius $r=\ell/2$ along an initial radial unit vector $u$ and launch both perpendicular probes along the parallel-transported $\pm v$. The spatial geodesic in each factor has constant speed equal to the length of the projected velocity. The emitter uses the cosine law in the surface, adds the squared distance in the line, and stops at the first tick outside. It receives neither a plane curvature nor an expected exit time. For orthonormal $u,v$, expanding that squared distance gives
$$R(t)^2=r^2+t^2-\tfrac13 K r^2t^2\,[(u\times v)\cdot e_z]^2+O(K^2\ell^6),$$
where $e_z$ points along the flat line. Thus the perpendicular product reads the curvature of the selected plane. The sum over the three planes of an orthonormal frame is the time-symmetric rest contraction $G_{uu}=K$ on this product, where the subscript denotes the observer's timelike direction. An aligned frame has sectional curvatures $(K,0,0)$. A tilted frame has $(4,36,9)K/49$.

A separate calculation integrates the disk cross-sections to get the proper three-volume, varies the radius until the volume matches the Euclidean reference and reads the derivative of volume for its area. At radius $1/16$ the tilted positive chart yields sectional readings $(0.08164,0.73504,0.18372)$ and their sum $1.000405$. The equal-volume area gives $1.000132$. In the negative chart the two estimates are $-0.999595$ and $-0.999868$. Across both signs and both frames the largest mismatch at this radius is $0.000456$. These areas are numerical metric integrals on supplied anisotropic geometries. The five space-form area comparisons above use independently sampled marker counts. [Product geodesics and volume integration](product_chart.py)

## Independent response on the freely formed collar

The 80-cell freely formed two-memory source has its own material detectors. For the six clock surveys of a tetrahedron rescaled to unit volume, let $b'_o,b'_d$ be the area-shape derivatives for an ordinary reversal-energy addition and strengthened pairing. Let $z_o,z_d$ be the same changes in a distinct probe's bridge-completion probability. At reference mean sign $sfz_d$, a balanced binary entropy has derivative $-f z_i\operatorname{atanh}(sfz_d)$ to leading order for direction $i$. Here all measured $z_o,z_d$ are negative, and all area slopes are positive. Ordinary calibration therefore gives a positive area coefficient. Transferring that coefficient to pairing tests $b'_d/b'_o=z_d/z_o$. This is a test on a dilute family of copies of a freely formed source and the specified binary entropy and tetrahedral instrument. The two-coin von Neumann entropy and the area of a proper ball are different observables. [Area and calibration definition](../../source_implication/README.md)

Reversing the bridge launch, switching between vacancy and excess, and testing all three existing leg appointments (0.5, 1, 2) gives twelve completion records. Their $z_d/z_o$ range is $0.17913$ to $32.17270$. The six area ratios range from $0.00089435$ to $0.00547804$. The original exposure-2 excess calibration transfers at most 0.602% of the required response. Across the wider menu, any fixed nonnegative mixture of the completion probabilities gives a ratio in their interval, since the ordinary changes all have the same sign. A nondecreasing differentiable combination of the shapes has a ratio in their interval for the same reason. Thus at most **3.059%** transfers even with both freedoms, at nonzero ordinary gain. The two spatial source profiles differ. An equation relating the full stress tensor to independent proper-ball geometry requires additional measurements. [Extended bridge comparison](collar.py)

From the repository root:

```sh
PYTHONPATH=code:. python3 -m matched_einstein_limit.fixed_volume.chord_kernel.measure
python3 -m pytest -q code/matched_einstein_limit/fixed_volume/chord_kernel
```
