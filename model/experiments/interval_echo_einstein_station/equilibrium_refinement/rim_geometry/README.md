# Can a clock aperture identify the entropy-balancing area?

The central observer times exchanges to five boundary receivers on the one-cell, two-load source. All ten ways to place both residents on that rim have ten boundary completions and the same interior entropy score. A geometric test asks how much of the rim a resident covers, calibrates this area on one resident, and predicts what two residents cover together. Every positive aperture with a single angular width misses at least one placement by **0.193858 nats**, even with adjustable width and directional sensitivities. A separate test on the full clock table finds five instances in which a smooth minimizing geodesic would have to branch. [Native reports, all aperture counts and certificate](results.json) · [Parent source and entropy measurements](../CLOCK_BALL.md)

## What the observer times

Use the 20-cell host and its enrolled one-cell region. Park the source for the authenticated read and reverse each message route. A separate material exchange hand calibrates the appointment duration. The center is port 0. Its six interior addresses and five rim addresses give 55 two-resident words after heralding both residents inside radius two. The six center/rim origins each address all twelve ports and return the borrowed registers to their initial state. Restricted to the rim, their measured distances are those of the cycle $2,4,9,8,6,2$, with one or two exchanges between different receivers. Calling a pair near or far refers to these readings.

The ordinary training transfer moves a uniform interior resident onto one uniformly chosen rim site while keeping two residents in the ball. The volume in the parent's marker instrument stays fixed. Its boundary completion score changes from zero to $\log 5$. The held words with two occupied rim sites each score $\log 10$. These are natural logarithms. Their expectation on mixtures gives the boundary term that cancels the interior Shannon differential at the uniform source. Neither the entropy measurement nor a specific smooth angular area enters the definition of the clock distances.

## All fixed-width angular apertures

Interpolate the rim distances by a circle of circumference five. Its five marked sites are equally spaced. Give each detector position any nonnegative weight. Its panel lights when at least one occupied site lies in its closed arc of length $L$, fixed across the detectors. A rigid angular offset or a positive measure on offset positions is allowed. Its lit area is the positive weighted measure of this union. This includes point-sampled panels and continuous rim coverage.

Write $s$ for the average one-resident lit area. Ordinary calibration gives $\eta=\log(5)/s$ for $s>0$. In the average over the five near pairs or the five far pairs, the sum of the two separate footprints is $2s$, including for unequal site gains. If $O_j$ is the mean intersection in pair class $j$, the predicted score is $\log(5)(2-O_j/s)$. Agreement requires both overlap ratios to equal

$$t=\frac{\log(5/2)}{\log 5}=0.56932344\ldots.$$

An arc contains $k$ consecutive marked sites. At one fixed width its possible occupancies differ by at most one, including atoms on arc endpoints. The complete counts per detector are

| $k$ | Contained near pairs | Contained far pairs |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 0 | 0 |
| 2 | 1 | 0 |
| 3 | 2 | 1 |
| 4 | 3 | 3 |
| 5 | 5 | 5 |

For each positive $k$, divide the last two columns by $k$. Positive averaging puts $(O_1/s,O_2/s)$ on one of the segments between adjacent rows, weighted by each row's contribution to $s$. A zero row contributes neither overlaps nor $s$. On the segment from $k=3$ to $4$ the two endpoints are $(2/3,1/3)$ and $(3/4,3/4)$. The smallest maximum absolute deviation from $t$ occurs where the near excess equals the far shortfall. Its size is $1/2-2t/3$. At smaller occupancies the far ratio is at most $1/3$. At larger occupancies both ratios are at least $3/4$. Either deviation exceeds the segment optimum. At least one individual held word therefore misses by

$$\frac{\log 5}{2}-\frac{2\log(5/2)}{3}=0.1938584683\ldots\quad\text{nats}.$$

Uniform angular measure and $L=3/(3/2-t)=3.22346144\ldots$ attain the bound. Their near and far overlap lengths are $L-1$ and $2L-5$. Any additional error on individual single-site gains can be measured separately. The bound uses the independent ordinary *mean* gain and holds without rotational symmetry of the detectors.

Allowing two different widths changes the conclusion. Poll the whole rim with probability $p$ and a uniformly selected individual site otherwise. A single lights the detector with probability $s=p+(1-p)/5$, and every pair lights it with probability $u=p+2(1-p)/5$. The choice $p=t/(5-4t)=0.20910205\ldots$ gives $(\log 5)u/s=\log 10$ on all ten pairs. Every $p$ admits its own exact one-resident calibration. Here the selected $p$ uses the two-resident completion count. This positive construction displays the optical freedom left after ordinary calibration. Finite clock data could identify the mixture independently by measuring the fraction of full-rim exposures.

## Independent retention by a material clock

Each rim receiver requests the local comparison from every rim station. Each routed report passes one independent coin test per exchange on the measured shortest route. The test reads a separate material exchange hand after half an appointment. With a fresh single load the native two-state evolution gives $\cos^2(\pi/4)=1/2$. Thus a report at distance $d$ passes with probability $2^{-d}$. Receiver hits require an occupied returned report. An expired or empty report contributes zero. Distinct reports use independent coins. All original source flags and routes are read and undone. The clock choice uses the measured exchange duration alone. The deposited stochastic counts sample the exact conditional hit law on each of the 55 admitted words, with 262,144 independent simulated starts per word and receiver. [Acquisition and minimax certificate](measure.py)

A single resident lights $5/2$ panels on average. Two residents one exchange apart light $59/16$, compared with $4$ for two exchanges apart. Independent single-resident calibration and held panel counts give **[2.37015, 2.37736]** and **[2.57123, 2.57867]** nats respectively, at simultaneous 99% coverage. Both exclude $\log 10=2.302585\ldots$. For a class with $N$ single-panel draws, Hoeffding's two-sided bound for independent bounded draws gives a half-width $5\sqrt{\log(600)/(2N)}$ on the mean five-panel total. A union bound covers the ordinary, near and far means. Interval division by the positive ordinary reading propagates calibration error.

For any independent, monotone radial retention on these measured sites, let the probabilities at distances zero, one and two be $a\ge b\ge c\ge0$. Across the cycle a single lights $a+2b+2c$ panels. Its near and far overlap sums are $2ab+2bc+c^2$ and $2ac+b^2+2bc$. Their difference is $(b-c)(2a-b-c)$. Equal held two-resident scores force $b=c$. This also holds for positive mixtures and nonnegative pixel weights after averaging over the five cyclic placements. Each contributing profile with $b>c$ gives a strictly positive average difference. A flat noncentral response is possible, and the measured fair-coin reader has $b=1/2,c=1/4$. The proof uses independent report retention. A common random footprint follows the intersection law of the preceding section instead.

## Proper distances and proper area

The full registered clock menu has five exact branching tests. For example, port 1 is one exchange from the center. Ports 8 and 9 both lie one exchange from 1 and two from the center. The two endpoints are themselves one exchange apart. Suppose a smooth Riemannian spatial metric realized all these numbers as attained minimizing distances. Fix any minimizing path from 0 to 1. Appending a minimizer from 1 to either endpoint gives a minimizing path by the measured equality of lengths. At 1 the continuation must have the same unit tangent in either case. Uniqueness of geodesic continuation makes the two endpoints at unit continuation time identical, contradicting their measured positive separation. The other forks occur at 5, 7, 10 and 11. The assertion concerns exact realization of these scheduled exchanges as proper lengths. A finite approximation with a controlled continuum error poses a different test.

Even the six-site subset that *can* be realized smoothly fails to select a proper area. Here is an explicit three-dimensional counterexample with the same central distances, the same ten rim-pair distances, equal central ray directions and the same **true proper volume**. Use polar variables $(r,\omega)$ with angular unit sphere $S^2$ and its round metric $g_\Omega$. Put $R=5/(2\pi)$ and choose a smooth positive $F(r)$ equal to $r$ near the origin and to $R$ for $r\ge 1$. Mark the center and five equally spaced equatorial directions at $r=2$. Choose a smooth $\psi$ supported in $(7/4,9/4)$ with $\psi(2)=0,\ \psi'(2)=1$, for example $(r-2)$ times a smooth bump equal to one near 2. Choose $B\ge0$ supported in $(1,3)$ and equal to one on the support of $\psi$. Put $h(\omega)=\omega_z^2$. For $|u|\|\psi'\|_\infty<1$ define

$$g_u=dr^2+F(r)^2[1+h(\omega)(B(r)+u\psi'(r))]g_\Omega.$$

These smooth, positive metrics have radial distance $r$ from the center, by the lower bound $ds\ge|dr|$ and the radial paths. All five central rays and their mutual directions agree. At $r\ge1$ the angular line element is bounded below by $R\,ds_\Omega$. The shorter equatorial arcs at $r=2$ attain lengths one and two. A competing path that reaches $r\le1$ has radial variation at least two. Thus all ten ambient rim distances agree throughout the family. On a two-dimensional angular sheet the positive factor multiplying $g_\Omega$ multiplies the area density by the same factor. Using $\int_{S^2}h\,d\omega=4\pi/3$ gives

$$V_u(2)-V_0(2)=u R^2(4\pi/3)\psi(2)=0,\qquad
 A_u(2)-A_0(2)=u R^2(4\pi/3)\psi'(2)=\frac{25u}{3\pi}.$$

Both signs of the area change are available. This auxiliary continuum family tests identifiability; the six marked sites alone admit it. It also preserves an ordinary calibration. Let a small smooth source response $u$ depend on $\langle k\rangle$ and on $P(k=2)$, with the uniform reference subtracted. The $0\to1$ ordinary transfer fixes the coefficient of $\langle k\rangle$. The held split-to-spread transfer fixes $\langle k\rangle$ and increases $P(k=2)$ by $1/2$. Its geometric response can take either sign with the same ordinary calibration and fixed true volume. Nothing here assigns this auxiliary metric to Interval Echo. A separate proper-volume measurement at neighboring radii distinguishes these metrics through $V'_u(2)=A_u(2)$. If $|A'|\le M$ on a measured shell of proper thickness $\delta$, its mean proper area differs from $A(2)$ by at most $M\delta/2$ before measurement error.

Run from the repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_einstein_station.equilibrium_refinement.rim_geometry
python -m "$pkg.measure"
python -m pytest -q --pyargs "$pkg"
~~~
