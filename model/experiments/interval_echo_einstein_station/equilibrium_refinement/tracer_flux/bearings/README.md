# Locating displaced markers on the two-resident source

An independently launched material tag identifies the same twelve bearings in the presence of either held boundary correlation. Its short-flight arrival records recover the angular information needed to average out an unknown marker-density Hessian. Two residents can also supply actual displaced marker records on one finite shell. Applying the curvature estimator directly to the available discrete word counts gives a measurable bias even in the flat record-coordinate description. [Flight counts](results.json) · [Marker placements and counts](placements.json)

## A reading independent of the resident preparation

Use the [five-rim tracer apparatus](../README.md), including its three-symbol terminal test, calibrated dimer hand, two mobile load-two residents, independent load-one tag and parked exterior contacts. On each fresh copy attempt to inject the tag at a named station. Record a blocked launch or its terminal station after the supplied shutter dwell. The two transfers exchange correlations between residents. Every station admits launches on both held preparations and the uniform 55-word ball. The reader receives anonymous transmitter and receiver labels throughout.

Write the native 660-word Hamiltonian as $H=D+T$, where $D$ keeps the tag's station and $T$ changes it. Between different stations $a$ and $b$, $T_{ba}$ is either zero or $-1/2$ times a permutation isometry on resident records. Exchanging the tag with a vacancy or a resident has the same amplitude. In the interaction picture $U_I(t)=e^{itD}e^{-itH}$, integration of its evolution equation gives

$$\|U_I(t)-I+itT\|\le\tfrac12t^2 L,\qquad L=\|[D,T]\|+\|T\|^2.$$

Indeed $\|U_I(s)-I\|\le s\|T\|$ and the rotated hopping differs from $T$ by at most $s\|[D,T]\|$. Rotating by $D$ at the end preserves every tag terminal probability. Integer row and column sums of the two doubled operators bound $L\le441/4$. For $r=t^2(441/8)$ the conditional probability on a non-neighbor is at most $r^2$, and on a neighbor it is at least $(t/2-r)^2$ when $t/2>r$. These inequalities extend to arbitrary mixed or coherent resident records on the admitted launch sector by linearity of the measurement probability and the isometry bound. A launch must have nonzero admission to furnish its conditional reading. At both registered dwells, $1/512$ and $1/256$, the threshold $t^2/16$ strictly separates the two classes. Switching off contrast preserves this angular calibration. Switching off tag exchange makes every offsite hit impossible.

On each of the three stocks, independent simulated cohorts at both dwells recovered all 30 links and the same six antipodal pairs. The counts retain blocked trials. An independently read isolated exchange at each shutter measures its transfer probability $q=\sin^2(t/2)$. One quarter of this measured probability sets the detector threshold, and its binomial interval sets the timing uncertainty on the registered short branch. Simultaneous 99% intervals cover both standard readings and all receiver entries using Bonferroni with the two shutters and three stocks. At the shorter dwell the uniform-stock upper bound over non-neighbors was $1.06\times10^{-9}$; the lower bound over neighbors was $8.95\times10^{-7}$. This first-hop angular test uses a supplied shutter. It makes no inference of continuum geodesic distance from the duration of a scheduled relay. [Hamiltonian reduction, bounds and receiver reconstruction](measure.py)

## Which bearings the measurements fix

Let $A_j$ be the indicator of $j$ steps in the recovered graph. Put $F=I-A_3$, $S=A_1-A_2$. The reader verifies by integer multiplication that $G=F+S/\sqrt5$ has unit diagonal, $G^2=4G$, trace 12, and $A_1G=\sqrt5G$. This is the rank-three normalized slow Gram of the angular calibration. Only its orthogonal frame is free. For one unit vector from each antipodal pair the six symmetric dyads $u_i u_i^T$ have Gram diagonal 1 and off-diagonal $1/5$. Its eigenvalues are 2 once and $4/5$ five times. The dyads span all symmetric three-dimensional forms. The unique isotropic weights are $1/6$ each, as follows by pairing $I/3$ with the six dyads. Thus the native directions suffice for the angular part of the [displaced-marker estimator](../../../../../../code/matched_einstein_limit/fixed_volume/displaced_markers/README.md) if the corresponding proper displacements and apertures are realized. [Original one-load angular calibration](../../../../../evaluation/patch_spatial_readout.json)

## Testing the displacements themselves

Use the declared local position reading of cumulative signed load in this Gram. A balanced pair at opposite stations has position zero. Relative to that reference a sharp pair at $i,j$ gives $2(u_i+u_j)$. These are records measurable by the [routed two-resident survey](../../CLOCK_BALL.md). Of the 55 admitted words, five give the common origin. Each remaining shell has 25 distinct endpoints, radius squared $8(1\pm1/\sqrt5)$ and ten available antipodal pairs. The missing receiver is the original exterior station. These radii are in local Gram units.

For the ten pairs on each shell, five bearings are equatorial to the exterior axis and five have squared axial cosine $x_\pm=(1\mp1/\sqrt5)/2$. Their isotropic average needs total tilted weight $1/(3x_\pm)$. On the inner shell this is $(5-\sqrt5)/6$, strictly between zero and one. Giving equal weight within each of the two orbits attains $I/3$ exactly. On the outer shell the largest squared axial cosine among the ten paired axes is below $1/3$, excluding every nonnegative isotropic average of those pairs. Five additional endpoints lack a partner; the paired-center proof of cancellation does not apply to them. The integer certificate checks the full dyad identity in $\mathbb Q(\sqrt5)$ and the separation of the marker positions. The inner marker placements can therefore average **every quadratic directional response** isotropically at this radius, with the same source and labeling. [Word-pair table and certificate](displacements.py)

Use a fresh uniform monitored source as a marker, enrolling its word independently of each nominated center. The center admits 30 words at the inner radius $\ell$, five at $\ell/2$, and each of the twenty displaced half-radius apertures admits one. Indeed distinct Gram positions on this source are separated by more than $\ell/2$. The reading can use conditional sampling of the center word and a native comparison of the decoded displacement; no continuum sprinkling enters the counts. In units of the common full-radius reference denominator, the normalized central counts are $C=30,H=40$ and the displaced average is $B=8$. The smooth marker formula applied to these exact counts gives

$$\ell^2\widehat R=3(18B+22H-40C)/(4H-C)=-264/65.$$

The native simulated counts on 22 independent banks give $-4.06117$, with simultaneous 99% interval $[-4.07724,-4.04512]$. The flat Gram completion has scalar curvature zero. This identifies the failed assumption when ordinary resident words furnish the markers: their intensity has atoms at resolved positions. The finite source supplies one usable angular stencil; a physical curvature measurement still needs a common proper radius, suitable marker density and control along a shrinking family. More copies at this single radius improve precision around the displayed biased value.

A mixture of the origin and distant word records can move the *mean* to a smaller radius. Its event spread can be bounded directly. If the first nonzero radius is $r_-$ and a single word has displacement $X$ with $|E X|=s<r_-$, then $E|X|^2\ge r_- E|X|\ge r_- s$. Thus $E|X-E X|^2\ge s(r_--s)$. A binary mixture of the origin and one inner endpoint attains equality along that ray. The mean of $n$ independent copies of that binary protocol has variance $s(r_--s)/n$. Relative RMS localization $o(s)$ consequently needs $n\gg r_-/s$, and absolute RMS localization $o(s^2)$ needs $n\gg r_-/s^3$ in the fixed Gram units. A collective position readout is a concrete alternative. Using it for gravity requires an operational connection between that aggregate and a common physical ball in the same matter experiment.

From the repository root run the two modules and tests with a single numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_einstein_station.equilibrium_refinement.tracer_flux.bearings
python -m "$pkg.measure"
python -m "$pkg.displacements"
python -m pytest -q --pyargs "$pkg"
~~~
