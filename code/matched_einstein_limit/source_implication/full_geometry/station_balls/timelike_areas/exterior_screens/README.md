# What a second screen at the same source can measure

On the active 80-cell Interval Echo board, marking whether the probe reaches port 1 halfway through the previous angular exposure changes its later arrival probability at that port by **0.012706**, with all starts retained. Three independent simulated banks measure how the change depends on how much path information the marker stores. A separate calculation for the two proposed smooth geometries gives an identifiability limit even for ideal passive ray screens. It also identifies one interior, finite-aperture measurement that would distinguish the geometries. [Measurement protocol](PROTOCOL.md) · [Native data](results.json) · [Two candidate geometries](../README.md)

## Marking a live path changes its late image

Use the same selected two-head state, emitting port, moving charges, reversal sectors and terminal time as the earlier 80-cell angular panel. At time 0.75, a two-state pointer marks occupation of port 1. It preserves the full head state conditional on that path. At time 1.5 read the pointer and the twelve terminal positions, retaining all other positions in a thirteenth bin. The marker is an added query during otherwise native propagation. The prepared state and generator are unchanged between the three arms. [First angular acquisition](../../port_apertures.py) · [Native Hamiltonian and reference](../../../../full_algebra/PROOF.md)

For each charge and reversal sector, write $U=\exp(-0.75 iH)$ and let $P$ project on port 1 with all head states retained. Write $A=UPU$ and $B=U(I-P)U$, applied to the same initial port and head distribution. The two pointer-bin amplitudes are $B+kA$ and $\sqrt{1-k^2}A$ for pointer overlap $k$. The final populations summed over both pointer outcomes are affine in $k$ because the interference term between the two paths carries exactly that factor. The no-marker setting has $k=1$, the full path marker has $k=0$, and a third independent bank at $k=1/2$ tests their predicted midpoint. This identity holds even when the head preparation is mixed.

Each acquisition simulates $2^{26}$ independent starts by a multinomial draw from the native Born probabilities. The blank minus fully marked terminal difference at port 1 has a simultaneous 99% interval **[0.0121777, 0.0132351]**. Ports 0 and 1 give resolved changes among all thirteen comparisons. Across all bins the largest deviation of the half-overlap observations from the midpoint of the other two banks is 0.0000436 against a joint allowance of 0.000529. The new blank agrees with the earlier independent angular bank to within 0.000033 in every bin. Hoeffding's bounded-event inequality and a union bound cover 52 frequencies on four banks at once. The intervals also allow the certified displacement of the selected reference. The archived path/terminal data contain all 78 bins, including impossible outcomes. The check tests the deposited counts against the native joint marginal law without requiring identical random draws on different machines. [Counts and simultaneous intervals](results.json) · [Acquisition](measure.py) · [Checks](test_screens.py)

A transparent monitor of a pre-existing trajectory preserves the late marginal when the early marker is ignored. The measured native pointer disturbs this marginal. Sequential port detections therefore need an instrument model if they are used to infer passive trajectories or count a common cross section. Other instruments, sufficiently weak monitoring, wave optics in curved geometry and alternative observer reconstructions can be investigated independently. The terminal ports here have addresses. They have no independently measured area element or affine distance from the launch.

## Exact ambiguity with all exterior centered ray observations

Use spatial normal coordinates $x$ at the marked center of the station chamber, the two ultrastatic metrics of the parent computation, and the same orthonormal frame at the center. Their deformations have the form $\chi(|x|)C_x^T B_\pm C_x$, where $C_x v=x\times v$, $B_+=\operatorname{diag}(3,1,2)$ and $B_-=\operatorname{diag}(1,3,2)$. The cutoff vanishes at and past 0.03. The metrics agree on an open exterior of its support.

In both metrics $h_g(x)x=x$. In polar coordinates there is an exact $dr^2$ term and no radial-angular cross terms. Every line $x=sv$ from the center is therefore an affinely parametrized spatial geodesic. The complete spacetime launch map in this chamber is $(v^0,v)\mapsto(t_0+s v^0,sv)$. For a fixed centered source event this map agrees for **all** initial directions, angular profiles, departure energies and speeds, including null, timelike and boosted launches. At every screen outside the support, all three of the following consequently agree across both metrics.

- The intersection events, arrival times and locally measured arrival directions agree. Both endpoint metrics agree at the outside screen and at the center.
- Differentiating the identical launch map with respect to the initial direction gives identical connecting fields of each radial beam. At every exterior intersection, their pairwise inner products and their derivatives in the affine parameter agree, since the metrics agree throughout an exterior neighborhood. This gives the full angular Jacobi map in common endpoint coordinates, including every infinitesimal pixel area, shear and every finite solid-angle integral. Any common exterior screen and common calibration with transverse intersections may be used.
- Subsequent evolution through any common exterior segment agrees by the geodesic initial-value equation. Returning to the source along the same reversed radial ray preserves the travel time. An adaptive sequence of registered launches or passive, nondisturbing readings with shared observation noise keeps the same transcript law, by induction on its conditional transition probabilities.

The two candidate manifolds therefore have identical data for this entire centered-ray protocol, as well as the same six native-calibrated station lengths. The construction supplies a class of geometries, ray launches and passive screens. It gives no equality for an experiment scattering from an interior pixel, moving the emission *event* off-center, measuring endpoints in the perturbed chamber, or solving a covariant wave equation in each candidate spacetime. Equating the quantum native path/terminal table above to a passive screen requires a separate prediction of its detector dynamics in the chosen metric.

There are only two alternatives in this comparison. For squared-range unit area $A_g(\ell)$ of the common-volume ball in the tagged axial frame $(5/3,4/3,0,0)$, any function of these exterior ray observations has worst absolute error at least $|A_+-A_-|/2$. The midpoint attains this error for the stated pair. The same lower bound holds for worst *expected* absolute error with randomized estimators, by the triangle inequality and the identical input laws. At $\ell=1/1024$ the bound is $4.4692388\,\ell^4$ from direct proper-area integration. The exact limiting normalized bound is $64\pi/45$. Arbitrarily fine angular resolution on the exterior screens leaves a fourth-order ambiguity. [Proper-volume integrals and analytical curvature limit](../README.md)

## One interior finite cap separates the two possibilities

An ideal measurement at a known proper radial distance $0<r\le0.02$ in the undeformed central frame gives a positive way to discriminate the two candidates. Open the circular beam cap of half-angle $30^\circ$ about the tagged $x$ axis and measure its physical cross section at that affine distance. Set $t=\operatorname{sinc}^2(r/a)$ with round chamber radius $a$. For a unit launch direction $n$, the angular area density per solid angle is $r^2\sqrt{D_\pm}$ with
\[
D_\pm=t^2+t r^2[\operatorname{tr}B_\pm-n^TB_\pm n]+ r^4 n^T\operatorname{cof}(B_\pm)n .
\]
Indeed $C_n$ is a quarter turn on the tangent plane, and the determinant of the tangential restriction of $B$ is its cofactor along the normal. Thus $D_+-D_-=-2r^2(t+2r^2)(n_x^2-n_y^2)$. The last angular factor is positive throughout the chosen cap. Each $\sqrt{D_\pm}$ lies between $t+r^2$ and $t+3r^2$, since the eigenvalues of $B_\pm$ lie in $[1,3]$. The cap integral of $n_x^2-n_y^2$ is $\pi\cos(30^\circ)\sin^2(30^\circ)=\pi\sqrt3/8$. Finally $t>0.97$ throughout this interval by $\sin u/u\ge1-u^2/6$. The exact signed cap-area difference therefore satisfies
\[
\frac{\pi\sqrt3}{8}(1-r^2/0.97)
\ \le\ \frac{A_{\rm cap,-}-A_{\rm cap,+}}{r^4}\ \le\ \frac{\pi\sqrt3}{8}(1+r^2/0.97).
\]
An absolute cap-area error below $0.339\,r^4$ suffices to assign the correct geometry by the midpoint for any stated radius, if the true metric is one of this pair. This target uses one real finite cap. Even a successful comparison of two candidates would still require a verified directional, affine and proper-area calibration on the same native source before using its result for the general Einstein premises. [Direct tangential determinant quadrature](test_screens.py)

Reproduce the archived records and independent checks from the repository root.

~~~sh
export PYTHONPATH=code:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m matched_einstein_limit.source_implication.full_geometry.station_balls.timelike_areas.exterior_screens.measure
python -m pytest -q code/matched_einstein_limit/source_implication/full_geometry/station_balls/timelike_areas/exterior_screens
~~~
