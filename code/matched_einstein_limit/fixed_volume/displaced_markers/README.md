# Measuring equal-volume area from displaced interior counts

A single observer can count every angular sector of every nested ball and still miss its curvature. A positive, spatially varying marker intensity can exactly cancel the volume distortion. Six nearby observers resolve this ambiguity if their clock rulers realize geodesic distances on the same spatial slice. They can recover the leading fixed-proper-volume area response without surface markers, an entropy-to-area calibration, or knowledge of the local marker intensity. On five supplied geometries the largest error in the mean scalar-curvature estimate falls from **0.34065 to 0.005082** over four refinements. A separate surface census checks the result. [Simulated counts and comparisons](results.json) · [Counter and decoder](measure.py)

## An exact obstruction even with all centered records

In a three-dimensional space of constant sectional curvature $K$, the proper volume in geodesic polar coordinates is $S_K(r)^2\,dr\,d\omega$. Here $S_K(r)$ is the proper circumference radius, equal to $r$ in flat space. Enroll local markers of intensity $m_K(r)=(r/S_K(r))^2$ per unit proper volume. This smooth, positive function extends to one at the center. Then the intensity measure is exactly $r^2\,dr\,d\omega$. All centered radial and angular bins, their subdivisions, intersections, count sums and unions have identical means at every radius below the cut locus. A Poisson process with this measure gives identical full joint count laws, including coincidences in overlapping bins. The common central marker calibration and all the additivity relations therefore coexist with either sign of curvature. Independently calibrated marker density throughout the ball would break this ambiguity.

The discrete spherical capacity calculation counts $20\,4^L$ screen cells and assigns a record capacity per cell. Relating this to an area in proper ruler units supplies further information beyond the count and its additive subdivision. The counterexample concerns the inference from counts with unknown local density. It leaves open a derivation from specified native instruments. [Discrete sphere calculation](../../../capacity_readback/discrete_spherical_refinement.py)

## The eight readings

Enroll a smooth Riemannian spatial slice, center $x$, an orthonormal frame $v_1,v_2,v_3$ at $x$, and an independent proper-length scale. Each center uses its own geodesic ball in this common slice. Every bank samples the same positive, smooth marker intensity $m$ up to a common multiplier. Normalize the mean count at radius $s$ by the known number of repetitions, instrument gain and flat reference volume $4\pi s^3/3$, giving $C(y,s)$. At $x$ take $C=C(x,\ell)$ and $H=C(x,\ell/2)$. At each of the six points $\exp_x(\pm\ell v_i)$ take $C(y,\ell/2)$. Call their average $B$. Then

$$
\widehat m=(4H-C)/3,\qquad
\widehat R=\frac{18B+22H-40C}{\widehat m\,\ell^2}
           =\mathrm{Scal}(x)+O(\ell^2).
$$

In normal coordinates the proper volume density is $1-\mathrm{Ric}_{ij}z^iz^j/6+O(|z|^3)$. Flat ball averages have second moments $s^2\delta^{ij}/5$. Odd Taylor terms integrate to zero, giving, locally uniformly in $y$,

$$
C(y,s)=m(y)+s^2q(y)+O(s^4),\qquad
q=\Delta m/10-m\,\mathrm{Scal}/30 .
$$

It follows that $(C-H)/\ell^2=3q/4+O(\ell^2)$. Along each pair of opposite geodesics the linear term in $m$ cancels. Their mean gives $(B-H)/\ell^2=\Delta m/6+O(\ell^2)$ because the change of $q$ contributes only at fourth order. Taking $18(B-H)-40(C-H)$ eliminates the unknown density Hessian. Richardson extrapolation gives $\widehat m=m+O(\ell^4)$. The factor common to all marker banks cancels in the quotient. [Geodesic volume expansion](https://doi.org/10.1307/mmj/1029001150)

**Which bearings suffice?** Allow any finite collection of opposite unit directions $v_j$. To second order the paired displacements identify the quadratic readings $\nabla^2m(v_j,v_j)$ and the two central radii identify $q$. Scalar curvature is identifiable for arbitrary local marker intensity exactly when $I/3$ belongs to the linear span of $v_jv_j^\top$. If $\sum w_jv_jv_j^\top=I/3$, trace gives $\sum w_j=1$ and the weighted mean of the paired readings gives the decoder above. Conversely, if $I/3$ is outside the span, there is a symmetric $E$ orthogonal to all these dyads with $\operatorname{tr}E\ne0$. Replace the local density Hessian by $\nabla^2m+E$ and the scalar curvature by $\mathrm{Scal}+3\operatorname{tr}(E)/m$. Every measured second-order coefficient agrees. Constant-curvature metrics and smooth densities with the prescribed jets realize both choices locally. This necessity concerns the leading data; higher-order measurements on additional source assumptions can carry more information. For nonnegative weights the rotational calibration error adds at most $9|\operatorname{tr}[(\sum w_jv_jv_j^\top-I/3)\nabla^2m]|$ to the inferred product $m\,\mathrm{Scal}$ at leading order.

## Accuracy needed for the Einstein premises

Suppose the normalized-count expansion has remainder at most $L s^4$. Along the surveyed geodesics bound the fourth directional derivative of $m$ by $M$ and the second directional derivative of $q$ by $Q$. With the six axial sites, if every normalized reading has absolute error at most $\epsilon$, the numerator divided by $\ell^2$ differs from $m\,\mathrm{Scal}$ by at most

$$
\left(\frac{3M}{4}+\frac{9Q}{4}+\frac{179L}{4}\right)\ell^2
     +\frac{80\epsilon}{\ell^2}.
$$

Indeed the paired displacement remainder in $B-H$ is at most $(M/24+Q/8+L/8)\ell^4$ and the centered remainder in $C-H$ is at most $17L\ell^4/16$. The reference-density error is at most $5L\ell^4/12+5\epsilon/3$. A known positive lower bound on $m$ transfers both estimates through the quotient. Mis-setting an individual ball radius by $\delta r$ changes its normalized volume by $O(\delta r/\ell)$, so independent range errors $o(\ell^3)$ suffice. Center displacements must also respect the paired geodesic arrangement to the required accuracy; a center-position error $o(\ell^2)$ in proper distance suffices for bounded density gradient and smooth ball means. A fixed-size ruler error does not survive arbitrary refinement.

At a common target proper volume $4\pi\ell^3/3$, the smooth area is

$$A|_V=4\pi\ell^2-\frac{2\pi}{15}\,\mathrm{Scal}(x)\,\ell^4+O(\ell^6).$$

This follows by expanding the volume, differentiating the smooth normal-coordinate integral for its area, and solving for its radius at the prescribed volume. Using $\widehat R$ makes the finite count contribution to the inferred area $O(\epsilon\ell^2)$ for bounded curvature and positive density. Thus $\epsilon=o(\ell^2)$ gives the $o(\ell^4)$ area error required by the matched theorem. With independent Poisson counts and bounded positive intensities, an expected count budget per setting growing faster than $\ell^{-4}$ gives this precision in probability. With fixed physical density, pooling independent copies requires a repetition budget growing faster than $\ell^{-7}$. The elementary counting intervals below cover sampling error about the finite-radius mean. Curvature bias is separate. [Finite tower error requirements](../../PROOF.md) · [Matched-volume geometry](../README.md)

For a source tangent, either obtain these error bounds also for the differentiated readings or divide a controlled preparation contrast by its amplitude. For a symmetric source secant of amplitude $a\to0$ with uniformly bounded third scalar-response derivative, independently reconstructing the two area coefficients has normalized derivative error $O((\ell^2+\epsilon/\ell^2)/a+a^2)$ for an exact common frame. The geometric Taylor remainders at the two endpoints contribute the same $O(\ell^2/a)$. Taking $\ell^2+\epsilon/\ell^2=o(a)$ suffices. A source-specific unknown $m$ is allowed on either arm. Applying this formula to the entropy assay requires an independently identified common physical variation.

On a slice with zero extrinsic curvature at the central event, the Gauss relation gives $\mathrm{Scal}=2G(u,u)$. Recovering different stress contractions requires further spacetime frames, or independent constraints on the missing components. Ordinary matter fixes the physical coupling on its own transfer. For the held collar tangent the measured $\delta A|_V$ must then cancel the independently read bulk entropy at that frozen coefficient. Existing native forked relay lengths fail the exact smooth geodesic test on their full table, and the registered two-load experiment supplies no sequence of admissible displaced shrinking geodesic balls. The smooth comparison certifies an inference that would work with such observer records. [Native range tests](../../../../model/experiments/interval_echo_einstein_station/equilibrium_refinement/rim_geometry/README.md) · [Distinct timelike interval construction](../../source_implication/bisected_diamonds/README.md)

## Counts and controls

Use the five constant-curvature spatial metrics and four reference radii in the independent equal-volume surface experiment. Plant the centered-blind marker $m_K$ above on each geometry. Each central mean is identically one. For the displaced mean use polar coordinates about the original center. At radial distance $r$, include the solid-angle cap satisfying the space-form geodesic cosine law for distance at most $\ell/2$ from the displaced observer. Integration of $r^2\,dr\,d\omega$ over this cap supplies the expected count. The forward simulator knows $K$. The decoder gets the eight counts, repetition factors and physical scales. For the additional flat control use $m=1+r^2/4$. Interpreting its central two-radius difference as curvature at uniform density yields $-4.5$. The displaced decoder gives zero in expectation.

At $\ell=1/16$ the sampled scalar estimates for $K=-1,+1$ are **−5.99771 and +6.00579**. Independent archived matched-area counts give **−5.99813 and +6.00230**. At the same radius the largest deterministic displacement bias across the five metrics is **0.005082**. At the four radii the detector budgets are $2^{17}\ell^{-6}$ per setting, achieved by aggregation. Fresh Poisson draws on all eight banks give two-sided Garwood intervals for each normalized mean. The union bound covers all 192 means at 99%, and interval arithmetic propagates these through the quotient. The independent surface counts never enter the displacement decoder. All comparisons concern supplied smooth metrics. [Archived surface census](../measurements.json) · [Finite-frame and flat-density checks](test_displaced.py)

From the repository root run the Python module with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -m matched_einstein_limit.fixed_volume.displaced_markers.measure
python -m pytest -q code/matched_einstein_limit/fixed_volume/displaced_markers
~~~
