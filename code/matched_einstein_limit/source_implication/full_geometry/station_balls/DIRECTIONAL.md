# Which additional local measurements recover curvature?

A single station can recover all six components of spatial Ricci from thirteen regional proper volumes. In spacetime, twelve signed finite-aperture screen-area comparisons can recover all nine independent null contractions, including energy-flux directions. These are sufficient geometric inputs with explicit uniform error bounds. The six clocks on the selected Interval Echo equilibrium currently give tagged station-pair detection times. They do not count these proper volumes or measure beam areas at calibrated affine distances. The geometry comparison preserves all fifteen calibrated lengths while changing the proposed local readings. [What the measured source fixes](README.md) · [Earlier static infinitesimal-rim and clock calculation](../../../rate_geometry/README.md)

## Spatial volume protocol

At a measured event select an orthonormal observer frame. Let \(n_i\) represent one unit direction on each of the six antipodal axes of a regular icosahedron. This menu reuses the existing twelve-port angular design. The six spatial positions of the radar survey are a different set. In particular the physical identification of these angular labels with tangent directions at the selected station would have to be measured. Rotate each launch label with the observer frame.

Take \(a=1/2, b=1/4\). For a normal radius \(r\) let \(U_n=\exp_p\{r(a n+b z):|z|<1\}\). Also take \(U_0=\exp_p\{br z:|z|<1\}\). Measuring membership requires a proper radial ruler and tagged launch directions. The shifted regions are images of Euclidean subballs in the central tangent space. They may differ from metric balls centered at their displaced endpoints. Every region lies inside the normal ball of radius \(r\). Write \(V_n=\operatorname{Vol}(U_n)\), and \(V_E=4\pi(br)^3/3\). Then
\[
 q_i=-\frac{6[(V_{n_i}+V_{-n_i})/2 - V_0]}{a^2r^2V_E}
   = n_i^\mathsf{T}\operatorname{Ric}_h n_i + E_i.
\]
In normal coordinates the quadratic term of proper-volume density is \(-\operatorname{Ric}_h(x,x)/6\). A tangent subball has second moment \(r^2(a^2nn^\mathsf T+b^2 I/5)\). Subtracting the central subball cancels the trace contribution. Averaging the antipodal regions cancels the cubic term.

Let \(P_i=n_i n_i^\mathsf T\). The six Gram entries per row are one on the diagonal and \(1/5\) elsewhere, with \(\sum P_i=2I\). Thus
\[
 D(q)=\tfrac54\sum_i q_i P_i-\tfrac14 I\sum_i q_i
\]
inverts the six quadratic readings. The smallest singular value of the measurement map in Frobenius norm is \(\sqrt{4/5}\). Under a rigid rotation of the launch frame the decoded tensor transforms by conjugation. The Frobenius error is at most \(\sqrt{15/2}\max_i|E_i|\).

**Explicit scale bound.** On each relevant segment assume the parallel radial tidal matrix has operator norm at most \(K\) and its second covariant derivative along the segment at most \(H\). Use \(Kr^2\le1\) within a normal region. Set \(C=13K^2/90+H/20\). If each proper volume has absolute error at most \(\epsilon_V\), then
\[
 |E_i|\le\frac6{a^2}C[(a+b)^4+b^4]r^2+
               \frac{12\epsilon_V}{a^2r^2 V_E}.
\]
Bounded \(K,H\) and \(\epsilon_V=o(r^5)\) give convergence. A source-variation or linear-response calculation must control the variation of these errors at its preparation spacing.

To verify the constant, in a parallel screen basis the two-by-two radial Jacobi matrix satisfies \(Y''=-\mathcal RY\), \(Y(0)=0\), \(Y'(0)=I\). Iterating its Volterra integral equation gives \(\|Y(t)\|\le \sinh(\sqrt K t)/\sqrt K\le2t\) and \(\|Y(t)/t-I\|\le Kt^2/3\). Inserting these bounds once more, the Taylor expansion of \(\mathcal R(t)\) with error \(Ht^2/2\) yields
\[
 \left\|Y(t)/t-I+\mathcal R(0)t^2/6+\mathcal R'(0)t^3/12\right\|
 \le (K^2/60+H/40)t^4 .
\]
In dimension two, \(\det(I+T)=1+\operatorname{tr}T+\det T\). Its trace and determinant errors sum to at most \(Ct^4\). In spatial normal coordinates \(\det(Y(t)/t)\) is the volume density at radial coordinate \(t\). The cubic trace reverses sign under \(x\mapsto-x\). The same calculation holds on the two-dimensional spacetime null screen using affine distance.

## All null contractions from finite solid-angle caps

At a common event choose a physical orthonormal tetrad, calibrate affine null coordinates to \(-u\cdot k=1\) and let \(k=u+n\). Write \(M_{ij}=R_{ij}+R_{00}\delta_{ij}\), \(p_i=R_{0i}\), so \(R_{kk}=n^\mathsf T M n+2p\cdot n\). These nine numbers determine every null contraction. Adding any multiple of the spacetime metric to Ricci leaves them invariant.

Measure the proper screen area swept out by each signed axis's solid-angle cap of half-angle \(60^\circ\) on the light cone at affine parameter \(r\), before a caustic. The cap solid angle is \(\pi\). To cancel the cubic remainder, also use the past cone at the same affine magnitude along the **reverse spacetime vector** \(-k\). This pairing differs from a future launch in the opposite spatial direction. In a time-asymmetric spacetime this requires up to twenty-four raw areas to make twelve paired readings. Define \(w_n=-6[(A^{\rm future}_n+A^{\rm past}_{-k})/(2\pi r^2)-1]/r^2\). With area error \(\epsilon_A\) per cap the preceding tidal bound gives a uniform error \(F=6Cr^2+6\epsilon_A/(\pi r^4)\). A single future-only acquisition is also consistent with a bounded first radial derivative, retaining a first-order truncation term. To achieve second-order accuracy with future-only launches, evaluate the twelve caps at two radii, using $2w_n(r/2)-w_n(r)$ where $w$ here uses future data only. This cancels the cubic density contribution and bounds the error by $9Cr^2+198\epsilon_A/(\pi r^4)$. In the explicit quadratic comparison below each measured cap agrees with its past partner by symmetry.

For a unit direction \(m\) uniformly distributed in this cap, azimuthal integration and the uniform distribution of \(n\cdot m\) on \([1/2,1]\) give
\[
 \langle m\rangle=m_1 n,\quad
 \langle mm^\mathsf T\rangle=\alpha I+\beta n n^\mathsf T,\quad
 (m_1,\alpha,\beta)=(3/4,5/24,3/8).
\]
Hence \(w_n=\alpha\operatorname{tr}M+\beta n^\mathsf T M n+2m_1p\cdot n\) up to \(F\). Form \(S_i=(w_{n_i}+w_{-n_i})/2\), and \(L=D(S)\). Decode
\[
 \widehat M=(L-\alpha\operatorname{tr}(L)I)/\beta,\qquad
 \widehat p=\tfrac1{8m_1}\sum_i n_i (w_{n_i}-w_{-n_i}).
\]
Indeed \(3\alpha+\beta=1\), and the inverse angular smoothing expands a tensor by at most \(1/\beta\) in Frobenius norm. The errors obey \(\|\widehat M-M\|_F\le \sqrt{15/2}\,F/\beta\) and \(\|\widehat p-p\|_2\le\sqrt3 F/(2m_1)\). Thus \(\epsilon_A=o(r^4)\) suffices under the stated uniform bounds. Apertures can overlap because each cap may be measured on a fresh preparation.

**A signed test of null balance.** Suppose a second instrument on the same prepared events independently reads the corresponding cap averages of a single symmetric stress tensor. Keep the matter coupling $\kappa$ fixed at its ordinary-source calibration. If all twelve paired cap residuals are at most $\delta$ in absolute value including geometric truncation, stress and measurement errors, the preceding inverse gives
\[
 \sup_{|n|=1}|R_{kk}-\kappa T_{kk}|
 \le \big(\sqrt{15/2}/\beta+\sqrt3/m_1\big)\delta < 9.62\delta .
\]
This follows from $|n^\mathsf{T}\Delta M n+2\Delta p\cdot n|\le\|\Delta M\|_F+2\|\Delta p\|_2$. Using the twelve signed residuals jointly improves the constant. For $z_i=n\cdot n_i$ set $b_i=(10/3)z_i^2-17/18$. The coefficients multiplying the two signed residuals on axis $i$ in the recovered contraction are $b_i/2+z_i/3$ and $b_i/2-z_i/3$. The frame identities $\sum z_i^2=2$ and $\sum z_i^4=6/5$ hold for every unit $n$. The second follows by applying the given tensor inverse to the rank-one tensor $nn^\mathsf T$ and contracting again with $n$. Hence the sum of all twelve squared coefficients is $377/108$. Cauchy–Schwarz gives the uniform, signed end-to-end bound
\[
 \sup_{|n|=1}|R_{kk}-\kappa T_{kk}|\le (\sqrt{377}/3)\,\delta < 6.48\delta .
\]
The bound accounts for arbitrary correlated signs across all twelve residuals. Thus vanishing independently bounded residuals on one refining family would force the local signed balance in every null direction, with no refit of the coupling. The finite theorem's other conditions, including a common covariant metric, stress and conservation across refinements, remain necessary.

## Angular records on the actual interval-echo source

There is a finite angular measurement with the required nine independent labels on the **same selected 80-cell preparation**. Launch a native mobile at port zero and inspect its location once after dwell 1.5. Leave the existing head state and fair load-charge mixture in place. Enroll the twelve ports in its first hardware cell and keep the outside outcome as a thirteenth bin. A weighted directional panel scores one at its pole, one half on each of its five neighbors and zero elsewhere. All twelve panels can be read from this common port record.

On the twelve native labels the neighbor matrix has eigenvalues $5$, $\sqrt5$ and $-1$ on the constant, the three dipoles and the five traceless quadratics respectively. The middle result follows by summing the five neighboring directions of a pole. For the last, their five azimuths are equally spaced with axial cosine $1/\sqrt5$, so their sum of traceless angular projectors is the negative of the pole's traceless projector. The unweighted pole-and-neighbors menu consequently has rank four across these nine components. Doubling the pole weight gives transfer eigenvalues $7$, $2+\sqrt5$ and $1$ before division by the score normalization two. The measured panel therefore has **rank nine** on those labels. The coefficients of this discrete angular instrument differ from the continuous solid-angle cap moments above.

With $2^{26}$ synthetic independent native starts the largest antipodal odd amplitude is **0.4723675** and the even panel spread is **0.2285594**. The union-bounded simultaneous 99% score error, including state preparation enclosure, is below **0.000243**. A homogeneous angular population has zero odd signal and zero even spread. The actual counts measure mobile location probability. Converting those counts into proper screen area requires independent affine range, ray launch normalization, spreading and detection calibration. An anisotropic point-source detector distribution by itself cannot set any value of physical Ricci. [Fresh 13-bin counts](port_apertures.json) · [Native angular experiment](port_apertures.py)

## Independent smooth controls

The spatial control sets \(g_x(v,w)=v\cdot w+(x\times v)^\mathsf TB(x\times w)\) with a non-diagonal positive matrix \(B\). Straight radii remain unit-speed geodesics, so ordinary Cartesian coordinates are the normal coordinates at zero. The producer integrates \(\sqrt{\det g}\) over the thirteen domains. Their target Ricci is \(3(B-\operatorname{tr}(B)I)\), verified independently by the coordinate derivative calculation. For the null control prescribe an anisotropic conformal spacetime \(e^{2\phi}\eta\) with \(\phi(t,x)=-x^\mathsf TMx/4-tp\cdot x/2\). The conformal Ricci formula gives the target at the origin directly. Along the null path \((t,x)=(v,vn)\) let \(q=n^\mathsf T M n+2p\cdot n\). Its physical affine distance is \(\int_0^v e^{-qz^2/2}dz\) and the physical screen area per solid angle is \(v^2 e^{-qv^2/2}\). The producer inverts this affine relation separately for every incoming direction and integrates the measured area over each cap. These controls do not assign a continuum geometry to the native clocks. [Implementation](directional.py) · [Volumes, areas and convergence errors](directional.json)

From the repository root, run the modules ending in `station_balls.directional` and `station_balls.port_apertures` in the same environment as the [six-clock replay](README.md). Run `pytest -q -W error code/matched_einstein_limit/source_implication/full_geometry/station_balls` for the complete test suite.
