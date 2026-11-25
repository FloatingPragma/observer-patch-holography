# Clock rates and fixed-volume areas

A clock can run at a different rate beside an energy source while local radar assigns the same spatial lengths to every pair of adjacent markers. In a static isotropic local model this happens whenever all material frequencies and reciprocal signal speeds acquire the same position-dependent multiplier. The fixed-volume rim then has exactly zero area response. This excludes a positive Einstein energy-density response for that completion of a whole-energy clock coupling.

## The spatial response required at the center

Normalize the central clock and a reference signal speed to one. Write the static material clock rate as $N(r)=1+a r^2+O(r^4)$ and the isotropic coordinate signal speed as $v(r)=1+s r^2+O(r^4)$. These expressions assume a regular spherical center, fixed comoving markers, reciprocal propagation and zero shift. Local out-and-back timing over a short interval $dr$ gives proper length $B(r)dr$, where $B=N/v$. A ring guide at the same radius gives circumference $2\pi rB$. The inferred metric on this comparison is

\[
 ds^2=-N^2dt^2+B^2(dr^2+r^2d\Omega^2),\qquad B=1+b r^2+O(r^4),\quad b=a-s .
\]

At the center the spatial Ricci tensor is $-8b\,\delta_{ij}$, the Hessian of $N$ is $2a\,\delta_{ij}$, and the static slice has zero extrinsic curvature. The Gauss constraint and static Ricci equations give the orthonormal components relative to flat vacuum,

\[
 G_{00}=-12b=12(s-a),\qquad G_{ii}=4(a+b)=4(2a-s),\qquad G_{0i}=0 .
\]

For density contrast $\rho$, isotropic pressure contrast $p$ and independently calibrated $\kappa=8\pi G$, a necessary test at this center is $12(s-a)=\kappa\rho$ and $4(2a-s)=\kappa p$. For a pressureless perturbation the required signal-speed coefficient is twice the clock-rate coefficient. A common multiplier $v=N$ gives $B=1$ for every profile and a Euclidean volume and area at every aperture. A central clock can detect a nonzero redshift in this very experiment.

The factor of two refers to the coordinate speed in the isotropic coordinate $r$. For example, the measured circumference defines the areal radius $R=rB(r)=r+b r^3+O(r^5)$. In this coordinate the radial signal speed is $dR/dt=1+(s+3b)R^2+O(R^4)$, whose quadratic coefficient is $-a$ for the pressureless solution. The rim circumference at matched volume and the redshift between material clocks have the same values in both coordinates.

The area test has its own normalization. For $V=4\pi\ell^3/3$, direct integration of $4\pi r^2 B^3 dr$ and evaluation of $4\pi r^2 B^2$ at the matched aperture give

\[
 \Delta A=\frac{16\pi}{5}b\,\ell^4+O(\ell^6)
         =-\frac{4\pi}{15}G_{00}\,\ell^4+O(\ell^6).
\]

Consequently the common multiplier misses the entire fourth-order signal at nonzero calibrated density. This conclusion uses the common rescaling of flight speed as a separate hypothesis. A native fluctuation of the seam generator could affect propagation, material rest frames or ruler formation in additional ways. [Finite tower and required area accuracy](../PROOF.md#geometric-normalization-and-the-three-errors)

## Finite comparisons

The ordinary fixture uses density one in its reference volume unit. The held fixture assigns the protected correlation-energy difference to the same volume unit. The run rechecks the four native energy words and their $1/4$ energy difference. The spatial embedding, static rates and density conversion belong to the comparison apparatus. [Native protected energies](../seam_overlap/README.md)

Both fixtures supply $N=1+\rho r^2/8$. The spatial control supplies $v=N/(1-\rho r^2/8)$. The rate-only control supplies $v=N$. At fixed markers both controls have the same clock-redshift profile. Each of the 32 to 256 radial shells carries midpoint rates, local round-trip paths and interior ring guides at two Gauss nodes. An independent central clock times consecutive emissions one proper-time unit apart at the rim. The volume decoder multiplies each clocked radial length by the mean of its two ring areas and sums the shells. It chooses the aperture using that interior sum alone. A separate full lap at the final rim measures its area. The decoder uses only integer ticks and their unit.

| Aperture radius | Inferred ordinary coupling with spatial response | Inferred held coupling with spatial response | Held coupling with common rate only |
| ---: | ---: | ---: | ---: |
| 0.40 | 1.5571510 | 1.5149597 | 0 |
| 0.20 | 1.5138166 | 1.5037125 | 0 |
| 0.10 | 1.5034261 | 1.5009262 | 0 |
| 0.05 | 1.5008548 | 1.5002314 | 0 |

The coupling estimator is $-15\Delta A/(4\pi\ell^4\rho)$. At the finest aperture the rate-only held clock reads a ratio of 1.000078125 and its rim area agrees with vacuum. In the spatial control the limiting ordinary and held couplings of $3/2$ follow from the supplied constitutive law. The held transfer tests the decoder and the changed speed response. A native closure would require clock and speed records on one preparation, an independently reconstructed stress and admission of the same proper balls in the regional theorem.

The rim redshift also tests the diagonal spatial component. At matched volume, $(N_{\rm rim}-1)/\ell^2=a+O(\ell^2)$. The clock and area therefore estimate $G_{ii}=4a-G_{00}/3$ independently of a pressure measurement. The spatial control gives $G_{ii}/\rho$ tending to zero. The common multiplier gives $G_{ii}/\rho=1/2$ despite its zero time component. Both components enter the local pressureless test.

The two-node ring rule integrates the coordinate $r^2$ weight exactly within each constant shell. For smooth $B=1+O(r^2)$, replacing $B^3$ by its midpoint value gives volume error $O(\ell^5/m^2)$ for $m$ shells. Rounding the radial and ring readings to a tick of size $t$ contributes $O(m\ell^2t)$ to volume. The corresponding area errors at matched volume are $O(\ell^4/m^2+m\ell t)$. The replay uses $m\propto1/\ell$, $t=\ell^5/2^{24}$ and an aperture search tolerance of one tick. The analytic errors are smaller than fourth order along this refinement. The executable reports the four tabulated finite panels.

```sh
PYTHONPATH=code .venv/bin/python -m matched_einstein_limit.rate_geometry.measure
PYTHONPATH=code .venv/bin/python -m pytest -q code/matched_einstein_limit/rate_geometry
```

## A reciprocal strain test on the native bridge

On a bridge with a dynamical strain register, ordinary rail energy and a change in two-head pairing both change subsequent material signals. At the first 80-cell long-wave appointment, an ordinary calibration predicts a faster crossing on decorrelation, while the independent material flight gives a slower crossing. A certified rank test identifies a fixed-marginal change of modular charge invisible to every probe of the final register. Two directional launch settings across four appointments allow the charge to be measured on repeated preparations. These experiments distinguish the information available to an observer from the response of a physical clock and a proper-volume ball. [Experiment, proofs and simulated counts](seam_strain/README.md)

## Which clock information fixes null curvature?

Even exact geometric light-return times at finitely many stationary clocks and every proper ball on one rest slice leave a component of null curvature undetermined. Work in four dimensions with central event \(p=(0,0)\), flat metric \(g_0\), and clock tracks at fixed spatial coordinates. Let \(b\) be any three-vector, and choose a smooth cutoff \(\psi\) supported away from every other station and equal to one near the center. Set \(f(t,x)=t(b\cdot x)\psi(x)\) and \(g_f=e^{2f}g_0\). On every clock track, \(f=0\) for every appointment. Conformal rescaling preserves the unparametrized null trajectories and their endpoint events. Every ideal geometrical light-flight time between those clocks is identical, including sequences with reflections. On the whole slice \(t=0\) the induced spatial metrics coincide. Every intrinsic rim area and rest volume there agrees. These statements hold for each value of \(b\), so varying \(b\) also preserves their first variations. The stationary clocks are held on their tracks. An independently observed acceleration, an interior moving clock or a calibrated affine distance carries different information.

At \(p\), the gradient of \(f\) vanishes and its only potentially nonzero Hessian entries are \(f_{0i}=f_{i0}=b_i\). For \(k=(1,n)\), \(|n|=1\), the conformal Ricci formula gives \(G_f(k,k)=R_f(k,k)=-4 b\cdot n\). The two signs of \(b\) give opposite null curvatures with identical clock and slice data. The independently specified conserved vacuum stress \(T=0\) satisfies null balance in the flat metric and violates it in a direction with \(b\cdot n\ne0\) in the perturbed metric. This is an ambiguity in a continuum interpretation of finite observations. Detector probabilities and material accelerations can carry information beyond ideal light-path times. A family with uniformly controlled derivatives and interior observations requires its own refinement estimate. For recovery of conformal structure from light observations see [Kurylev, Lassas and Uhlmann](https://arxiv.org/abs/1405.3386).

Suppose path data determine a smooth reference metric \(\bar g\) up to the positive scale \(g=w^{-2}\bar g\). Normalize \(w(p)=1\) with a material clock. Along a \(\bar g\)-affine null ray through \(p\), normalized to unit energy in the central frame, the four-dimensional conformal Ricci formula gives
\[
 R_g(k,k)=R_{\bar g}(k,k)+2w''(0).
\]
An independently calibrated local clock followed in a known timelike direction gives \(w=d\bar\tau/d\tau_g\). The clocks here sample the three events at ray parameters \(-h,0,h\). If the three scale readings have absolute errors at most \(\epsilon\) and \(|w^{(4)}|\le M\) on the segment, twice their central second difference estimates the curvature correction with error at most \(Mh^2/6+8\epsilon/h^2\). Taylor expansion bounds the truncation error. Absolute summation of the stencil weights bounds the clock error. Errors in the reference curvature, ray parameter or comparison frame require separate bounds. Clock errors \(o(h^2)\) and a bounded fourth derivative suffice for this scale contribution to converge. In the example \(w(s)=\exp[-(b\cdot n)s^2]\) near the origin, which gives the signed curvature above. The nine-direction Einstein theorem can consume these geometric readings once an independent stress instrument gives the matching nine contractions at the same coupling. [Tensor theorem](../PROOF.md) · [Derivative and error checks](test_null_scale.py)

## Static null curvature from a directional rim and clock array

Take a smooth static metric \(g=-N^2dt^2+h\) with \(N>0\). At a common center \(p\), in its static orthonormal frame, let \(k=N^{-1}\partial_t+n\) with \(h(n,n)=1\). The warped-product connection gives
\[
 R_g(k,k)=\operatorname{Ric}_h(n,n)+
  [\Delta_h N-\nabla_h^2N(n,n)]/N .
\]
Indeed the additional nonzero Christoffel symbols are \(\Gamma^0_{0i}=\partial_i\log N\) and \(\Gamma^i_{00}=N\nabla^i N\). Substitution gives \(R_{00}=N\Delta_h N\), \(R_{0i}=0\) and \(R_{ij}=R^h_{ij}-N^{-1}\nabla_i\nabla_jN\), independently of a field equation.

Let \(j(r,n)\) be the proper shell area per unit launch solid angle at spatial geodesic distance \(r\) on the static slice. Two transverse Jacobi fields start at zero with orthonormal unit derivatives. Each is \(r e_a-r^3 R(e_a,n)n/6+O(r^4)\) in its parallel frame. Their area gives \(j/r^2=1-r^2\operatorname{Ric}_h(n,n)/6+O(r^3)\). Read the clock rate \(N(\exp_p rn)/N(p)\) independently by a static clock comparison. Average each reading at opposite directions before forming a quadratic coefficient. Then
\[
 D_n=-6(j_{\rm even}/r^2-1)/r^2,\qquad
 L_n=2(N_{\rm even}/N(p)-1)/r^2
\]
converge to spatial Ricci and the normalized lapse Hessian in direction \(n\). On the six icosahedral axes, let \(\bar L\) be the arithmetic mean. Because \(\sum_i n_i n_i^T=2I\), the estimate of null Ricci is \(D_n+3\bar L-L_n\). The six projectors have diagonal Gram entries 1 and off-diagonal entries 1/5. They span all even quadratic forms. Under staticity a complete matter comparison also requires independently vanishing mixed time-space stress. Adequate readouts have absolute error \(o(r^4)\) in area per unit solid angle and \(o(r^2)\) in clock ratio, including angular calibration and proper-radius errors. A finite-width aperture needs its actual angular moments or a vanishing aperture bias. A general time-dependent spacetime has additional components. [Small geodesic ball expansions](https://archive.ymsc.tsinghua.edu.cn/pacm_download/117/6254-11511_2006_Article_BF02395060.pdf)

The native pole and its five neighbors contain exactly one port from each antipodal pair. Equal weighting on this menu has rank one on the six even directional readings. Doubling the weight of the pole gives the averaged reading \(m_i=(6\bar q+q_i)/7\) and hence \(q_i=7m_i-6\bar m\). This recovers all six components on the finite angular design. If the six normalized aperture estimates of the spatial Ricci coefficient have independent absolute error bounds $e_D$ and the six lapse-Hessian estimates bounds $e_L$, each recovered null contraction has error at most $11e_D+13e_L$, apart from the geometric expansion errors. The spatial weights are 6 at the selected pole and −1 at the other five. The lapse weights are −11/2 and 3/2 respectively. These absolute-sum bounds also apply to correlated errors, and each is sharp for freely signed component errors. On one cell the native one-load Hamiltonian preserves the three coordinate waves of the first spatial triplet. Prepare the axial member or an equiprobable mixture of its two transverse members, let it propagate and read every terminal address. Every ordinary star has probability exactly 1/2 for both preparations. Scoring the pole one, each neighboring hit one half and all other hits zero gives 3/8 and 1/4 at the selected pole. The simulated scores are 0.37535 and 0.25013 in 1,048,576 trials per preparation. This reads a material angular distribution. Calling it the shell area of a proper ball requires an additional physical relation. [Native counts, angular ranks and static geometric controls](static_directions.json)

An independent signed control uses the prescribed spatial product of a two-surface of sectional curvature \(K\) and a line, with lapse \(1+c z^2/2\) near the center. Along geodesics the radial Jacobi sectional curvatures are zero and \(K(1-n_z^2)\). Solving the scalar Jacobi equation gives the area and clock data at finite radius. Direct static Ricci contraction gives \((K+c)(1-n_z^2)\). The angular decoder recovers both signs with maximum absolute error below 0.00021 at radius 1/16 for \((K,c)=(1,-2),(-1,2)\) and exactly zero in the flat unit-lapse control. [Calculation and tests](test_static_directions.py)

There is also a static ambiguity in whole-shell measurements. In Minkowski coordinates fix finitely many spatial clock positions including the origin and choose a smooth radial cutoff supported away from every other clock, equal to one near the origin. Let \(f(x)=a(x_1^2-x_2^2)\chi(|x|)/2\). The metrics \(g_\pm=e^{\pm 2f}\eta\) give all these clocks identical proper time and zero proper acceleration at all times. They have the same unparametrized light paths and clock-to-clock optical times. Interchanging the first two spatial axes is an isometry between the two induced metrics fixing the origin. Every centered geodesic ball therefore has exactly the same total proper volume and rim area in the two metrics, for all normal radii and at every appointment. The clock and ray labels stay fixed in the comparison. At the origin \(R_\pm(k,k)=\mp2a(n_1^2-n_2^2)\), since \(df=0\) there. Thus even with exact static scalar data, a labeled ray has two possible curvatures. The two models have different angular shell and interior lapse data. This comparison concerns geometric optics and centered scalar shells. It establishes no degeneracy of the entire microscopic law.

**First-order stationarity for every rest ball in an open region.** Let \(g_a=e^{2au}\eta\) in the same example, with \(u=(x_1^2-x_2^2)/2\) on the neighborhood where the cutoff equals one. The common source can be vacuum, \(T=0\) and constant zero matter entropy. Consider any center \(x\) and radius \(r\) with its closed Euclidean ball inside that neighborhood. At \(a=0\) the proper-distance variation at the endpoint \(x+rn\) is \(\int_0^r u(x+sn)\,ds\). This follows by first variation of length at the unperturbed straight geodesic, holding its endpoints fixed. At fixed proper radius the coordinate endpoint radius thus changes by the negative of that integral. The spatial volume form contributes \(3u\,d^3x\). Combining both variations gives
\[
 \partial_a V_a(x,r)|_0=\int_{S^2}\left[3\int_0^r s^2u(x+sn)\,ds
       -r^2\int_0^r u(x+sn)\,ds\right]d\Omega=0 .
\]
The spherical mean equals \(u(x)\) because \(u\) is harmonic. Equivalently expand this quadratic at \(x\): its linear term averages to zero and its Hessian is traceless. Geodesic coarea gives \(A_a=\partial_r V_a\). Smoothness in a normal neighborhood then gives zero area variation at fixed radius and at fixed volume. Every such ball on every static slice satisfies first-order generalized-entropy stationarity in its rest frame. The stress Ward identity is exact. But \(\partial_a R_{g_a}(k,k)|_0=-2(n_1^2-n_2^2)\) in the background null frame. This gives a flat reference obeying balance and a vacuum perturbation violating linearized balance while all of those rest-frame tests agree to first order. Boosted observer balls, directional areas and interior lapse measurements have additional information. These continuum models test the logical sufficiency of the proposed scalar inputs. They are not claimed as microscopic Interval Echo solutions.
