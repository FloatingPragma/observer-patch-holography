# A fixed-volume area that distinguishes two moving-frame geometries

Six radar stations on the 80-cell equilibrium supply fifteen calibrated distances. Two smooth spatial geometries realize them all, and also give identical volumes and areas for every sufficiently small ball centered on the marked station at rest. Take the same tagged axial observer with velocity $u=(5/3,4/3,0,0)$ in both ultrastatic spacetimes. Direct volume integration on its geodesic rest slice followed by a separate boundary-area integration gives distinct readings. This calculation exhibits a specific observable missing from the station records. The moving aperture and its metric are supplied in the two constructions. [Source clocks and the two smooth completions](../README.md) · [Integrated radii, volumes and areas](results.json)

| Volume radius $\ell$ | Target proper volume | Proper area in $+$ | Proper area in $-$ | $(A_+-A_-)/(2s\ell^4)$ |
| ---: | ---: | ---: | ---: | ---: |
| 1/128 | $1.99737082\times10^{-6}$ | 0.000762993401320 | 0.000763027268680 | −4.54560002 |
| 1/256 | $2.49671352\times10^{-7}$ | 0.000191498450621 | 0.000191500540154 | −4.48723802 |
| 1/512 | $3.12089190\times10^{-8}$ | 0.0000479213381823 | 0.0000479214683587 | −4.47282968 |
| 1/1024 | $3.90111488\times10^{-9}$ | 0.0000119832524761 | 0.0000119832606056 | −4.46923882 |

## How the observer's ball is constructed

Take the round chamber in the existing station model with curvature radius $a=0.1/(\pi/2-\pi/20)$ in the calibrated station range unit. For a tangent vector with spatial normal coordinate $x$, the round metric $h_0$ has radial coefficient one and transverse coefficient $\operatorname{sinc}^2(|x|/a)$, where $\operatorname{sinc}(z)=\sin(z)/z$. Write $C_xv=x\times v$. The two spatial metrics add $\chi(|x|)C_x^\mathsf T B_\pm C_x$, with $B_\pm=2I\pm\operatorname{diag}(1,-1,0)$ in inverse squared range units. Thus $s=1$ in those units. Choose a smooth nonnegative cutoff equal to one for $|x|\le0.02$ and zero by $0.03$, inside the undeformed round part of the chamber. Both additions are nonnegative and annihilate the radial tangent, preserving the minimizing station paths. A quarter turn interchanges the two metrics locally and preserves all centered rest-ball measurements. The first axis remains tagged to a measured station path when the observer moves. [Exact general Ricci calculation](../smooth.py)

Let $y$ denote coordinates in the observer's orthonormal rest space and $\gamma=5/3$. The spacetime exponential sends $y$ to spatial normal coordinate $x=Ly$, $L=\operatorname{diag}(\gamma,1,1)$, and coordinate time $t=(4/3)y_1$. Indeed the radial lines are affinely parametrized geodesics for both spatial metrics and the time coordinate advances linearly in the ultrastatic product. The pulled-back rest metric is
\[
 H(y)=L^\mathsf T h_\pm(Ly)L-(\gamma^2-1)e_1e_1^\mathsf T .
\]
Its radial tangent satisfies $H(rn)n=n$. The radial coordinate is therefore intrinsic proper distance throughout the sampled normal neighborhood. The slice is spacelike throughout this neighborhood because $h_\pm\ge\operatorname{sinc}^2(0.02/a)I$, and $\gamma^2[1-\operatorname{sinc}^2(0.02/a)]<1$. The largest sampled ball and its volume solver bracket lie wholly inside the cutoff plateau.

For target volume $4\pi\ell^3/3$, solve separately in each metric for $r_\pm$ using
\[
 V(r)=\int_0^r \rho^2 d\rho\int_{S^2}\sqrt{\det H(\rho n)}\,d\Omega .
\]
Then compute $A(r_\pm)=r_\pm^2\int_{S^2}\sqrt{\det(E_n^\mathsf TH(r_\pm n)E_n)}\,d\Omega$ with $E_n$ an orthonormal tangent basis on the unit sphere. The solver never uses the small-ball curvature coefficient. For each native pair with calibrated distance $d_{ij}$ the construction joins its chambers by a cylinder of length $d_{ij}-0.2$. The two station-to-outlet paths each have length 0.1. The receipt carries all fifteen calibrated distance intervals and their cylinder lengths; those lengths determine the external gluing and the same interior patch works throughout their certified range. The table uses order-ten Gauss–Legendre quadrature in radius and polar cosine with twenty equally spaced azimuths. Raising both quadrature orders, checking the proper-volume residual and differentiating the integrated volume independently reproduce the boundary calculation. These are numerical quadrature checks. Areas and volumes are in squared and cubed range units. [Integral and root solver](measure.py) · [Independent resolution and geometry checks](test_measure.py)

## Why the signed limit is fixed

At the center the spatial Ricci tensors are $(2/a^2-6b)I\pm3s\operatorname{diag}(1,-1,0)$ with $b=2$. Scalar curvature agrees between the two completions. In the product spacetime the boosted contraction differs by $G_+(u,u)-G_-(u,u)=6s(\gamma^2-1)=32s/3$. The geodesic rest slice has zero second fundamental form at its center, hence intrinsic scalar curvature $2G(u,u)$. The fixed-volume area formula now predicts
\[
 \lim_{\ell\downarrow0}\frac{A_+-A_-}{2s\ell^4}
 =-\frac{4\pi}{15}\,3(\gamma^2-1)=-\frac{64\pi}{45}=-4.46804289\ldots .
\]
Each halving reduces the measured error by approximately four. A native measurement of these tilted areas would distinguish the two compatible completions, provided the instruments identify that boost, the same central event and the common target proper volume.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -B -W error -m matched_einstein_limit.source_implication.full_geometry.station_balls.timelike_areas.measure
python -m pytest -q code/matched_einstein_limit/source_implication/full_geometry/station_balls/timelike_areas
~~~
