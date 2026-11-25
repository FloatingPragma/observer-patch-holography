# Where two elliptical cluster components can put their peaks

The two major peaks of the Bullet Cluster lensing reconstruction lie outside opposite ends of the interval between the X-ray plasma centres along a measured sky coordinate. No sum of two monotone elliptical projected profiles centred on the plasma can place both strict peaks there. Each contribution may have its own orientation, axis ratio, normalization and radial profile. The test uses the published peak and plasma positions. [Clowe et al., §3 and Table 2](https://arxiv.org/pdf/astro-ph/0608407) · [Coordinates and calculation](results.json)

## Two-centre theorem in the sky plane

Let the gas centres be $0$ and $v$. Let the projected convergence be a constant plus $f(x^TAx)+g((x-v)^TB(x-v))$, with $A,B$ positive definite symmetric $2\times2$ matrices and $f,g$ nonincreasing. This permits sharp edges and plateaux. At any point off both centres the two outward normals are $Ax$ and $B(x-v)$. Unless they point in exactly opposite directions, a small step opposite the sum of their unit vectors strictly reduces both quadratic radii. The convergence cannot fall along that step. Hence every off-centre strict maximum satisfies

$$
Ax+tB(x-v)=0,\qquad t>0.
$$

Write $M=A^{-1}B$. Its eigenvalues are positive since it is similar to the positive symmetric matrix $A^{-1/2}BA^{-1/2}$. Put $T=\operatorname{tr}M>0$ and $D=\det M>0$. Solving the two-dimensional linear equation with the adjugate identity gives

$$
x(t)=(I+tM)^{-1}tMv=\frac{tMv+t^2Dv}{1+tT+t^2D}.
$$

Choose any linear sky coordinate $u$ oriented so that $C=u(v)\ge0$. Put $b=u(Mv)$. A projected peak below zero implies $b<0$, because the numerator is $tb+t^2DC$. A projected peak above $C$ implies $tb>C+tTC$, hence $b>0$. These conditions are incompatible, even at different parameter values. Peaks at the centres are interval endpoints. This proof concerns two dimensions and two fixed shape matrices. The scalar projection inequality is machine-checked as *twoCentrePeakProjection*. [Lean proof](../../../../Lean/ObserverPatchHolography/EinsteinBranch/HaloStress.lean) · [Example of a strict peak beyond one centre](test_peaks.py)

## The measured displacements

The western peak is 7.1 arcsec east and 6.5 arcsec north of the subcluster BCG (brightest cluster galaxy). The eastern peak is 21.2 arcsec west and 17.7 arcsec north of the southern main BCG, or 2.5 arcsec east and 11.5 arcsec south of the northern one. Table 2 lists one main BCG position. Applying each offset to that position brackets the assignment. On the local east/north sky plane, take the unit coordinate vector $(0.992308,-0.123793)$ along the smaller eastern displacement. [Positions in §3 and Table 2](https://arxiv.org/pdf/astro-ph/0608407) · [Archived sky coordinate conversion](../bullet_cluster.py)

| Projected distance | Arcsec |
|---|---:|
| Western peak beyond western plasma | 36.475 |
| Separation of the two plasma centres | 75.754 |
| Eastern peak beyond eastern plasma, southern reference | 21.811 |
| Eastern peak beyond eastern plasma, northern reference | 48.943 |

Suppose each peak position can move by at most $p$ and each source centre by at most $s$ in the sky plane. If $s+p<21.81057$ arcsec, both peak projections retain their strict ordering beyond the two centres. Crossing the opposite gas centre also requires a larger displacement. With peak positions held fixed, at least one effective centre must therefore move by **96.25 kpc** at the published plate scale. With both gas centres held fixed the same lower bound applies to the largest peak displacement. This is a necessary deterministic displacement under the model. The authors report bootstrap offsets for individual peaks. A sampling probability for the two-peak projection requires the joint centroid distribution. The authors report positions of reconstructed peaks and estimate their centroid uncertainties by bootstrap. Any displacement between the reported centroid and a strict maximum of the modelled convergence at that resolution belongs in $p$.

The western peak alone bounds the ellipticity. The centre directions from that peak make an angle of $\delta=1.53136^\circ$. For an ellipse of minor-to-major axis ratio $q$, the angle between its radius and quadratic normal is at most $\pi/2-2\arctan q$. The tangent of the angle in principal coordinates is $(q^{-2}-1)|w|/(1+q^{-2}w^2)$ when $w$ is the slope of the radius. The maximum occurs at $|w|=q$. To turn the two centre directions into opposite normals, the smaller axis ratio must be at most $\tan(\delta/4)=0.006682$. Equality is attainable for that one peak. Orient the two ellipses to attain the limiting turns in opposite directions and scale the two quadratics to cancel their normals there. Subtract each quadratic from a sufficiently large constant and clip at zero. Near the target the sum is strictly concave with a local maximum at the target. The two-peak theorem applies even to those extreme shapes. [Independent limiting-shape calculation](test_peaks.py)

The western inequality also holds for any finite positive sum of elliptical profiles centred on either plasma clump, with each term assigned its own matrix. If every axis ratio exceeds the bound, all outward normals at the target lie in a common open semicircle. One short displacement then reduces every elliptical radius. With a combined peak and centre displacement allowance $e$ below both centre distances $d_W,d_E$, the angle can increase by at most $\arcsin(e/d_W)+\arcsin(e/d_E)$, by tangency to the two error disks. For allowances of 5, 10 and 20 arcsec the resulting upper bounds on the smallest axis ratio are 0.0521, 0.0984 and 0.1984. These allowances are sensitivity parameters. [Computed angle bounds](results.json)

A collar source that follows the galaxies or their overlap histories could have a displaced or irregular projected distribution. The static galaxy profile supplies no merger transport law. Identifying such a source requires predicting its positions from merger records independently of the lensing map. These data constrain the location of the gravitating correction. A microscopic identification would require a separate measurement. [JWST strong and weak lensing reconstruction](https://arxiv.org/abs/2503.21870)

From the repository root, replay the comparison and geometric controls with

~~~sh
PYTHONPATH=code python -m cosmology.rar_deep_regime.merger_peak.measure --check
PYTHONPATH=code python -m pytest -q code/cosmology/rar_deep_regime/merger_peak
~~~
