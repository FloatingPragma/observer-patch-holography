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
