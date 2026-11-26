# Stress conservation on the two sightlines through J0946+1006

On the spherical Chabrier force continuation, the two Einstein rings require the outer line of sight to have the larger radial-pressure column. In static attractive gravity, stress conservation then requires a positive weighted mean of tangential minus radial pressure. For the aligned optical interpretation the mean must reach **0.007570 solar masses per cubic parsec** after dividing pressure by the speed of light squared. A spherical source with isotropic total pressure fails this test even if its pressure is negative, provided its radial energy density plus pressure is nonnegative. The result applies to the total stress in this stellar and optical model. [Measured aperture contrasts and assumptions](../README.md) · [Numerical bounds](results.json)

## Which measurements enter

Gavazzi et al. measured two ring radii of $1.43\pm0.01$ and $2.07\pm0.02$ arcsec in the same galaxy. In the adopted flat distance calibration their central values are 5.115 and 7.404 kpc. The archived orbital comparison uses its two stellar light components, a Chabrier stellar mass range and the square-root acceleration continuation with scale calibrated on SPARC. The stellar velocity dispersion supplies an aperture check. The orbital force at every radius in the lens is a model inference. Varying the stellar calibration matters: on the nominal heavier Salpeter branch a positive secondary bend of 0.10775 arcsec fits the aligned rings with zero radial pressure. [Discovery imaging](https://arxiv.org/abs/0801.1555) · [Two stellar components](https://arxiv.org/abs/1111.4215) · [Outer-source spectroscopy](https://doi.org/10.1093/mnras/stab1399) · [Orbital inputs and control](../results.json)

The aligned calculation cancels a common external surface-density sheet and allows a positive secondary deflection. A second calculation uses the two-source elliptical imaging fit in circular apertures compared with the same spherical orbital continuation. Interpreting that comparison as pressure requires the spherical gravitational model to reproduce the circular aperture contrasts. The two reductions use the same images. Both allow any nonnegative exterior orbital mass measure above the second ring under the stated interior force law. Published marginal-parameter boxes have no assigned joint confidence level. [Two-plane optical reconstruction](https://arxiv.org/abs/1403.5278) · [Aperture comparison from the elliptical fit](../elliptical.json)

## Why an outward-increasing pressure column needs tangential stress

Let $R_1<R_2$ be the projected radii, $r$ the spherical areal radius, and $P_i=\int_{-\infty}^{\infty}p_r(\sqrt{R_i^2+z^2})\,dz$. The spherical weak-field identity of the aperture experiment gives $(P_2-P_1)/c^2\ge L$ from the measured optical contrast and the orbital bound. At a single radius $P/c^2=\bar\Sigma_o-\bar\Sigma_w$ in that approximation. An allowance $E$ for omitted terms or an adverse systematic in the **difference**, in surface-density units, replaces $L$ by $L-E$. The positive bounds below apply when this remainder is positive. [Weak-field optical and orbital masses](https://arxiv.org/abs/astro-ph/0512213) · [Projected pressure derivation](../../aperture_virial/README.md)

For an absolutely continuous radial pressure with convergent columns and integrable weighted derivative, set

$$
W(r)=2\left[\sqrt{(r^2-R_1^2)_+}-\sqrt{(r^2-R_2^2)_+}\right].
$$

Integrate $p'_r$ between the two radii seen at each fixed $|z|$, then change order. The available length in $z$ at each $r$ is $W$, hence $P_2-P_1=\int Wp'_r\,dr$. This weight is positive for all $r>R_1$. Its integral obeys $\int W(r)\,dr/r=\pi(R_2-R_1)$. Indeed, express the square-root difference as an integral of $R/\sqrt{r^2-R^2}$ over $R$ from $R_1$ to $\min(r,R_2)$ and interchange the positive integrals. The inner integral is $\int_R^\infty R\,dr/[r\sqrt{r^2-R^2}]=\pi/2$.

In a static spherical metric with lapse $e^{\nu(r)}$, radial conservation of the diagonal stress tensor gives $p'_r=-(\epsilon+p_r)\nu'+2(p_t-p_r)/r$, where $p_t$ is either of the two equal tangential pressures and $\epsilon$ is energy density. Take $\nu'\ge0$ and $\epsilon+p_r\ge0$. Integrating this equation yields

$$
\frac{\int (W/r)(p_t-p_r)\,dr}{c^2\,\pi(R_2-R_1)}
\ \ge\ \frac{L}{2\pi(R_2-R_1)}.
$$

The shells contributing to this mean extend beyond the outer sightline. It need not attain the bound between the two rings, and neither individual line-integrated pressure needs to be negative. An extra component source or flux exchanging momentum with the collar must enter conservation of the **combined** stress. Subtracting measured ordinary stresses from the left side would be needed to attribute the excess specifically to a collar. [Static spherical stress conservation](https://doi.org/10.1140/epjc/s10052-024-13442-w)

There is also a lower bound on the total amount of positive stress difference. Maximizing $2\sqrt{r^2-R_1^2}/r^3$ gives $W/r^3\le4/(3\sqrt3 R_1^2)$, with the maximum at $r=\sqrt{3/2}R_1$ when that radius lies below $R_2$. With $[x]_+=\max(x,0)$, it follows that

$$
Y=4\pi\int_{R_1}^{\infty}r^2[(p_t-p_r)/c^2]_+\,dr
\ \ge\ \frac{3\sqrt3\,\pi}{2}R_1^2 L.
$$

$Y$ is the positive pressure difference integrated with areal radial measure, expressed in mass units. An independent model with spherical collisionless matter at speeds $|v|\le v_{\max}$ satisfies $2p_t\le(\epsilon/c^2)v_{\max}^2$ and $p_r\ge0$, by adding the two transverse kinetic stresses in the local static frame. Its total mass-energy integral $4\pi\int_{R_1}^\infty r^2\epsilon\,dr/c^2$ must then exceed or equal $2c^2Y/v_{\max}^2$. With an illustrative speed cap of 1,000 km/s the required mass is at least $4.314\times10^{15}\,M_\odot$ in the aligned quoted-error case. The rings do not measure that speed cap or an upper limit on the exterior halo mass. A non-kinetic effective collar has no such particle speed bound.

| Optical interpretation | Angular box | Minimum column difference ($M_\odot$/pc²) | Minimum weighted $(p_t-p_r)/c^2$ ($M_\odot$/pc³) | Minimum $Y$ ($M_\odot$) |
| --- | --- | ---: | ---: | ---: |
| Aligned | Quoted | 114 | 0.007570 | $2.400\times10^{10}$ |
| Aligned | Doubled | 76 | 0.004830 | $1.577\times10^{10}$ |
| Aligned | Extended inner | 10 | 0.000609 | $1.980\times10^{9}$ |
| Elliptical fit, circular apertures | Quoted | 149 | 0.009895 | $3.137\times10^{10}$ |
| Elliptical fit, circular apertures | Doubled or extended | 61 | 0.003877 | $1.266\times10^{10}$ |

## Can an arbitrarily massive exterior supply this stress using only slow particles?

There is a stronger test when the **entire** static, spherical energy and stress are kinetic. Write $\beta=\sup |v|/c$. Positivity of radial kinetic pressure gives $m(r)\le M_o(r)$ from the radial field equation $M_o=m+4\pi r^3p_r/c^2$ at the weak-field order of the optical experiment. Here $dm=4\pi r^2\epsilon\,dr/c^2$. Conservation and $2p_t\le\beta^2\epsilon$ then give
$$
 (P_2-P_1)/c^2\le {\beta^2\over4\pi}\int h\,dm,
 \qquad h=W/r^3 .
$$
The nonincreasing function $H$ equal to $H_*=4/(3\sqrt3 R_1^2)$ below $R_2$ and to $h$ above it dominates $h$. The outer expression decreases since it is proportional to $1/[r^3(\sqrt{r^2-R_1^2}+\sqrt{r^2-R_2^2})]$. Both mass measures are nonnegative by hypothesis. Therefore $\int h\,dm\le\int H\,dm\le\int H\,dM_o$. The last inequality compares cumulative masses on every interval starting at zero, then integrates over the level sets of $H$.

For $r\ge R_2$ write $c_i=\sqrt{1-R_i^2/r^2}$. The projected orbital mass contrast per unit spherical shell mass is $w_o=-(c_1-c_2)/[\pi r^2(1+c_1)(1+c_2)]$. On the same shell $H=2(c_1-c_2)/r^2$. Thus $H/(4\pi|w_o|)=(1+c_1)(1+c_2)/2\le2$. When $\beta\le1/\sqrt2$, the maximum allowed kinetic contribution from each exterior shell is no greater than the optical contrast that its orbital mass removes. This controls arbitrary nonnegative exteriors yielding finite aperture integrals.

More explicitly, the observed optical contrast is the orbital contrast plus $(P_2-P_1)/c^2$. Use its lower bound from the aligned arcs or the elliptical fit. Subtract the archived *upper bound on the truncated interior orbital contrast*, before adding the negative exterior term. The remainder is precisely the guaranteed margin $L$ in the table. Dropping the combined exterior term gives
$$
 L\le {\beta^2 H_* M_o(R_2)\over4\pi},\qquad
 \beta^2\ge {3\sqrt3\,\pi R_1^2 L\over M_o(R_2)} .
$$
To bound the denominator conservatively, put the upper limits of both Chabrier stellar components inside $R_2$ and use the largest training acceleration scale in $M_o=M_b+R_2\sqrt{M_b a_0/G}$. The full upper angular radius gives $M_o(R_2)\le4.071\times10^{11}\,M_\odot$ for the quoted errors. The resulting necessary velocity ceilings are **at least $0.343c$** for the aligned branch and **$0.392c$** for the elliptical branch. Doubling both radius errors gives $0.278c$ and $0.249c$. The widest inner range gives $0.098c$ and $0.249c$. All these bounds are below $1/\sqrt2$, so they hold trivially for higher ceilings as well. These values test whether slow collisionless matter can provide the pressure required by this *particular* force and optical model. They are not measured particle velocities. The projected inequalities do not require that any observed star moves at those speeds or that a non-kinetic effective collar has particles.

Every lower endpoint is rounded down after propagating the archived angular ranges and distance enclosure. Bounds use the maximum possible shell separation for the mean and the minimum inner radius for $Y$ in each box. The independent regression integrates two sightlines through a smooth conserved anisotropic profile with positive radial pressure and compares the actual column difference with the conservation prediction. [Bound calculation](measure.py) · [Independent projection and source checks](test_stress.py)

From the repository root, run

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code
python -m cosmology.rar_deep_regime.double_ring_pressure.stress_conservation.measure
python -m pytest -q code/cosmology/rar_deep_regime/double_ring_pressure/stress_conservation
~~~
