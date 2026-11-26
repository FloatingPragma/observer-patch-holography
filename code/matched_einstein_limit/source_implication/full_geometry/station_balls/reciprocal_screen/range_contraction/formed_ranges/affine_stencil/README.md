# A five-station spatial second difference on the formed collar

Five observers can remove a constant reading and every common linear spatial trend from a material detector signal. The observer at port 0 uses radar distances to place four neighboring stations in a comparison chart. All four averaging weights are positive. Applying those weights to independent local return counts gives a **positive second difference** when correlations form between the two parked heads. Vacancies give **0.0034133 to 0.0034152** per squared radar unit, and excesses give **0.0036798 to 0.0036818**, at simultaneous 99% coverage. [Native interval bounds and all terminal counts](results.json)

## Where the five stations are placed

The source is the registered four-unit free encounter on the 80-cell board, followed by removal of the forming excess without reading its final position. Its two marked heads have joint state $\rho$. The comparison stock $\tau$ takes its left and right heads from independent copies of the same encounter, preserving the marginals and parked energy. Use the original four stations and the first registered additional receiver, in the order $(0,1,38,39,5)$. Every directed radar experiment retains its misses at the six-unit deadline. Equal-probability vacancy and excess probes are pooled **before** fitting the five clock gains. Neither embedding coordinates nor the circuit layout enter the distance decoder. [Certified directed records for both stocks and species](../results.json)

Write $d_{ij}$ for the gain-corrected reciprocal distances and $z$ for the inverse radius of a spherical comparison space. In the registered interval $0.620<z<0.623$, the five-by-five cosine Gram $C_{ij}=\cos(zd_{ij})$ has one determinant zero. Outward intervals bracket the root, verify a strictly negative derivative throughout the interval and positive leading minors through order four. All arcs have $zd_{ij}<\pi$. This gives one common three-dimensional spherical realization of the five distances on this branch, unique up to an isometry. The fitted $z$ is 0.6216127666 on the formed pooled stock. The same checks succeed for the product stock and when the two probe species are used separately. [Branch and interval calculation](chart.py)

At the central receiver let $v_i$ be the normal-coordinate displacement to each neighbor along its shortest arc. Their tangent Gram is recovered directly from the ranges,

$$G_{ij}=d_{0i}d_{0j}
\frac{C_{ij}-C_{0i}C_{0j}}{\sin(zd_{0i})\sin(zd_{0j})}.$$

The first three tangent vectors are independent. The last has their coordinates $a=G_{1:3,1:3}^{-1}G_{1:3,4}$. Normalizing $(-a_1,-a_2,-a_3,1)$ to sum to one gives the unique weights with $\sum_iw_i v_i=0$. In the formed pooled chart they are 0.26460, 0.23145, 0.23500 and 0.26895. The mean squared displacement $m_2=\sum w_i d_{0i}^2$ is 8.51788018 in squared radar units. For any scalar field sampled at these stations the measured combination

$$D(f)=\frac{2}{m_2}\left(\sum_iw_i f_i-f_0\right)$$

vanishes on constant and linear fields. On the quadratic $f(v)=|v|^2/2$ it gives exactly one. More generally its value on $v^THv/2$ is $\operatorname{tr}(QH)$, where $Q=\sum w_i v_iv_i^T/m_2$. The certified inequalities $0.29I<Q<0.40I$ bound the angular sensitivity in every direction. For a smooth field with third directional derivative bounded by $B$ along these rays, the error against $\operatorname{tr}(Q\,\mathrm{Hess}\,f(0))$ is at most $B\sum w_i d_{0i}^3/(3m_2)$. This finite experiment supplies neither such a bound for physical detector fields nor a vanishing-radius family.

## Independent local return and energy controls

At each of the five ports an observer prepares one sharp vacancy or excess. Following four units of unmonitored evolution, it records whether the probe is back at that same port. There are no intermediate relays, coherent rail preparations or receiver-side phase references. Each bank retains both outcomes. The five-place field is the change in this return probability from $\tau$ to $\rho$, using the fixed geometry weights of the formed pooled chart for both species and all stocks. An independent ordinary standard dephases the left reversal of $\tau$ and raises the parked energy by half a native exchange unit. Thirty fresh binomial banks each use $2^{42}$ simulated starts. All native probabilities are independently enclosed by four-tick integer propagators and formed-state bounds; centered ticket rounding and Hoeffding's inequality give joint 99% count intervals. The scheduled four-unit shutter is part of the instrument. [Local detector calculation](returns.py)

| Readout | Vacancy | Excess |
|---|---:|---:|
| Pairing-induced weighted change before normalization | 0.0145371 to 0.0145451 | 0.0156722 to 0.0156803 |
| Ordinary weighted change per exchange unit | 0.0050957 to 0.0051118 | 0.0964039 to 0.0964199 |
| Equivalent ordinary exchange units | 2.8438 to 2.8544 | 0.16254 to 0.16266 |

Both detectors resolve the sign of the five-point signal. Calibrating each detector on the same independent ordinary half-unit standard gives disjoint inferred energies. The field being sampled is a local return probability; assigning these coherent returns a universal coupling to regional stress would require further evidence.

## What changes in the registered positions?

Freeze the reference weights learned on the product stock and apply them to the formed chart's positions. Their weighted displacement has length **0.0019764540 to 0.0019764591** radar units. Every exact common affine map would have preserved the zero displacement. An affine map fitted to port 0 and the three independent neighbors misses the fifth position by **0.0073433984 to 0.0073434172** radar units. The fifth receiver is held out of this affine-map fit. All five receivers enter the spherical-radius fit for each stock. Common translations, rotations, anisotropic strains and changes of ruler scale do not alter the zero criterion. The response is also nonzero in the independently fitted vacancy-only and excess-only charts. These results concern the measured five-station spherical branch. Extending the chart to additional stations and identifying proper distance and spacetime curvature are independent physical questions. [Six-station and alternative-curvature tests](../../../../../README.md)

From the repository root run

```sh
export PYTHONPATH=code:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=matched_einstein_limit.source_implication.full_geometry.station_balls.reciprocal_screen.range_contraction.formed_ranges.affine_stencil
nice -n 15 python -m "$p.returns"
nice -n 15 python -m pytest -q --pyargs "$p"
```
