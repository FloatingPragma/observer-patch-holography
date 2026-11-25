# Reading curvature from the boundary of an equal-volume ball

A rest observer in a curved three-dimensional space counts markers inside a ball, then measures the area of its boundary with a second marker population. A separate triangle of mirrors measures curvature through round-trip light times. The two estimates of the rest-frame Einstein tensor agreed within 0.00063 curvature units at the smallest radius. Keeping the radius fixed gave about five thirds of the correct value. [Observer counts](measurements.json) · [Reproduction](measure.py) · [Experiment and thresholds](PROTOCOL.md)

## Geometry and instruments

The comparison supplies an ultrastatic spacetime with constant spatial sectional curvature $K$. On its time-symmetric rest slice the spatial scalar curvature is $6K$ and the Einstein contraction $G_{uu}$ is $3K$. Prepare curved charts with $K=-1,-1/2,0,1/2,1$ and a separately counted flat reference. Proper radial distance is $r$. The isotropic tangential radius is $S_K(r)$, equal to $\sin(\sqrt K r)/\sqrt K$ for positive curvature, $r$ at zero curvature and $\sinh(\sqrt{-K}r)/\sqrt{-K}$ for negative curvature. The proper volume is $4\pi\int_0^r S_K(s)^2\,ds$ and the boundary area is $4\pi S_K(r)^2$.

For each flat reference radius $\ell$, an observer stops the radial aperture upon capturing its prescribed number of interior marks. A fresh set of interior marks checks the volume, and separately sampled surface marks supply the area. The forward simulator samples the volume of the threshold mark from the gamma law for a homogeneous Poisson process. Repeated preparations give a simulated flat count budget ranging from $2^{28}$ at the coarsest radius to $2^{52}$ at the finest. A right-angle pair of equal radial mirrors independently measures the radial distance and the mirror-to-mirror geodesic using ideal clocks. In local proper units, if that latter length is $b$, then $3(2\ell^2-b^2)/\ell^4=K+O(\ell^2)$. The observer uses the clocks and signs of the area differences. The measured differences and curvature estimates can have either sign.

## Coefficient and refinement

At fixed flat-reference volume $4\pi\ell^3/3$, expanding the volume and solving for the geodesic radius gives $r/\ell=1+K\ell^2/15+2K^2\ell^4/175+O(\ell^6)$. The boundary then satisfies

$$
\frac{A}{4\pi\ell^2}=1-\frac{K\ell^2}{5}-\frac{3K^2\ell^4}{175}+O(\ell^6).
$$

The leading curved-minus-flat area is $-(4\pi\ell^4/15)G_{uu}$ with a geometric remainder of order $\ell^6$. At a common radial coordinate $r=\ell$ the corresponding area is instead $4\pi\ell^2[1-K\ell^2/3+O(\ell^4)]$. This change of experimental condition gives a factor of $5/3$ in the inferred contraction as the radius shrinks. These formulas use the time-symmetric rest slice. [Small-ball proof and its hypotheses](../PROOF.md)

| Flat radius | Largest difference between area and triangle estimates of $G_{uu}$ | Largest difference between equal-volume area result and $3K$ |
| ---: | ---: | ---: |
| 1/2 | 0.04121 | 0.06272 |
| 1/4 | 0.01110 | 0.01451 |
| 1/8 | 0.00285 | 0.00436 |
| 1/16 | 0.00063 | 0.00115 |

At the finest radius, fixed-radius readings on the four nonflat charts were 1.66579–1.66755 times $3K$. The independently sampled interior check on every matched chart and reference stayed within the six-standard-deviation threshold. In this simulation the marker budget scales as $\ell^{-8}$. Both the area uncertainty from a counting fluctuation and that from the volume stopping rule scale as $O(\ell^6)$ in standard deviation. A fixed count budget would eventually lose the curvature signal in the subtraction.

This experiment supplies curved geometry, light propagation, marker densities and proper clocks. It measures one contraction of the geometric area term. The entropy of the bulk and boundary, its area coefficient, generalized-entropy stationarity and a material source that produces the geometry require independent measurements. The registered Interval Echo carrier has its own clock and region measurements. [Regional material experiment](../../../model/experiments/interval_echo_einstein_station/README.md)

Run from the repository root.

Each enrolled seed drives a gamma stopping volume followed by independent Poisson reads. The deposited and replayed realizations have distinct counts at the largest budgets. Each complete set of clock ticks and integer counts has its own digest. Evaluating the same decoder and the original gates on each gives largest finest-grid area–radar discrepancies of 0.00062724 in the deposit and 0.00055380 on both replay platforms. The registered ceiling is 0.015. [Count verification and physical gates](measure.py)

```sh
python code/matched_einstein_limit/fixed_volume/measure.py
python -m pytest -q code/matched_einstein_limit/fixed_volume
```

## Proper volume from boundary observations

For an observer with timelike velocity \(u\), take the local geodesic rest slice \(\exp_p(u^\perp)\) and its induced spatial metric \(h\). Before the cut locus the ball of proper radius \(r\) has \(V(r)=\int_{d_h(p,x)\le r}\sqrt{\det h}\,d^3x=\int_0^r A(s)\,ds\). Thus measured areas of nested rims at calibrated proper separations determine the enclosed volume. The area of one rim generally does not. A local calibrated triad gives another prescription. The squared lengths of three independent short displacements and their three pair differences give their Gram matrix by polarization. The volume of their parallelepiped is the square root of its determinant. In the translation-invariant rank-three repair-record completion these rules give \(4\pi r^3/3\) exactly in its Gram units. Identifying the positions and length scale across physical patches as the source changes is a further requirement. [Rank-three theorem and position hypothesis](../../../paper/tex_fragments/DERIVATION_TECHNICAL_SUPPLEMENT_PORT.tex)

A single boundary can determine volume if it also measures paths through the interior. On a smooth nontrapping three-dimensional rest region with strictly convex boundary, launch unit-speed spatial geodesics at entry point \(x\) and inward direction \(v\). Record the length \(L(x,v)\) to first exit. Sample entries with density proportional to \(\cos\theta\,d\Omega\,dA\), where \(\theta\) is the angle with the inward normal. Then
\[
              V=\frac{A}{4}\,\overline L .
\]
Indeed the Liouville volume of the unit tangent bundle is \(4\pi V\). Geodesic flow preserves that volume. Following the last boundary entry backwards parametrizes it, up to a measure-zero set, by the incoming flux element times path length. Integrating \(\cos\theta\) over the inward hemisphere gives \(\pi\), so the total flux is \(\pi A\). This is Santaló's formula with constant integrand in three spatial dimensions. In an ultrastatic spacetime a reflected light ray measures twice this length in proper time with \(c=1\). Time-dependent geometry or a varying lapse requires an independent conversion of light records to lengths on the chosen rest slice. Missed exits and incorrect angular sampling change the estimator. [Paternain, Salo and Uhlmann, *Geometric inverse problems*](https://users.jyu.fi/~salomi/lecturenotes/GIP2D_driver.pdf)

In the five ultrastatic comparison charts above, a geodesic and the centre span a totally geodesic two-plane. Bisecting the chord gives a right triangle. Its adjacent half chord equals \(\arctan(\mu\tan(kr))/k\) for \(K=k^2>0\), \(\mu r\) for flat space and \(\operatorname{artanh}(\mu\tanh(kr))/k\) for \(K=-k^2\), with \(\mu=\cos\theta\). The clock emitter takes these formulas and the measured radial aperture. The decoder uses only the sixteen angle-resolved round-trip ticks, their cosine flux weights, the proper time unit and the separately counted rim area. A held interior count checks its predicted volume. Over all five charts, the largest absolute disagreement divided by the target volume is \(1.80\times10^{-4},9.41\times10^{-6},5.45\times10^{-7},4.26\times10^{-8}\) in order of refinement. These are simulated records on supplied geometries. In a flat ball a radial beam alone gives \(L=2r\) and overestimates volume by 50% if used as the flux mean. An unweighted mean over inward solid angle underestimates it by 25%. [Clock generation and decoder](chords.py) · [Matched boundary and interior measurements](chords.json)

At fixed volume \(V\), let \(\overline L_i\) be the mean proper chord in source preparation \(i\). Exact subtraction gives \(\Delta A=-4V\Delta\overline L/(\overline L_0\overline L_1)\). Consequently a positive change of collar charge with matched ordinary stress at positive coupling predicts a *longer* mean chord. In the four-dimensional Einstein normalization the leading prediction is \(\Delta\overline L=\kappa\Delta K/(6\pi^2\ell)+o(\ell^3)\) when \(V=4\pi\ell^3/3\) and \(\Delta K=O(\ell^4)\). It can be tested without using stress to infer a path length. The independently calibrated volume match must be sufficiently accurate. Near a smooth point \(dA/dV=2/\ell+O(\ell)\), so a volume-matching error \(o(\ell^5)\) is sufficient for an area error \(o(\ell^4)\). Absolute area errors \(o(\ell^4)\) and mean-chord errors \(o(\ell^3)\) suffice in the boundary formula. Shared calibration errors can instead be propagated on the actual contrast. A fixed fractional precision alone gives no such control. [Finite algebra](../../../Lean/ObserverPatchHolography/EinsteinBranch/OverlapCoupling.lean) · [Native boundary-exit clocks](../../../model/experiments/interval_echo_einstein_station/collar_scaling/README.md)

From the repository root, run `PYTHONPATH=code python -m matched_einstein_limit.fixed_volume.chords` to compare both boundary and interior readings with the deposited experiment.

## When a smaller area implies a smaller volume

For nested proper geodesic balls centered on the same observer before the cut locus, radial shells satisfy \(dV=A\,dr\). Label the shells by their enclosed volume \(v\). Then \(r(v)=\int_0^v dw/A(w)\). If measured areas are smaller at each matched volume, reaching the same volume requires a longer radial trip. Inverting the increasing radius functions gives less enclosed volume at the same proper radius. This comparison needs the nested area profile and an independent radial clock. The entropy measurement on a single boundary supplies neither.

For finite matched shells of volume increments \(w_i\), let \(a_i,b_i>0\) be the reference and perturbed shell areas, \(s_i\) their measured regional entropy contrasts, and \(e_i=s_i+\eta(b_i-a_i)\) the entropy-balance errors. Let \(c\) be the discrepancy between the independent radial-clock contrast and the two sums of shell thicknesses reconstructed from inverse areas. Then exactly
\[
  r_1-r_0=c+\sum_i\frac{w_i(s_i-e_i)}{\eta a_i b_i}.
\]
The finite proof bounds this radial shift from below when the instrument supplies bounds on the signed errors. Areas equal to volume increments divided by proper shell thickness give zero chart discrepancy by construction. Independently measured rim areas require testing that discrepancy. A differential stationarity result only controls the first variation, and finite state changes must include their remainders. [Kernel-checked finite shell identity](../../../Lean/ObserverPatchHolography/EinsteinBranch/NestedVolumes.lean)

For a smooth three-dimensional geodesic ball the small-radius expansions independently imply \(\delta V|_r=(r/3)\delta A|_V+o(r^5)\) for the leading scalar-curvature contrast. The interior counts at fixed radius and the surface counts stopped at matched volume test the two sides on separate preparations of the supplied curved charts. At radius 1/16 the measured ratios \(3\Delta V/(r\Delta A)\) range from 0.99921 to 1.00069 for the four nonzero curvatures. The flat control has two noise-level contrasts, so their ratio does not estimate a curvature coefficient. [Direct count decoder](nested.py) · [All detector counts](measurements.json)
