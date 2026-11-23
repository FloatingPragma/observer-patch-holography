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

```sh
python code/matched_einstein_limit/fixed_volume/measure.py
python -m pytest -q code/matched_einstein_limit/fixed_volume
```
