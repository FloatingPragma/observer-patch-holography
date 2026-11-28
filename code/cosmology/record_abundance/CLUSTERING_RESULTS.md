# How strongly can the mobile-double record field cluster?

Planting one ordinary double at the inherited equilibrium occupation of 0.16 produces an average of **2.04 additional mismatching seams**. Random native exchanges carry this change outward. A uniform response on large scales has 0.80952 times the fractional amplitude of the corresponding ordinary-occupation fluctuation. Assigning the cuts fixed mass in a fixed proper volume per cell imposes a stronger constraint: even the most favorable occupation arrangement cannot produce the measured densities of galaxy clusters. [Registered comparison](CLUSTERING.md) · [Full results](clustering.json)

## Same-cell and distinct-cell responses

Let $A_i$ be the sharp mobile-double occupation, independently Bernoulli with probability $p$ in the stationary preparation. This describes the source's actual independent initialization and its preservation by random permutations. A fixed-total equilibrium has finite-population corrections; the formulas here use the original distribution over totals. Let $D_{ij}=A_i+A_j-2A_iA_j$ indicate disagreement. For this test each cell receives half of its three incident seam readings, $S_i=\tfrac12\sum_{j\sim i}D_{ij}$ in seam units. Its mean is $3p(1-p)$. This local assignment defines the proposed density. The pre-existing transported halo profile uses its own mark routing on a comparison chart. [Source and stationary law](README.md) · [Original occupation instrument](../../../model/candidate_laws/interval_echo/observations.py)

Write $P_{ik}(t)$ for the probability that an initially tagged load at $k$ is at $i$ after $t$ services. It follows from the ordinary single-tag exchange law. The exact conditional difference between planting and removing the initial load, with all other initial occupations unchanged in distribution, is

$$
\mathbb E[S_i(t)\mid A_k(0)=1]-\mathbb E[S_i(t)\mid A_k(0)=0]
=\frac{1-2p}{2}\left(3P_{ik}(t)+\sum_{j\sim i}P_{jk}(t)\right).
$$

Conditional on a complete swap history, the two preparations differ at one terminal site. The three edges incident on it are the only ones whose status can change. A vacant neighbour gains a disagreement, and an occupied one loses it. Averaging the signed difference gives $1-2p$ per such edge. Multiplying the displayed formula by $p(1-p)$ gives the unequal-time covariance with the initial sharp occupancy. At release half the integrated response is at the original cell and the other half is shared by its three neighbours. At positive lag both move according to the same ordinary transfer.

For a uniform change of $p$, $d\ln\mathbb E S/d\ln p=(1-2p)/(1-p)=17/21$ at $p=0.16$. For an output eigenmode of the degree-three adjacency with eigenvalue $\lambda$, the fractional cross-response to an ordinary occupation perturbation acquires the factor $(1-2p)(3+\lambda)/[6(1-p)]$ relative to the ordinary propagated response on that mode. Any residual quadratic fluctuation has zero cross-covariance with the initial linear occupation. Distinct edges sharing an endpoint have equal-time covariance $p(1-p)(1-2p)^2$, disjoint edges have exactly zero, and an edge paired with itself has variance $2p(1-p)[1-2p(1-p)]$. Thus the same-cell and shared-endpoint contributions are fully accounted for. These identities use independent Bernoulli inputs and permutations; Gaussian initial conditions, a sound speed and collisionless evolution are separate physical claims. [Exact small-graph checks](test_record_clustering.py)

At service eleven on 320 cells the held mean response of 8,192 matched simulated preparations was:

| Distance from the initially changed cell | Model response in seams | Held response ± one standard error |
| --- | ---: | ---: |
| 0 | 0.1275 | 0.1259 ± 0.0039 |
| 1 | 0.3267 | 0.3297 ± 0.0060 |
| 2–3 | 0.9895 | 0.9800 ± 0.0112 |
| 4 and above | 0.5963 | 0.6021 ± 0.0109 |

Across all distances the forecast at occupancies 0.16, 0.50 and 0.84 is respectively +2.04, 0 and −2.04 seams. The measurements at the same appointment are +2.0376 ± 0.0138, +0.0352 ± 0.0192 and −2.0498 ± 0.0139. Increasing the ordinary load in a nearly full board therefore reduces the charge in this assignment. Sharing exchange words between the planted and empty copies defines paired sampling. These counts measure the conditional finite material source under its own clock. [Dynamics and all appointments](clustering.py)

## Fixed cell volume and the cluster data

Every cell obeys $0\le A_i\le1$ and $0\le S_i\le3/2$ pointwise, even far from equilibrium. Relative to the specified cosmic background at $p=0.16$, the largest ordinary density ratio is $1/p=6.25$, and the largest cut-density ratio is $1/[2p(1-p)]=3.72024$. Positive conversions of the two signals into mass cancel from each ratio. Their sum, if it supplies all gravitating matter, is bounded by 6.25 times its cosmic mean for every choice of their positive mass factors. Spatial averaging with positive volume weights preserves both caps. These statements require fixed source capacity and common physical cell density in the cluster and background. They also require the regional mass to be the sum of these local readings. The most favorable three-seam cut already allows more clustering than locally independent occupations at any $p$.

Lemze et al. combined strong and weak lensing and X-ray data and tabulated the enclosed three-dimensional mass of Abell 1689 at six fixed radii. At the largest radius, $1.554\,h^{-1}$ Mpc, they find $(9.56\pm1.80)\times10^{14}\,h^{-1}M_\odot$. At redshift 0.183 use the conservative upper reference $\rho_{\rm crit,0}(1+z)^3$ for the cosmic matter mean, valid for $\Omega_{m,0}\le1$ with conserved nonrelativistic matter. The $h^2$ in the measured mass divided by proper sphere volume cancels the $h^2$ in the reference. Even at the quoted lower mass the ratio is **107.43**, over **17 times the combined model cap**. All five smaller apertures exceed it more strongly. The reconstruction uses spherical deprojection; the authors discuss line-of-sight elongation and departures from hydrostatic equilibrium. The mass uncertainties at nested radii are correlated. This comparison uses the authors' distances and redshift and requires no numerical value of $\Omega_{m,0}$ within the stated bound. [Lemze et al., Table 3 and §6](https://doi.org/10.1111/j.1365-2966.2008.13116.x)

An independent check uses the dark component alone. Gonzalez et al. inventoried hot gas and stars, including diffuse intracluster stars. For systems above $2\times10^{14}M_\odot$ their weighted mean baryon fraction within $r_{500}$ is $0.136\pm0.005$ in their reference cosmology. At the upper fraction the inferred mean nonbaryonic density across the sample is **429.5** times each cluster’s critical-density reference, using the definition of $r_{500}$. In a flat cosmology with positive components the cosmic mean of the putative cut matter is at most $\rho_{\rm crit}(z)$. At the inherited cap it can therefore supply **at most 0.87% of the sample's nonbaryonic mass fraction at these apertures**, evaluated in the same cluster weights. Gas/stellar conversions and the total-mass calibration enter the measured fraction; this comparison uses the quoted uncertainty and does not assign it a formal rejection significance. The fixed-radius lensing test above avoids using $r_{500}$ to define the overdensity being tested. [Gonzalez et al., §§V.3.3–V.4](https://arxiv.org/abs/1309.3565)

These observations exclude the fixed-volume, equal-mass, local-cut assignment as the dominant source of cluster mass. A different physical density of observer cells, long-lived transported charges arriving from outside the region, changing record gain, or gravitating material beyond these components changes the prediction. The conditional finite exchange dynamics supplies no ruler, stress or cosmic perturbation transfer that determines those possibilities. [Exact finite field and observational arithmetic](clustering.json)

~~~sh
PYTHONPATH=code:. python -Bm cosmology.record_abundance.clustering
PYTHONPATH=code:. python -Bm cosmology.record_abundance.clustering --fetch
PYTHONPATH=code:. python -Bm pytest -q code/cosmology/record_abundance
~~~
