# A shared clock, sharp optical focusing and material flux

At the volume-matched stationary head of the 80-cell radar, no common clock and scalar conversion from extensive energy current to flux density satisfy all three signed sharp-screen focusing equations. The three images use a vacancy probe, an excess probe and a fair mixture whose raw arrivals are pooled before the image is reconstructed. The result permits arbitrary smooth dependence of the clock and conversion on the complete source state. [Interval certificate](certificate.json) · [Verifier](measure.py)

## A test at a single event

Fix the common native exchange dwell $t=4$ and the launch ensemble of ten emitters about the bridge with ports 1 and 38. Write $J$ for the centered two-dimensional Jacobi image matrix, $G=-\operatorname{tr}(J_{tt}J^{-1})$ and $\theta=\operatorname{tr}(J_tJ^{-1})$. For a smooth common clock $\lambda(t,s)$ on source preparation $s$, with nonzero $\lambda_t$, ordinary differentiation gives
$$
-\operatorname{tr}(J_{\lambda\lambda}J^{-1})=
aG+b\theta,\qquad a=\lambda_t^{-2}>0,\quad b=a\lambda_{tt}/\lambda_t.
$$
Both clock jets may depend on the whole preparation. At fixed $s$ they are the same numbers for every image at this event.

Let $g_i=G_{i,\mathrm f}-G_{i,\mathrm r}$, $e_i=\theta_{i,\mathrm f}-\theta_{i,\mathrm r}$ and $j_i$ be the independently metered extensive cell energy current, for $i=-,+,p$. A common conversion to signed, equally normalized null-stress differences requires $ag_i+be_i=h j_i$. The free scalar $h$ absorbs the assigned flux area, orientation sign and coupling. This is a necessary condition even if the two tagged preparations have different densities and longitudinal pressures, since those terms cancel in the opposite-null difference. [Bridge energy continuity, phase meter and launch definition](../README.md#the-bridge-energy-current-on-the-same-untagged-source)

| Image | $e_i$ | $j_i$ | $g_i$ |
|---|---:|---:|---:|
| Vacancy | [0.025015, 0.025453] | [0.0094789, 0.0095275] | [0.51200, 0.51712] |
| Excess | [0.006120, 0.006242] | [−0.031258, −0.031114] | [0.17652, 0.17889] |
| Fair mixed | [0.014841, 0.015267] | [−0.010890, −0.010793] | [0.07219, 0.07729] |

The determinant with rows $(e_i,j_i,g_i)$ lies in $[0.00020771,0.00022816]$. The three homogeneous equations have only $a=b=h=0$, contradicting $a>0$. The measured discrepancy also has a positive lower bound in exchange-dwell units. For every real $u,k$,
$$
\max_i |g_i+u e_i-kj_i|\ \ge\
\frac{|\det[e,j,g]|}{\|e\times j\|_1}>0.1229 .
$$
This follows by pairing the residual with $e\times j$. The residual is divided by $a$ here, so an arbitrary change of clock units cannot make this discrepancy arbitrarily small. [Interval replay](test_clock_area.py)

## Signed variations with an unfixed baseline

Allow the three odd residuals to have arbitrary source-independent offsets. Alone, their differences in preparation give three constraints on three changing coefficients. Indeed, the nonzero determinant above implies that those three odd residuals can locally be held at arbitrary constants by solving for $(a,b,h)$ near an initial choice with $a>0$. Thus the absolute test alone does not rule out a response-only model with unconstrained baselines.

One additional requirement comes from the known fair matter preparation. At the same event and normalized null direction, its material stress expectation is the mean of the two tagged expectations whenever stress is an affine source observable and the flux area is common to the three readings. This requirement permits distinct vacancy and excess stresses, both of which can vary freely. For a forward quantity write $m(Z)=Z_{p,\mathrm f}-(Z_{-,\mathrm f}+Z_{+,\mathrm f})/2$. In addition to the three odd residuals, $a\,m(G)+b\,m(\theta)$ must then have zero variation. This condition can also have an unfixed constant background.

Form four-vectors $g=(g_-,g_+,g_p,m(G))$, $e=(e_-,e_+,e_p,m(\theta))$, and $j=(j_-,j_+,j_p,0)$. For each exact preparation chord $D$ the material controls fix all three energy currents, so $\partial_D j=0$. Differentiating the four residuals at the reference gives
$$
g_D+u e_D+v_D g+w_D e-k_Dj=0,\qquad
u=b/a,\quad(v_D,w_D,k_D)=(\partial_Da,\partial_Db,\partial_Dh)/a .
$$
Here $v_D,w_D,k_D$ may differ without constraint between chords. Only $u$ must agree at the same reference. Taking determinants gives the necessary value
$$
u=-\det[g,e,j,g_D]/\det[g,e,j,e_D].
$$
The equal-work chord requires $u\in[-10.269,-2.873]$. A chord fixing both works and varying the even $I\otimes X/2$ component requires $u\in[-286.520,-93.370]$. A third chord fixing both works and varying $I\otimes Y/2$ requires $u\in[-70.105,-41.099]$. Thus arbitrary constant residuals and source-dependent changes of all three coefficients do not recover the full signed equations under affine stress pooling. [Material preparation constraints and positivity](../energy_odd.py) · [Required clock intervals](certificate.json)

A stronger tracer interpretation requires every image in a given direction to measure the same curvature, so its two tracer differences must vanish. Absolute focusing in the forward direction would require $u\in[-35.032,-34.591]$ from vacancy versus excess and $[-0.630,-0.279]$ from vacancy versus the actual fair-pooled image. There is also an independent response-only test. On the equal-work chord in the original optical certificate, form four rows from the two tracer differences for each of the two directions. The columns are $\partial_D\theta_i-\partial_D\theta_j,G_i-G_j,\theta_i-\theta_j,\partial_DG_i-\partial_DG_j$. Their determinant exceeds $5.97\times10^6$. Thus allowing arbitrary changing clock coefficients and independent target variations for the two directions leaves these four equations inconsistent. [Baseline and source-slope enclosures](../optical_certificate.json)

## Interval reconstruction

The tagged and pooled image moments and their first two time derivatives come from the native integer propagator. Each contracted moment tensor is affine on the conditional even and odd head blocks. The verifier encloses the derivatives of the centered and normalized image, including both moving centroids, in the reference uncertainty ball. It divides the allowed parity interval $[0.1340,0.1346]$ into ten pieces. On each piece it evaluates the forward-minus-reverse quantity at the common center. For its range it bounds the **difference** of the two gradients, using the Hilbert–Schmidt radius in both conditional states and the width in parity weight. This uses the same state in both directions. Directly subtracting two broad enclosures of focusing loses that correlation. The control is affine on unnormalized traceless blocks; dividing the conditional-state gradients by their block probabilities converts to these coordinates. The code includes the independently bounded correction from the rational tangent proposal to the exact native preparation chord. [Certified stationary head](../../../full_algebra/reuse/INVISIBLE.md) · [Native screen moments](../optical_certificate.py)

These conclusions test a specific event, material flux, pair of opposite directions and observed screen. A source-dependent shift of which physical event is compared, tracer-specific affine clocks or a different transported transverse screen changes those identifications. The calculations do not decide the general continuum null-balance premise.

Run from the repository root.

    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.directional_stress.regional_modular.clock_area.measure
    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/directional_stress/regional_modular/clock_area
