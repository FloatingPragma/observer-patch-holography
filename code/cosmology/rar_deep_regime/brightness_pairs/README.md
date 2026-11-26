# Stellar surface brightness and the orbital excess at matched baryonic force

Does the fainter part of a galaxy consistently require more apparent extra gravity? In the SPARC rings selected below, the answer depends on the galaxy. At low baryonic accelerations, 70 galaxies supplied two annuli differing by at least a factor of two in local 3.6 µm stellar surface brightness, with baryonic accelerations matched within 10%. At twice the individual quoted speed errors, two fainter annuli had a larger relative acceleration than every constant-amplitude square-root model permits, eight had a smaller one, and 60 overlapped its prediction. All 70 fainter annuli were farther out in the same galaxy. Their common telescope distance cancels from the comparison. A smaller visible surface brightness therefore cannot by itself enforce the proposed monotone increase of the square-root amplitude on these nominal mass models. [Frozen selection](PROTOCOL.md) · [All enrolled rings and error boxes](results.json)

## What is being observed?

In OPH, local record density $\rho$ counts **records made and retained by the source system per physical volume**. At fixed physical depth $\delta$ the dense-recovery premise bounds the defect by $\kappa B\exp[-(\delta\rho^{1/3}-r_0)/\zeta]$. On the sparse erasure model, the retention probability $p$ yields defect $(1-p)/2$. The bound permits different actual stresses at the same density; it is an upper bound and does not imply that every sparse patch gravitates more. Establishing a galaxy-level law requires an independent measurement of local retention, a calibration between recovery failure and gravitating stress, and the ordinary gravity response on the matching physical region. The paper gives no measured mapping from a stellar luminosity surface density to $\rho$ or $p$. Nor do the deposited galaxy orbital/lensing records contain the physical observer retention counts. A repeat observation from another telescope can refine the orbital or lensing measurement but leaves the records formed in the galaxy intact. [Dark-sector paper, galactic collar and dense-recovery sections](../../../../cosmology/oph_dark_matter_paper.tex) · [PR-100, PR-135–137](../../../../docs/registers/OBSERVABLES.json)

Here the operational proxy is $\Sigma_{3.6}={\rm SBdisk}+{\rm SBbulge}$, the local stellar light per projected area. It tracks starlight, while the catalogue separately models the gas force. Changing the assumed distance scales disk area and luminosity by the same factor in the nearby-galaxy catalogue convention. Measuring the luminosity field is distinct from acquiring neighbors' retained physical records. The authors of SPARC have already tested acceleration-relation residuals against galaxy properties. This additional retrospective comparison holds the galaxy fixed, matches two nominal baryonic accelerations, and uses a bound valid for **every** square-root acceleration scale. It supplies neither a second galaxy population nor a new survey. [SPARC catalogue and original mass-model paper](https://doi.org/10.3847/0004-6256/152/6/157) · [Primary residual analysis](https://arxiv.org/abs/1610.08981)

## Selection and falsifiable sign

SPARC's parent sample has 175 galaxies. Quality 1 or 2 and inclination at least 30° leave 153. For each acceleration cut the protocol picks the greatest positive brightness ratio from qualifying pairs using only the photometric/baryonic columns. Each galaxy then contributes at most one contrast. Speeds and speed uncertainties are first opened after both selections. The catalogue was used in earlier project studies, so this is a retrospective, velocity-blind selection, not a prospective new observation. The SPARC gas contribution retains its sign and its published helium correction. Disk and bulge mass-to-light ratios are fixed at 0.5 and 0.7. Published near-zero positive surface brightness is accepted; rounding can influence the largest contrast. Both selected annuli must individually satisfy the cut. [Data-column definitions](../data/ReadMe) · [Exact-rational selection and comparisons](measure.py)

Write $t=g_{b,f}/g_{b,h}$ for baryonic acceleration in the faint and bright annuli. At angular radii $\theta_f,\theta_h$ and observed rotation speeds $v_f,v_h$, the orbital acceleration ratio is
$$
 R=(v_f/v_h)^2\,\theta_h/\theta_f.
$$
A common inclination correction cancels from the speed ratio. Distance cancels from the radius ratio. At fixed angular mass and light profiles, $v_b^2$ and physical radius both scale linearly with the assigned distance, so the catalogue baryonic accelerations entering the cuts also stay fixed under this common distance rescaling. A local warp, pressure support or a radius-dependent change to the mass model can change the comparison.

For $g_{\rm orb}=g_b+\sqrt{a_0g_b}$ let $u=\sqrt{a_0/g_{b,h}}$. The predicted ratio is
$$
 R_0=\frac{t+u\sqrt t}{1+u},\qquad
 \min(t,\sqrt t)\le R_0\le\max(t,\sqrt t).
$$
This interval covers any $a_0\ge0$ and uses no fit. For the stronger brightness hypothesis allow a local multiplier of the excess, $g_{\rm orb}=g_b+A(\Sigma_{3.6})\sqrt{g_b}$, with $A_f\ge A_h\ge0$ when the annulus is fainter. Replacing the faint amplitude by the bright one can only decrease the predicted ratio. Hence $R\ge\min(t,\sqrt t)$ is a necessary condition even when each galaxy has its own normalization. An upper measured ratio below that threshold is a conditional counterexample to the monotone proxy. An excess above the upper shape bound is compatible with such a proxy but also with other changes to the orbital model. A ratio inside the interval does not demonstrate a shared global amplitude.

| Sample | Speed-error box | Larger in faint annulus | Overlaps shape interval | Smaller in faint annulus |
|---|---:|---:|---:|---:|
| Both $g_b<1.2\times10^{-11}$ m/s² | Catalogue central values | 38 | 8 | 24 |
| Same 70 galaxies | ±1 quoted error | 8 | 50 | 12 |
| Same 70 galaxies | ±2 quoted errors | 2 | 60 | 8 |
| Both $g_b<1.2\times10^{-10}$ m/s² | Catalogue central values | 52 | 7 | 38 |
| Same 97 galaxies | ±1 quoted error | 16 | 62 | 19 |
| Same 97 galaxies | ±2 quoted errors | 4 | 82 | 11 |

There were 83 eligible-quality/inclination galaxies with no qualifying deep pair, and 56 with no qualifying wider pair. All 97 fainter rings in the wider sample are at greater galactocentric radius; 65 of the wide pairs are literally the same pairs as in the deep sample. In each error box the two reported speeds vary independently (positive lower bound, or zero), and the most favorable endpoints of their ratio are compared to the full model interval. Exactly rational comparisons avoid rounding near a sign boundary. The boxes carry no joint or population confidence level. Adjacent radii share instrumental and galaxy systematics. Inclination warps, shared gas/stellar mass modeling, noncircular motion, the use of a local disk-force relation away from a baryon-enclosing exterior, and photometric profile uncertainties remain untested. All three outcomes are reported; no independence test is assigned to the wider control.

Two examples have nearly the same catalogue distance. In D631-7 (7.72 Mpc), stellar brightness falls from 26.06 to 2.10 solar luminosities per square parsec. The entire two-error observed ratio exceeds 1.44, while the constant-amplitude shape ratio is at most 0.994. In UGC05764 (7.47 Mpc) brightness falls from 11.43 to 0.44 in the same units. Its observed ratio is below 0.783 throughout the two-error box, while the minimum ratio for *any* nondecreasing excess amplitude towards the fainter ring exceeds 1.008. These illustrative galaxies have distinct gas profiles and are members of the full 70-galaxy count. They are each seen from one distance. No record follows the same galaxy moved to a second distance. This test cannot separate low stellar light from large galactocentric radius in this enrollment. [Measured source rows](../data/table2.dat)

## How much can stellar calibration change the eight lower ratios?

Hold the eight original pairs and the two-error speed boxes fixed. Allow separate common disk and bulge mass-to-light ratios in $[0.4,0.6]$ and $[0.56,0.84]$ for each galaxy, with its signed gas contribution fixed. These ranges give a sensitivity test with no assigned probability. For a given annulus, the baryonic force is affine in the two mass-to-light ratios. Both annuli remain positive at all four corners. At an intermediate point, write the two affine forces as convex combinations of the corner values with the same weights. Their ratio is a weighted mean of the corner ratios, with weights proportional to the bright-ring forces. The extrema therefore occur at corners. The deep cut and both force-matching inequalities are also linear constraints on the forces. [Sensitivity protocol](PROTOCOL.md) · [Exact-rational calculation](sensitivity.py)

Five pairs remain below the monotone bound throughout this stellar range. Independently, five pairs have a two-error orbital upper endpoint below $10/11$, the smallest possible lower shape bound for *any* positive forces satisfying the original 10% match. Three pairs belong to both groups. Thus seven conflicts persist at every allowed stellar calibration that also keeps their two forces matched. In UGC07125 a disk ratio of 0.6 permits the lower shape inequality while retaining the original cuts. Reaching that inequality establishes a necessary condition only. UGC06667 alone stays deep and matched and remains discrepant throughout the entire stellar range. This calculation keeps the original annuli and does not ask whether another pair would become the largest brightness contrast at a changed calibration. [Results for all eight galaxies](sensitivity.json)

Even without a stellar model one can bound the force change needed. Let $H$ be the largest orbital ratio permitted by the two speed boxes. The lower shape bound can be reached only if the true force ratio is at most $q=\max(H,H^2)$. Independent relative changes of magnitude at most $e$ in both nominal baryonic accelerations can lower their ratio from $t$ to $t(1-e)/(1+e)$. Consequently the required tolerance is at least $(t-q)/(t+q)$, ranging from 0.807% in UGC12732 to 13.061% in UGC02259. This is the most favourable opposing adjustment at the two radii. Published speed errors alone supply no error budget for the inferred baryonic forces.

## Replay

From the repository root, run:

```sh
python code/cosmology/rar_deep_regime/brightness_pairs/measure.py
python code/cosmology/rar_deep_regime/brightness_pairs/sensitivity.py
python -m pytest -q code/cosmology/rar_deep_regime/brightness_pairs/
```

The full receipt includes both catalogue file hashes, original table line numbers, photometric contrasts, signed gas fractions, speeds, error intervals and all selected outcomes. The rule was frozen in commit df3949ddb before reading the selected velocities.
