# Can each trail stock supply the outer rotation of individual galaxies?

At the catalogue values, **40 of 118** eligible bulgeless SPARC galaxies require more spherical gravitating mass than the 32-seam stock can supply. We give that stock a ratio of **5.386** times the entire observed galaxy baryon mass, above its paid yield of **5.134**, and put all that extra mass inside the measured orbit. Fifty percent more gas, stellar mass-to-light ratios up to 0.8 and simultaneous favorable shifts of two quoted errors reduce the count from 40 to **one, NGC 3741**. Its remaining shortfall depends on the pressure correction and on unmeasured gas. [Counts and rotation curves](orbits.pdf) · [Individual inputs, bounds and results](results.json)

The 14-seam stock fails for 90 catalogue centers and 12 galaxies in the two-error test. Four central measurements exceed the 64-seam stock, but none does after the one-error shifts and mass allowances. Each comparison is between one galaxy's measured rotation and its own baryon inventory. There are no weights from lens samples or stacks. These are upper-bound tests. Falling below a cap means that the stock can be large enough at that radius. A spatial profile is still needed to predict the speed.

## Sample and mass calculation

We started with all 175 galaxies in the deposited SPARC tables. Quality 1 or 2 and inclination at least 30° leave 153. Requiring zero bulge force at *every* published radius removes 31. Four others lack a usable outer point. Each of the remaining 118 has positive tabulated total stellar luminosity and HI mass. We select the largest radius with a positive measured speed and error, a relative speed error at most 10%, and positive nominal resolved baryonic force. For 100 galaxies that is the final tabulated radius. The force-budget inequality holds at every radius under the tested spherical model. We therefore use the full sample, and separately report the 97 galaxies with outer baryonic acceleration below the project's existing deep cut of \(1.2\times10^{-11}\,\mathrm{m/s^2}\). [SPARC catalogue and original rotation-curve references](https://astroweb.case.edu/SPARC/) · [Column definitions](../data/ReadMe) · [Frozen selection](PROTOCOL.md)

At stellar mass-to-light ratio \(\Upsilon\) and gas multiplier \(f_g\), the matched catalogue quantities are

\[
M_b=\Upsilon L_{3.6}+1.33f_g M_{\rm HI},\qquad
V_b^2=\Upsilon V_d|V_d|+f_g V_g|V_g|.
\]

The nominal values are 0.5 solar units and 1. The catalogue disk force assumes unit stellar M/L. Its gas force already includes helium and can have a negative signed contribution. The mass inventory uses the **total** luminosity and atomic gas, including any material beyond the radius. The force calculation uses the measured baryonic geometry at the chosen radius. SPARC generally omits molecular gas. A gas multiplier assumes that additional gas follows the released gas profile. Diffuse or molecular material with a different distribution would require its own inventory and force. [SPARC mass models, section 4](https://arxiv.org/abs/1606.09251)

For positive spherical stock of total mass at most \(qM_b\), the measured orbit must satisfy

\[
V_{\rm obs}^2\le V_b^2+Gq M_b/R,
\qquad
q_{\rm required}=\frac{(V_{\rm obs}^2-V_b^2)R}{GM_b}.
\]

The bound concentrates the entire extra stock inside \(R\). A spherical profile with some material outside would contribute less. The native trail actually receives total work \((34N-51)/4\) for \(N=14,32,64\), relative to an initial ordinary energy of 50.5. We give each stock the slightly larger ratio \(8.5N/50.5\). Applying it to a galaxy assumes the same gravitational mass-energy conversion, availability per *current* ordinary baryon and a spherical additional component. Baryon loss, a reservoir imported from outside the present galaxy, different service frequency, or substantial disk-like forces could change this test. [Native energy receipt](../../../../model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary/stock/trail.json) · [Previous cosmological stock comparison](../../cold_acoustic_polarization/TRAILS.md)

## How many galaxies exceed each bound?

| Allowed inputs | 14 seams, cap 2.356 | 32 seams, cap 5.386 | 64 seams, cap 10.772 |
|---|---:|---:|---:|
| Catalogue values | 90 | 40 | 4 |
| One-error box, M/L ≤0.8, gas ≤1.5× | 23 | 2 | 0 |
| Two-error box, M/L ≤0.8, gas ≤1.5× | 12 | 1 | 0 |
| Three-error box, M/L ≤1.0, gas ≤2× | 4 | 0 | 0 |

At the exact paid ratios of 2.104, 5.134 and 10.520, the central counts are **96, 43, 4** and the two-error counts are **14, 1, 0**. The four central failures even at the loose 64-seam cap are NGC 2915, NGC 3741, ESO 444-G084 and UGC 00128. They require ratios of 16.28, 15.74, 13.30 and 10.93, respectively. All four are in the deep subset. Within that subset the central counts for the three loose caps are 84, 40 and 4, and the two-error counts are 12, 1 and 0.

Distance changes scale radii linearly, flux-based masses quadratically and resolved squared speeds linearly. At fixed projected velocity the inferred speed scales as \(\sin i_{\rm cat}/\sin i\). We raise distance and inclination up to the stated catalogue errors, cap inclination at 90°, and lower the observed speed by the stated number of ring errors. We use the quoted luminosity errors at the same time. The gas multiplier extends down to 0.7 and M/L down to 0.3. We search all normalization endpoints to account for signed component forces. At a fixed radius and distance the ratio is a linear-fractional function with positive mass denominator, so its extrema over that rectangle occur at corners. These boxes retain the shapes of the published gas and stellar profiles. They are deliberately favorable sensitivity tests, with no probability assigned to their joint corners.

The data had appeared in earlier project studies. The native trail ratio and the present cut were set without fitting these galaxy measurements. The existing deterministic five-fold split labels are recorded for every galaxy. Previously designated fold zero contains 24 eligible galaxies, with central exceedances of 16, 6 and 1, and two-error exceedances of 1, 0 and 0. NGC 3741 belongs to fold one and had explicitly been studied before. This comparison tests a new prediction on public data. [Earlier five-fold split](../held_galaxies.py) · [Earlier NGC 3741 aperture comparison](../../record_abundance/MASS_RATIO_RESULTS.md)

## NGC 3741 sets the two-error limit

SPARC reports distance 3.21 ± 0.17 Mpc from the red-giant-branch tip, inclination 70° ± 4°, luminosity \((2.8\pm0.1)\times10^7L_\odot\), and total HI mass \(1.82\times10^8M_\odot\). The stellar-plus-helium-corrected-gas mass is \(2.561\times10^8M_\odot\) at catalogue center. The outer measured ring at 7.00 kpc has \(V=51.60\pm3.24\) km/s, \(V_g=13.34\) km/s and unit-M/L \(V_d=4.14\) km/s. The maximum total speed from the middle stock is **32.15 km/s** at these values.

With distance increased to 3.55 Mpc, inclination to 78°, stellar M/L to 0.8, luminosity raised two errors and 1.5 times the helium-corrected gas, the allowed baryon mass reaches \(4.734\times10^8M_\odot\). The measured speed floor becomes **43.35 km/s** at an adjusted radius of 7.74 kpc. The middle stock supplies at most **41.57 km/s**, requiring \(q\ge5.96\). The previous ring at catalogue radius 6.54 kpc requires 5.48 under the same box. These neighboring points share observational systematics.

Either of two specific changes closes the remaining outer-ring deficit. At the favorable two-error settings the required gas multiplier is **1.639**, equivalent to another \(4.12\times10^7M_\odot\) above the 1.5× gas allowance if it follows the same disk profile and supplies its share of additional stock. Alternatively, with the mass normalizations and angles held at those settings, the galaxy would need to lie at **3.86 Mpc**. For an additional centrally concentrated spherical seed and its proportionate stock the corresponding required seed is \(4.25\times10^7M_\odot\). The SPARC tables contain no direct upper bound on all diffuse/ionized gas.

The original HI data-cube study traces a symmetric warp and models radial flows of 5–13 km/s. The authors applied positive pressure-support corrections of **2.7 and 4.8 km/s** at the last two radii, and included half the corrections in their velocity uncertainties. Removing the full corrections as an additional favorable sensitivity gives required ratios of **4.71 and 4.52** in the two-error box. Both are below 5.386. With the one-error box they remain above it at **6.04 and 6.07**. The radial flows are distinct from this pressure correction and were already modeled. The data challenge a universal 32-seam budget at the corrected values, while the small residual in the generous box depends on how those outer velocities and gas inventories are treated. [Gentile et al., sections 3.5–3.6](https://arxiv.org/abs/astro-ph/0611355) · [Earlier GMRT observations](https://arxiv.org/abs/astro-ph/0502307)

## Replay

Run from the repository root using the project's scientific Python environment.

~~~sh
export PYTHONPATH=code:.
python -m cosmology.rar_deep_regime.individual_trail_budget.measure
python -m cosmology.rar_deep_regime.individual_trail_budget.figure
python -m pytest -q code/cosmology/rar_deep_regime/individual_trail_budget
~~~

The receipt records hashes of both SPARC tables and the independently supplied native energy file. The tests recalculate the central ratios for every galaxy, check the nuisance envelopes against sampled interior points including a negative gas force, and reproduce the two outer-ring sensitivity thresholds. [Per-galaxy receipt](results.json) · [Source data provenance](../data/PROVENANCE.md)
