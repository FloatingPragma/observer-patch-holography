# M10 · Row 26 · Gravity in the Newton limit

**Scope**: The pinned native experiment concerns direct contact with a compact pair. The finite exterior and weak-field derivations use separately specified spherical or stationary matter, common metric dynamics and an ordinary calibration.

**Local reviewer draft. The official mark is 🚧.** OpenWave asks for inverse-square attraction between masses through gravitoelectromagnetism (GEM). The strength earns credit when the model produces it. On a spherical OPH construction, five sets of freely released clocks recover the attraction due to ordinary and joint record stocks. The finite Einstein theorem explains what forces the sign and range once the measured source and geometry meet its premises. On the adopted Interval Echo pin an exact calculation instead identifies a contact interaction for the compact pair and a third load. [Criterion](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L130) · [M10 cell](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L462) · [Adopted source and law](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/openwave/xperiments/m10_oph/research/m10_theory_canonical.md)

**Draft cell, unchanged mark**

> 🚧 Not yet tested. The earlier electron-anchored conversion gives a Newton constant 1.18 × 10⁴⁰ times too large. Five release radii in a separate spherical construction recover inverse-square attraction at a supplied coupling. On the pinned board the compact pair contributes a contact potential. A native long-range gravitational source and its physical strength remain open.

## From measured records to the Einstein equation

In OPH observers compare the material records inside a small ball and survey its boundary at fixed proper volume. Suppose the independently weighed positive source raises the bulk entropy. Stationarity then requires a boundary-area decrease of exactly the corresponding magnitude. This statement fixes a comparison of two observed areas at the same enclosed volume. At fixed proper radius, positive local spatial curvature instead reduces enclosed volume relative to the reference. The comparison specifies which quantity the observer holds fixed. [Fixed-volume area expansion](../../../code/matched_einstein_limit/PROOF.md#geometric-normalization-and-the-three-errors)

Write $s,a$ for the measured bulk-entropy and boundary-area changes. Let $C=4\pi\ell^4/15$ and $\eta>0$ be the area coefficient. In each timelike direction $u$ the three measured discrepancies are $s-2\pi C T_{uu}$, $a+C(G_{uu}-\Lambda)$ and $s+\eta a$. Combining their signed values leaves exactly $\eta C(G_{uu}-\Lambda-2\pi T_{uu}/\eta)$. For the theorem's ten calibrated directions, finite inversion bounds every tensor component by

$$
\max_{\mu,\nu}|G_{\mu\nu}+\Lambda g_{\mu\nu}-8\pi G_N T_{\mu\nu}|
\ \le\ \frac{255}{32\pi\eta\ell^4}
(e_{\rm bulk}+\eta e_{\rm area}+e_{\rm stationarity}),
\qquad G_N=\frac1{4\eta}.
$$

Here the three $e$'s bound the corresponding absolute discrepancies, using the same preparations, balls, clocks and stress normalization. If the combined errors shrink faster than $\ell^4$ along the admitted physical refinement, the exact Einstein equation follows in the limit under the theorem's common-tetrad convergence hypotheses. A local experiment that establishes these premises can establish a macroscopic gravitational consequence without resolving a microscopic gravitational pull. [Finite tensor bound and limit](../../../code/matched_einstein_limit/PROOF.md#quantitative-einstein-residual) · [Lean theorem](../../../Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean)

There is an affirmative instance. In the regular homogeneous interior of the published spherical construction, the permanent clock-pair instrument has bulk, area and stationarity errors $O(\ell^5),O(\ell^6),O(\ell^5)$, respectively. Independently counted stocks enter the measured shell action as two pressureless masses, and the cofinal construction preserves enrolled readers. The bound gives a component residual $O(\ell)$ and the limiting equation for their sum. Its exterior and surface matching use the specified shell law. The interior proof alone does not measure the exterior of the pinned board. [Source, shell matching and theorem](../../../model/candidate_laws/interval_echo_published_records/spherical_dust/PROOF.md) · [Finite clock pairs and refinement](../../../model/candidate_laws/interval_echo_published_records/spherical_dust/clock_pairs/REFINEMENT.md) · [Shell identities in Lean](../../../Lean/ObserverPatchHolography/EinsteinBranch/SphericalDust.lean)

## The gravitoelectric limit and attraction between extended bodies

Use $c=1$, signature $(-+++)$, a common subtracted vacuum and the isolated, asymptotically flat weak-field limit. Write $g=\eta_{\rm flat}+h$, $\bar h_{\mu\nu}=h_{\mu\nu}-\eta_{{\rm flat},\mu\nu}h/2$, and choose harmonic coordinates $\partial^\mu\bar h_{\mu\nu}=0$. At leading order the reconstructed Einstein equation gives $\Box\bar h_{\mu\nu}=-16\pi G_N T_{\mu\nu}$. For a stationary, slow compact source with $T^{00}\simeq\rho$, $T^{0i}\simeq j_i$ and negligible material stresses at retained order, set $\bar h_{00}=4U$, $\bar h_{0i}=-4A_i$. The stationary equations and isolated boundary data give

$$
\nabla^2 U=-4\pi G_N\rho,\quad
\nabla^2\mathbf A=-4\pi G_N\mathbf j,\quad
(U,\mathbf A)(\mathbf x)=G_N\int
\frac{(\rho,\mathbf j)(\mathbf y)}{|\mathbf x-\mathbf y|}\,d^3y .
$$

Regularity through the source and decay at infinity select the Green solution. An imposed harmonic field or external tide would change the readings. For matter with appreciable pressure, write $\Phi=-h_{00}/2$ for the static lapse potential. Trace reversal gives $\nabla^2\Phi=4\pi G_N(\rho+\sum_i p_i)-\Lambda$ to this order, with $\Phi=-U$ in the slow pressureless limit. Use a common vacuum and include boundary forces when comparing masses. [Static trace and pressure boundary terms](../../../code/matched_einstein_limit/PROOF.md#what-the-finite-bound-guarantees) · [Finite force-balance identities](../../../Lean/ObserverPatchHolography/EinsteinBranch/StaticMassFlux.lean)

Since $g_{00}=-1+2U$, a marked clock released from rest has initial acceleration $\mathbf a=\nabla U$ at Newtonian order. It points toward positive mass. For a spherical source, Gauss's law fixes its total inward flux on every enclosing sphere while the area grows as $4\pi r^2$. Hence $U=G_N M/r$ and $a_r=-G_N M/r^2$ outside. For **two disjoint spherical bodies** at separation $d$, average the potential of either body over the other. It is harmonic in that body's volume, so every spherical layer gives its central value by the mean-value property. The interaction energy is therefore $-G_NM_1M_2/d$ for arbitrary radial densities. Taking its separation derivative gives equal and opposite attractive forces of magnitude $G_NM_1M_2/d^2$. Dividing by each body's independently defined inertial mass gives the two recoils under universal stress coupling. This is the instantaneous weak-field, slow-motion two-body consequence. The simulated releases below measure the test-body case. For an extended nonspherical source, multipoles can change individual directions even when the enclosing flux gives the same mass. [Static shell averages](../../../code/matched_einstein_limit/source_implication/material_clocks/PROOF.md)

The spherical model also gives an exact finite-shell check before surface crossing. Let $R$ be circumferential radius and $w=v^2$ the squared radial drift on its flat spatial slices. The extrinsic-curvature eigenvalues are $v',v/R,v/R$ up to a common sign. Integrating the normal Einstein constraint gives

$$
\frac{R_2w(R_2)-R_1w(R_1)}2
=G_N\int_{R_1}^{R_2}4\pi R^2\rho(R)\,dR .
$$

A regular center fixes the integration constant. With no stress outside the body, every exterior shell must have $Rw/2=G_NM$. The stationary radial clock Hamiltonian $H=vp+\sqrt{m^2+p^2}$ gives $\ddot R=w'/2=-G_N M/R^2$ in the probe's proper time. This explains why attraction is unavoidable on this branch. An additional measured joint stock gravitates at the ordinary gain when its action carries the same pressureless stress. The identification of a native modular charge with that stock is a separate source premise. [Finite-shell derivation and sensor definitions](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/README.md#what-the-einstein-constraint-forces)

## Five release radii and a changed geometry

The ordinary source has individual frequency $1/3$ and zero joint stock. The held source keeps that frequency and adds $8/33$ in coincidences. The spherical law assigns $\mu=G_N M$ with $G_N=1/4$ in its model units. A fleet released from rest at each of $R=2,4,6,8,16$ reports its circumference divided by $2\pi$ before release and at proper time $s=1/8$, in ticks of $2^{-44}$. Two probe masses, 1 and 3, use the same rulers and clock rule. A separate reader takes only the timestamp and the two rounded radius tickets. For any inverse-square release, integrating its acceleration twice bounds $\mu$ between $kR_1^2$ and $kR_0^2$ with $k=2(R_0-R_1)/s^2$. It widens the bounds for quantization and intersects all releases. [Hamiltonian integration and finite enclosures](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/radial_release.py) · [Raw tickets](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/radial_release.json) · [Independent ticket reader](../gravity/headline.py)

| Preparation | Largest lower bound | Smallest upper bound | Result |
|---|---:|---:|---|
| Blank | 0 | $2.92\times10^{-11}$ | Admits zero |
| Ordinary | 0.0833333093 | 0.0833333396 | Admits $1/12$ |
| Held | 0.1439393266 | 0.1439394095 | Admits $19/132$ |
| Held, changed outer geometry | 0.1591204184 | 0.1297233663 | Empty intersection |
| Held, doubled coupling | 0.2878785231 | 0.2878788431 | Admits $19/66$ |

The changed carrier agrees throughout $R\le4$, including the complete inner clock, echo and release experiment. For $R>4$ it multiplies $w$ by $1+q(1-q)^3$, $q=4/R$. The source, asymptotic mass coefficient, apparatus and local clock action are the same. The outer tickets rule out every single inverse-square strength. The geometric boundary term increases from $19/132$ at 4 to $5377/33792$ at 16 across an independently empty annulus. A second control doubles $G_N$ and passes the shape test with double the strength. Equal ordinary and joint response predicts an increase in $\mu$ of $2/33$ after ordinary calibration. The assigned gain and the equal source weights retain their provenance. [Control proofs and numerical bounds](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/README.md#spherical-alternative-and-five-independent-release-locations)

## What a native detector can establish

On the adopted M10.4 eleven-written-seam board, the accepted pair on ports 156, 157 and 166 has squared overlap $729/730$ with an exactly stationary antisymmetric component. Attach a freely moving third binary load with the isometry $J$ on the other 237 ports. For the full three-load generator $H_3$ the integer column calculation gives $H_3 J=J[h_E+(1+3\kappa)I-(1+2\kappa)|163\rangle\langle163|]$, where $h_E$ compresses the one-load generator to the exterior and retains its diagonal. The cancellation holds at every exposure. At native strength 8 the direct potential is $-17$ at the contact port 163 and exactly zero at the other 236 ports. Exchange-only gives $-1$ at the contact. [Native integer source identity and certified transition probabilities](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/written_source.py) · [Receipt](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/written_source.json)

A sharp terminal detector at the contact distinguishes the full bounded preparation. With a probe launched there, the probability of both source acceptance and a detector click after unit exposure is at least 0.95330159. With the source removed, the same full-board detector clicks with probability 0.60048154. Launching at neighbor 162 gives at most 0.03957053 for the prepared source versus 0.18245716 with the source removed. The comparison accounts for exhaustion and the accepted component outside the exact sector. This detects contact scattering. It rules out interpreting the sector's direct potential as a distributed Newton potential. Collective excitations and emergent clock geometry have different source hypotheses.

## An independent current test for the GEM alternative

A separate weak-field calculation distinguishes tensor gravity from a scalar law that duplicates the static force. It adds a source current and two counter-propagating light pulses to the continuum GEM branch. This is an optional discriminant. The published Newton attraction criterion can be tested with sources at rest.

An isolated thin spherical shell of radius $a=1$, mass $M=1/8$ and angular speed $\Omega=0,\pm1/32$ has $J=2Ma^2\Omega/3$. The calculation retains the current to first order in angular speed. Rotational deformation and supporting stress of order $\Omega^2$ lie beyond that approximation. A receiver on a fixed equatorial light guide dates both returns on its own clock, with the guide held stationary relative to the asymptotic reference. Source reversal holds the receiver, guide and scalar field fixed. The Green integral of the mass current gives tangential $A=G_NJ/(2r^2)$ outside the shell. Indeed the azimuthal average leaves $\frac12G_N M\Omega a\int_{-1}^1 u\,du/\sqrt{r^2+a^2-2rau}$. Orthogonality of the Legendre polynomials selects degree one. The density integral selects degree zero and gives $U=G_N M/r$.

In the specified stationary metric on the guide take $g_{tt}=-1+2U$, $g_{t\phi}=-4rA$, $g_{\phi\phi}=r^2(1+2U)$, with $r$ fixed by the empty reference. Integrating the two null roots yields a signed proper return-time difference $-8\pi G_N J/[r\sqrt{1-2G_N M/r}]$ for a positive spin in the chosen orientation. The physical linear GEM prediction is its leading term $-8\pi G_NJ/r$; the square-root factor is the exact clock conversion in the specified approximate metric. The acquisition integrates the shell numerically, solves both light directions and quantizes both arrivals. An independent analytic enclosure checks every integer ticket. The part odd under source reversal is $-0.004122952818$, $-0.002053344321$ and $-0.001024657101$ at reference radii 4, 8 and 16. The same sources in the scalar-only control give zero at each radius. These are synthetic return times in model units. Applying this result to M10.4 requires a derivation of the mass-current response from its records. [Loop generator and independent ticket checks](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/current_loop.py) · [Both directions and all controls](../../../model/candidate_laws/interval_echo/openwave/row_26_newton/current_loop.json)

## Filing decision and replay

**HOLD.** OpenWave revision `c9b5707ac03285eb4c885115245ac2cd5e90273d` identifies the M10.4 source pin `09f062f77f7f540f2783f128ae6314c5d8b818ea`. To claim the positive spherical result for that pin, show that the stress inferred from its observer records sources the reconstructed lapse or curvature with the ordinary coefficient and appropriate boundary data. Establish the units independently for a physical-strength claim. OpenWave records an electron-anchored unit attachment giving $G$ too large by $1.18\times10^{40}$ and light speed $0.00196\,c$. Neither the dimensionless spherical gain nor the optional loop calculation repairs that attachment. The finite dust proof, controlled shape test and native contact theorem retain their stated domains. These are local exploratory results. [Current mark](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L462) · [Pin decision](https://github.com/openwave-labs/openwave/pull/754) · [Filing gates](../../../docs/OPENWAVE_SUBMISSION_GUIDELINES.md#9-filing-gates)

Run from the research repository root. The first command reads the stored radius tickets independently. The next three regenerate the native contact result, the five-radius clocks and the mass-current comparison. The final command verifies the underlying finite-source and spherical-reader regressions. The two acquisition programs compare their results to local receipts by default, as does the current-loop check.

~~~sh
python openwave/interval_echo/gravity/headline.py
python -W error -m model.candidate_laws.interval_echo.openwave.row_26_newton.written_source
python -W error -m model.candidate_laws.interval_echo.openwave.row_26_newton.radial_release
python -W error -m model.candidate_laws.interval_echo.openwave.row_26_newton.current_loop
python -m pytest -q -W error model/candidate_laws/interval_echo_published_records/spherical_dust/test_exterior_reader.py model/candidate_laws/interval_echo/openwave/row_04/test_open_pair.py
~~~

For a longer source-to-theorem replay, the first command verifies the finite measured clock-pair banks and boost directions against their archive. The second checks the quantitative theorem in the Lean kernel using the repository toolchain.

~~~sh
python -m model.candidate_laws.interval_echo_published_records.spherical_dust.clock_pairs.experiment
PYTHONPATH=code:. python -m matched_einstein_limit --lean-check
~~~
