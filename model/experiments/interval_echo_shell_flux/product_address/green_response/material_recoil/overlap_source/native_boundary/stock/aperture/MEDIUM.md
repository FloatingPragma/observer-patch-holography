# A measured ordinary medium carries the effect of a written contact

A laboratory observer in the executable model can infer an exterior release from source observations on separate copies. On the 240-port chart one thermal ordinary unit spends 2.415% of its time on the written seam's two endpoints, compared with 0.794% on the blank chart. A two-unit mobile clock one hop away slows. The independently measured change in its initial frequency relative to the blank arm is between −0.2611 and −0.2524 per native time. Reading the ordinary locations and calibrating with an ordinary neighbor predicts −0.2562 to −0.2557. The local spectral reader assigns the marked contact 4.340 nats relative to the calibrated reference, with interval [4.2826, 4.4012]. This scalar stays fixed when the supply is reversed. The predicted motion also uses the observed distribution of ordinary material. The first exterior cut encloses the seam and its nearest neighbors. Its boundary edges are entirely unmarked. [Preparation and clock read](MEDIATED_RELEASE.md) · [Records](mediated_release.json)

## What the source readings prove

Let $p_x$ be the probability of recording the ordinary unit at port $x$, and $s_x$ the sum over its neighbors on unmarked edges. Inserting a two-unit probe in the equal real wave over the other $N-1$ ports preserves the ordinary address. It can also be done coherently, leaving the ordinary address unread. For each exterior live edge oriented from $v$ to $u$, write $J_{v\to u}$ for inward probe probability per native time. The actual exchange gives

$$
(N-1)\dot J_{v\to u}(0)=4\kappa[(s_u+p_u)-(s_v+p_v)]
-\partial_t\operatorname{Im}\rho_{uv}(0).
$$

The last term is measured before insertion by a two-port interference analyzer. At equilibrium it vanishes. The equilibrium exterior field driving the release is thus $4\kappa(s_x+p_x)$ in native units. At any exterior cut, summing the oriented bond result gives the acceleration of the enclosed probability. The common coefficient is fixed by the ordinary contact. An orbital interpretation of these release rates is a separate hypothesis about spatial scale and persistent motion.

| Ordinary preparation | Prediction from independent source reads for the first exterior cut | Native acceleration change |
| --- | ---: | ---: |
| Equilibrium, address unread | [0.023011, 0.023064] | +0.023038 |
| Position read before cold release | [0.023020, 0.023056] | +0.023038 |
| Same probabilities with aligned amplitudes | [0.022878, 0.022931] | +0.022905 |

The entries are changes relative to the corresponding blank board, in probability per squared native time. The intervals are rounded outwards at joint 99% counting coverage for the respective experiments. The undisturbed equilibrium position banks at time 1/32 resolve an additional enclosed probability of $8.74\times10^{-6}$ with error at most $2.57\times10^{-6}$. At the next exterior cut the inferred equilibrium acceleration changes by between $-5.69\times10^{-5}$ and $-2.94\times10^{-5}$. The experiment measures changes in enclosed probability. Individual orbital accelerations require an additional trajectory measurement. [Undisturbed protocol](UNDISTURBED_MEDIUM.md) · [Time, position and quadrature counts](undisturbed_medium.json)

### Finite derivation

Write a material configuration as $(i,j)$ for the double and single addresses, with $i\ne j$. The columns $W_z(i,j)=1_{j=z}$ describe the conditional launches up to normalization. On a magnitude-two written edge the native memory factor evolves identically on all seven possible local material words in this sector. On the other edges the loads hop and can exchange places. In the resulting joint generator $H$, with the common isolated energy subtracted, the diagonal contrast term is

$$V(i,j)=-\kappa[4m_i+m_j+4\,1_{i\sim_U j}],$$

where $m$ counts incident writes and $U$ is the live graph. Because every exchange moves the single address exactly as on the live one-unit graph,

$$(H W_z)(i,j)=\tfrac12 (L_U)_{jz}+V(i,j)1_{j=z}.$$

For an initial ordinary density $\rho$, the launched ensemble is $W\rho W^*/(N-1)$. Differentiate the native current on one live seam. For every spectator address $z$ outside $u,v$, the contribution in units of $1/(N-1)$ is $p_z[V(v,z)-V(u,z)]$. The one transition exchanging the loads supplies $\operatorname{Re}[L_U/2,\rho]_{uv}+\operatorname{Re}\rho_{uv}[V(v,u)-V(u,v)]$. On an exterior edge both $m$ vanish, the exchange contrast difference vanishes and the single-load generator has the same commutator there as $L_U/2$. The ordinary evolution gives $\operatorname{Re}[H_b,\rho]_{uv}=-\partial_t\operatorname{Im}\rho_{uv}$. These identities prove the exterior formula for an arbitrary ordinary preparation on this sector.

If the ordinary address is sharply read first, the current simplifies on every live edge to

$$(N-1)\dot J_{v\to u}=4\kappa(m_u-m_v)(1-p_u-p_v)+F_u-F_v,\qquad
F_x=4\kappa s_x+(4\kappa+1/2)p_x.$$

Integer comparisons with the complete native current verify 151,200 blank edge/address pairs and 150,960 written pairs. The source reader uses only the contact chart, marks, ordinary calibration and position frequencies to predict the exterior acceleration. Its independent binary and multinomial readings retain the unoccupied outcomes. [Native edge checks and source decoder](medium_field.py) · [Finite results](medium_field.json)

### Which clock changes, and over what range?

A charge-conserving comparison transfers the double from a remote uniform reference mode to port $x$. Trials in which the ordinary unit occupies $x$ are flagged in an admission bank. The inserted and reference arms evolve on the open board. Projecting the joint generator onto the insertion port, subtracting the ordinary generator in the reference arm and accounting for their shared exchange gives the initial frequency relative to the remote hand,

$$\omega_x=d_x/2-(4\kappa+1/2)\left(m_x+\frac{s_x}{1-p_x}\right),$$

where $d_x$ is its original degree. One ordinary neighbor on a blank contact calibrates the decrement to [32.4951, 32.5037]. The same coefficient on a written contact follows from the load values two and one. This is an insertion clock of a specified material species. The supplied thermometer sets $\beta=\log(8/7)$, the interferometer prepares the comparisons, and the common contact chart uses hop distance.

A local change can modify the diagonal thermal weight at distance $d$ only through walks that reach and return from the altered endpoints. Moments through order $2d$ agree. For one marked seam, the shifted bath matrices $3I-H_b$ have nonnegative entries and row sums at most 11. The ratio of the local unnormalized Gibbs weights therefore differs from one by at most $2e^{11\beta}(11\beta)^{2d+1}/(2d+1)!$. At eight hops the bound is $1.69\times10^{-11}$ on any such degree-completed carrier. Normalizing to one ordinary unit also produces a small global depletion; it accounts for the frequency shift of opposite sign farther out. The bound applies to the unnormalized local weight at the supplied temperature.

To test a common potential, subtract each measured clock's diagonal kinetic contribution, keep the independently specified live hopping, and calibrate one gain on an ordinary-only release. After a sharp address read this construction predicts a first exterior acceleration increase of 0.022932. Propagating the source counts through the joint and scalar equations bounds the difference between them to [0.0001019, 0.0001101]. The exact joint motion includes an exchange that leaves a record of which ordinary address supplied the missing probe port. Both the inward sign and the small disagreement survive counting error. [Scalar comparison](CLOCK_TO_FLIGHT.md) · [Uncertainty calculation and held times](clock_to_flight.json)

### A driven medium can transmit the change farther

Use the same marked seam in an externally supplied ordinary stream. Each appointment offers one of eight tags. The free sharp single load transfers on a blank half contact with probability one half and stays put on the written contact. Collect the unit at seven surveyed hops and return it to a chosen endpoint. Compare both source locations with their respective blank controls. The same passive record changes the stationary occupation per unit throughput at all 61 interior addresses four to six hops away. At a four-hop station, independent position and injection meters resolve a change between +2.313 and +2.339 occupied appointments per injection when fed from one end, and between −2.235 and −2.225 when fed from the other. The actual read occurs at appointment 2,049. The stationary calculation includes all idle offers and the drain service.

Let $L$ be the original grounded graph Laplacian, let $b$ have entries $+1$ and $-1$ at the two marked endpoints, and let $z'$ be the written-board profile per unit injection divided by 16. Deleting the passive edge gives $L'=L-bb^T$. At a fixed injected current the original profile $z$ and $z'$ have the same right-hand side, so $L(z'-z)=b(b^T z')$. Its two source entries cancel exactly. Grounding fixes the solution uniquely. Reversing the feed reverses the voltage drop and the sign of every nonzero change on the entire interior. There is an extended redistribution with no added monopole at fixed throughput. The throughput itself changes when a single circulating unit is held fixed.

The measured ordinary occupancies also predict opposite changes of released-probe flow for the two feeds. At the first exterior cut, the count intervals for acceleration are [−0.00960, −0.00956] and [+0.00914, +0.00918] per squared native time. Opposite signs persist in the predictions across cuts two through six. Independent position banks at time 1/32 resolve the first exterior cut for each feed. The reservoirs create the background stream. This demonstrates how a measured overlap can have a remote effect whose sign depends on the independently measured supply. [Protocol](DRIVEN_MEDIUM.md) · [Exact profiles and independent counts](driven_medium.json)

## What persists in isolated material?

Both boards receive the same uniform in-phase ordinary launch. They then evolve freely. The probe arrives at an age drawn uniformly from zero to $T$. Separate observers count positions and read two-port interference on other copies of that preparation. A randomized fringe setting measures the correction to the position term directly. On every exterior cut, an integer comparison of the full 57,360-configuration release with the separate source observable verifies the operator identity for **every** one-unit preparation. The probe observable is the change of initial acceleration of enclosed probability relative to the blank board.

| Source or preparation | First exterior cut | Second exterior cut |
| --- | ---: | ---: |
| In-phase launch, recorded source, $T=65{,}536$ | [+0.0110218, +0.0110289] | [−0.0059604, −0.0059566] |
| Same launch, certified infinite uniform-age limit | [+0.0109717, +0.0110785] | [−0.0059931, −0.0059251] |
| From this window size onward, the indicated sign is guaranteed | 19,433,775 | 22,900,210 |

A positive entry means written material increases the initial inward transfer of probe probability. Both intervals in the first line use source measurements only and joint 99% counting coverage on the executable model. The independent long-time bounds use the same preparation, map and native dynamics. At the measured longest window ten of the twelve boundaries show an outward change. The result establishes a lasting, spatially varying response in age-averaged records even with no continuing external supply. [Free evolution protocol](ISOLATED_MEDIUM.md) · [Source counts](isolated_medium.json) · [Long-time enclosure](rest_limit.json)

The uniform random-port launch stays invariant on both boards and gives exactly zero contrast at all ages. Two other preparations agree at every single port and every nearest-neighbor analyzer before free evolution. They use the same 120 disjoint nonneighbor pairs and opposite interference phases. The equal mixture is the invariant launch, so linearity forces opposite responses for every age window. Pairs were selected for sensitivity at $T=64$. At the held $T=65{,}536$, the first-cut intervals are [+0.0065578, +0.0065638] and [−0.0065638, −0.0065578]. Comparing the two paths in a prepared pair distinguishes the launches. The independently measured modular reading is common to all of them. A source law for the sign therefore has to specify the ordinary preparation or its accessible correlations. [Recorded comparisons](isolated_medium.json) · [Linearity theorem](../../../../../../../../../../Lean/ObserverPatchHolography/EinsteinBranch/SourceApertures.lean)

The stationary local field also has a measurable budget. A source observer counts whether an ordinary unit appears in each closed live neighborhood. Summing these probabilities counts each occupied port once plus its live degree. On a degree-six board one unit therefore supplies at most $4\kappa(1+6)=224$ in native units at $\kappa=8$. The positive and negative parts of a written-minus-blank field each satisfy this cap. In the observed thermal preparation their summed magnitudes are [2.7292, 2.7339] and [3.0683, 3.1442], respectively. These outward-rounded intervals propagate the common ordinary clock gain and the neighborhood count uncertainty. They use the position component of the source observable. [Budget and calibration](isolated_medium.json)

### Certifying the sign at arbitrarily long uniform ages

For the uniform rest launch $\psi$, the blank evolution preserves every initial source probability. Let $O$ be the common, unnormalized cut observable obtained by the exact integer contraction, and $A=2H$ the written integer generator. Orthogonal spectral subspaces contribute independently in the infinite uniform-age limit. The verifier encloses this limit using rational integer residuals. It also bounds all transients at finite $T$. [Reproducible verifier](rest_limit.py)

1. Round the complete numerical eigenbasis to a rational matrix $V$, with rational proposed eigenvalues in $D$. Exact row-sum bounds give $\|V^*V-I\|_2<3/4$ and $2\|AV-VD\|_2<r$, $r=10^{-9}$. The resolvent of $D+tV^{-1}(AV-VD)$ outside the union of radius-$r$ disks exists by the Neumann bound for all $0\le t\le1$. Each isolated disk therefore contains exactly one true eigenvalue. There are 107 singleton disks. The calculation grants arbitrary spectra and signs to the remaining subspace.
2. For a singleton let $g$ be the separation from its proposed center to the nearest other center. Its normalized trial column has residual below $r$. Expanding it in the true orthogonal eigenspaces bounds the sine of its angle to the true singleton by $\eta=r/(g-r)$. Its launch weight differs from the trial weight $p$ by at most $\eta$. Its diagonal observable differs by at most $2\eta\|O\|$. The error in its weighted contribution is at most $\eta(1+2p)\|O\|$. The entire unassigned subspace contributes at most $\|O\|(1-\sum(p-\eta))$ in absolute value, including its internal coherences at finite times.
3. For a singleton and a different true eigenvalue, the averaging factor has modulus at most $4/((g-2r)T)$ since $H=A/2$. Cauchy–Schwarz bounds the sum of absolute orbit-read cross terms in that row by $\|O\|$. Counting the conjugates as well gives a sufficient total transient allowance of $\sum 8\|O\|/((g-2r)T)$. Subtract the exact blank reading and divide by the $N-1$ allowed probe addresses. All products controlling the reported intervals use integer or rational arithmetic with outward rounding. The age thresholds are sufficient bounds and need not be sharp.

## A Poisson response measured with one ordinary calibration

A source observer can predict a change in where a second observer finds a freely released probe. On the 240-port board reserve a port as an ordinary source station. Start the mobile ordinary unit in the blank lowest band. On the calibration arm a double load occupies the clamped station during one part in 10,000 of each cycle. On the held arm a measured contact elsewhere on the board is written for that same fraction. The blank comparison uses the same empty clamp. The two sources act on the same 239 accessible sites. The kinetic exchange and the ordinary-neighbor clock set the coupling before the held comparison. Alternating the contacts changes the unit's coherent evolution during each cycle. Mixing blank and fully written trials would give a different response. [Construction and controls](POISSON_RESPONSE.md) · [Finite source experiment](WEAK_SOURCE.md)

### Why an inverse Laplacian appears

Write the blank real one-unit generator as $H_0=D-a$, with $a_{xy}>0$ on a live link, and its positive normalized lowest mode as $\phi$. The site probability $p_x=\phi_x^2$ is measured on separate copies. For a changed generator $H_0+\epsilon V$, let $\psi$ be a real eigenmode with energy $E$. The same blank coefficients $c_{xy}=a_{xy}\phi_x\phi_y$ determine all the following fluxes. With $L_c u(x)=\sum_y c_{xy}(u_x-u_y)$ the two eigen-equations give the exact identity

$$
L_c(\psi/\phi)_x = -\epsilon\phi_x(V\psi)_x+(E-E_0)\phi_x\psi_x .
$$

At zero exposure set $s_x=-\phi_x(V\phi)_x$. Differentiation at the simple blank lowest eigenvalue gives $L_c u=s-(\sum_xs_x)p$ with $\sum_xp_xu_x=0$. The measurable occupation derivative is $2p_xu_x$. The background subtraction follows from conservation of one ordinary unit. In the unreserved uniform case, with blank Laplacian $L$ and written incidence $m$, this reduces to $L(N\partial_\epsilon p)=4\kappa(m-\bar m)$. Its flux through any enclosing cut is exactly the enclosed source less the uniform compensation. At $\kappa=0$ the uniform launch stays invariant even when the written edge no longer exchanges. Boundary edge counts are archived for all three native boards. Dividing the flux by a demonstrated three-dimensional physical area and identifying its gradient with orbital acceleration require further measurements.

An independent spatial survey tests that last step on the native chart. One $\pi$ swap transports an ordinary signal across a chosen blank contact with unit probability. Count the minimum number of such ticks from the measured seam and count crossing contacts at each integer radius. On the 3,840-port board, the mean area for radii 16–31 divided by that for 8–15 is $1167/668=1.747$ in common link units. Even granting arbitrary positive areas to within-owner and between-owner links puts this ratio between 1.683 and 1.753. For area proportional to squared swap distance the ratio must exceed $1147/323=3.551$ even with up to one hop of freedom in the origin. Only 1,835 of 3,840 ports are inside the largest enrolled radius. Thus the calibrated graph flux has a different radial law on this survey from a three-dimensional spherical field. The comparison tests this finite hop-distance and link-area reconstruction. [All chart sizes and exact area ratios](area_survey.json)

The ordinary clamped source gives $V_{\rm ord}=-(4\kappa+1/2)\operatorname{diag}1_{\text{neighbor of station}}$ after a common energy subtraction. Every matrix entry agrees with the projection of the native joint generator. Writing the held seam gives $V_{\rm w}=-L_{\rm edge}/2-\kappa\operatorname{diag}m$. An independent ordinary insertion clock gives $4\kappa\in[31.9951,32.0037]$. At duty $10^{-4}$, all-site relative errors between the predicted changes and the simulated position counts are **0.495% for ordinary loading and 0.192% for the held write**. All 239 ports enter each norm. The finite duty, finite timing window, input calibration and counting error contribute to these discrepancies. [Transfer data and input intervals](weak_source.json) · [Weighted-row identity in Lean](../../../../../../../../../../Lean/ObserverPatchHolography/EinsteinBranch/SourceApertures.lean)

A two-stage interferometric filter prepares the blank band from the common uniform launch. Its simulated pass probability is 0.940205; the rejected preparations are recorded. Rational spectral enclosures bound the unwanted excited amplitude by $1.636\times10^{-7}$. Each accepted copy evolves for a uniformly chosen integer number of completed symmetric cycles, from zero through $2^{21}-1$. The decoder uses independent blank counts for the spatial coefficients and interval perturbation of the grounded Laplacian for the source-input bands. The calculation uses projective clamping. For a finite realization, divide each fragment into $2^{54}$ monitored pieces. Only exchanges on the six incident seams can move the clamp, so $\|(I-P)HP\|\le3$. Two Duhamel integrations bound the difference of the retained and projected propagators on a fragment of duration $t$ by $9t^2/(2\,2^{54})$. Summing the filter, evolution and release intervals gives at most $2.46\times10^{-11}$ in unnormalized propagation norm. This bound describes a possible monitored implementation; the reported counts simulate the projected dynamics. [Filter, exposure accounting and certificate](weak_source.py)

### A held motion record

After the final cycle, park the clamp empty and return the internal seams to blank. The prepared ordinary unit remains. Insert a second, distinguishable mobile double and count its position at release and at native time $1/16$. The blank comparator undergoes the same procedure. A different observer counts ordinary positions and reads a two-port fringe on separate copies just before insertion. The exact exterior source identity, checked by integer contraction on all **56,882** allowed two-load configurations, gives the change in the acceleration of enclosed probe probability. These source measurements determine its inward sign on all eight tested cuts, four for the held write and four for the ordinary standard.

| Distance from held seam in hops | Finite source-only interval for extra initial acceleration | Predicted extra enclosed probability at time $1/16$ | Simulated terminal interval |
| ---: | ---: | ---: | ---: |
| 2 | $[3.9637,4.1052]\times10^{-5}$ | $5.643\times10^{-8}$ | $[4.9918,6.4109]\times10^{-8}$ |
| 3 | $[3.3537,3.6202]\times10^{-5}$ | $4.918\times10^{-8}$ | $[4.2849,5.7040]\times10^{-8}$ |
| 4 | $[1.9358,2.2192]\times10^{-5}$ | $2.857\times10^{-8}$ | $[2.2632,3.6823]\times10^{-8}$ |
| 6 | $[0.6959,1.0070]\times10^{-5}$ | $1.192\times10^{-8}$ | $[0.6569,2.0760]\times10^{-8}$ |

The third column applies the *first-order* source prediction, including the blank two-load detector propagation, without using the terminal data. The second column uses measured source frequencies at the actual finite duty, including the ordinary current derivative. The release changes include all positions and launches on the prepared band. The initial $2^{54}$-trial terminal bank included zero on all eight cuts. The precision repeat used 2,048 independent batches of $2^{50}$ draws per setting. It resolves positive motion changes at all four written cuts and three ordinary cuts at joint 99% counting coverage. These are executable-model frequencies. The initial acceleration and finite-time position change are distinct observables. [Initial test](poisson_flight.json) · [Source-only and precision terminal records](poisson_precision.json)

### Which timed source reading entails inward motion?

On a four-owner repair path, observers can read when a retained contact has magnitude two. Across all 81 initially blank load placements, the smallest positive stationary fraction on any seam is $49/408$. The middle seam for placement 2200 has fraction $424/1055$, while 1100 has zero. Thus a single such stationary seam does not supply the calibrated weak duty $10^{-4}$ by its occupancy alone. In one 512-cycle 2200 record the middle seam was active on 309 ticks. A supplied relay exposed a separate 239-site chart whenever the recorded seam was active. Full-tick exposures decreased the initial enclosed probe acceleration at every enrolled cut relative to blank. [Enrollment and source banks](CONTACT_TIMING.md) · [Timetables, classes and results](contact_timing.json)

A laboratory observer can distinguish two ways of delivering the same weak dose. Center a short pulse on each of those 309 ticks on every copy, or apply the full active ticks on a known fraction $q=512/(309\times10^4)$ of copies. At every completed tick both arms have exactly the same expected integrated exposure. In the second arm the law of total probability gives $a_{\rm mix}-a_0=q(a_{\rm on}-a_0)$ for any release contrast. In the first arm each copy propagates through ordered contacts. For uniformly sampled ages 1 through 512, the source-only two-hop bands are $[3.4120,3.5533]\times10^{-5}$ and $[-1.6110,-1.6102]\times10^{-6}$, respectively. A positive change means an increase in the initial acceleration of the enclosed probability. Negative changes are relative to the blank preparation. The short-pulse arm also increases the three- and four-hop changes and decreases the six-hop change. [Settings and independent source counts](WEAK_TIMING.md) · [All four cuts](weak_timing.json)

For a piecewise-constant observed contact $b(t)$ and a proposed mean $p$, let $F(t)=\int_0^t (b(s)-p)\,ds$ and $\bar H=H_0+pV$. Write $U(t)=e^{-iF(t)V}W(t)$. Then $W$ evolves under $e^{iFV}\bar H e^{-iFV}$. The derivative of this conjugation with respect to $F$ has norm $\|[H_0,V]\|$. Unitary variation of constants and $\|e^{-iFV}-I\|\le |F|\|V\|$ give the finite-matrix bound
$$
\|U(t)-e^{-it\bar H}\|\le |F(t)|\|V\|+\|[H_0,V]\|\int_0^t |F(s)|\,ds .
$$
The endpoint term vanishes for matched total time. The integral retains timing information even then. For a unit-norm input and a bounded detector $O$, the absolute error in its expectation is at most $2\|O\|$ times the displayed bound. This controls when a recorded timetable admits a mean-generator approximation. The actual sign intervals above use the native finite detector identity and measured counts, independently of whether this bound is tight.

Four new native service runs used identical short pulses of fractional tick width $10^{-4}/(424/1055)$ on active readings. No withheld events or renormalization were used. For ages sampled uniformly through 8,192 and 32,768 ticks, **every run had a positive source-only release change at two, three, four and six hops**. In the longer window the actual mean duties ranged from $9.3482\times10^{-5}$ to $1.0309\times10^{-4}$. The six-hop acceleration intervals lie within $[6.1174,10.1101]\times10^{-6}$ across the four runs. The decoder uses the ordinary calibration, site occupancies and real two-port interference recorded on separate copies. Integer contraction of the native two-load dynamics fixes the initial probe acceleration from those records. This discharges the premise for these age ensembles. A physical reason for the short pulse width, a spatial population on galaxy scales and an ordinary-matter conversion are still needed for the astronomical inference. [Prospective settings](REPAIR_EXPOSURE.md) · [Ordered pulses and source-only results](repair_exposure.json)

Leaving the same contact connected for each entire active tick tests the autonomous availability signal at full strength. Every one of the four runs then gives a decrease relative to blank at both four and six hops in both long age windows. At 32,768 ticks all four two-hop signals increase. The source readings thus predict a spatial change of sign even after long operation. Condensing the written exposure into short pulses on each active tick is a material condition of the tested extended inward response. [Full-duration control and source banks](full_repair.json)

### A connected modular charge measured by ordinary pairs

There is an uninterrupted interaction in the native law. A written zero and a written magnitude two remove the same hopping contact. For a single mobile ordinary unit every one of their local columns is the vacuum head generator on that material word, so all material motion agrees for any history of passive evolution at fixed written locations. The same statement holds for a distinguishable single and double together. With two ordinary singles there is one exception. When both endpoints of the written seam are occupied, the magnitude-two head sees an additional term. For coin parity even it is $-4\kappa$ or $-(4\kappa+1)$ according to sign-reversal parity; for coin parity odd and for a written zero it vanishes. The term commutes with the vacuum head generator. Integer checks cover all 99 output rows on every enrolled input and both generator coefficients. Ordinary phase comparisons at empty and doubly occupied endpoints enclose $-32$, $-33$ and zero at $\kappa=8$. Each location pair in the uniform two-unit launch has probability $1/28680$. This supplies a small initial exposure without imposing a pulse duration. Subsequent pair probabilities evolve under the native Hamiltonian. [Frozen setting and complete columns](PAIR_CONTACT.md) · [Ordinary frequency banks](pair_contact.json)

Two observers see a connected modular increment on the original 2200 repair source. Read the adjacent-owner seams 10–20 and 14–24. Accept when both are written and each has magnitude zero or two. The sharp stationary service gives the joint counts $((48,176),(112,112))$ over 448 accepted words and 1,662 rejected words. Let $P$ be this reference and $Q=P_A\otimes P_B$ the independently re-paired source at the same marginals. The change from $P$ to $Q$, read against $\log Q-\log P$, is
$$
 q=\sum (Q-P)\log(Q/P)=\delta\log\!\frac{P_{00}P_{11}}{P_{01}P_{10}}
   =\frac{\log(11/3)}{14},\qquad \delta=P_{11}-Q_{11}=-\frac1{14}.
$$
There are six admitted sharp head labels at each location. All 36 joint words occur and their likelihood ratio $P/Q$ is exactly constant within each magnitude block. The magnitude statistic therefore gives the **full connected logarithm on this prepared support**. A separate 36-outcome bank gives $[0.09228,0.09333]$ nats. Retaining the sufficient four bins gives $[0.09274,0.09287]$ nats. These are conditional reference readings on a sharp service with both margins preserved. [Full reference audit](FULL_PAIR_REFERENCE.md) · [Both sets of raw counts](pair_contact.json)

Park two ordinary units coherently across the two endpoints of the first or the second seam. Recombine those alternatives after a dwell $\tau$. For a sharp magnitude-two head the calibrated phase relative to a written zero is $g=1/2+(e^{32i\tau}+e^{33i\tau})/4$. For any binary conditional read $Y_{ab}$ the difference of independently paired and native outcomes is $(Q_{11}-P_{11})(Y_{00}-Y_{01}-Y_{10}+Y_{11})$. For this fringe the mixed difference is $|1-g|^2/2$. At $\tau=1/32$ the source predicts an increase of $[0.0084412,0.0084461]$ in the bright output probability. Independent ordinary-pair counts give $[0.0084333,0.0084540]$. A single ordinary unit is exactly blind to the same charge. The laboratory implication is specific to the observable and its ordinary calibration.

The same local instrument can read an unknown charge. Withhold the held head table, count the two bright ports, and invert the calibrated linear relation. Using one million simulated trials per setting, including new reference and clock banks, the null gives [−0.0531, 0.0567] nats, the half-strength source [−0.0012, 0.1086], and the full-strength source [+0.0335, +0.1508]. Only the full-strength interval excludes zero. These counts establish shot-noise sensitivity for the ideal instrument. Its physical implementation, environmental noise and time per trial are unknown. [Withheld-source meter protocol](PAIR_METER.md) · [All three ordinary-only reconstructions](pair_meter.json)

For the unparked pair write $S$ for the two occupied sites, $m_i$ for incidence on a written contact and $E_U(S)$ for the number of live edges with both ends occupied. Up to the common isolated energy its generator at fixed head parities is
$$
 H=L_{2,U}/2-\kappa\Big(\sum_{i\in S}m_i+2E_U(S)\Big)
       +\sum_{e}\lambda_e\,1_{S=e}.
$$
Here $L_{2,U}$ exchanges a unit with an empty neighboring site and $\lambda_e$ is the occupied-pair shift read above. Starting with the incoherent uniform distribution makes every normalized identity invariant and all source changes zero. The in-phase launch has the same initial position counts. After 64 time units the increased connected logarithm decreases the fraction of ordinary material inside the two-, three- and four-hop boundaries. All times 1, 4, 16 and 64 and the single-seam controls are retained. [Uninterrupted motion and independent position counts](pair_motion.json)

#### Release of a separate, freely moving probe

Restore both seams to blank after the ordinary pair evolves. Let $W$ insert a double in an equal-amplitude superposition of the $N-2$ empty sites, before normalization. Exchange commutes with this insertion, $T_3 W=W T_2$. Let $N_R$ count the two singles inside the receiver boundary. In the launch subspace the enclosed double-position operator has pullback $|R|I-N_R$. The kinetic contribution to $(N-2)$ times initial enclosed acceleration is therefore $[T_2,[T_2,N_R]]$. The contrast potential for a double at $x$ and singles at $S$ is $-2\kappa E(S)-4\kappa s_x(S)$, where $s_x$ counts its single neighbors. In the acceleration, only its difference across a boundary edge $u\in R$, $v\notin R$ enters. When both are vacant the diagonal contribution is $4\kappa(s_u-s_v)$. When a single and the double swap, the symmetric off-diagonal entry between ordinary pairs $B+u,B+v$ is $\kappa(s_u(B)-s_v(B))$. These identities follow from $-[H,[H,P_R]]$ and include the interference terms of undisturbed material. Sharp position reading removes the off-diagonal entries. Independent construction of the complete three-load matrix on an irregular seven-site chart verifies every entry for three boundaries at both stiffnesses.

For each of the four conditional source preparations, read ordinary position or ordinary pair interference on separate copies and fix $4\kappa$ from the ordinary standard. Multiply the resulting mixed response by $-q/\log(3/11)$. This predicts how a different observer's freely released double initially moves relative to blanking the same contacts on the native joint source. The held columns use independent ordinary-material readings on the two actual source arms. No release outcome is needed. The entries below are changes in enclosed probe probability per squared native time, in units of $10^{-9}$, after 64 native time units of source evolution.

| Boundary (hops) | From connected charge and conditional calibration | From held ordinary-material records |
| ---: | ---: | ---: |
| 2 | [−473.19, −467.41] | [−487.40, −453.42] |
| 3 | [−305.78, −297.26] | [−329.87, −274.34] |
| 4 | [−163.99, −155.37] | [−188.53, −131.40] |
| 6 | [+90.85, +97.42] | [+72.46, +117.14] |

At age 16 the source predicts decreases at two and six hops and increases at three and four. At ages 1 and 4 the source intervals include zero on every cut. The full native source measurement thus closes a signed release implication on this chart. A positive connected logarithm alone does not force increased inward acceleration at each radius. The readout, prepared material and its age specify the implication. [Probe preparations](PAIR_PROBE.md) · [All 32 comparisons, instruments and counts](pair_probe.json)

### What a permanent write changes

At full exposure the exact one-unit operator is $H=L_U/2-\kappa\operatorname{diag}m$. On the two-dimensional subspace of the two written endpoints its Rayleigh quotient is at most $5/2-8=-11/2$. Off those endpoints the quadratic form is nonnegative. Min–max and $H\ge-8I$ put exactly two modes in $[-8,-11/2]$; all remaining energies are nonnegative. For each such bound mode the exterior equation implies that the maximum absolute amplitude beyond hop $r$ is at most $(6/17)^r$ times its global maximum, by induction on the rows with degree at most six. Consequently the occupation beyond distance $r\ge1$ is at most $2(180/289)^r/(1-180/289)$, capped at one. This uses a bound of $2\cdot5^j$ sites at distance $j$ from the two endpoints and holds uniformly in carrier size.

For any real bound eigenmode let $p_x=\psi_x^2$. Outside the source,
$$(L_U p)_x=4E p_x-\sum_{y\sim_U x}(\psi_x-\psi_y)^2.$$
Thus the summed outward drop of the **measured occupation probability** across a cut enclosing the source is nonnegative. For the stationary probe field $F=4\kappa(I+A_U)p$, the only extra term when commuting $L_U$ past $A_U$ is $\sum_y(d_x-d_y)p_y$. The degree spread is at most six while $-4E\ge22$. This gives $L_U F\le0$ outside a one-hop collar. A strictly occupied exterior then forces an inward initial transfer of the probe even on an irregular chart. These are *absolute* release and source statements. The appropriate blank transfer can itself be positive on an irregular chart and must be subtracted to claim extra attraction.

Starting one ordinary unit at either measured endpoint with equal probability puts 97.52%, 97.49% and 97.49% of its spectral weight in the two bound modes on the 240-, 960- and 3,840-port boards. The general energy bound for this preparation is at least $11/16$. All launches enter the full comparison. At the 65,536-native-time uniform window on the smallest board, its first two exterior changes from blank are $[1.3125,1.3130]$ and $[-0.029779,-0.029767]$ in enclosed probability per native time squared. Strong binding and long-range weak susceptibility are distinct predictions of the same native equation. [All cuts, exposures, source counts and kinetic calibration](poisson_response.json)

### Do insertion clocks and coherent propagation define the whole motion?

We also allowed a separate kinetic exchange amplitude on each live link in the 240-port thermal experiment. Let $p_x$ be the ordinary occupancy and $C_{xy}$ the measured real two-port coherence. The off-diagonal coherent generator inferred between successful insertions at neighboring ports is
$$G_{xy}=-\frac{1-p_x-p_y+C_{xy}}{2\sqrt{(1-p_x)(1-p_y)}}.$$
Together with the diagonal insertion phases this agrees element by element with the compression of the full native relative generator. For the uniform successful launch $f_x=\sqrt{(1-p_x)/(N-1)}$, the exterior site quantity $(Gf)_x/f_x=-4\kappa s_x/(1-p_x)$ follows, where $s_x$ is the ordinary occupancy of the live neighbors. The exclusive coherent prediction for the first-cut change is [0.0229316, 0.0229453]. The inclusive release exceeds it by [0.00009820, 0.00009954]. Ten of twelve cuts resolve a residual change. All use the same ordinary calibration.

A pulse of duration $\pi$ on one blank native contact swaps every one of its nine load words exactly at $\kappa=8$. Every successfully inserted double arrives at the opposite port. Interfering that route with a direct insertion measures the ordinary records that accompany the swap. On the first registered exterior link its conditional visibility rises from [0.9961740, 0.9961760] to [0.9962346, 0.9962367] on writing the source. A closed reversible spatial propagator on the probe coordinates alone would give unit visibility for a certain transfer. A common gravitational description would have to account separately for the ordinary material records carried along with the exchange. Repeating the release at zero contrast stiffness on the same thermal chart gives vanishing initial acceleration at all twelve exterior cuts on both arms in the native calculation. [Calibrated exchange experiment](exchange_geometry.json)

### Which records can the repair use?

Each of four owners prepares two ordinary units with probability $p$ and zero otherwise. Every placement of a given inventory reaches the same stationary class. The five class sizes in order of increasing inventory are 27, 884, 2,110, 884 and 27. Observers mix these populations before selecting the trials where both of the first two recorded seam magnitudes are zero or two. Each trial retains its acceptance flag. Rejected trials have the same receiver outcome on both arms.

On accepted copies, exchange the second head between independent partners that share an available witness record. The repaired pair retains both one-head distributions with that witness. For one witness value let $a,b,c,d$ be the four native cell probabilities per attempted patch, and let $n$ be their sum. Independent re-pairing changes the doubly positive cell by $(bc-ad)/n$. Sum over witness values to obtain $D$. At fixed one-head marginals, the connected logarithm on the receiver's two heads changes by $q=-D\ell$, where $\ell=\log(P_{00}P_{11}/P_{01}P_{10})$ is measured on the original pooled pair. A detector with measured mixed response $m$ changes by $Dm=-qm/\ell$ whenever $\ell\ne0$. The weight of accepted copies is already included. Passive omission leaves the physical pairing intact and its measured fringe change is consistent with zero. [Enrollment](RECORD_DENSITY.md) · [Exact class counts, odds and measured charge](record_density.json)

At $p=1/16$, independent simulated source counts give the following signed results per attempted patch. The third-head reading includes its blank outcome. The conditional ordinary calibration was measured on a separate 240-port board. Positive release means increased initial acceleration of enclosed probe probability at the stated cut, following 64 units of ordinary preparation time.

| Records available to the repair | Connected charge, nats (exact law, rounded) | Bright probability change on independent receiver | Sign of two-hop / six-hop release from source |
| --- | ---: | ---: | --- |
| First two heads only | $+2.728\times10^{-6}$ | $[-6.331,-5.688]\times10^{-6}$ | increased / decreased |
| Third head also available | $-2.315\times10^{-4}$ | $[5.1220,5.1285]\times10^{-4}$ | decreased / increased |
| Four final ordinary loads available | $-9.984\times10^{-6}$ | $[2.1845,2.2488]\times10^{-5}$ | decreased / increased |
| Third head and four loads available | $-1.251\times10^{-5}$ | $[2.7393,2.8036]\times10^{-5}$ | decreased / increased |

At low filling a repair without a witness has $D=-(3041616/51527255)p^2+O(p^3)$. With access to the third head, $D=(56/663)p+O(p^2)$. Both share limiting odds $341887/278520$. On the two-head receiver the reference is faithful at every interior filling. A receiver also retaining the third head would see a repaired outcome outside its original support with probability 0.004519 per attempt at $p=1/16$. Its reference logarithm requires a separate treatment.

The measured sign also depends on the reference. At $p=1/4$ a shutter accepting a fraction $q_{1/16}/q_{1/4}$ of blind repairs gives exactly the same positive connected-log change as the unthinned $1/16$ source. The changes in the doubly positive probability are respectively $+5.3473\times10^{-6}$ and $-5.1084\times10^{-5}$ per attempt. The same calibrated fringe and two-hop release therefore change in opposite directions. Keeping the reference odds with the charge resolves the sign on both preparations. The finite result identifies what must accompany a source measurement in any proposed common gravitational law. [Source-to-receiver prediction and held fringe counts](record_receiver.json)

### Reading the source with a longer dwell

Observers can increase the source signal using the two measured nonzero ordinary frequencies. At a dwell $2\pi/65$ their phases lie on either side of a half turn. The predicted mixed bright response is $(1+\cos(\pi/65))^2/8=0.4994$, close to its upper bound of one half. For the two-doubles, two-empty native reference, this gives a **3.564–3.571 percentage point** increase per accepted pair. An ideal test of that signed alternative against zero requires at most **14,504 accepted receiver trials per arm** for a 1% error bound in each direction, conditional on the source and frequency calibration intervals. Fresh model counts gave 7,953 and 8,668 bright outcomes. Source selection accepts 448 of 2,110 native words. A re-paired preparation draws heads from independent accepted copies.

At filling $1/16$ retain every attempted patch. The conservative sufficient counts per arm are 34.67 billion with blind repair and 4.35 million with the third seam available. For a predicted minimum shift magnitude $d$ we used the preassigned direction, a threshold $d/2$ and Hoeffding's bound $\exp(-Nd^2/4)$ for each of the two error probabilities. These numbers describe independent ideal detector shots. A physical analogue requires measured preparation efficiency, coherence over the dwell, and a mapping of its readouts to the OPH source. [Source-readout protocol](SOURCE_DWELL.md) · [Counts, detection bounds and phase calibration](source_dwell.json)

### Does the measured source carry energy at equilibrium?

Re-pairing the accepted 2200 stock increases the energy of a separately equilibrated ordinary pair by about 0.000494 native units. Source counts predict the held energy reading. At the same temperature the positional work is short-ranged, and a clamped ordinary comparison gives opposite relative strengths at the near work and two-hop release meters. The source reference also gives a quantitative test of stationarity at fixed one-head marginals. [Thermal energy, reversible work and finite derivation](THERMAL_RESULTS.md)

### Reproduction

From the repository root, with SciPy installed, run the commands below. The count-only modes reconstruct the source predictions without propagating the probe.

```sh
export PYTHONPATH=code:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture
python -Bm "$pkg.mediated_release" --counts-only
python -Bm "$pkg.mediated_release"
python -Bm "$pkg.medium_field"
python -Bm "$pkg.undisturbed_medium" --counts-only
python -Bm "$pkg.undisturbed_medium"
python -Bm "$pkg.clock_to_flight"
python -Bm "$pkg.driven_medium"
python -Bm cosmology.paid_trail_bao.mediated_clock
```

The isolated experiments replay with `python -Bm "$pkg.isolated_medium" --counts-only`, `python -Bm "$pkg.isolated_medium"` and `python -Bm "$pkg.rest_limit"` in this same environment.

For the Poisson response use the same package and thread settings. Replay the modules in this order.

```sh
python -Bm "$pkg.poisson_response"
python -Bm "$pkg.weak_source" --counts-only
python -Bm "$pkg.weak_source"
python -Bm "$pkg.poisson_flight"
python -Bm "$pkg.poisson_precision" --counts-only
python -Bm "$pkg.poisson_precision"
python -Bm "$pkg.contact_timing"
python -Bm "$pkg.weak_timing"
python -Bm "$pkg.repair_exposure"
python -Bm "$pkg.full_repair"
python -Bm "$pkg.pair_contact" --counts-only
python -Bm "$pkg.record_density"
python -Bm "$pkg.record_receiver"
python -Bm "$pkg.source_dwell"
python -Bm "$pkg.pair_meter"
python -Bm "$pkg.pair_motion"
python -Bm "$pkg.pair_probe" --counts-only
python -Bm "$pkg.exchange_geometry"
python -Bm "$pkg.area_survey"
PYTHONPATH=code:. python -Bm cosmology.paid_trail_bao.poisson_transfer
PYTHONPATH=code:. python -Bm cosmology.paid_trail_bao.exterior_marks
```
