# Regional modular energy, mechanical work and the sharp optical screen

On the central bridge of the 80-cell Interval Echo experiment, the complete regional load state has three possible modes. Its modular Hamiltonian is exactly measurable using four native material observables. One source control raises the bridge work by the same amount in both probe species while the twelve radar means and both complete marked forces remain fixed. The resulting regional entropy increases in both species. The forward sharp optical responses have opposite signs. An independent, single-appointment measurement of the local currents and their rates recovers the full sharp optical response. [Stationary head source](../../full_algebra/reuse/INVISIBLE.md) · [Clock-calibrated work](../traction/README.md) · [Paper passage](../../../../../paper/tex_fragments/INTERVAL_ECHO_REGIONAL_MODULAR_WARD.tex)

## Which regional entropy is measured?

Prepare the volume-matched stationary state of the uniformly reused radar. Launch a sharp mobile vacancy or excess at one of the ten equally weighted ports of the two registered rings and let it propagate for four native ticks. The regional observer retains the two load registers at ports 1 and 38 and traces over all other loads, both parked heads, the launch tag and reversal parity. In the sector with one mobile the outside load configuration is orthogonal to an excitation inside. The exact local state is therefore a direct sum of a scalar outside probability and a two-port occupied density. A missed load contributes to the outside probability.

Let $D$ select the two occupied modes. Let $X,Y,Z$ be their Pauli matrices, all extended by zero on the outside outcome, with $Z=P_1-P_{38}$. Read occupation $d$, real and imaginary interference $x,y$, and imbalance $z$. Then
$$
\rho_R=(1-d)\oplus(dI+xX+yY+zZ)/2,\qquad
\lambda_\pm=(d\pm u)/2,\qquad u=\sqrt{x^2+y^2+z^2}.
$$
All three eigenvalues are positive on the reference enclosure. Subtract the outside eigenvalue of $-\log\rho_R$ to obtain
$$
K_R=A D+B(xX+yY+zZ),\qquad A=\log(1-d)-\tfrac12\log(\lambda_+\lambda_-),\qquad B=-\log(\lambda_+/\lambda_-)/(2u).
$$
The native strain actuator supplies work $W=(17D-X)/2$ when its independently measured material-frequency ratio is held fixed. The oriented mobile current from 1 to 38 is $j=Y/2$. The reference state fixes the four coefficients in the exact first law for **every** admissible source tangent,
$$
\delta S_R=(A+17Bx)\,\delta d -2Bx\,\delta\langle W\rangle+2By\,\delta\langle j\rangle+Bz\,\delta z .
$$
These measurements can also be predicted on the fixed hardware from a complete timestamp reconstruction of the incoming head state. The local work and the port interferometer can test that reconstruction independently. [All eigenvalue and coefficient bounds](regional.json) · [Source and propagation certificate](measure.py) · [Full incoming-state timestamp analysis](../../full_algebra/reuse/clock_history.py)

## A finite boundary law

For the actual freely flying state the entropy rate is $\dot S_R=\langle i[H,K_R]\rangle$. The internal bridge Hamiltonian acts only on the local pair. Its contribution vanishes exactly because the reduced state commutes with its logarithm. Hamiltonian terms supported entirely on the complement commute with the lifted regional observable. The remaining contributions are the ten incident boundary hoppings and the endpoint diagonal actions on the head memories. Both endpoint degrees equal six, so the common blank diagonal action is scalar on the occupied pair and commutes with $K_R$. The active marked diagonal terms can change interference with no direct hop across the boundary. The verifier splits the complete native commutators of each local Pauli observable at the two endpoints before contraction with the measured modular coefficients. The results use nats per exchange tick.

| Source | Boundary hopping | Endpoint diagonal actions | Total entropy rate |
|---|---:|---:|---:|
| Vacancy | [−0.206758, −0.206600] | [−0.00005541, −0.00005519] | [−0.206813, −0.206655] |
| Excess | [+0.214863, +0.215153] | [+0.00051585, +0.00051802] | [+0.215380, +0.215670] |

The earlier source pair with equal ranges and equal complete marked forces also illustrates the distinction between the stored and regional entropy. Its dimer entropy slopes are in $[2.9466,2.9521]\times10^{-5}$ for vacancy and $[-1.2724,-1.2669]\times10^{-5}$ for excess. [Earlier exact source direction](../connection/EQUAL_FORCE_ENTROPY.md) · [Regional contractions](regional.json)

## Equal work at exactly the same survey area and volume

On each parity block allow all fifteen traceless Hermitian head coordinates, with both parity weights fixed. Form the 22-row physical response matrix of twelve charge-pooled range means, two independently metered full marked forces, and four local dimer means for each charge. A rational right inverse $K$ satisfies $\|I-MK\|<9.74\times10^{-8}$ for the exact matrix $M$. Choose a target vector $v$ which vanishes except for $-2$ on each species' $X$ row. The exact direction $K(MK)^{-1}v/16788$ realizes this target. Its Hilbert–Schmidt distance from the recorded rational-center proposal is below $9.77\times10^{-8}$. Its norm is less than one. The native Pauli preparation slack exceeds $0.20$ at the stationary state. For either native preparation lottery, a displacement of size $|s|$ costs at most $16|s|\|D_{\rm src}\|_2$ in this slack, by Cauchy–Schwarz on the fifteen coordinates of each block. This proves feasibility throughout $|s|\le0.01$. Individual head marginals may change.

By affineness all twelve ranges, both marked forces, and both species' $d,z,y$ are exactly constant on this family. The tetrahedral area and volume reconstructed from the twelve ranges are therefore exactly constant. Each work mean increases by precisely $s/16788$. The regional entropy first laws reduce to their independently calibrated work components. For a finite change the occupied block derivative is $-X/16788$ and its squared Hilbert–Schmidt norm is $2/16788^2$. Along the segment the smallest occupied eigenvalue decreases by at most $|s|/16788$. The resolvent bound on the logarithm then bounds the entropy second derivative below by minus this squared norm divided by the eigenvalue lower bound.

| Source | Regional entropy per work at fixed $d,y,z$ | Entropy increase from 0 to 0.01 |
|---|---:|---:|
| Vacancy | [0.67421, 0.67593] | above $4.0153\times10^{-7}$ nats |
| Excess | [0.98443, 0.98601] | above $5.8625\times10^{-7}$ nats |

Independent simulated trials on positive preparations around the published conditional fixed-state centers at odd weight 0.1343 compare the two endpoints $s=\pm0.01$. The work instrument resolves its eigenvalues $0,8,9$. A separate regional density measurement resolves the four local means and infers the entropy from its three eigenvalues. With $2^{59}$ simulated trials per setting and simultaneous 99% coverage, the observed entropy increase divided by the work increase lies in $[0.5633,0.7950]$ for vacancy and $[0.8417,1.1415]$ for excess. The population proof above includes uncertainty of the true stationary state and the exact source control. [Exact source control](control.py) · [Rank and finite-change bounds](control.json) · [Complete finite counts](counts.json)

## One fair source with its preparation charge tag discarded

The carrier starts at unit filling. In any occupied two-port region, the local load words are $(0,1),(1,0)$ for a vacancy or $(2,1),(1,2)$ for an excess. Away from the mobile the observed ports have the single common word $(1,1)$. A fair classical laboratory choice between the two charges, with its choice erased, therefore leaves one five-dimensional local state
$$
\rho_{\rm fair}=(1-\bar d)\oplus M_-/2\oplus M_+/2,\qquad \bar d=(d_-+d_+)/2,
$$
where $M_q$ is the occupied two-by-two block of the conditional result. The two occupied blocks are distinguished by total **local** load. No historical tag is needed. Reserve one outside mode shared by both branches before evaluating the entropy. Writing $f(p)=-p\log p$ gives
$$
S(\rho_{\rm fair})-\tfrac12(S_-+S_+)
 =f(1-\bar d)-\tfrac12\sum_q f(1-d_q)+\bar d\log2.
$$
Both $d_q$ are constant on the entire equal-work control. Consequently every *finite entropy change* is precisely the mean of the two conditional changes. The signed slope lies in $[4.94001,4.94977]\times10^{-5}$ nats per unit displacement. The finite gain from $s=0$ to $s=+0.01$ lies in $[4.93893,4.94977]\times10^{-7}$ nats. Its upper bound uses concavity at the reference and its lower bound averages the two finite curvature bounds. [Primitive native load words](../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/ranging/measure.py) · [Five eigenvalues and exact bounds](control.json)

An independent simulated untagged reader performs tomography on the two observed local-total-load sectors and includes every miss. The outside probability is one minus the two measured occupied probabilities. A separate work reader adds symmetric and antisymmetric outcomes across charges before computing work. Both use the same fair preparation with no input charge report. [Untagged five-mode and work counts](counts.json)

## Reading optical focusing from local material observations

For a sharp receiver $y$ define $j_{z\to y}=i(P_zHP_y-P_yHP_z)$. This is minus the derivative with respect to a local port phase twist of the directed bond's Hamiltonian, at zero phase. Native continuity at every receiver and for each launch independently gives
$$
\dot p_y=\sum_{z}\langle j_{z\to y}\rangle,\qquad
\ddot p_y=\sum_z\langle i[H,j_{z\to y}]\rangle .
$$
The incident current probes touch one edge. Their acceleration operators involve at most two hops from the receiver. The native meter expands each operator in real or imaginary two-port interference, single-port population, and joint head Pauli measurements. Misses contribute zero. At a single appointment, at most one population setting, six current settings and 49 acceleration settings per receiver and parity resolve its three time derivatives including order zero. With both tagged five-pixel images, these material readings determine the throughput, moving centroids and full normalized $J,\dot J,\ddot J$. A separate raw-arrival replay at the appointment and three pairs of adjacent clock times reproduces the terminal material reconstruction by extrapolated finite differences. The terminal comparisons are simulated population readings. [Measurement expansion and independent propagation](screen_ward.py) · [Raw arrival stencils and terminal reconstruction](screen_ward.json)

On the *same exact equal-work control*, the signed derivative of $G=-\operatorname{tr}(\ddot J J^{-1})$ with exchange dwell as the time parameter has the following enclosures. Both sources give strictly increasing regional entropy.

| Sharp angular image | Vacancy $\delta G$ | Excess $\delta G$ |
|---|---:|---:|
| Forward | [−23.402, −20.902] | [+9.8046, +9.9755] |
| Reverse | [+37.1259, +38.9444] | [−3.1999, −3.0365] |

The fair untagged optical receiver pools both charges' hit probabilities **before** normalization, centering, differentiation or inversion of the angular image. Its signed focusing derivative lies in $[-120.437,-119.515]$ forward and $[84.717,85.370]$ reverse on the same source direction as the five-mode entropy measurement. Calculating the mean of the two per-charge contractions gives a different observable. The directly simulated five-mode tomography and untagged work measurement compare $s=-0.01$ to $s=+0.01$ in a symmetric two-arm contrast, with simultaneous 99% bounds $[7.2800,12.4679]\times10^{-7}$ nats and $[1.14088,1.23683]\times10^{-6}$ native work units. [Unresolved-charge counts](counts.json) · [Pooled image calculation](optical_certificate.json)

Suppose that the same work change on the two tracers produces equal scaled focusing for a given direction. Even allowing an arbitrary gain per direction, a common source-independent reparametrization with $u=\ddot\lambda/\dot\lambda$ requires $\delta G_-+u\,\delta\theta_-=\delta G_++u\,\delta\theta_+$ in that direction. Here $\theta=\operatorname{tr}(\dot J J^{-1})$. Forward propagation requires $u\in[2.3685,2.5904]$ and reverse propagation requires $u\in[6.0157,6.3701]$. Both demands cannot hold at one clock appointment. The enclosures include uncertainties in the true conditional stationary states, parity mixture and exact source tangent. [Independent fixed-point certification](optical_certificate.py) · [Full numeric intervals](optical_certificate.json)

The terminal Ward instrument measures native load currents and accelerations. Interpreting the screened statistic as proper null curvature needs identified common events, a geodesic with an independently calibrated affine clock, and proper transverse screen transport. Translating the dimer's extensive strain work or three-mode modular energy into $T_{kk}$ needs an independently measured spatial stress, energy flux and proper element on the same ray. The finite counterexample tests a specific identification of those instruments. The independently metered linear response and local boundary law remain available as source inputs to other observer constructions.

From the repository root, run `OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/directional_stress/regional_modular`. The modules `measure`, `control`, `screen_ward`, `optical_certificate` and `counts` under `matched_einstein_limit.source_implication.directional_stress.regional_modular` each accept `--write`.


## The bridge energy current on the same untagged source

Assign each cell its internal Hamiltonian $h_a$ and half of each incident bridge action. For the oriented bridge $e=(1,38)$, with physical action $b_e$, the continuity equation for the receiving cell gives
$$
J_E=\tfrac i2[b_e,h_b-h_a],\qquad \dot E_b=\sum_{e\to b}J_{E,e}.
$$
The finite generator satisfies this identity entry by entry in each charge and parity block. The current has an independent phase-kick measurement. Following the common launch and four-tick flight, apply the isolated bridge action for phase $\alpha$ and read the cell energy difference $h_b-h_a$. Half the derivative of that phase scan at zero equals $J_E$. The bridge has spectrum $\{0,8,9\}$ in native units, so 19 equispaced phases recover the derivative exactly in population. This extensive energy flux includes the material action terms. The two-port load quadrature alone measures transported occupation. [Regional continuity and primitive current meter](../README.md) · [Phase-kick settings](../phase_read.py)

Append one native energy-current effect per charge to the 22-row source map. The two extra constraints keep the charge-unresolved current fixed even when the incoming charge choice is erased. The first control changes both works equally. Two other controls fix all 24 readings and vary, respectively, the even-parity head coordinate $I\otimes X/2$ or $I\otimes Y/2$. Their conditional and untagged regional load density matrices are constant. The controls use matter and radar calibration only, so the released angular images are held predictions. Full interval inversion supplies reachable chords parameterized by $-0.01\le s\le0.01$ about the same stationary preparation, with all prescribed affine equalities exact along their entire length. The source perturbations are preparations about the stationary reference; stationarity of each perturbed head under repeated radar reuse is a separate condition. [Source controls and fixed-point propagation](energy_odd.py)

In an orthonormal frame with $k_\pm=(1,\pm n)$, the null-stress difference is $4T_{0n}=-4T^{0n}$. This algebra cancels the density and normal pressure for every state. Interpreting the difference of the two native screen readings as its response assumes that the two images sample the same event with equally normalized opposite directions, and that the bridge current has a fixed conversion into proper flux density. With those identifications the current derivative vanishes on either control. A common, source-independent clock $\lambda(t)$ changes the directional odd focusing derivative to $(\delta G_{\rm odd}+u\,\delta\theta_{\rm odd})/\dot\lambda^2$, where $u=\ddot\lambda/\dot\lambda$ and $G_{\rm odd}=G_{\rm f}-G_{\rm r}$. The three source controls independently test its possible value. A different transported proper screen, or a source-dependent proper area, event or clock, would need its own measurement.

| Material source control | Exact $\delta W$ per charge | $\delta G_{\rm f}-\delta G_{\rm r}$, charges pooled | Required $u$ |
|---|---:|---:|---:|
| Equal work | $1/57309$ | [105.145, 106.457] | [154.1904, 165.5747] |
| All 24 readings fixed, even $IX/2$ driven | 0 | [−104.355, −103.445] | [67.2528, 68.9376] |
| All 24 readings fixed, even $IY/2$ driven | 0 | [−329.399, −327.216] | [48.1397, 48.9169] |

The reference fair current lies in $[-0.010890,-0.010793]$ native energy units per tick. A vanishing current *derivative* is certified along each chord. The two all-fixed controls have disjoint clock intervals despite identical values of every one of the 24 calibration records along their whole chords. Thus even a flux conversion and clock jet supplied as arbitrary single-valued functions of these contemporaneous records cannot connect the pooled sharp-screen difference to the extensive current. This conclusion uses equally normalized opposite directions at a common event as the hypothesis being tested. A separately measured dependence on additional local records or temporal history lies outside the test. [Rank, current and directional optical enclosures](energy_odd.json)

An independently simulated reader pools both charges before counting the signs of native cell energy outcomes at 19 phases. At $2^{60}$ shots per setting the simultaneous 99% intervals for the symmetric $s=-0.01$ to $s=+0.01$ current contrasts all include zero. The deterministic pooled phase means at the rounded control centers independently reproduce the zero slopes; finite shots establish a precision interval. Removing the quadratic-contrast term in flight while holding preparations and meters fixed also gives a separate changed-law control. [Phase read, sampled counts and changed-law readouts](energy_phase.json) · [Sampling implementation](energy_phase.py)

The complete nine-mode reduction on the four radar stations uses the two additional port loads and their actual coherences with the bridge. It tests the entropy of the same four registered sites as the range comparison, including an independently selected control that fixes the populations of all three number sectors. [Four-station entropy, conditional contributions and measurement protocol](STATION_REGION.md)

### An untagged optical read with a running material clock

Accumulate a tick at the mean instantaneous occupation of the two-site bridge in the vacancy and excess ensembles. At the common appointment the rate of this clock lies between 0.0119761 and 0.0119785 and its acceleration between 0.0006640 and 0.0006902. The model computes the rate of each occupation by the Hamiltonian commutator. On two controls holding both occupations and both energy currents fixed, the untagged odd focusing still has a negative derivative on this clock. The corresponding upper bounds are −720346 and −2277806 in its time units. [Running clock and exact control verifier](occupation_clock.py) · [Interval enclosures](occupation_clock.json)

A stronger test allows the clock rate to be any smooth positive function of the native dwell and both occupations. Four independently selected head changes hold the existing 24 material and radar readings together with both occupation time derivatives exactly fixed. The clock speed and the ratio of its acceleration to its speed must then have zero source derivative at this appointment, by the chain rule. This leaves only a single free reference ratio, common to every change. Write h and g for the signed source derivatives of the expansion and focusing of the untagged image. A fixed flux-density conversion requires g + u h = 0 for every change. The determinant of the rows (h,g) on the last two changes lies between 47.51 and 69.85. Thus clocks of this class cannot restore the signed equation for this image and the metered current. This test uses a single pooled image. A clock using additional material history, a conversion with additional source-dependent area information, or a different transported screen would require its own independent measurements. [Archived source coordinates and determinants](occupation_clock.json)
