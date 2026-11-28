# Row 26 · Gravity between masses

**Hold for filing.** OpenWave asks for attraction proportional to $1/r^2$ between masses through the gravitoelectromagnetic limit and credits the strength only when the model produces it. Its M10 cell records an untested attraction and a failed physical unit attachment. The pinned M10.4 native calculation below identifies a contact interaction. A separate spherical model produces the desired radial dependence and a quantitative change from overlap records at a supplied coupling. [Official criterion](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L130) · [M10 assessment](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L462) · [Canonical M10.4 specification](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/openwave/xperiments/m10_oph/research/m10_theory_canonical.md)

## Native written source and a freely moving probe

Use the 240-port board, its eleven source-written seams and the free-flight Hamiltonian of M10.4, source revision `09f062f77f7f540f2783f128ae6314c5d8b818ea`. Add one binary load at an exterior port and ask whether a sharp occupation detector at the common contact, port 163, clicks after one unit of exposure. The third load and detector are additional test apparatus. The two-load source has an exactly stationary antisymmetric component with squared overlap $729/730$ in the successful preparation. Its ports are 156 and either 157 or 166. Bounded retries also give an exhaustion outcome. [Pin preparation](../row_04/open_pair_release.py) · [Independent raw-board assembly and source amplitudes](../row_04/open_pair_check.py)

Let $J$ attach the normalized compact pair to a load on any of the other 237 ports. These 237 initial words are orthonormal after attachment. Let $h_E$ be the compression of the one-load generator to those ports, retaining its diagonal, and $H_3$ the full three-load free generator. On this entire space,

$$
 H_3 J=J\big[h_E+(1+3\kappa)I-(1+2\kappa)|163\rangle\langle163|\big].
$$

For a binary occupation word the diagonal of twice the generator is the number of occupied–vacant live seams times $1+2\kappa$, plus $2\kappa\sum_{i\in{\rm occupied}}(6-\deg_{\rm original}i)$. Every live occupied–vacant seam has exchange entry $-1$. The written seams act as zero. The program applies these local integer entries to both sides of each of the 237 columns. At the common contact the two amplitudes for entry onto the occupied tips cancel. Elsewhere the source energy adds the same phase. This proves invariance under unrestricted native evolution for arbitrary exposure at the two evaluated strengths, $\kappa=8,0$, and by affine dependence for all real $\kappa$. At strength eight the direct interaction energy is $-17$ at the contact and exactly zero at all 236 other ports. Free waves can scatter from the contact. Graph steps to it span 1 through 18 and serve only as addresses. [Integer certificate and terminal laws](written_source.py) · [Receipt](written_source.json)

| Probe launch | Contact click, vacant board | Contact click, vacant exterior compression | Contact click, compact source | Accepted **and** click, full bounded preparation |
|---|---:|---:|---:|---:|
| At the contact | 0.600482 | 0.599192 | 0.994011 | 0.95330159–0.99613644 |
| Neighbor 162 | 0.182457 | 0.183096 | 0.002712 | 0–0.03957053 |

These are synthetic outcome probabilities. The source-free full board permits the probe to enter the three source ports. The vacant exterior compression gives a separate comparison with the same accessible ports as the compact-source sector. A rational Taylor sum at unit exposure bounds the transition-amplitude error below $10^{-13}$. For the full preparation the trace-distance bound on any click probability conditional on acceptance is $1/\sqrt{730}$. Multiplying by the success probability $1-(364/729)^8$ retains the failed preparations as their own outcome. This gives a resolvable *local* change even for the full prepared source. A contact potential and wave propagation supply no distributed static $1/r$ potential on the protected branch. Other source sectors, collective behavior and a possible effective clock geometry require their own evidence. The separate clamped thermal experiments also distinguish a short-range attraction from a Newton potential. [Thermal range](mobile_mediator/thermal_range/README.md) · [Finite-density guide](mobile_mediator/finite_density/README.md)

## Spherical alternative and five independent release locations

On the published spherical dust construction the ordinary source stock is $M_o=1/3$. The held source adds joint stock $8/33$, for $M_h=19/33$. These are independent counts in the posted record law. The candidate assigns $\mu=M/4$ and exterior drift $v=-\sqrt{2\mu/R}$. The observers' distance record is a spherical circumference divided by $2\pi$. A fleet of identical test clocks released from rest on a sphere supplies its final circumference after proper duration $1/8$. The reader uses quantized initial and final radius tickets with unit $2^{-44}$ and the probe clock appointment. These synthetic areal radii use the spherical survey's circumference convention. [Source counts and clock Hamiltonian](../../../interval_echo_published_records/spherical_dust/EXTERIOR_READER.md) · [Areal-radius survey](../../../interval_echo_published_records/spherical_dust/survey.py) · [Einstein tensor and regional stationarity on this branch](../../../interval_echo_published_records/spherical_dust/PROOF.md)

The generator advances $H=vp+\sqrt{m^2+p^2}$ for masses 1 and 3, starting with zero radial velocity. Write $w=v^2$. Conservation of energy and proper-time normalization give $\dot R^2=E^2-1+w(R)$ and $\ddot R=w'(R)/2$, extended to release by continuity of Hamilton's equations. Both probe masses therefore obey $\ddot R=-\mu/R^2$ throughout the spherical exterior before surface contact. The observable is a change in areal radius per squared wrist-clock time. Each freely falling accelerometer reads zero. This is the test-body, spherically symmetric gravitoelectric limit. It covers neither mutual recoil of two comparable sources nor a moving-source gravitomagnetic measurement.

An independent reader asks whether **any one** $GM$ could explain the finite release tickets. For a monotone inward release from $R_0$ to $R_1$ in time $s$, put $k=2(R_0-R_1)/s^2$. Integrating the hypothesized acceleration against $2(s-u)/s^2$ gives $kR_1^2\le GM\le kR_0^2$. The reader expands both radii and the displacement by their quantization errors and intersects the intervals from launches at 2, 4, 6, 8 and 16. It uses neither input source mass nor drift. The unchanged carrier also has an independently enclosed cycloid endpoint. [Reader, controls and bounds](radial_release.py) · [All radius tickets](radial_release.json)

| Preparation | Largest lower bound from radius tickets | Smallest upper bound from radius tickets | Common $GM$ |
|---|---:|---:|---|
| Blank | 0 | 0.000000000030 | Admits zero |
| Ordinary | 0.0833333093 | 0.0833333396 | Admits $1/12$ |
| Held | 0.1439393266 | 0.1439394095 | Admits $19/132$ |
| Held, changed outer geometry | 0.1591204184 | 0.1297233663 | Empty intersection |
| Held, doubled coupling | 0.2878785231 | 0.2878788431 | Admits $19/66$ |

The changed carrier agrees identically with the spherical carrier for $R\le4$, including the two-clock redshift, echo and free releases inside that range. For $R>4$, multiply $w=2\mu/R$ by $1+q(1-q)^3$, $q=4/R$. This keeps the clocks, local Hamiltonian, spherical symmetry, asymptotic $1/R$ coefficient and both types of test body. The change is twice continuously differentiable at 4. Its inward acceleration at $R>4$, in units of $\mu/R^2$, is $B=1+q(1-q)^2(2-5q)$. Each source has the same finite stock in both receivers. The altered geometry fails the common-strength test even after allowing every ticket consistent with an a priori trajectory enclosure. The largest certified lower bound is 0.15910711 while the smallest certified upper bound is 0.12997565.

For this enclosure, take $f=R_0-1/300$. Until first exit from the annulus the acceleration magnitude is at most $4\mu/f^2$, since $|B|\le4$. The release displacement is at most $L=2\mu s^2/f^2<1/300$ for every enrolled launch and source, excluding first exit by continuity (the acceleration is inward since $B\ge1-3\max_{[0,1]}q(1-q)^2=5/9$). On the annulus the acceleration can differ from its initial value by at most $84\mu L/f^3$, from $|2B+qB'|\le84$. Multiplication by $s^2/2$ bounds the endpoint error around the initial-acceleration prediction. The program evaluates these enclosures with rational arithmetic, includes both roundings and checks both probe masses.

Doubling the chosen coupling gives a different strength with the same inverse-square shape. The extra held response is $G(8/33)=2/33$ on this branch, fixed once the ordinary coupling and the weight of joint stock in the source law are specified. To claim this as a prediction of the pin, the missing implication must take measured native source and clock records into this geometry with the same normalization. A physical-strength claim also requires independent length, time and mass standards. The spherical calculation supplies neither this source map nor an observed cosmological profile.

## Can the cubic energy–lapse construction reproduce the ordinary strength?

For a finite grounded cubic lattice with $q=1/8$ and any finite profile with $F(0)=1$ and grounded exterior, the unit-flux bound gives energy $D\ge 8/3$ and coupling $g\ge 1/(48\pi)$. Matching the specified gap and signal speed to ordinary units leaves a dimensionless coupling ratio greater than $3.78\times10^{42}$ relative to ordinary gravity. Changing the time and length units together preserves that ratio. The separate spherical experiment takes the ordinary coupling as input. [Unit map, profile and bound](cubic_strength/strength.py) · [Computed strengths](cubic_strength/results.json) · [Finite flow inequality](../../../../../Lean/ObserverPatchHolography/EinsteinBranch/GravityCalibration.lean)

## What the Einstein constraint forces

On the flat spatial slices of this spherical reconstruction, the eigenvalues of the extrinsic curvature are $v',v/R,v/R$ up to a common sign. The normal-normal Einstein constraint is therefore $2(Rw)'/R^2=16\pi G\rho$, with $\rho$ the source density in the slice normal's frame. Integrating over a finite annulus gives

$$
\frac{R_2w(R_2)-R_1w(R_1)}{2}
=G\int_{R_1}^{R_2}4\pi R^2\rho(R)\,dR .
$$

At a regular center the inner boundary term vanishes. Every exterior sphere enclosing the measured stock and no additional stress must then have $Rw/2=GM$. For a stationary exterior the wrist-clock release theorem above gives $\ddot R=-GM/R^2$. Thus the radial dependence is compulsory **if** the Einstein constraint with independently identified stress, the spherical slicing, regular center and vacuum exterior hold. This derivation works at finite shell separation. Its logical direction matters: the spherical executable supplies its drift and checks the constraint, while the native written-source experiment supplies a different Hamiltonian. The homogeneous-interior entropy theorem has a stated domain. Applying its equation on arbitrary native exterior shells requires an exterior premise on those same records. [Einstein tensor, source action and boundary matching](../../../interval_echo_published_records/spherical_dust/PROOF.md)

For the held spherical carrier the computed geometric boundary term is exactly $19/132$ at each of the five radii. At radii 4 and 16 the changed carrier gives $19/132$ and $5377/33792$. It would require nonzero gravitating stock between those spheres, where the source ledger is empty. These rational values are a constraint audit of the assigned geometry. The wrist-clock tickets give a separate observational rejection of the changed profile.

Both commands run from the research repository root with its model dependencies.

~~~sh
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.written_source
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.radial_release
~~~
