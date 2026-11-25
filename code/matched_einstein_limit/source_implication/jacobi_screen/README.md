# Two transverse image axes on the stationary source

At native dwell 4 the vacancy and excess images require incompatible clock accelerations to describe the same tidal matrix. At dwell 4.5 the vacancy image reverses orientation relative to the excess image. A second experiment interchanges the two ends of the apparatus. At fixed exchange time, the ordinary-calibrated energy current fails to predict the collar response of the direction-dependent focusing statistic. These tests constrain the interpretation of this particular beam in the effective metric.

## Preparation and readings

The stationary two-parity source, its trace-norm enclosure of radius $6\times 10^{-7}$, and its previously calibrated ordinary-matter and unpairing variations come from the same 80-cell Interval Echo preparation. Each variation preserves the surveyed radar volume to first order. The forward launch pixel is chosen uniformly from $(0,5,9,8,7)$. The five receiver pixels are $(39,40,47,46,42)$. A sixth outcome records a miss. The charge sign selects a vacancy or an excess. Treating the two carriers as tracers of one background assumes that their insertion leaves the inferred geometry common to both. The signed reads are also compatible with distinct propagation or backreaction. The receiver marks its appointment with the opposite-end probability $\sin^2(t/2)$ of a blank exchange dimer on the known branch $\pi<t<2\pi$. This identifies exchange dwell. A proper null affine parameter requires a separate physical identification. [Source and gates](../beam_response/README.md) · [Station and volume response](../full_algebra/README.md) · [Dimer clock](../material_clocks/README.md)

Let $P_{yx}$ be the unconditional hit probability for sender $x$ and receiver $y$. Give each pentagon its two marked bearing coordinates. Write $W=\sum P_{yx}$ and $\langle\cdot\rangle=\sum P_{yx}(\cdot)/W$. The response of the detected centroid to the two independent first-harmonic launch biases is

$$
J=2(\langle y x^T\rangle-\langle y\rangle\langle x\rangle^T).
$$

All 25 probabilities enter before normalization. The misses enter the unconditional counts and fix the throughput $W/5$. We evaluate $U=e^{-itH}$, $-iHU$ and $-H^2U$ on four internal-state columns at each of the five launch ports and differentiate the full normalized image. Both moving centroids and the throughput contribute to $\dot J,\ddot J$. Define $B=\dot J J^{-1}$, $C=\ddot J J^{-1}$, $\theta=\operatorname{tr}B$ and $F=-\operatorname{tr}C$ when $J$ is invertible. Writing $\sigma=(B+B^T)/2-\theta I/2$ and $\omega=(B_{21}-B_{12})/2$ gives the identity

$$
F=-\dot\theta-\theta^2/2-\operatorname{tr}(\sigma^2)+2\omega^2.
$$

A geometric Jacobi matrix maps angular separation at a vertex to proper transverse displacement along a particular null ray. Its second derivative with respect to an affine parameter in a parallel screen measures the optical tidal tensor. Here both hardware screens stay at fixed sites while the receiver varies dwell. The finite experiment measures transmission, angular distortion and time dependence. Interpreting the result as curvature requires establishing that the observations identify corresponding events along a ray and measuring the proper screen and its transport. [Fleury, Pitrou and Uzan, Eqs. 3.5–3.10](https://arxiv.org/pdf/1410.8473)

## Two tests of a common ray

If two point-vertex images represent one geometrical Jacobi map $D$ up to smooth, invertible changes of launch and detection basis, then $J_q=L_q D R_q$. The sign of $\det J_-/\det J_+$ at nonsingular samples equals the sign of $\det(L_-R_-)/\det(L_+R_+)$. That sign is constant throughout any connected interval of regular calibrations, even when each basis varies with dwell. Both probes can share any monotone reparametrization of the appointments.

| Native dwell | Vacancy $\operatorname{sgn}\det J$ | Excess $\operatorname{sgn}\det J$ |
| ---: | ---: | ---: |
| 4 | + | + |
| 4.5 | − | + |
| 5 | + | + |

The strict sign reversals exclude a single point-vertex map in this detector interpretation. They make no assumption about intermediate visits or caustic locations. A seeded multinomial replay with $2^{62}$ virtual shots per setting reproduces the six signs with simultaneous coverage at least 99% under the nominal forward sampling model. The deterministic enclosures propagate the uncertainty of the stationary state separately. These data are simulation-generated. [Deterministic enclosures](certificate.json) · [Synthetic counts and timing stencils](counts.json)

A second test allows different initial images and derivatives. Suppose only that the two image matrices satisfy $J_{q,\lambda\lambda}=-K J_q$ for the same optical tidal matrix $K$ in the supplied common, parallel receiver basis. A common constant change of receiver basis and independent constant input calibrations preserve this condition. Set $u=\ddot\lambda/\dot\lambda$. The chain rule gives $C_q-uB_q=-\dot\lambda^2 K$. At dwell 4 the $(1,1)$ component of $C_- - C_+=u(B_- - B_+)$ requires $-9.9060<u<-9.9000$, while the $(1,2)$ component requires $19.7725<u<19.7785$. Thus no common scalar clock can identify the two full matrices with one tidal tensor in that screen. The weaker trace equation alone allows $-13.729<u<-13.718$ and gives $-222.45<\dot\lambda^2\operatorname{tr}K<-222.29$. Imposing exchange dwell as affine instead gives $F_-\simeq-319.016$ and $F_+\simeq-121.560$. Finite central differences at three decreasing step sizes approach the analytic values. Switching off seam contrast gives $F_-\simeq3.454$ and $F_+\simeq3.505$ on the same screens. [Complete time derivatives](results.json) · [Clock and changed-law comparison](counts.json)

## Direction-odd focusing and physical energy current

Choose either end for launch with probability one half and retain the direction tag. Use the two tag-conditional images to form $f=F_{\rm forward}-F_{\rm reverse}$ and $s=\theta_{\rm forward}-\theta_{\rm reverse}$. On that same tagged mixture, read the oriented native Hamiltonian energy current $j$ across the $1\to38$ bridge without conditioning on the tag. The two directions use five equally likely launches and the same internal preparation. Transposition of the propagation kernel gives the reverse amplitudes because the native generator is real and symmetric. Independent propagation from the opposite launch sites checks this identity, including the two time derivatives. [Source variations](source_response.json) · [Bridge current](../directional_stress/bridge.json) · [Current enclosures](../directional_stress/certificate.json)

In an orthonormal observer frame with signature $(-,+,+,+)$ and $c=1$, the normalized directions $k_\pm=(1,\pm n)$ give $T(k_+,k_+)-T(k_-,k_-)=-4T^{0n}$. Density and longitudinal pressure cancel. To apply this relation to the tagged experiment one would have to identify opposite rays at a common spacetime event and the corresponding proper-area/time normalization of $j$. If those identifications hold with a fixed positive conversion, the ordinary-matter response calibrates the factor in the odd null balance and predicts the held-out collar response.

| Exchange-time derivative at dwell 4 | Vacancy | Excess |
| --- | ---: | ---: |
| $\delta_o f$ | [9742, 10088] | [−13606, −13239] |
| $\delta_p f$ | [2883, 3006] | [−4053, −3918] |
| $\delta_o j$ | [−3.1100, −3.0651] | [−2.0737, −2.0445] |
| $\delta_p j$ | [−0.8794, −0.8603] | [−0.5983, −0.5860] |

Here $o$ denotes ordinary matter and $p$ unpairing. The signs of the two ordinary focusing responses differ, while both ordinary bridge-current responses are negative. The certificate also compares $(\delta_p f/\delta_o f)/(\delta_p j/\delta_o j)$ within each charge separately, retaining the correlation between the surveyed-volume corrections in both source directions. This held-out ratio is in $[1.05125,1.05694]$ for the vacancy and $[1.03003,1.03492]$ for the excess. Equal response at a coupling calibrated on ordinary matter requires unity.

For a common nonaffine exchange clock, normalized focusing is proportional to $f+u s$. With fixed source-independent clock and units, the two ordinary responses must satisfy $\delta_o f_q+u\delta_o s_q=c\,\delta_o j_q$ for the same signed factor $c$. These two equations imply $c>0$, contrary to the sign $c=-4\kappa$ times a positive conversion for this orientation. The conditional geometric interpretation cannot be recovered by merely choosing another fixed clock acceleration.

For a source-dependent clock, differentiating $f+u s=cj$ at fixed appointment gives

$$
\delta f+u\,\delta s+s\,\delta u=c\,\delta j+j\,\delta c .
$$

Using both unperturbed reference currents and both unperturbed directional images gives $-2.211<u<-0.602$ and $-25.502<c<-22.221$. These reference observations admit a common calibration of the expected sign. With fixed flux units ($\delta c=0$), the ordinary response requires $89767<\delta_o u<131334$ for the vacancy and $262529<\delta_o u<296074$ for the excess. The collar requires $26505<\delta_p u<39167$ and $77696<\delta_p u<88186$, respectively. Both pairs are disjoint. Thus even allowing the same clock law to depend smoothly on the source fails in this interpretation.

A common affine rescaling of both sides scales their null contractions identically. A freely adjustable $\delta c$ would represent an additional, independently unmeasured relation between the two instruments, for example the conversion of energy current per bridge to energy flux per proper area. The determinant of the two-charge matrix with rows $(s_q,-j_q)$ lies between $-0.004587$ and $-0.004373$. Consequently any two optical derivatives can be fitted locally if both $\delta u$ and $\delta c$ may be chosen freely for each source direction. The measured derivatives would require $102484<\delta_o c<166055$ and $30098<\delta_p c<49646$ in the declared units. This freedom explains why a fit using the optical readings as inputs carries no independent evidence for null balance. One proposed conversion takes the bridge's proper area to be a fixed fraction of the surveyed tetrahedron's total boundary area $A$, giving $c=-C/A$ with $C>0$. The selected station obeys $\eta\delta A=-\delta S$ at fixed surveyed volume with positive gain and strict positive entropy variation in each arm. Hence $\delta A<0$ and $\delta c=C\delta A/A^2<0$, the opposite of the required sign in both arms. This excludes that particular area conversion. An independently identified local transverse area and lapse could behave differently. [Finite entropy and gain bounds](../full_algebra/certificate.json)

The sign and matrix tests above apply to the specified instrument and frame. They allow a different neutral mode, a different observer reconstruction, or an effective limit with a proven error estimate.

## Forward plus reverse and a common collar correction

The material stretching comparison selects whole-seam work $b$ on the central contact. If both its longitudinal work and lapse energy are identified with the stress-density components on that contact at fixed positive conversion, the even null-source combination $T(k_+,k_+)+T(k_-,k_-)=2(T^{00}+T^{nn})$ is proportional to $b$. This is a separate material identification. At the fixed affine exchange clock the ratio of measured collar focusing to that predicted from ordinary work lies in $[0.9849101,0.9948901]$ for vacancy and $[0.7790982,0.7884013]$ for excess. Both ratios include the uncertainty of the stationary source and the shared volume adjustment. [Material stretch and independent work readings](../directional_stress/traction/README.md) · [Even and odd responses](source_response.json)

Both odd ratios exceed unity and both even ratios lie below unity. Let $H$ be an outward lower bound on the larger observed-to-predicted ratio and $L$ an upper bound on the smaller. For any fitted positive scalar $a$ multiplying the predicted collar source, the worst fractional mismatch $\max_i|r_i/a-1|$ is at least $(H-L)/(H+L)$. For two known endpoints the optimum is $a=(H+L)/2$, where the two errors are equal in magnitude. The certified lower bounds on the irreducible disagreement are 2.7547% for vacancy and 13.2880% for excess. If the same scalar applies to both charges, the bound is 14.2882%. These are bounds for the native clock and the specified stress and curvature assignments.

A stronger absolute-reference test would additionally identify $b$ with all of $T^{00}+T^{nn}$ on the reference state. Both reference values of $b$ are strictly positive. Positive even focusing on both carriers would require $u<-45.5455$ for the vacancy and $u>16.6865$ for the excess. Allowing a common gain of either sign gives $-8.1782<u<-8.1699$ and a negative gain. Additional baseline stress would change this absolute test. [Absolute reference check](source_response.json)

Three reciprocal coherent apertures also test a common clock and a common arbitrary regional source within each charge at one appointment. Their inferred focusing and the unnormalized pixel probabilities depend on a measured phase. A delayed receiver-only test preserves every emission and every free history. Independent phase averaging and inversion of the measured intensity matrix exactly recover the sharp detector table. [Coherent gates, delayed choice and optical calibration](coherent_ports/README.md)

## Numerical enclosures and reproduction

The fixed-point Chebyshev recurrence evaluates only the amplitude columns and receiver neighbourhoods needed for $U$, $HU$ and $H^2U$. It uses the degree, numerator scaling, Bessel enclosures and recurrence error of the stationary source's propagation certificate. The bound on the absolute row sum of $H$ is 32. Each order-$p$ amplitude entry has error below $32^p10^{-19}$, so each four-dimensional block has error below $4\,32^p10^{-19}$. Binomial products give the first two probability derivatives. Each translated integer-bearing coefficient has total absolute weight at most 100, verified for every moment. The resulting moment bounds include the kernel error, the stationary trace-norm radius and, for matched variations, the common compensating source rotation. Each division is evaluated using outward interval arithmetic. [Kernel enclosure](../full_algebra/native_precision.py) · [Two-derivative certificate](certificate.py)

From the repository root:

```sh
PYTHONPATH=code python -m matched_einstein_limit.source_implication.jacobi_screen.measure --write
PYTHONPATH=code python -m matched_einstein_limit.source_implication.jacobi_screen.certificate --write
PYTHONPATH=code python -m matched_einstein_limit.source_implication.jacobi_screen.counts --write
PYTHONPATH=code python -m matched_einstein_limit.source_implication.jacobi_screen.source_response --write
python -m pytest -q code/matched_einstein_limit/source_implication/jacobi_screen
```
