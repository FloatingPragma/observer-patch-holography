# What the native measurements determine about the marked force

The two stationary marks exert a directly measurable force on the regional translation observable of the twelve-cell chart. Eight additional four-tick, single-port experiments recover its full expectation, including its response to an ordinary source and to changes in head pairing. At the original four-tick appointment, two preparations have identical registered stress and identical optical focusing, yet distinct marked force. These statements use the same rooms, marked native interactions, source variations and mobile species. The derived quantity is the force in the finite momentum balance. A geometric reading of that force requires additional metric information. [Regional Ward identity](../translation/README.md) · [Optical focusing](../../jacobi_screen/README.md)

## An observable the four-tick optical screen misses

Write $H=H_0+V$, where $V$ is the native action of the two written contacts relative to blank. The chart gives a partial translation $S$ and $Q=(S-S^\dagger)/(2i)$. With the convention that momentum flux enters the chart, the exact rate is $\dot Q=\sum_e J^Q_e+F$, with $F=i[V,Q]$. Each cell separately obeys native energy continuity. The propagation uses the degree reference and the active squared contrast at its registered strength. [Twelve-cell construction and force certificate](../translation/certificate.json)

Launch the original ten tagged inputs on the two rings with equal weights. After four native ticks, make one of two fair, recorded phase interventions. The first applies either no pulse or a $\pi$ pulse to both port 0 and port 229. The second applies a $\pi$ pulse to one of those ports, chosen fairly. Both controls have the same one-port pulse statistics. For every incoming density operator $\rho$, the difference of the two averaged output states is
$$
2(P_0\rho P_{229}+P_{229}\rho P_0).
$$
The projectors act on the mobile port and leave the heads intact. Every port remains enrolled, so this identity applies to each input tag, head state, charge and reversal parity.

The two ports occupy cells two steps apart, and the shortest native path connecting them has nine seams. Consequently every effect supported in one cell or two adjacent cells has the same mean for the two interventions. This includes each registered cell energy, bridge work and energy flux. Port 0 is a nongateway and port 229 a gateway. The registered translation preserves gateway type, and each translation flux connects gateways. Every such regional or cellular shift and each boundary flux therefore also has the same mean.

For a local arrival effect $E_y$ at any port $y$, including a head analyzer there, expand $\operatorname{ad}_H^k(E_y)$ in terms of $H^\ell E_y H^{k-\ell}$. An entry connecting 0 to 229 requires $\ell+(k-\ell)\ge d(0,y)+d(y,229)\ge9$. All arrival probabilities and their first eight right-time derivatives agree across the controls. The forward and reverse normalized screen images are nonsingular at the reference for both species. Thus their Jacobi focusing, expansion and source derivatives agree wherever the images remain nonsingular. The verifier also propagates every source and receiver pair and compares the hit matrices and their first two time derivatives before contraction with either head. [Support proof and propagation](measure.py) · [Singular values and complete intervals](results.json)

In contrast, $F$ connects ports 0 and 229. The following outward-rounded bounds are differences between the two phase interventions. The source derivatives include the registered right-head rotation that holds volume fixed and the certified uncertainty of the reference source.

| Mobile charge | Reference force contrast | Ordinary derivative | Pairing derivative |
|---|---:|---:|---:|
| Vacancy | [4.00507, 4.00512] × 10⁻⁷ | [−8.24017, −8.12286] × 10⁻⁵ | [−2.46933, −2.41926] × 10⁻⁵ |
| Excess | [1.84296, 1.84298] × 10⁻⁶ | [5.43401, 5.51248] × 10⁻⁴ | [1.60731, 1.64077] × 10⁻⁴ |

The direct force check uses pair interference and head Pauli measurements. Independent simulated binomial counts on the two arms, with $2^{76}$ trials per setting and a simultaneous 99% Hoeffding bound, give vacancy slopes $(-8.187\pm0.103,-2.447\pm0.031)\times10^{-5}$ and excess slopes $(5.473\pm0.011,1.625\pm0.004)\times10^{-4}$ at the nominal preparations. These are simulated readings of the finite model. A simpler diagnostic prepares a fresh real superposition on the two ports with either sign. Its force contrast lies between −0.444290 and −0.444285 on the same stationary heads. [Count records, bounds and fresh-preparation comparison](results.json)

Any source-universal numerical connection determined by those screen records takes identical values under both interventions. Every registered instantaneous stress expectation entering its proposed contraction also agrees. Such a contraction must assign the same force, contradicting the table. This identifies exactly which records are insufficient. A connection measured by another instrument, stress components with longer spatial support, joint correlations between observables, and a limit theorem restricted to a specified class of smooth preparations are separate questions.

## Translating both written contacts

In the 99-state primitive Hamiltonian the blank/written off-diagonal blocks vanish. A seam term acting on a different memory commutes with the first memory's blank/written projector, including when the two seams share a load. Each address projector $P_A$ therefore commutes with the free Hamiltonian. For an arbitrary candidate momentum contribution $K$ on all address sectors,
$$
P_A\,i[H,K]\,P_A=i[H_A,P_A K P_A].
$$
In each sector $h=2H_A+cI$ and $f=4F$ have integral entries. The left and right marked force contributions separately have zero trace against powers $h^j$ for $0\le j<7$ and each gives $+32$ against $h^7$. Thus $\operatorname{Tr}(H_A^7 F)=+1/8$, split equally between the marks. Every commutator on the right of the compression identity has zero trace against $H_A^7$. Adding a momentum on the already available marked-address sectors cannot absorb this force while keeping the measured boundary flux on every original state. [Integer tests for both individual marks](placements.py) · [Sector results](placements.json)

There is a distinct exact translation comparison. Prepare the original address pair $A=(0,1;38,39)$ on one copy and the translated pair $B=(228,229;218,219)$ on another. Assemble the full material Hamiltonian of each independently at the same energy origin. With $P=S^\dagger S$ and $P'=SS^\dagger$, both marks moving gives $P'H_BP'S=SPH_AP$. The squared Hilbert–Schmidt residuals for the doubled Hamiltonians are 7840 and 8112 if both addresses stay in place, 3920 in both parities if only the left moves, and 3920 and 4192 if only the right moves. Moving both gives zero in all four charge/parity sectors. Thus
$$
H_B S-S H_A=(I-P')H_BS-SH_A(I-P).
$$
The corresponding off-diagonal address operator has a boundary-only Ward identity. It has zero expectation on every address-diagonal input. A coherent comparison would need to prepare and recombine the address branches. For paired launch amplitudes $w_B(0)=S w_A(0)$, the simulated cross amplitude $\langle w_B(t),S w_A(t)\rangle$, averaged over the original launch and head ensemble, starts at one. At four ticks its real part is 0.999806082 (vacancy) and 0.997812135 (excess). Independent neighboring flight appointments reproduce the real boundary slopes $+0.000138376$ and $-0.000948152$. This calculation establishes joint translation covariance between the two native preparations. A recoil of the stationary addresses during one free flight would require dynamics that change their addresses. [Cross-branch measurements and boundary rates](placements.json)

## Recovering the complete material force from later local records

All eight unordered off-diagonal port pairs in $F$ join a marked port $x$ to the chart image $f(y)$ of a port on the same marked seam. The two seams are $(0,1)$ and $(38,39)$, with images $(228,229)$ and $(218,219)$. The integer verifier sums all eight blocks and their adjoints and compares the result to the full operator in each charge/parity sector. For each pair let $z=f(y)$. At the four-tick appointment enroll the two-port subspace without distinguishing its ports, recording all outside outcomes as zero. Apply the correlated and anticorrelated phase controls with a relative phase $\theta$ on $z$. Allow an additional four ticks of native free flight. Read a hit at the single local port $x$ and measure a product of left and right head Pauli operators, including the identity. Misses count as zero. This preserves unconditional normalization across trials.

For a density block $X=\rho_{xz}$ at the original appointment, set $A=U(4)_{xx}$ and $B=U(4)_{xz}$. After the two phase controls, the difference of the unnormalized detected head-density matrices is
$$
D_\theta=2\big(e^{-i\theta}AXB^\dagger+e^{i\theta}BX^\dagger A^\dagger\big),\qquad
X=\tfrac14 A^{-1}(D_0+iD_{\pi/2})(B^\dagger)^{-1}.
$$
For both species and parities, rational fixed-point propagation and an inverse-norm bound certify $A,B$ invertible on every measured pair. On the primary pair the earliest native transfer block $(2H)^9_{0,229}$ already has determinant $2^{20}$ in every sector. A separate calibration uses each of the four preparations $|0\rangle,|1\rangle,|+\rangle,|{+i}\rangle$ for each of the port rail and two head factors. Product Pauli measurement at the receiving port on these 64 product inputs identifies the map $[A\ B]$ up to one common phase. That phase cancels in the formula for $X$. This calibration uses source-independent transfer probabilities and head records. The two-phase read then determines $2\operatorname{Re}\operatorname{Tr}(F_{zx}X)$, the contribution of this pair to the original force. Summing it over the eight pairs recovers the full $F$ on any joint state in the sector. [Transfer calibration and local inversion](late_read.py)

| Mobile charge | Recovered original force | Ordinary derivative | Pairing derivative |
|---|---:|---:|---:|
| Vacancy | −1.7246852 × 10⁻⁶ | −1.8144538 × 10⁻⁴ | −5.2923299 × 10⁻⁵ |
| Excess | −1.0168160 × 10⁻⁵ | +2.0663503 × 10⁻⁴ | +6.1800493 × 10⁻⁵ |

The reconstructed values lie inside the independently certified intervals in the regional Ward experiment. This delayed read uses simulated population probabilities. The largest Frobenius norms of the reconstruction weights are $7.80\times10^7$ and $1.32\times10^8$, so estimating the force from finite late-hit samples needs a separate precision analysis. Independent direct force counts are available above. The delayed protocol recovers material coherence that the original optical experiment misses. Identifying a proper metric and a connection that couples to that material on shrinking common regions remains a physical premise of the Einstein tower. [All local readouts, weight norms and reconstructed forces](late_read.json) · [Independent regional certificate](../translation/certificate.json)

## Testing the force against sharp optical focusing

Return to the original source and ten equally weighted four-tick launches. All propagation, marks, charge and reversal conventions are identical in the force and image experiments. The force $M=\langle i[V,Q]\rangle$ is unconditional in the launch tag. The forward and reverse optical images sort the same ensemble by launch tag at the receiver. Each $G=-\operatorname{tr}(\ddot J J^{-1})$ differentiates the transmission covariance at a *fixed* pair of screens with respect to exchange dwell. Its units are inverse exchange-time squared. The extensive $M$ is a directed shift derivative per exchange time. The shift moves central cell 0 to 19, whereas the central bridge joins cells 0 and 3. The chart supplies no identity of their directions, transverse area, or normalization. The scalar rule $\delta G=c\,\delta M$ therefore tests an explicit identification. The native Ward identity for $\dot Q$ makes no such prediction.

The two derivatives use the original stationary heads. Ordinary variation equalizes the reversal populations. The held variation removes part of the head pairing. Each includes its independently surveyed right-head rotation. The twelve-range tetrahedron measured with the six-tick survey has the same volume to first order in either variation. Four ticks specifies the separate force and image appointment. No step identifies the tetrahedron with a proper ball at an arrival event. [Original source and volume matching](../../full_algebra/README.md) · [Four-tick images and time conventions](../../jacobi_screen/README.md)

For a read $Y$, write $\delta_oY=Y_o^0+\alpha_oY_r$ and $\delta_pY=Y_p^0+\alpha_pY_r$, where $Y_r$ is the common rotation response. Both instruments use the same $\rho=\alpha_o/\alpha_p$. The combination $D_Y=\delta_oY-\rho\delta_pY=Y_o^0-\rho Y_p^0$ therefore removes the large rotation exactly. At native affine exchange dwell the ordinary data fix $c=\delta_o G/\delta_o M$. The held comparison is
$$
\mathcal H=\frac{\delta_pG}{c\,\delta_pM}
  =1+\frac{D_M/\delta_pM-D_G/\delta_pG}{\rho+D_G/\delta_pG}.
$$
To enclose it, contract fixed-point native propagation with the certified source, retain the common angle across both records, and differentiate the measured centroid, throughput and both image axes. All intervals include the stationary-source and propagation errors. The signed force and its two derivatives also agree with the independent delayed measurement above. A constant rescaling of either meter leaves $\mathcal H$ unchanged. [Calculation](sharp_focus.py) · [Full intervals and calibration tests](sharp_focus.json)

| Sharp focusing statistic | Vacancy, held / predicted | Excess, held / predicted |
|---|---:|---:|
| Forward | [0.94908, 0.95362] | [0.99888, 0.99910] |
| Reverse | [1.03946, 1.04084] | [0.98063, 0.98151] |
| Forward minus reverse | [1.01828, 1.01897] | [0.99317, 0.99346] |
| Forward plus reverse | [1.07954, 1.08446] | [1.01374, 1.01518] |

A separate numerical check propagates fresh source columns from each end by sparse matrix exponentiation. It measures the normalized image at both signs of a $2^{-22}$ source perturbation and differences the two outputs. All tested head states are positive. All eight resulting ratios, using the independently reconstructed delayed material force, fall inside the fixed-point intervals in the table. [Independent signed-flight regression](test_connection.py)

Every interval excludes unity, even with an ordinary gain fitted separately to each charge and direction. The closest, forward excess, differs from the force prediction by at least 0.0910%. To resolve that case, the enclosure covers the common pairing rotation and the ratio of the two rotations with 16 intervals apiece, retains the same rectangle for the optical and material records, then takes the union of all 256 bounds. The ordinary readings also require different gains across the two species. For example, odd focusing per unit marked force is between $-56.01$ and $-53.29$ million for vacancy and between $-66.40$ and $-63.53$ million for excess. These are conversion factors between the two specified native instruments. [Unrounded bounds and signed predictions](sharp_focus.json)

The two coefficients in a separate test *are* identifiable. Fit $\delta G_{\rm odd}=a\,\delta j+b\,\delta M$ to the two species' ordinary responses, with $j$ the independently metered central energy current. For even focusing replace $j$ by whole-seam work. The determinant of each ordinary two-by-two meter matrix excludes zero. Retain both ordinary coefficients on the held preparation. In units of each charge's signed ordinary focusing, the measured-minus-predicted collar residuals are $[0.00246,0.00599]$ and $[-0.00261,-0.0000255]$ for current plus force, and $[-0.02396,-0.00772]$ and $[-0.04392,-0.02940]$ for work plus force. The earlier phase interventions give an independent structural check. Their material force changes while the optical response and both contact records agree, so a constant additive force coefficient valid on every controlled preparation would have to vanish.

Under a source-independent change of clock, let $u=\ddot\lambda/\dot\lambda$ at the appointment. Clock-normalized focusing is $(G+u\theta)/\dot\lambda^2$, where $\theta=\operatorname{tr}(\dot J J^{-1})$. An arbitrary signed gain absorbs the common scale. For any one charge, ordinary and held proportionality then require
$$
(D_G-m\,\delta_pG)+u(D_\theta-m\,\delta_p\theta)=0,\qquad m=D_M/\delta_pM .
$$
Even allowing the gain to differ by charge, odd focusing requires $u\in[-8.706,-5.843]$ for vacancy and $u\in[28.737,139.940]$ for excess. Even focusing requires $[30.102,40.905]$ and $[4.625,9.810]$. The single forward and reverse directions also give disjoint pairs. A single clock acceleration across the two tracers cannot make the four-tick material force proportional to these optical changes. Source-dependent clocks or units introduce additional derivatives and would require independent measurement.

In a metric theory a coordinate force from the connection and optical tidal focusing are physically different measurements. The first can involve the first derivatives of the metric, while the second involves curvature. Local inertial coordinates can remove the connection at an event with nonzero curvature. Einstein's null equation compares proper curvature to a null contraction of a stress tensor. The present force experiment establishes a native material contribution to translation balance. A physical comparison to Einstein's equation still requires independently identifying a proper optical frame and the stress contraction on the same region. [Regional translation comparison](../translation/README.md) · [Jacobi screen and clock interpretation](../../jacobi_screen/README.md)

## Equal radar ranges with a different measurable material force

The twelve ranges and both marked-force measurements together have rank fourteen on the eighteen correlation coordinates. A second, independently constructed source pair keeps all fourteen means exactly fixed while entropy increases. [Equal ranges and equal force with increasing entropy](EQUAL_FORCE_ENTROPY.md)

Start the same 80-cell hardware at the volume-matched stationary state of repeated, uniform use of all 24 radar settings. The twelve registered six-unit mean ranges average the two probe charges equally. On separate copies use the same four-tick, ten-launch force experiment described above, resolving the force probe's charge. The four-tick force is a linear function of the incoming pair of unnormalized head states. This allows an exact test along the affine source path $\chi(t)=\chi_*+tD$. The stationary state is the center of this path. Neighboring states are physically preparable perturbations; stationarity of the neighbors under the reuse channel is a separate property. [Stationary source and exact range kernel](../../full_algebra/reuse/INVISIBLE.md)

The source experiment selected $D$ by projecting the entropy gradient into the kernel of the twelve directed range effects and correcting for the certified propagation error. It preserves both parity weights and every one-head marginal in both blocks. No force measurements enter its selection. Every state on the certified segment has exactly the same twelve pooled mean ranges, hence the same tetrahedral area and volume computed from them. The force can still change because it reads joint head correlations in the transported probe.

Each parity block has nine orthonormal correlation coordinates $B_i=\sigma_a\otimes\sigma_b/2$. Compress the full operator $i[V,Q]$ back through the four-tick flight and the ten equiprobable launches to obtain a four-by-four force effect $C_b$ on each incoming head pair. Form the 18-vector $m$ with entries $\operatorname{Tr}(C_b B_i)$. The signed force slope in direction $D$ is $m\cdot D$. Let $d$ be the serialized, force-independent center of this direction and $\epsilon_D$ its Hilbert–Schmidt error, including serialization. With per-block force-effect error $\epsilon_C=2\times10^{-15}$, an outward error bound for the computed slope is
$$
\|m\|_2\,\epsilon_D+\sqrt8\,\epsilon_C\,(\|d\|_2+\epsilon_D).
$$
The implementation uses the slightly larger factor 3 in place of $\sqrt8$. The positive slopes below hold for the exact range-kernel vector, with the original sign convention that the marked force adds to incoming momentum flux. [Force-effect contraction and native local check](equal_radar.py) · [Signed intervals](equal_radar.json)

| Force read | Exact force slope per unit displacement | Force increase from $0$ to $10^{-7}$ |
|---|---:|---:|
| Vacancy | [1.2543014, 1.2546353] × 10⁻⁸ | [1.2543014, 1.2546353] × 10⁻¹⁵ |
| Excess | [4.8911021, 4.8912701] × 10⁻⁷ | [4.8911021, 4.8912701] × 10⁻¹⁴ |

At that smaller step, the source certificate also proves a strictly positive increase in entropy. Between the two native endpoints $t=\pm0.01$, all twelve ranges still agree and the force difference is exactly 0.02 times the tabulated slope. This establishes a finite measurement that distinguishes the sources despite identical radar summaries. The conclusions apply to the registered charge-pooled means. Charge-resolved arrival histories, a proper observer metric, null energy and astrophysical observations require their own measurements.

An independent local simulation prepares both published conditional fixed-point centers at midpoint odd weight 0.1343, perturbed by $\pm0.01d$. Both arms satisfy the native preparation inequalities. The four-tick flight followed by the calibrated four-tick local read reconstructs all eight force contributions from two phase quadratures at single-port detectors. The measured difference divided by 0.02 lies in each certified slope interval. These positive calibration arms use the numerical centers. The conclusion for the *exact* stationary chord follows from linearity of the force expectation and the explicit correction bound above. The receipt also gives an enclosure of the absolute force at the exact stationary state over the entire matched parity-weight bracket. [Local arm and parity records](equal_radar.json)

The direct binary force meter returns a signed outcome of magnitude 35. On the nominal finite endpoint pair, $2^{82}$ simulated trials per arm yield differences of $2.8082\times10^{-10}$ and $9.7796\times10^{-9}$ for vacancy and excess. The simultaneous 99% Hoeffding half-width is below $1.1660\times10^{-10}$ for both signed differences. Twice that half-width lies below even the smaller certified exact endpoint separation, so the same shot budget also suffices to resolve the sign on the exact endpoint pair with at least 99% simultaneous probability. This trial count is a property of the specified binary meter and estimator. The calibrated delayed local readings above are simulated population probabilities. [Trial counts and instrument scales](equal_radar.json)

Run the regression checks from the repository root with `OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/directional_stress/connection`. The modules `measure`, `placements`, `late_read`, `sharp_focus` and `equal_radar` in `matched_einstein_limit.source_implication.directional_stress.connection` each accept `--write` to regenerate their corresponding JSON results.
