# What the stationary heads do to a passing pilot

On the connected 80-cell board the stored pairing changes how a pilot propagates. Two effects are exactly common to the vacancy and excess. One is the magnitude of a delayed coherent coupling between the two marked contacts. The other is a loss of spatial interference caused by fluctuations in the head state. Their signed responses remain different after calibrating each on ordinary excitation. Integrating out the heads leaves the direct propagation coefficients on **all 956 exterior ports** unchanged at every retarded frequency. Signals observed there can still visit the marked region and return. [Source certificate](results.json) · [Connected-board frequency replay](retarded.json) · [Native pilot counts](counts.json)

## Where the source enters the coherent propagation

Let $H_b$ be the active one-carrier generator on the connected board in left reversal sector $b$. Comparing with the parked source removes its stationary energy in that sector. An overall energy choice for a given carrier can be calibrated on a fixed blank port. For the prepared head state $\chi$ write $M_n=\mathbb E_\chi H^n$, tracing over both coins and summing parity. The coherent retarded propagation for $\operatorname{Im}z>0$ is $G(z)=\mathbb E_\chi(z-H)^{-1}$. An unobserved head also contributes incoherent transmission power, which is measured separately below.

For any port vector $u$, $-\operatorname{Im}\langle u,G(z)u\rangle>0$ if $u\ne0$, by the resolvent identity and $\operatorname{Tr}\chi=1$. Thus $G$ is invertible. Define the optical self-energy $\Sigma(z)=zI-G(z)^{-1}-M_1$. This is an effective frequency-dependent operator on the pilot. Using the isometry $Vu=u\otimes|\sqrt\chi\rangle$ into a purification, $G=V^\dagger(z-H\otimes I)^{-1}V$. Eliminating the orthogonal complement of $V$ gives the usual Schur complement. In particular the anti-Hermitian part obeys $\operatorname{Im}\Sigma(z)\preceq0$ for $\operatorname{Im}z>0$. This construction uses the prepared heads and active propagation. The thermal operator $-\beta^{-1}\log\operatorname{Tr}_{\rm exterior} e^{-\beta H}$ concerns a different equilibrium experiment. [Stationary preparation](../../full_algebra/PROOF.md) · [Thermal regional comparison](../../../../../model/experiments/interval_echo_einstein_station/PAIR.md)

For the four marked endpoints $M=\{0,1,38,39\}$ and exterior $B$, the exact active generator splits into the marked block, the common exterior $H_B\otimes I$ and mark-to-exterior couplings $V_B\otimes I$. This follows from the native law: only the two marked contacts act on the stored heads, and the only dependence on left parity is on the left marked contact. The code verifies every exterior and coupling matrix entry in both carrier and parity sectors. Put $R_B=(z-H_B)^{-1}$ and $\Gamma=V_B R_B V_B^\dagger$. A sixteen-dimensional marked inverse for each parity gives

$$
Q(z)=\mathbb E_\chi [(zI-H_{MM}-\Gamma\otimes I)^{-1}],\qquad
\Sigma_{MM}=zI-\Gamma-Q^{-1}-(M_1)_{MM}.
$$

The remaining blocks of $\Sigma$ vanish. Indeed the mixed blocks of $G$ are $Q V_B R_B$ and $R_B V_B^\dagger Q$; its exterior block is $R_B+R_BV_B^\dagger QV_BR_B$. Inverting reproduces the original exterior and both interface couplings exactly. Thus even a low-frequency pilot can scatter from a source-dependent effective interaction, but the optical potential of this preparation has no tail acting directly on external ports. A pilot on a route dynamically confined to that exterior has no head dependence. A distance inferred from unrestricted global records is a separate question. [Full-board and Schur solves](retarded.py)

At large $|z|$ the first two memory coefficients are

$$
\Sigma(z)=\frac{M_2-M_1^2}{z}
+\frac{M_3-M_1M_2-M_2M_1+M_1^3}{z^2}+O(z^{-3}).
$$

The instantaneous $M_1$ depends only on the parity-resolved left marginal and the right marginal summed over parity. The same holds for $M_2$: the two head-dependent contacts have disjoint position support, so their direct product is zero. Both marginals are preserved by the entire unpairing lottery and its compensating right rotation. Consequently these matrices, and the coefficient of $1/z$, are exactly collar-independent over the full board. This statement also holds after any source-independent change of the scalar energy zero. A phase reference carrying a different bare energy in each parity would constitute a different interferometer.

The earliest path from 0 to 39 is $0\to1\to38\to39$. Write $L=H_{1,0}$ and $R=H_{39,38}$; the connecting bridge hopping is $-I/2$. The third moment on this connection is $-\langle RL\rangle/2$, whereas the $1/z^2$ coefficient is the connected contribution

$$
(\Sigma_2)_{39,0}=-\tfrac12\,[\langle RL\rangle-\langle R\rangle\langle L\rangle].
$$

In each of $M_1M_2,M_2M_1,M_1^3$ the route contributes $-\langle R\rangle\langle L\rangle/2$, which gives the subtraction displayed above. The mark matrices act on different coins, while the parity-dependent left gap is retained. The native third-moment block for excess is $-8|00\rangle\langle11|$ for even parity and $-9|00\rangle\langle11|$ for odd. Vacancy gives the adjoints, also for each of the single-contact blocks. Thus the entire connected coefficient is conjugated between charges for **every** head density. Its magnitude and real part in the calibrated blank-port phase convention are exactly common. The magnitude of the coefficient is invariant even under a source-dependent change of the port phases. [Integer path calculation and stationary intervals](measure.py)

At the blank port frequency with native linewidth one there is also a finite-frequency certificate. For vacancy the real collar derivative of $\Sigma_{39,0}$ is in $[0.1875191,0.1924228]$ and its collar/ordinary ratio is in $[0.2968417,0.2974534]$. For excess the corresponding intervals are $[0.00660359,0.00689865]$ and $[0.2965246,0.2977807]$. Both responses are nonzero and the two ratio intervals overlap. The certificate rounds the sixteen full-board resolvent columns to dyadic integers, computes their residuals exactly, and uses the inverse-norm bound $1/\operatorname{Im}z$. The reduced four-port inverse is enclosed using the stationary trace-norm radius, an exactly inverted rational center and the inverse perturbation identity. Its derivative is $Q^{-1}(dQ)Q^{-1}$ on the off-contact 39–0 entry. The tangent enclosure includes the full survey-volume compensation. [Finite-frequency proof and certificate](finite.py) · [Outward slopes](finite.json)

The exterior terminal exposes a different response even with the same calibration. Launch and receive at blank port 5. At the same unit linewidth the real ordinary return susceptibility increases for both charges. Under unpairing it **increases for vacancy**, by $[7.1973,7.3581]\times10^{-6}$ per unit source change, and **decreases for excess**, by a magnitude in $[0.8753,0.8859]\times 10^{-8}$. Both signs follow from the exact connected resolvent residual and the source enclosure. The compact optical potential explains how an unchanged local operator can have a changed return signal: the exterior propagator takes the pilot to the marked region and back. The real return at one and the same blank port is unaffected by the choice of that port's phase. The real collar/ordinary ratios are $[0.0370916,0.0401422]$ for vacancy and $[-0.000196861,-0.000178438]$ for excess. These real returns cannot describe a common signed scalar after ordinary calibration. The two ratio intervals for the imaginary part overlap. This finite-frequency susceptibility is a model prediction; a previously tested sharp-return tomography identifies its source, while the new direct count experiment below addresses the early-time cubic and dark probability. [Certified exterior transfer](finite.json) · [Independent blank-port state tomography](../../remote_heads/README.md)

The five finite-frequency settings in the replay reference the onsite energy of blank port 5. The two carriers' nominal magnitude ratios range from about 0.2971 to 0.2973 across widths $1/4,1,4$ and detunings $0,\pm1$. A direct 3,840-dimensional solve independently checks the marked inverse. An independent sharp pilot launched and read at blank port 5 also changes its complex return signal through propagation to the marked block and back. The support statement follows exactly from the matrix identities. [Finite settings and residuals](retarded.json)

## How spatial interference separates coherent motion from head noise

On fresh copies alternate equally between pairs $(0,39)$ and $(1,38)$. Prepare the spatial plus superposition and, after free connected evolution, read the orthogonal minus superposition. Blank swap routes connect both pairs without touching either marked seam. A balanced blank splitter, local phase calibration and reversed preparation at readout implement the two basis states. Counts include every miss. On pair $(a,b)$ the leading dark probability is $t^2\langle D_{ab}^\dagger D_{ab}\rangle/4$ with $D_{ab}=H_{aa}-H_{bb}+H_{ab}-H_{ba}$. On the direct bridge the last terms cancel, and the other pair has no direct hopping.

The two differences are $c_b I\pm q[(2+b/4)Z_L+2Z_R]$, with $b=0,1$ and $c_b=-b/4$. Averaging the two squared expressions erases every charge sign. The complete quadratic coefficient has the two head matrices $\mathrm{diag}(4,0,0,4)$ and $\mathrm{diag}(145,1,1,145)/32$. Its reference value is $[3.0025090,3.0025145]$. Subtract the squared means of the two differences, with the same weights, to isolate the variance contribution. Those means are invariant on the collar path, so the **entire collar response of this dark port** is due to head fluctuations. Removing the mean makes the reference coefficient $[3.0017008,3.0017064]$.

This distinguishes a spatial Hamiltonian from the reduced channel. For any smooth unitary on the positions the quadratic probability into an initially orthogonal output is the squared transition matrix element of the instantaneous Hamiltonian. Its time derivative or a later time-dependent phase cannot contribute at that order. The first reduced generator on every spatial state is $-i[M_1,\,\cdot\,]$ and is unchanged by unpairing. A source change implemented solely by a smooth deterministic Hamiltonian on positions with this measured generator would therefore have zero collar derivative of the quadratic dark probability. Here the derivative lies strictly between $-2.227$ and $-2.198$. A metric-only description through that unitary mechanism misses a physical fluctuation channel. This leaves open descriptions that retain the additional degrees of freedom or account for matter-dependent noise.

A sharp local launch supplies a related check. For the earliest 0-to-39 amplitude its leading inclusive power is $\langle K^\dagger K\rangle$, and the coherent power visible against a source-preserving reference is $|\langle K\rangle|^2$, where $K=(H^3)_{39,0}$. Their ratio lies $[6.12,6.19]\times10^{-8}$ for vacancy and $[6.01,6.08]\times10^{-8}$ for excess. Tracing over the heads has observable consequences even with arbitrarily many repetitions of a sharp launch. [Inclusive and coherent powers](results.json)

## Ordinary calibration and the held collar

| Observable at the same stationary source | Held collar slope divided by ordinary slope |
|---|---:|
| Magnitude of connected delayed hopping | [0.2968265, 0.2974299] |
| Incoherent dark-port quadratic coefficient | [0.5314033, 0.5358198] |
| Independently read head entropy | [0.5643237, 0.5645286] |

Each ratio uses matched volume on the original six-unit survey. The common volume-compensating rotation is eliminated in the ratio enclosure. Dividing the delayed coupling by the cube of the measured blank bridge hopping and the dark coefficient by its square eliminates any common choice of time unit. The bridge magnitude $1/2$ is source-independent and attained by both charge probes. Smooth separate recalibrations with nonzero derivative cancel in the slope ratios. Even choosing the held scalar after seeing both pilot readings leaves at least 28.229% relative error on one of the two. This bound follows from $(b-a)/(b+a)$ at the nearest endpoints $a=0.29742990$ and $b=0.5314033$. Including the independently measured entropy raises the minimum to 30.970%. The reduced optical potential is a real physical effect of coupling to the stored degrees of freedom. Identifying it with a universal proper metric or either null contraction of the Einstein tower would require additional agreement. [Exact clock calibration](../velocity/README.md) · [Source-enclosed transfer comparisons](results.json)

The coherent cubic term has a material read. Prepare spatial superpositions of ports 0 and 39 with four relative phases, then record a sharp click at 39 after dwell $t$. The arm starting at the receiver returns with matrix $I+O(t)$; the other arm first arrives at order $t^3$. Subtract opposite launch phases to read the cross term $\langle K_{39,39}^\dagger K_{39,0}\rangle$. Its third coefficient is $i\langle(H^3)_{39,0}\rangle/6$. The same single-hop procedure on each marked contact measures the factors in the connected subtraction. The receiving reference amplitude is chosen proportional to $t^3$ on the remote link and to $t$ on a single-hop measurement.

The replay uses native dwells $1/128,1/64,1/32$ and $2^{74}$ virtual trials per setting, with exact integer counts from binomial batches. The three-dwell estimate removes first and second powers of $t$ in the divided signals. Its simultaneous 99% sampling intervals concern that finite statistic. The zero-dwell coefficients are certified separately from the integer generator. In both species the three-dwell counts give a positive held change for the real connected memory and a negative change for the dark probability. [Complete counts](counts.json) · [Acquisition and reduction](counts.py)

From the repository root,

~~~sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.influence.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.influence.retarded --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.influence.finite --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.influence.counts --write
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/matched_einstein_limit/source_implication/common_front/influence
~~~
