# Three phase-sensitive readings of one native source

Three native probes with the same incident density, spatial intensity and receiving aperture give incompatible values of the clock acceleration required for a common directional focusing response. This holds for vacancy alone and excess alone at one clock appointment. A separate matrix test rules out a single optical tidal equation for the three images even if the common receiving frame can shear, stretch and rotate during the observation. The phase response is visible in raw click probabilities. [Native observations](results.json) · [Full source enclosures](certificate.json)

## Equal incident state and equal local current

Use the five launch and five receiving ports of the parent experiment. At each end the first and second port share a blank edge, $0\leftrightarrow5$ on the left, $39\leftrightarrow40$ on the right. The local two-port blank Hamiltonian gives the gate

$$
G_s=\frac15\begin{pmatrix}4&3is\\3is&4\end{pmatrix},\qquad s=\pm1.
$$

A blank dimer calibrated to transfer probability $9/25$ gives the positive phase at dwell $2\arctan(3/4)$. The complementary appointment on the other branch gives the negative phase. The scalar phase common to the two-port block drops out of every sharp input and detector probability here. The other three ports are untouched. This setting was chosen with the sharp-aperture results in view. Plus and minus denote pulse settings and are independent of the carrier's charge. [Preparation and dimer](measure.py)

With a uniform, recorded choice of launch port, the incident spatial density is $G_s(I_5/5)G_s^\dagger=I_5/5$. The gate acts as identity on both stationary head blocks. After summing over the external launch-pixel label, the complete propagating density is exactly the same for every setting, including every source perturbation. The correlation with that external label can differ. The free Hamiltonian, clocks, routes and source stay fixed. Every measurement of regional stress made prior to the final gate, using repetitions from the same uniformly tagged ensemble, therefore has the same probability law for every setting. This identity holds for any regional stress observable. The complete receiver aperture remains the same projection, so total transmission is identical at every dwell. With the same pulse magnitude, each tagged launch also has the same position probabilities for both pulse phases. The first-harmonic spatial calibration $2X^T|G_s|^2X/5$ has positive eigenvalues for both signs, where $X$ is the matrix of five marked pentagon bearings.

Both ends use the measured gates. For a particular forward setting the amplitude is $G_R U G_L$. For the reverse trial the apparatus applies the same gates at the same physical ends and the amplitude is its transpose, since the native Hamiltonian and each gate are symmetric. The experiment retains separate forward and reverse directions. Four separately observed settings use the sign pairs $(++),(--),(+-),(-+)$. The first three determine the exclusion. A release with the exchange contrast removed supplies a comparison on the same apertures. [All phase settings](results.json)

## Clicks before optical inference

For a vacancy released at the first forward pixel, the probability of a click at the first receiver pixel is 0.00245382 at $(++)$ and 0.00210062 at $(--)$. For an excess the corresponding values are 0.00275539 and 0.00267969. At $2^{38}$ virtual shots the vacancy settings give 674,498,147 and 577,389,867 clicks. The complementary counts, the other settings, and the reverse direction are retained. Uniform 99% binomial bounds distinguish the two phases in both directions and both charge sectors. These are simulated instrument records.

For each charge, every unnormalized hit or miss probability is linear in the incident source density. Randomizing the two gate signs independently and discarding them yields the classical spatial mixing $|G|^2 P |G|^2$ on the sharp probability table $P$. The code verifies the identity for each individual hit and both time derivatives. The coherent settings retain their phase labels and are evaluated separately. Their distinct clicks therefore precede any centroid normalization, trace extraction or averaging across species. [Phase and dephasing control](measure.py) · [Synthetic counts](results.json)

There is an exact hardware correction for the phase settings. Each measured $2\times 2$ intensity block of $Q=|G|^2$ has eigenvalues $1$ and $7/25$. Thus $P=Q^{-1}\bar P Q^{-1}$ reconstructs the sharp probability table exactly from the four independently phase-randomized settings. This also holds after differentiating in dwell or in the state. The inversion amplifies spectral-norm errors by at most $(25/7)^2$. Applied to the native records it restores the original sharp image, including its established held-out mismatch with the bridge current at the native affine-clock interpretation. This correction removes an actual instrument effect without fitting a gravitational coupling. [Reconstructed sharp response](results.json) · [Comparison at the sharp aperture](../source_response.json)

## An area and clock cannot calibrate these three images

For each setting, measure its angular response matrix and its first two exchange-dwell derivatives at 4 on the selected stationary source. Write $s$ for forward minus reverse expansion and $f$ for forward minus reverse trace focusing in this clock. At any one member of the source family a common clock acceleration $u=\ddot\lambda/\dot\lambda$ and a common physical conversion of stress require $f_i+u s_i=h$ for every setting $i$. The right side $h$ is arbitrary. It can include any chosen coupling, area, proper time or unmeasured regional stress, since the same unconditional source is used for the three settings. The points $(s_i,f_i)$ must therefore be collinear.

| Source | $u$ required by $(++)$ and $(--)$ | $u$ required by $(++)$ and $(+-)$ | Three-point determinant |
| --- | ---: | ---: | ---: |
| Vacancy | [−47.816, −43.459] | [−78.845, −67.884] | [1016.44, 1705.03] |
| Excess | [29.302, 31.743] | [36.936, 37.128] | [9.620, 15.072] |

The determinant is $(s_2-s_1)(f_3-f_1)-(s_3-s_1)(f_2-f_1)$. Its strict sign excludes a common value even if $h$ and $u$ differ between sources or between the two charges. The test makes no use of pressure, a stress decomposition, an ordinary-matter fit, or even the value of the bridge current. Its physical premise is that these apertures measure the same pair of opposite null directions at the same effective event. When that interpretation is applied to the native images it fails. The result allows probe-dependent paths, diffraction and charge-dependent propagation. [Outward determinant and slope enclosures](certificate.json)

## A delayed choice at the receiver

In a separate experiment every emission is sharp. No gate acts at the sending end. After the free-flight interval, the receiver applies two independently phased splitters on disjoint blank edges of its pentagon before reading the pixel. The left receiver pairs are $(0,5),(9,8)$ and the right receiver pairs are $(39,40),(47,46)$. Timing here refers to the end of free transmission. The final read can be performed in a separate isolated receiver operation. The choice leaves the whole propagation history unchanged even conditional on the emitter's direction and its saved launch pixel. For each individual launch, the total hit probability and its complement are also unchanged.

For the three receiver settings $(--),(+-),(-+)$ the direction-odd determinant is in $[0.8112,1.4250]$ for the vacancy and $[538.86,558.51]$ for the excess. The four possible minors are reported. The other three vacancy enclosures include zero. The independent moving-screen rank test below has a nonsingular minor for each of the two charges, two directions and two acceleration components. Its largest preconditioned error is below 0.285. This localizes the effect to phase-sensitive analysis of the outgoing beam with a fixed source. The exact intensity correction for the delayed experiment acts only on the left of the hit table, $P=Q^{-1}\bar P$. [Delayed-choice data and enclosures](delayed.json) · [Reproduction](delayed.py)

## Removing the choice of a moving screen

If true transverse displacements satisfy a common Jacobi equation, a common change of affine parameter and any twice-differentiable invertible recalibration of the receiver screen give, at each appointment,

$$
\ddot J_i=A\dot J_i+B J_i
$$

for common $2\times2$ matrices $A,B$. This follows by differentiating $J_i=L(t)D_i(\lambda(t))R_i$ twice and using $D_{i,\lambda\lambda}=-K D_i$. The independent invertible $R_i$ are fixed at launch. This necessary condition also allows arbitrary initial displacements and opening rates. It can be tested for each charge and propagation direction separately before inferring any curvature trace or dividing by an image determinant.

For each of two transverse output components and each of the six columns (three physical settings, two bias axes) form a row $(\dot J_{1j},\dot J_{2j},J_{1j},J_{2j},\ddot J_{kj})$. A common equation requires rank at most four. The certificate checks five-row minors. A rounded dyadic preconditioner $P$ computed from an approximate inverse obeys $\|I-PM\|_\infty<1$ in outward arithmetic on the indicated minor $M$. The Neumann series then proves $M$ invertible on the full source enclosure. All eight charge, direction and component combinations have a nonsingular minor. The largest certified upper bound on $\|I-PM\|_\infty$ is below 0.166. A planted common equation passes the consistency condition and changing one acceleration row gives a nonzero minor. [Rank certificate](certificate.py) · [Positive control and independent gate checks](test_coherent.py)

These measurements use fixed hardware positions and vary receiver dwell. Identifying those samples with proper events along a null ray requires a separate propagation, scale and screen calibration. If the phase settings select different effective rays or probe phase-sensitive wave structure, the common-ray hypothesis does not apply. The conflict concerns this proposed readout of curvature. A theorem about gravity in an effective macroscopic limit requires a separately specified coarse-graining or neutral probe and quantitative control of its errors. [Geometric optics and data for the sharp screens](../README.md)

Run from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.jacobi_screen.coherent_ports.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.jacobi_screen.coherent_ports.certificate --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.jacobi_screen.coherent_ports.delayed --write
python -m pytest -q code/matched_einstein_limit/source_implication/jacobi_screen/coherent_ports
```
