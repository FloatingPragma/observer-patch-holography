# Does the collar delay the phase of a remote signal?

At the stationary 80-cell source a coherent vacancy or excess released from ports 0 and 5 gives a source-sensitive interference signal at ports 38 and 39. The phase of this signal changes with both ordinary excitation and the removal of collar pairing. A receiver analyzing the same position counts at different frequencies infers incompatible signed delay changes after ordinary calibration. Allowing a common change of time origin and rate still excludes a single scalar response, separately for each charge. [Source-enclosed results](results.json) · [Synthetic acquisition counts](counts.json)

## What a spectral line can establish

The two interventions change the initial stationary head density. All flight generators, their spatial contacts and the free clock remain fixed. For fixed emitter and receiver operations and one scanned free dwell $t$, insert the spectral resolution $e^{-itH}=\sum_E e^{-itE}P_E$ into the probability of any complete outcome. Each signal has the form

$$
f_s(t)=\sum_\nu c_\nu(s)e^{-i\nu t},\qquad \nu\in\{E-E'\}.
$$

Multiple fixed-duration operations surrounding that interval merely change the coefficients. A finite sequence of scanned dwell lengths $a_jt$ with fixed operations between them has frequencies drawn from finite sums of the $a_j\nu_j$. Fixed finite banks of independent preparations admit the same argument on a product Hilbert space. The Laplace transform and its source derivative have at most simple poles in this same finite set. Moving an occupied pole at first order would produce a double pole in the source derivative. Thus a shifted occupied microscopic resonance or a dilation of a nonzero occupied frequency cannot be driven by either change of head preparation alone. A signed fringe may lose a line through cancellation. The apparent maximum of unresolved lines can also move as their relative weights change.

There is an exact statement about which energies transmit. If $K_E=P_yP_EP_x$ is the compressed spectral projector from launch $x$ to receiver $y$, the inclusive line power is $\operatorname{Tr}(\chi K_E^\dagger K_E)$. Faithfulness on the enrolled head space makes this positive precisely when $K_E\ne0$, independently of the ordinary or collar preparation. At a single marked dimer both material charges exhibit the native angular frequencies 8 and 9 in the two reversal sectors. The contact matrices obey $H_b^2=\nu_b H_b$ exactly; a sharp release has far-port population $\sin^2(\nu_b t/2)$ times a nonzero, rank-one head weight. Changing the squared-contrast coefficient from 8 to 7 shifts the two frequencies to 7 and 8. [Integer-projector and changed-generator control](test_fringes.py) · [Stationary-source faithfulness](../velocity/certificate.json)

For an exponentially windowed, nonzero, continuous-time finite signal $R_s(\omega)=\int_0^\infty e^{(i\omega-\eta)t}f_s(t)\,dt$, $\eta>0$, the transform is rational. Its group delay $\operatorname{Im} R_s'/R_s$ is real rational wherever the phase is defined and tends to zero at large real frequency. Two preparations therefore cannot differ by an exact nonzero constant group delay on an open band. For an infinitesimal change the same conclusion follows if the first nonzero coefficient of the large-frequency expansion stays nonzero with the same order nearby. Positive remote arrival intensities on this faithful source satisfy that condition because their shortest-path coefficients are strictly positive. The sampled apparatus below admits a separate tangent result even for signed fringes. These finite statements allow frequency-dependent changes and do not prescribe an emergent dispersion relation of an interacting thermodynamic limit. [Earliest nonzero arrivals](../README.md)

## How the receiver measures a phase delay

Use all native seams during free flight. The adjacent unmarked sending ports 0 and 5 form the launch pair. On fresh heads, release the balanced superpositions $(|0\rangle+e^{i\alpha}|5\rangle)/\sqrt2$ for $\alpha=0,\pi/2,\pi,3\pi/2$. Blank-seam exchange supplies the balanced splitter, and a relative native phase supplies the quadratures. At each exchange dwell $t_n=n/2$, $1\le n\le12$, read one of the three positions 38, 39 or elsewhere. Both receivers see the same experiment. No terminal head state is selected. [Native rail phases](../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/phase_geometry/measure.py) · [Free transfer kernels](../../full_algebra/native_precision.py)

Write $K_{yx,b}(t)=\langle y|e^{-itH_b}|x\rangle$ for the four-coin transmission block at fixed charge and left parity. Subtracting the opposite launch phases yields

$$
C_y(t)=\tfrac12(p_0-p_\pi)+\tfrac i2(p_{\pi/2}-p_{3\pi/2})
      =\sum_b\operatorname{Tr}\!\left(\chi_b K_{y5,b}^\dagger(t)K_{y0,b}(t)\right).
$$

Free evolution before and after an external timing offset must use the same convention in both arms. Here the two launch sites evolve simultaneously for exactly the same duration, so an unobservable uniform energy shift in a parity block cancels. The receiver can compute the following quantities from the counts alone,

$$
M_k(\omega)=\tfrac1{12}\sum_{n=1}^{12} t_n^k e^{i\omega t_n} C_y(t_n),
\qquad
\tau_y=\partial_\omega\arg M_0=\operatorname{Re}(M_1/M_0).
$$

This is the spectral group delay of a measured position-interference fringe over the specified observation window. Interference permits it to lie outside that window. Its interpretation as a null flight time would need an independent propagation argument. The three settings use phase advances $1,i,(3+4i)/5$ per half tick. Every registered reference visibility is bounded away from zero. The source derivatives follow by the ordinary quotient rule on the complex moments. In the following table one endpoint and one carrier preparation are fixed within each row. Each change of column reuses the same position counts. [Outward rational-kernel and interval calculations](measure.py)

| Receiver 38, collar delay slope divided by ordinary slope | Zero phase advance | Quarter turn | Phase advance $(3+4i)/5$ |
|---|---:|---:|---:|
| Vacancy | [0.3002951, 0.3049158] | [0.5212251, 0.5733160] | [0.3229990, 0.3243727] |
| Excess | [0.4376064, 0.7801875] | [0.4381409, 0.4607514] | [0.1415073, 0.1494889] |

At zero advance, the vacancy ordinary response in native clock units is in $[-0.538,-0.292]$ and its excess counterpart is in $[0.051,0.173]$. Both the frequency dependence and this sign difference occur on the same registered source. More generally a finite clock receiver has $M_0=P(e^{i\omega/2})$ for a complex polynomial $P$. Its first source-delay change is the frequency derivative of the single-valued periodic rational function $\operatorname{Im}[(\partial_sP)/P]$. If this derivative is constant on an open frequency interval with visibility, rational continuation makes it constant wherever defined. A periodic rational function with constant derivative has derivative zero: real poles of that function would also be poles of the derivative, and its integral over a period vanishes once any removable points are filled. This excludes an exactly achromatic, nonzero infinitesimal delay from reweighting these fixed bins.

## Can a common change of clock restore agreement?

Write $O_i,P_i$ for the ordinary and collar signed delay slopes in channel $i$. For a single common scalar source response, allow a separate static gain in every channel. Both the common clock origin and rate may change independently along each intervention. With $t_{\rm physical}=b+a t$ at fixed physical analysis frequency, the phase transforms by $M_{\rm physical}=e^{i\omega b}M(a\omega)$ and the reference delay by $b+a\tau(a\omega)$. At $a=1,b=0$ this adds an arbitrary common offset and an arbitrary multiple of $\tau_i+\omega_i\partial_\omega\tau_i$. A receiver holding phase advance per native tick instead adds the rate term proportional to $\tau_i$. To cover both choices, grant independent coefficients to $\tau_i$ and $W_i=\omega_i\partial_\omega\tau_i$. Eliminating each channel's static gain by its ordinary calibration requires **some** constants $r,A,B,C$ satisfying

$$
P_i=r O_i+A+B\tau_i+C W_i
$$

for all channels of a given charge, provided the latent ordinary scalar has nonzero response. If the latent ordinary slope vanishes, the ordinary column itself must already be a combination of the three clock columns; the same rank test excludes that case. The original volume control uses the same right rotation in both arms, with $\rho=\alpha_o/\alpha_p$ independently measured on that survey. The combination $H_i=P_i-O_i/\rho$ eliminates this rotation exactly. Column substitution preserves the rank of the above equation. Any five-channel determinant with rows $(1,\tau_i,W_i,O_i,H_i)$ must vanish. All channels with the two enrolled endpoints and three frequency settings were retained. [Volume adjustment and matched tangents](../../action_response/four_ticks.py)

For each charge, use the three settings at receiver 38 and the quarter and rational turns at receiver 39. In this same five-row order the determinant is in **[5180.19, 6431.61] for vacancy** and **[−329.01, −270.03] for excess**. Over all six possible five-channel subsets, the source enclosure certifies nonzero determinants for five vacancy subsets and four excess subsets. The proof checks the original fixed-volume directions and both species separately. It permits arbitrary shared origin and rate derivatives, arbitrary scalar coupling and separate static channel gains. Source-dependent gains independently assigned to the receivers or a genuinely dispersive law would describe additional physical observables. [Determinants and all channel intervals](results.json)

An independent sparse exponential calculation supplies synthetic multinomial counts at the nominal source. The five physical settings are equilibrium, ordinary, unpaired and the two opposite quarter-turn rotations used to estimate the exact rotation tangent. For $2^{58}$ trials per setting, all three position outcomes are retained. With $x$ hits out of $N$, inversion of the two count concentration bounds gives a radius $(\sqrt{2xL}+2L)/N$; the union bound includes both receivers, both charges, all five source settings, all phases and all twelve appointments. Propagating these simultaneous 99% intervals gives nonsingular determinants in six vacancy subsets and four excess subsets. A separate control that adds a clock delay $\Delta$ to all already-completed records shifts every spectral delay by exactly $\Delta$, as the phase transformation above predicts. These are simulated-model counts, with no astronomical measurement involved. [Count program](counts.py) · [Counts and nominal 99% rank witnesses](counts.json)

This experiment tests a particular finite route to a source-responsive common metric. Resolving ordinary and collar changes in remote phases supplies material evidence of propagation through the prepared source. A shared proper delay, tensor null contraction and shrinking-ball Einstein balance still require a physical construction that survives channel and frequency calibration.

From the repository root,

~~~sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.fringes.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.fringes.counts --write
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/matched_einstein_limit/source_implication/common_front/fringes
~~~

The archived quadrature counts have an exact SHA-256 over their canonical integer table. The test reconstructs the entire confidence and rank receipt from those original counts. Every bin in both the archived and fresh acquisition is checked against the live native probability using a simultaneous 99% Bernstein bound. Both acquisitions independently pass the original confidence, original-denominator and rank checks. At $2^{58}$ starts even valid changes of summation or sampling implementation can change individual integer counts. [Count acquisition and archive replay](counts.py) · [Frozen-count and physical gates](test_fringes.py)
