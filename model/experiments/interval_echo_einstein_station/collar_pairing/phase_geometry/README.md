# Reading the formed collar with material probes

A mobile excess on the 80-cell federation leaves correlations between two seam memories. A second excess can detect these correlations by starting at an unmarked port and asking a single question at a clock appointment. Did it return to its launch site? Its return probabilities over time determine the entire joint state of the two stored coins. A separate instrument transfers the coins to two material position registers and measures their regional modular response. These experiments provide a source reading for a comparison with an independently reconstructed metric. [Source formation](../../mobile_collar_source/README.md) · [Measurements and counts](results.json)

## One port determines the source state

The forming excess starts at port 0, evolves freely for four exposure units and is removed with its location unread. The two marked contacts occupy (0,1) and (38,39). Write $\rho$ for their joint coin density in the sign-even, magnitude-one sector and $\pi=\rho_A\otimes\rho_B$ for the independent preparation with exactly the same marginals. Both preparations have parked energy 16. A fresh excess starts at port 5. All 960 ports and all of the federation's native contacts participate during this free flight. A separately prepared exchange clock supplies timing in native exposure units. Each trial gives one sharp position reading at a single chosen time. An outward excursion, return path or intermediate hit need never be identified. All misses count. [Source and clock conventions](../ranging/README.md)

**Finite identifiability theorem.** For the specified native generator and launch port 5, the function giving the return probability on any nonempty open interval determines every density matrix on the four-dimensional two-coin space. The claim assumes known flight dynamics, a calibrated clock and fresh copies of the same source preparation.

Let $H_2=2H+cI$ be the integer flight matrix and $A_j=(H_2^j)_{55}$ its four-by-four return block. With $U(t)=(e^{-itH})_{55}$, the measured probability is $\operatorname{Tr}(U^\dagger U\rho)$. The $n$th derivative of $U^\dagger U$ at zero is $(-i/2)^n$ times

$$B_n=\sum_{j=0}^n (-1)^j {n\choose j} A_j^\mathsf{T}A_{n-j}.$$

The scalar $c$ cancels from the effect. The even $B_n$ are symmetric, and the odd $B_n$ are skew. Record the ten upper-triangle coordinates of the even matrices and the six strict upper-triangle coordinates of the odd ones in disjoint columns. The orders $0,4,6,7,8,9,10,11,12,13,14,15,16,17,18,20$ give a square integer matrix with determinant congruent to 429 modulo 1009 and 128 modulo 1013. The receipt contains the complete determinant recomputed by unbounded integer path sums. The effect derivatives therefore span every Hermitian observable. Equal probability curves on an open interval have equal derivatives at zero by analyticity. Their difference would have zero pairing with every Hermitian observable, so the two states coincide. Adding launches and readings at port 40 gives full rank through degree 10. If the squared-contrast coefficient is zero during interrogation of the same stored source, $H=H_{\rm position}\otimes I_4$ in this sector and the effect space has rank one. [Rank certificate and finite-time records](remote.py)

In a simulated sample of $2^{40}$ returns per preparation at each of four appointments, a single terminal distinguishes the paired and independent sources at exposure 4. The population difference is $-6.83643\times10^{-6}$. The interval from the counts is $[-1.180\times10^{-5},-1.933\times10^{-6}]$ at 99% simultaneous coverage across exposures 0.5, 1, 2 and 4. The other three intervals contain zero. Recovering all the source coordinates is much less stable. At the 32 quarter-unit appointments through exposure 8, the smallest singular value of the known-trace Pauli design is $1.1054\times10^{-11}$ with one port and $2.7935\times10^{-8}$ with two.

## Transferring the source to material rails

An alternative destructive read uses two disjoint marked contacts on a freshly formed and parked copy. Admit one material excess on the first port of each marked contact and isolate these contacts from the other seam and site services during interrogation. The choice of switches, input rails, clock appointments and sharp port readers is part of the apparatus. The full 99-state contact restricts exactly to

$$h_+=8(I-|\psi_+\rangle\langle\psi_+|),\qquad
 |\psi_+\rangle=(|01\rangle+|10\rangle)/\sqrt2,$$

in the basis of material position followed by coin. This follows by applying the integer contact generator to all columns of the sign-even subspace. At exposure $\pi/8$ its propagator is $-(Z\otimes Z)\mathrm{SWAP}$. Starting each material rail in $|0\rangle$ transfers the arbitrary two-coin state to the two material registers, locally conjugated by $Z$. The head outputs are fixed. A vacancy conjugates the contact by $I\otimes X$ and gives the local rail correction $ZX$. Linearity gives the same transfer with any coin entanglement or external reference. At zero contrast the native marked generator vanishes on these columns and the outgoing rails are source independent. [Contact projection and gates](measure.py)

Four disjoint blank seams swap the two rails onto output dimers (5,9) and (40,45) in exposure $\pi$. Direct position reads measure $Z$. A half exchange followed by a position read measures $Y$ up to sign. A preceding quarter turn from the degree-five site reference gives $X$ up to sign. Its energy gap is $8[(1+q)^2-1]$ for charge $q=\pm1$. The recorded signs include these fixed orientations and the transfer phase. Products of the two position records therefore recover all 15 nonconstant Pauli means of the input, using independent copies for different settings.

## Regional information observed on a held cohort

The fixed-reference modular observable is $K=-\log\rho$ for the faithful formed preparation. In the strengthening orientation the slope is
$\Delta K=\operatorname{Tr}[(\rho-\pi)K]$. The product-marginal modular term vanishes since both single-coin marginals agree. An independent cohort of excess probes calibrates $K$ through the two material rails. The smallest estimated eigenvalue is $7.0674\times10^{-9}$. Taking its logarithm is sensitive to the reference estimate. The specified simulation therefore aggregates $2^{74}$ Bernoulli trials per nonconstant Pauli word for reference calibration. This is a large ideal-ensemble resource requirement. The recorded integers come from binomial count sampling. They represent model records, with no observations of galaxies or laboratory matter used.

Each held cohort uses $2^{18}$ trials per word per preparation. All counts in the following intervals use independent simulated trials, with one common reference. The population slope is $-0.29384843$ nats per unit change of pairing.

| Held material | Estimated slope | Simultaneous 99% interval including reference uncertainty |
|---|---:|---:|
| Excess | −0.28900 | [−0.37913, −0.19887] |
| Vacancy | −0.30573 | [−0.39587, −0.21560] |

For reproducibility, let $N_r=2^{74}$, $N=2^{18}$, $\alpha_r=0.004$ and $\alpha_t=0.003$. Each reference Pauli mean has simultaneous error at most $e=\sqrt{2\log(30/\alpha_r)/N_r}$ by Hoeffding and the union bound. Pauli orthogonality gives $\|\widehat\rho-\rho\|_{\rm op}\leq\|\widehat\rho-\rho\|_2\leq\sqrt{15}e/2=d$. Operator monotonicity of log bounds the reference observable error by $L=-\log(1-d/\lambda_{\min}(\widehat\rho))$. For the held signed Pauli outcomes use the measured coefficients $k_j=\operatorname{Tr}(-P_j\log\widehat\rho)/4$. Conditional on calibration, weighted Hoeffding bounds the error of the paired-minus-independent estimate by $2\sqrt{\log(2/\alpha_t)\sum_jk_j^2/N}$. Adding $2L$ covers the unknown true reference since $\|\rho-\pi\|_1\leq2$. The union bound over the reference and the two held cohorts gives 99% coverage. [Raw counts and uncertainty budgets](results.json)

A material observer can separate any two stored coin states in this sector and can compare a proposed metric signal with an independently calibrated source response. Converting the clock or a linear combination of its readings into geometric distance requires an independent operational definition. The measurements give neither a proper small-ball volume nor a relation of its signed variation to physical stress.

Run from the repository root with Python, NumPy, SciPy and SymPy.

```sh
python3 -m model.experiments.interval_echo_einstein_station.collar_pairing.phase_geometry.measure
python3 -m pytest -q model/experiments/interval_echo_einstein_station/collar_pairing/phase_geometry
```
