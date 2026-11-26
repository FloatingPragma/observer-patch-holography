# Entropy, area and volume read from the same retained preparations

Two independent readouts resolve increasing entropy and decreasing clock-survey area on the prepared two-head source. The finite ordinary-gain residual includes zero, whereas doubling that independently calibrated gain gives a residual of the opposite sign from the entropy increase. These are **aggregated synthetic terminal records** for the deposited numerical controller and its Born probabilities. Every route timeout and every reversal-sector outcome is counted. The experiment uses 54 independent banks with $2^{66}=73,786,976,294,838,206,464$ synthetic trials per bank. A multinomial draws each cohort at once; individual microscopic histories were not simulated that many times. [Frozen protocol](PROTOCOL.md) · [All terminal counts and probability tickets](records.json)

| Finite difference per unit increase of the twirl probability | Central estimate | Simultaneous 99% interval |
|---|---:|---:|
| Two-head entropy, nats | +0.41249 | [+0.410813, +0.414157] |
| Survey tetrahedron area, squared clock units | −0.01840 | [−0.019575, −0.017216] |
| Survey tetrahedron volume, cubed clock units | +0.000004 | [−0.000342, +0.000349] |
| Entropy plus ordinary gain times area, nats | −0.01253 | [−0.041422, +0.016367] |
| Entropy plus twice the ordinary gain times area, nats | −0.43754 | [−0.493658, −0.381423] |

Both endpoint stocks are the attracting preparations at $1/8\pm10^{-5}$ twirl with the same reversal weights and original controller. The compensating native phase acts on the first head's $Y$ axis. The source preparation uses the modeled survey to match its reference volume. The independent clock banks then test that volume on the two prepared stocks. Every difference divides positive endpoint minus negative endpoint by $2\times10^{-5}$. The analytic derivative and the predicted secant are separate records. In particular, the model predicts a finite residual around −0.01279 at this step, even though the central derivative cancels at the ordinary gain. The count intervals include that finite residual and zero. Their width does not establish exact cancellation. [Physical preparation](../README.md) · [Volume-matched endpoint selector](../held.py) · [Separate ordinary-gain calibration](../../../PROOF.md)

**Head readout.** Independent copies of each stock supply the published 15 two-head Pauli settings. A bank records the parity and the sign of the Pauli outcome. Together these give four disjoint outcomes. For parity $b$ its unnormalized block is $X_b=(w_b I+\sum_{P\ne I}m_{bP}P)/4$. Count both signs of the first setting to estimate $w_b$ and subtract the two signed counts for each $m_{bP}$. The eight eigenvalues of the two blocks give $-\sum_i\lambda_i\log\lambda_i$, including the entropy of the recorded parity mixture. Every estimated coefficient uses the registered terminal outputs. The calibration supplies neither the entropy nor its slope. [Native head and return readers](../../stationary_baths/forecast/measure.py)

**Clock readout.** Each of the 12 directed routes has an independent count bank with equal probe-charge sampling. A successful round trip reports tick 2 through 12; a timeout reports its own outcome. Scores in hand units are $2/4,3/4,\ldots,12/4$ and $3$ for a timeout. Divide each total score by all starts. Refit the four clock gains on each endpoint, construct the six corrected squared ranges and use the four triangle areas and tetrahedral Gram volume. The measured object is this enrolled four-clock tetrahedron. Proper geodesic-ball area, spatial curvature, regional matter stress and a physical null-balance identity require separate measurements.

## Uncertainty and independent replay

Each row uses nonnegative dyadic tickets summing to $Q=2^{36}$. Collect 15 signed Pauli scores in each parity, two parity masses and twelve bounded clock scores at each of the two endpoints. There are at most 88 such expectations. For $N=2^{66}$ independent starts and $h=\sqrt{\log(17600)/(2N)}$, Hoeffding's inequality and a union bound give simultaneous coverage at least 99%. In units of the outcome range, the ticket allowance is $d/Q$ for $d=4$ or 12 terminal categories. The controller's finite settling allowance is $2(1-(1/8-10^{-5}))^{256}$ in trace distance norm for the exact twirl contraction with these fixed parity weights. Each score also carries an explicit $10^{-10}$ numerical allowance. Confidence concerns sampling the registered numerical laws. The parent bound on the approximation to the exact grouped native channel is a different bound and is not absorbed into a claim about an exact microscopic endpoint or proper metric. [Controller precision](../README.md)

The independent verifier computes all (2(15\times4+12\times12)=408) single-bin marginal deviations as integers. Every count differs from its archived expectation by less than $1.291\sqrt N$; the registered validation accepts at most $3\sqrt N$, whose simultaneous Hoeffding tail is below $816e^{-18}$. This check tests the sampled count law and is separate from the 99% inference above. For entropy, each block estimate is a Gaussian-integer Hermitian matrix divided by $4N$. At both endpoints of each of four disjoint dyadic brackets its integer characteristic determinant has opposite signs. Consequently all four eigenvalues are enclosed without numerical diagonalization in the verifier. Orthogonality of the 16 Pauli matrices bounds the Frobenius error in each block by twice the common coefficient error; Weyl's inequality transfers it to each eigenvalue. The entropy gradient is then bounded over those positive intervals. For the clocks the verifier computes each mean exactly from integer ticks, bounds the area and volume derivatives over their complete range box and applies the mean-value inequality. The ordinary gain uses its independently established interval $[23.10423617469,23.10423621780]$. All nonlinear interval operations in the independent certificate use 45-digit interval arithmetic. [Root brackets and portable response certificate](certificate.json) · [Verifier independent of NumPy and random-number generation](audit.py)

The seed protocol was committed before sampling in commit 4acbb609b. The sampling generator was NumPy 2.3.5 PCG64 with seed 362741903 on Darwin arm64 and Python 3.12.6. A platform-matched seed replay uses the deposited dyadic tickets. The mathematical verification uses the actual archived counts on every platform. Counts are encoded as decimal strings to preserve all bits in JSON readers.

From the repository root:

```sh
export PYTHONPATH=code:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=matched_einstein_limit.source_implication.full_algebra.reuse.retained_equilibrium.observer_banks
python -m "$p.audit"
python -m "$p.measure"
python -m "$p.measure" --replay-seed
python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/reuse/retained_equilibrium/observer_banks
```
