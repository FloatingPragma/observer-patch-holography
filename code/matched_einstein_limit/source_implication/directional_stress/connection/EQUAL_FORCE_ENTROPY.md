# Equal radar and equal force with increasing entropy

At the volume-matched stationary state of the uniformly reused radar survey on the active 80-cell board, two physically preparable sources have exactly the same twelve directed mean ranges and exactly the same complete marked forces for vacancy and excess. Each parity weight and both one-head marginals in each parity block also agree. The head entropy is strictly greater on the displaced source. The fourteen measured means therefore cannot determine a first-order compensation for the entropy change at this stationary source. [Stationary preparation and twelve-range proof](../../full_algebra/reuse/INVISIBLE.md) · [Native force apparatus](README.md#equal-radar-ranges-with-a-different-measurable-material-force)

## Four remaining directions

Use coordinates $x_{bi}$ on the two unnormalized head blocks, with parity $b$ even then odd, and normalized Pauli products $B_i=\sigma_a\otimes\sigma_c/2$ in the order $XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ$. All eighteen coordinates preserve the trace and both partial traces of each block. Their Euclidean norm equals the Hilbert–Schmidt norm of the joint perturbation.

For a four-tick force probe of charge $q$, write $L_\ell$ for embedding the heads at launch port $\ell$, $W_{qb\ell}=\exp(-4iH_{qb})L_\ell$ and $F_{qb}=i[V_{qb},Q]$. Define
$$
C_{qb}=\frac1{10}\sum_{\ell=1}^{10} W_{qb\ell}^\dagger F_{qb} W_{qb\ell},\qquad
m_{q,bi}=\operatorname{Tr}(C_{qb} B_i).
$$
The sum runs over both five-port launch rings with equal weight. The mean force is $\sum_b\operatorname{Tr}(C_{qb}\chi_b)$, including all launches and outcomes. Consequently its response to an incoming source change is exactly affine. The directly calibrated delayed eight-pair read implements the same full operator. The force gradient uses the declared chart, active native action, and fixed exchange-time appointment. [Four Hermitian effects and two gradient rows](force_gradients.py) · [Independent delayed read](late_read.py)

Stack the twelve charge-pooled range rows and the two force rows in the $14\times18$ matrix $A$. All are population expectations on fresh copies of the source. The exact range and force effects are approximated by rational matrices with certified errors. Each of the two force rows has Euclidean error below $6\times10^{-15}$. The operator norm error of the entire stack is below $2.494\times10^{-14}$.

The verifier builds a rational approximate right inverse $K$ and a proposal $d$ from the component of the entropy gradient in the approximate kernel. It corrects the proposal's nominal residual using exact rational arithmetic. The independent bound $\|I-AK\|<5.679\times10^{-7}$ proves that $A$ has rank fourteen. Its kernel on this slice has dimension four. For the *true* response map set
$$
D=d-K(AK)^{-1}Ad.
$$
Then $AD=0$ exactly and $\|D-d\|_2<5.676\times10^{-7}$. Each of the fourteen native means is affine along $\chi(t)=\chi_*+tD$. All fourteen coincide on the entire feasible chord, including at finite displacement. Since the tetrahedral area and volume are computed from the twelve ranges, both are constant there. [Right-inverse and exact-kernel certificate](force_entropy.py) · [Bounds, direction and complete force rows](force_entropy.json)

## Entropy and a direct head read

Let $\sigma_b$ be the normalized conditional stationary blocks. The coordinate entropy score is $g_{bi}=-\operatorname{Tr}(B_i\log\sigma_b)$. The scalar log parity weights vanish on these traceless changes. The independent conditional-state certificate encloses each $\sigma_b$ in a Hilbert–Schmidt ball of radius $2.526\times10^{-6}$. Exact positivity tests give lower eigenvalue bounds 0.20 and 0.21 throughout the respective balls. The integral formula for the logarithm bounds the Hilbert–Schmidt perturbation of $\log\sigma_b$ by the state radius divided by its lower eigenvalue. The verifier combines that error, the rational logarithm error and the exact-kernel correction. It gives
$$
1.1964104\times10^{-5}<\left.\frac{d}{dt}S(\chi_*+tD)\right|_{t=0}<4.7077286\times10^{-5}\quad\text{nats}.
$$

The native preparation inequalities have positive slack throughout $|t|\le0.01$. For $0\le t\le10^{-7}$ the full density matrix has minimum eigenvalue above 0.026. On this chord the entropy second derivative is bounded below by $-\|D\|_2^2/0.026$. Integrating proves a finite entropy gain above **$1.0041\times10^{-12}$ nats** at $t=10^{-7}$ with all fourteen means exactly unchanged. The center is stationary under reuse. Only native preparability is needed for the neighboring states. For any differentiable correction formed from the fourteen means and the preserved marginals, its derivative on this direction vanishes.

For any positive integer $n$, independent preparations on $n$ copies give $S(\chi(t)^{\otimes n})=nS(\chi(t))$. The fourteen means on every copy, their averages across copies and the reconstructed tetrahedra all have exactly the same values on the two arms. The positive entropy derivative per copy persists for every $n$. This is a family of repeated surveys at one physical scale. A spatial refinement requires its own proper-ball construction.

A parked, parity-blind $X\otimes X$ measurement distinguishes the sources. This product was chosen for its largest predicted contrast among the nine product axes on this witness. Its signed mean slope lies in $[-0.470564,-0.470559]$. The difference in plus probability between $t=-0.01$ and $t=+0.01$ is consequently about $-0.004706$. Positive endpoint calibration states use the independently published fixed-point centers and odd weight 0.1343. A seeded simulation measures both endpoint probabilities with $2^{22}$ trials per endpoint. Simultaneous 99% Hoeffding bounds resolve the difference. The short, one-sided chord above certifies the entropy increase. [Exact direction, measurement probabilities and counts](force_entropy.json)

Full arrival-time records contain more information than these fourteen summaries. Eighteen selected return-time bins have full-rank affine response $y=y_*+R(x-x_*)$ on this same slice. The native force is exactly forecast from its population means by
$$
M_q=M_{q,*}+(R^{-T}m_q)\cdot(y-y_*).
$$
This assumes fixed hardware and the same incoming preparation on fresh copies, and uses unconditional clock-bin frequencies. A known reference supplies the constant offset. The independent force meter can test that prediction. A rank condition on measurement records does not specify an observer's proper metric or the stress tensor entering a null equation. [Timestamp effects](../../full_algebra/reuse/clock_history.py)

From the repository root, regenerate the certificate with `OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.directional_stress.connection.force_entropy --write`. The focused test is `OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/directional_stress/connection -k changes_at_exactly_equal`.


Thirty-one charge-pooled bins also recover the parity weight and both marginals. Frozen weights predict the two absolute forces on the entire normalized enrolled algebra. [Clock calibration and held forecasts](../../full_algebra/reuse/CLOCK_HISTORY.md#absolute-force-from-clock-counts)
