# A long-wavelength density clock and a driven load pulse

On the unchanged blank Interval Echo federation, the spontaneous fluctuation of a long density wave determines the response to an independently released load pulse at every time in a homogeneous product medium. At an equilibrium with correlated loads, changing the filling still moves the two readings in the same direction. Changing the correlations while holding the one-port marginals fixed moves them in opposite directions. The sign reversal persists on the longest measured mode of the 20-, 80- and 320-cell boards. [Exact fractions, complete-cell time records and count tables](results.json)

## One stationary source for the ordinary calibration and the held test

Let `P(f)` assign an independent unit load to each port with probability $f$. Its density is constant within each conserved load-number sector, hence stationary under the native Hamiltonian. Prepare
$$
 \rho(f,e)=(1-e)P(f)+\tfrac e2[P(f-1/8)+P(f+1/8)] .
$$
Erase the lottery label. Both instruments use independent copies of the same resulting density. At $f=1/4,e=1/2$ all words have positive probability. Varying $f$ changes every one-port mean while keeping the connected covariance $e/64$ between distinct ports fixed. Varying $e$ holds all one-port densities fixed and changes the overlap correlations throughout the board. This stationary intervention also changes the variance of the total load.

At the common reference a sharp read on any seam gives probabilities $(73,23,23,9)/128$ for $00,01,10,11$. Their $e$ derivatives are $(1,-1,-1,1)/64$. The marginal modular terms have zero expectation derivative. The independently readable joint modular slope is
$$
 \operatorname{Tr}[(\partial_e\rho_{AB})(-\log\rho_{AB})]
   =\frac{1}{64}\log\frac{529}{657}=-0.00338587 .
$$
A separate synthetic sample of $2^{30}$ seam preparations gives a simultaneous 99% interval $[-0.00341044,-0.00336319]$. The native mean bond energy changes by $-(1+2\kappa)/64$ per unit $e$. The degree-reference terms have constant expectation along this variation.

## How the two records are acquired

Each board supplies a registered balanced pattern $g$ with $\sum_i g_i=0,\sum_i g_i^2=1$. Select the lowest nonconstant graph Laplacian eigenspace numerically and project the port-zero indicator into it. Round the unit vector in units of $2^{-24}$ and remove its mean exactly. The saved integer profile specifies every sign and weight. Its Rayleigh quotient and squared eigenvector residual are exact rationals. The generator uses the native port incidence. No spatial metric enters the enrollment. [Integer port profiles](profiles.json)

* **Passive clock.** Read a sharp occupation word at release and again at the appointment. From each pair compute $F(0)F(t)$, with $F=\sum_i g_i n_i$, retaining every word. Denote the resulting mean by $C(t)$. For the stationary diagonal preparations above there is a positive, symmetric spectral measure for this correlation. Its current-weighted squared frequency is $P=C^{(4)}(0)/[-C''(0)]$.
* **Driven transport.** In fresh copies, choose a port $i$, reset it to full or empty by a rail swap, and discard its prior load inclusively. Read the terminal occupations. Form the signed, mode-weighted sum $R(t)=\sum_i g_i[\langle F(t)\rangle_{i,1}-\langle F(t)\rangle_{i,0}]$. This reads how an independently released spatial density pattern disperses. Its squared-frequency statistic is $T=R^{(4)}(0)/[-R''(0)]$. The test also records $R(t)$ directly at finite time.

For a product filling the difference of the two reset states at $i$ equals $(n_i-f)P(f)/[f(1-f)]$. Summing with $g_i$ proves $R_f(t)=C_f(t)/[f(1-f)]$ for **every appointment** and therefore $T=P$. For a mixture, the driven record averages the conditional $R_f$ with the preparation weights. The passive $C(0)$-normalized record averages the same conditional signals with weights proportional to the preparation probability times $f(1-f)$. A component with larger fluctuations contributes more to a spontaneous density record. A controlled reset is issued equally often in each component of the lottery. This distinction is observable after the label has been erased.

## Transfer on the same equilibrium

Entries are derivatives at $f=1/4,e=1/2$, in native inverse exchange-time squared per unit of the indicated variation. The ordinary measurements give gain $(\partial_f T)/(\partial_f P)$. Multiplying the independently read $\partial_e P$ by this gain predicts the held transport response. The final two columns can be compared without choosing frequency units.

| Ports | Profile Rayleigh quotient | Ordinary passive | Ordinary driven | Held passive | Held driven, predicted | Held driven, measured | Measured / predicted |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 2.76393 | 443.601 | 580.854 | 8.235 | 10.783 | −18.152 | −1.68333 |
| 240 | 0.0441425 | 964.215 | 1253.010 | 20.322 | 26.409 | −36.594 | −1.38569 |
| 960 | 0.0103158 | 958.335 | 1245.000 | 20.292 | 26.362 | −36.260 | −1.37547 |
| 3840 | 0.00255460 | 956.800 | 1242.807 | 20.310 | 26.381 | −36.143 | −1.37003 |

On the 12-port cell direct diagonalization includes every one of the 4,096 words. The passive calculation sums spectral weights. The driven calculation propagates the independent full and empty releases and reads terminal probabilities. At time $1/8$ in native units, finite shots on either side of the common source yield the following differences in normalized signal, upper arm minus lower arm.

| Change | Passive difference | Driven difference |
|---|---:|---:|
| $f=1/4\pm1/32$, fixed $e=1/2$ | +0.00024884 | +0.00033650 |
| $e=1/2\pm1/4$, fixed $f=1/4$ | +0.00001972 | −0.00010000 |

The simultaneous 99% halfwidth of each difference is below 0.0000042. Each arm uses $2^{50}$ synthetic binary scores. Draw the sender and receiver labels independently with probabilities proportional to $|g_i|$. In the passive arm compress the complete pair of words to the product of the selected occupation signs and pattern signs. Its expectation is $4C/\|g\|_1^2$. In the driven arm use a fair full/empty reset and multiply its sign by the receiver occupation sign and the two pattern signs. Its expectation is $R/\|g\|_1^2$. Each binomial parameter is rounded to a multiple of $2^{-38}$. The finite-count bounds include that rounding and a two-sided Hoeffding union bound over all settings within each law. The cell probabilities come from binary64 diagonalization. Independently released one-load dimers read $\sin^2(t/2)$ at the registered appointments. The common release schedule is supplied.

## Finite calculation on all three federations

For a real pattern $F(s)$ diagonal in word $s$, the first two even moments are the squared norms of $[H,F]$ and $J=[H,[H,F]]$ in the stationary density. Writing $L$ for the combinatorial Laplacian and $\chi=f(1-f)$, the normalized second moment is exactly $g^TLg/2$ at every homogeneous filling. It has the same value for the driven measurement and the normalized passive measurement at every admitted mixture. This positive agreement has zero filling derivative and cannot calibrate a material-dependent rate by itself.

The diagonal of $J$ is $\tfrac12\sum_i(Lg)_i n_i$. For words differing by moving one load from $i$ to empty $j$, its off-diagonal entry is
$$
 J_{st}=\mathbf 1_{i\sim j}\frac{g_i-g_j}{2}(H_{ss}-H_{tt})
       +\frac14\sum_{k\sim i,j}(g_i+g_j-2g_k)(1-2n_k).
$$
On the binary sector $H_{ss}=6\kappa N+\tfrac12\sum_i d_i n_i-(1+2\kappa)\sum_{i\sim j}n_i n_j$, where each seam is counted once. Thus each off-diagonal entry is an affine function $a+\sum_k b_k n_k$ of spectator bits. The two orderings of disjoint swaps cancel because $F$ is additive. All remaining terms have been enumerated. After division by $\chi$, the diagonal contribution to $\langle J^2\rangle$ is $\|Lg\|^2/4$, and each unordered endpoint pair contributes
$$
 2\big[(a+f\sum_k b_k)^2+\chi\sum_k b_k^2\big].
$$
This yields a quadratic polynomial in $f$. Integer arithmetic on the registered profiles gives exact rational scatter coefficients, derivatives and transfer ratios. Direct native-Hamiltonian checks include a graph with unequal degrees and the complete 12-port spectrum. The changed law $\kappa=0$ preserves the graph, exchange, alphabet and preparations. On the 3,840-port board the driven change from the product reference to its equal mixture drops in magnitude from 36.14276 to 0.11699. [Implementation](measure.py) · [Independent checks](test_measure.py)

### Finite-appointment signs on the federations

The same four finite preparations used for the cell count test give rigorous signs for actual signals on each larger board. Each normalized product signal has a positive spectral measure of mass one. The two mixture protocols use convex sums of these measures. An upper bound on every native energy difference is $W=(1+\kappa)|E|+\kappa(6N-2|E|)$, obtained by summing the positive binary bond and reference norms. The sixth moment of each measured spectral distribution is at most $W^2 m_4$. At a common dwell $t$, the two arms of each comparison have identical constant and second-order terms. Taylor's theorem for cosine therefore gives
$$
 \left|\Delta(t)-\frac{(m_{4,+}-m_{4,-})t^4}{24}\right|
 \leq\frac{W^2(m_{4,+}+m_{4,-})t^6}{720}.
$$
The replay chooses a single dyadic appointment per board that bounds the remainder below half of each leading difference. Both finite ordinary signals increase. These are outward roundings of the exact bounds for the two held changes.

| Ports | Native dwell | Passive signal change | Driven signal change |
|---:|---:|---:|---:|
| 240 | 1/32768 | $[6.92,9.30]\times10^{-21}$ | $[-1.58,-1.34]\times10^{-20}$ |
| 960 | 1/131072 | $[6.32,8.47]\times10^{-24}$ | $[-1.43,-1.21]\times10^{-23}$ |
| 3840 | 1/524288 | $[6.12,8.20]\times10^{-27}$ | $[-1.38,-1.17]\times10^{-26}$ |

These appointments demonstrate finite signal separation mathematically. Resolving the displayed changes by sampling would require much larger exposures than the 12-port trial. The exact finite bounds justify the signs independently of the time resolution available in a particular macroscopic observer apparatus.

## A local coherent collar and the physical scope

On the registered seam $(0,5)$, transfer density from the symmetric single-load state to the antisymmetric one, tensored with $P(1/4)$ on all spectator ports. The signed operator is $-|10\rangle\langle01|-|01\rangle\langle10|$. Both full side marginals are fixed for every bipartition separating the two endpoints. Every sharp-word probability and the entire total-load distribution are fixed. The native bond energy rises by one per unit transfer. The product reference admits both signs for amplitudes below $3/16$. This intervention is localized and its perturbed states can evolve after release. Its modular derivative at the product reference vanishes.

The passive instrument reads a full sharp word at release, which erases this coherence exactly. Each driven release outside the seam preserves it. On the 240-, 960- and 3,840-port patterns the signed change of the driven second time derivative is respectively $-9.1953940\times10^{-4}$, $+5.1044240\times10^{-5}$, and $+3.0726059\times10^{-6}$. These are exact-rational signs for the registered profiles. On the 12-port cell the driven derivative at time $1/8$ is −0.01631473. The passive difference is exactly zero even at finite time. The large-board local signal becomes small, so this test supplies no nonzero limit for a single fixed seam.

The measurements sample long spatial patterns of the native load density. The squared-frequency statistic is a ratio of spectral moments and can receive contributions from fast modes. Its contrast dependence and the shrinking second moment do not identify an individual acoustic pole, a ray flight time or a proper clock. The data concern blank binary federations, a supplied grand-number source lottery and specified instruments. A proper metric and a geometric null-stress law on the parked-head source would require an independent identification. The comparisons here constrain proposals that identify these particular material readings with one common rate.

Run the measurement and tests from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -m model.experiments.seam_program.collective_rate.long_wave.measure
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -m pytest -q -W error model/experiments/seam_program/collective_rate/long_wave
```
