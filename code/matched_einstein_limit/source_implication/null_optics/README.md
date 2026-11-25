# Do two material probes measure the same directional travel time?

At the certified 80-cell entropy maximum, send either a vacancy or an excess between four material stations and read the first arrival. At half-unit ticks and a six-unit timeout, their mean times cannot represent a common directed travel table, even with independent positive clock rates for every sender and probe charge and arbitrary additive delays at both ends of each trip. This test uses the equilibrium itself and allows asymmetric travel in every direction. A separate result identifies a source-sensitive early transmission coefficient that both probes do share. [Certified equilibrium and its preparation](../full_algebra/README.md) · [One-way clock measurements](results.json) · [Rational sign bounds](certificate.json)

## A clock comparison with free rates and offsets

The four stations are $(0,1,38,39)$. At intervals $h=1/2$, ask the receiving port whether the tag has arrived. Stop on the first hit; a trial with no hit receives the full timeout $T$. There is no return leg. Every copy starts with the same freshly prepared collar and an isolated clock reference. Write $Y^q_{ij}$ for this capped mean with probe charge $q$. A shared directed table $D$ with positive sender gains and additive endpoint delays would require

$$g_i^q Y^q_{ij}=D_{ij}+a_i^q+b_j^q,\qquad g_i^q>0.$$

The twelve numbers $D_{ij}$ are entirely free, including both directions of every edge. The rates and offsets for a vacancy may differ from those for an excess. For a source-change experiment the same criterion can be applied to a signed response table, with gains held over the source change and arbitrary changes in both endpoint delays. A source-dependent change in the probe's propagation or detection away from the endpoints falls outside this model.

Subtract the species equations and take their weighted sum with any matrix $W$ of zero row and column sums. The endpoints and the common table cancel. On the **unperturbed exact equilibrium** at timeout six, take

$$
W=\begin{pmatrix}
0&1&-300&299\\
-171&0&299&-128\\
-118&289&0&-171\\
289&-290&1&0
\end{pmatrix}.
$$

The four excess row sums $\sum_j W_{ij}Y^+_{ij}$ are positive and the four vacancy row sums are negative. The smallest certified absolute margin is 0.19816 in weighted clock units. Multiplication by arbitrary positive station rates preserves every sign. The difference therefore cannot be zero. This certificate requires no condition on the gains across preparations and uses no geometric fitted distances.

## A four-edge test on a physical unpairing

The state has two unnormalized coin blocks $X_b$ for even and odd left reversal parity. Replace the right coin with an independently drawn right marginal while keeping the left side, including its parity. The signed tangent is $D_b=\operatorname{Tr}_R X_b\otimes\operatorname{Tr}_{b,L}(X)-X_b$. This physical change preserves both coin marginals and has certified positive intrinsic entropy response on the same preparation. Write $Y$ now for the change of each clock mean. Compare the two receivers at each of the two outer senders.

| Positive difference convention | Excess clock change | Vacancy clock change |
|---|---:|---:|
| $Y_{0,1}-Y_{0,38}$ | +0.02097755 | −0.00344490 |
| $Y_{39,38}-Y_{39,1}$ | +0.02342093 | −0.00330258 |

The two receiver orderings are opposite. The corresponding four-edge weight has zero row and column sums, so all endpoint offsets again cancel. The worst certified sign margin is 0.00328098 at timeout six, 0.00107117 at four and 0.00001481 at two. On the freely formed even source, strengthening pairing relative to independent coins gives the reverse sign pattern with margin above 0.00052177 at six. On the equilibrium's ordinary unread-sign intervention alone, positive fitted gains are feasible on all three timeouts for this unrestricted directed-table test. Feasibility for a single tangent does not establish a geometric distance. [Preparation and signed interventions](../full_algebra/PROOF.md)

## Holding the survey volume fixed during both physical changes

The physical operation also admits a volume-matched version. Rotate the right coin about its measured marginal Bloch axis. This unitary preserves the coin marginals, the total entropy and the parked energy. For raw direction $D$ and rotation tangent $R$, the velocity on the cap-six survey's constant-volume surface is $D+\alpha R$, with $\alpha=-V'(D)/V'(R)$. The tetrahedral survey determines $\alpha$ without consulting the one-way clocks. At the published centre the angles per unit change are $-294.1305$ for unpairing and $-989.9471$ for ordinary dephasing. [Physical rotation and finite-volume matches](../signed_transfers.py)

With **volume-matched unpairing**, the same four-edge excess splits are $+0.049908,+0.004383$, and the corresponding vacancy splits are $-0.012724,+0.005317$. The last sign has changed. That particular four-edge certificate no longer excludes a common travel response. Instead, use the zero-row, zero-column matrix
$$
W=\frac{1}{1000}\begin{pmatrix}
0&1000&-883&-117\\
50&0&-117&67\\
62&-112&0&50\\
-112&-888&1000&0
\end{pmatrix}.
$$
At the exact enclosed equilibrium every excess row sum is positive and every vacancy row sum negative, with a certified signed margin exceeding **0.0007037**. The enclosure includes $\alpha\in[-296.674,-291.616]$, the marginal rotation axis, and the errors of the distinct one-way and two-way effects. To obtain it, differentiate the tetrahedral volume on the enclosed two-way clocks. For each coordinate Pauli rotation, enclose a one-way observable's commutator on the full stationary state ball. Contract those three intervals with the enclosed unit marginal axis and propagate the volume correction into the eight weighted rows. Thus the full-directional exclusion also holds on the entropy-increasing, fixed-volume unpairing tangent. [Interval calculation and witness](matched.py)

For **volume-matched ordinary dephasing**, the centre's four-edge excess splits are $+0.031213,-0.063752$ and the vacancy splits are $-0.086149,+0.028303$. Their signs persist over the certified equilibrium enclosure. At the centre a numerical fit to this single tangent admits positive gains with their sum normalized to one and the smallest equal to 0.0071604. No positive separating margin exists for that centre's ordinary-only table at the reported numerical precision. The ordinary result and the certified failure for the matched pairing change have different scopes. The independently verified matched entropy and scalar-area balance continues to hold on the selected survey.

## One directional source measurement that both probes share

As a fixed part of the apparatus on each fresh equilibrium copy, sharply read the left contrast sign and discard its outcome. On the two reversal blocks this sends $X_0,X_1$ to $(X_0+X_1)/2$ in each parity. Launch the tag and make a single remote arrival check after a short exposure $t$, travelling from 0 to 39 or back. The first three hops connect these stations. For the twice-native integer generator $H_2$, the first nonzero block in the amplitude is $(-it/2)^3(H_2^3)_{ji}/3!$. The squared block before division by $2^6(3!)^2$ is 4096 or 5184 times a rank-one coin projector, according to parity. Its index is 00 or 11. Either reversing the route or reversing the probe charge exchanges the indices. These statements follow from integer multiplication on all 960 probe positions.

For the physical unpairing above the total coin tangent has zero single-coin marginals. Its 00 and 11 entries therefore coincide. Averaging the two reversal coefficients through the unread sign operation gives the same signed change of the leading $P(t)/t^6$ in both directions and with either probe. Its exact coefficient is $(145/72)(D_0+D_1)_{00,00}$, numerically $-0.232266$ on the certified equilibrium. The equality for the two species and the two directions is algebraic. Without any intervention the freely formed even source also has a common leading coefficient for its strengthening tangent, about $+0.030936$.

This is an intensity coefficient. The first nonzero time power remains six. The unread-sign operation erases the ordinary parity-dephasing tangent exactly because applying the same dephasing twice has the same effect as applying it once. Measuring capped mean times with this same preprocessing still fails the four-edge pairing test at timeout six. Shared attenuation can coexist with different material arrival distributions. Identifying freely propagating null rays and their optical expansion with a common proper metric requires an independent propagation observable.

## Error control and replay

Write $F_n$ for the first-hit amplitude after $n$ inspections. The positive mean-time effect satisfies $0\le E\le T I$ and equals
$$E=T I+\sum_{n<T/h}(nh-T)F_n^\dagger F_n.$$
The operator recurrence is $F_n=U_{ji}(nh)-\sum_{k=1}^{n-1}U_{jj}((n-k)h)F_k$. The certificate reuses the bounded, integer Chebyshev evolution for both charge and parity sectors, including the native source amplitude. It propagates rational error bounds through each first-hit sum. Each four-by-four multiplication has norm error at most $8/2^{100}$ from fixed-point rounding. The final operator error is below $10^{-13}$ per edge. Separately, the exact entropy maximum lies within $2\times10^{-7}$ in Hilbert–Schmidt norm of the supplied centre, or less than $6\times10^{-7}$ in trace norm. Partial trace and the product difference bound make the unpairing tangent's trace-norm uncertainty at most three times this. These bounds cover the signs on the exact, enclosed source. The linear program suggests feasible gain sets; the displayed obstructions use explicit integer weights and rigorous sign bounds. [Floating comparison](measure.py) · [Certificate](certificate.py)

Capped capture time includes detector coupling, repeated interrogation and missed events. The exclusion applies to interpreting these reported directional times as a common travel table. A geometric contribution to a more complete propagation model is compatible with the result.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python3 -m matched_einstein_limit.source_implication.null_optics.measure
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python3 -m matched_einstein_limit.source_implication.null_optics.certificate
OPENBLAS_NUM_THREADS=1 python3 -m pytest -q code/matched_einstein_limit/source_implication/null_optics
```
