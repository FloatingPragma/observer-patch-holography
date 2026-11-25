# Can resident return clocks measure the area in the two-load entropy balance?

The same monitored 80-cell source supports two independent measurements. Exterior records recover the boundary information that cancels the interior entropy variation. Four clocks on those exact exterior-facing ports also assign a tetrahedral area and volume. **The clock area fails the identification.** There is a finite change in the two-load source that leaves its mean energy and the entire distribution of each registered three-outcome clock read unchanged, yet decreases the interior entropy by at least **0.0000553 nats**. An independently prepared source with dyadic probabilities reproduces the separation in synthetic counts. A fifth measured station also gives a direct test of the assumed three-dimensional geometry. [Integer certificate](certificate.json) · [Finite source lottery and measured counts](counts.json)

## Two measurements on the identical region

Use the original adjacent two-cell patch, the same native contact Hamiltonian, two resident double loads, and the same half-unit monitoring of the complete sharp word. Park the exterior contacts during each evolution interval. The possible records are the 276 pairs of occupied ports among the 24 ports. Preparation of a sharp word, exterior reads, undo operations and mixing via independent source lotteries use the existing hardware. The unbiased monitored kernel $P_{xy}=|\langle y|e^{-iH/2}|x\rangle|^2$ preserves the uniform distribution. Its numerical spectra in the parent experiment support convergence from individual words. The uniform preparation can also be sampled directly. The local 80-cell patch comprises cells 0 and 3 with actual cut ports $3,10,36,44$. These are the same ports and word ordering used by the entropy meter. [Original experiment](../README.md) · [Patch, generator and clock construction](measure.py)

Let $k$ count occupied rim ports, and $i$ denote the interior word. Each $i$ has $d_k=\binom4k$ compatible exterior records. At the uniform reference the bulk entropy differential is $dH(I)=-d\langle\log d_k\rangle$. There is an independent way to read this boundary information. Let $h_{ab}$ be the probability that $a$ is occupied and $b$ empty before a launch, and let $\mu=\langle k\rangle$ come from a separate exterior count. Then $h=\sum_{a\ne b}h_{ab}=\langle k(4-k)\rangle$. On any preparation the probabilities of one and two occupied rim ports are $h-2\mu$ and $(3\mu-h)/2$. These give the entire expectation of $\log d_k$ using only initial exterior records. The two identities are checked on all 276 words.

The independent clock watches the same two residents. It starts a signal when $a$ is occupied and $b$ empty, records the first monitoring step at which $b$ is occupied and $a$ empty, then waits for the first return to the initial occupancy condition. The signal follows occupations; no tagged third particle enters the two-load Hamiltonian. Every interval lasts one half of a native time unit. At each interval an independent twelve-fair-bit draw has probability $1/4096$ to expire the clock, including simultaneous expiry and return. Failed initial condition, expiry after an eligible start, and successful return are all logged. For return step $T$ and $z=1-1/4096$, success conditional on a launch has probability $s=\mathbb E z^T$. Half the expected elapsed time up to the return or expiry is $1024(1-s)$ native time units. The clock reads this scalar from the two-outcome terminal distribution following an eligible start. Individual timestamps beyond that summary can carry further information.

Average the two measured directional distances on each rim pair. When the resulting Gram matrix is positive definite, its four marked vertices determine a Euclidean tetrahedron. We measure its surface by summing its four triangle areas and its volume by $\sqrt{\det G}/6$. As controls, we score every ineligible start with its expected timeout and repeat the eligible-start measurement with a fixed 24-interval cap. These give different calibration responses. The bit-to-nat conversion in the entropy measurement plays no role in constructing the clock lengths. [Three clock menus on 20-, 80- and 320-cell hosts](results.json)

## Ordinary calibration and independent held changes

For each host train on a probability transfer from the uniform $k=0$ class to the uniform $k=1$ class. On the 80-cell source its change in mean native energy is $-34.27368$ per unit transfer; reversing the transfer raises energy. To match volume, choose the pair of equal-energy words in the same $k$ class with the largest measured volume response. Moving probability within this pair has zero first-order interior entropy response. Neither the area nor the entropy enters its selection. The first-order calibration at fixed volume gives $\eta=+0.000451157$ nats per squared clock unit. On the 20-cell two-cell patch the same prescription gives a negative coefficient. The local 320-cell patch repeats the 80-cell answer. The host change by itself supplies no shrinking-ball limit.

Fresh acquisitions use the two signed stocks at displacement $\pm1/1024$ for each transfer. A scalar root solve matches the calculated finite volume on each arm, independent of the area and entropy reads. The certificate below supplies an additional exact matched-volume curve. Each acquisition uses $2^{55}$ aggregated **simulated** starts per setting, with a separate multinomial draw for the resolved interior word. All reported intervals have at least 99% simultaneous coverage under the enrolled sampling model. The union bound covers at most 2,000 count marginals and sums, using half-width $\sqrt{2XL}+2L$ at count $X$, $L=\log(4000/0.01)$. The clock intervals include the independent $5.3\times10^{-8}$ native-kernel error on the conditional success probability. The finite secant calibration lies in $[0.0003535,0.0006258]$.

| 80-cell held transfer | Interior entropy secant, nats per unit | Rim area secant, squared clock units per unit | Calibrated residual, nats per unit |
| --- | --- | --- | --- |
| Zero to two occupied rim ports | $[-1.794,-1.788]$ | $[-41830,-40134]$ | $[-27.97,-15.97]$ |
| Fixed mean native energy and mean rim load | $[0.4283,0.4333]$ | $[-21545,-19850]$ | $[-13.06,-6.58]$ |

In the first held row the two measured changes have the same sign for every positive area coefficient. The second row also preserves the independent mean rim-load reading. For both rows the measured volume-difference interval contains zero. Source preparation, declared compensating words and every early/late/failure count are archived. [Acquisition](acquire.py) · [Calibration, held controls and counts](counts.json)

## Exact failure for the entire registered clock record

On the 80-cell uniform source form a linear readout matrix on the 276 word probabilities. Give it rows for normalization, native mean energy, nine independent launch flags, twelve successful return flags and the indicator of $k=2$. The nine launch flags span all twelve since $h_{ab}-h_{ba}=\langle n_a-n_b\rangle$. Together these rows have rank **24**. This is a statement about the true native propagator with the enrolled timeout.

The verifier quantizes proposals for the eigenvectors and eigenvalues of the half-integer Hamiltonian. Exact integer matrix products bound their orthogonality defect and eigen-equation residual. Interval sines and cosines bound the phase entries. Duhamel's formula transfers these errors to $e^{-iH/2}$; the row $\ell^1$ error of the proposed probability matrix is below $6.45\times10^{-12}$. To compute a discounted first entrance, solve $x=zPx$ outside the target with prescribed entrance values. The inverse has maximum-row-sum norm at most 4096. Entrance uncertainty in the return leg propagates with gain at most one. Exact integer residuals therefore bound each complete-flight probability error by $5.3\times10^{-8}$. Finally a dyadic proposed right inverse $B$ satisfies $\|AB-I\|_\infty<0.001305$ for the unknown exact table $A$. [Independent integer and interval calculation](certificate.py)

Take $D=B(AB)^{-1}e_{24}$. Only the $k=2$ probability changes, by one unit, among the 24 rows. In particular all three terminal probabilities for each physical route stay exactly fixed along $p(t)=p_*+tD$. The bound $\|D\|_\infty<1.877$ makes both signs physically available for $|t|<0.0019306$. Since $3\,dP_1+4\,dP_2=0$, the bulk entropy slope along this route is

\[
\frac{dH(I)}{dt}=\frac{4\log 4-3\log6}{3}=0.0566330122\ldots .
\]

The sign follows already from $4^4>6^3$. Concavity of the interior entropy gives $H(I;p_*)-H(I;p(-1/1024))>0.0000553$ nats. All twelve full three-outcome terminal laws stay fixed, so every area and volume obtained solely from those laws stays exactly fixed, including the all-attempt version. The independent mean rim occupancy changes by $2t/3$ and helps explain how the exterior information meter distinguishes the stocks. A practical arm uses a separate dyadic lottery with $2^{50}$ tickets and all weights deposited. Against fresh uniform records, its measured entropy change is in $[-0.0001013,-0.0000917]$ nats. Each true native success-probability difference on this approximate arm is below $7.48\times10^{-10}$. The exact corrected curve and the deposited finite-ticket approximation are separate constructions. [Exact bounds](certificate.json) · [Finite-ticket arm](counts.json)

## Does this patch form a physical ball?

Four vertices alone fit their own six positive-definite ranges by construction. Challenge that fit with a fifth actual station, port zero, whose distances to all four rim stations come from the same resident clock. Its squared distance from the three-dimensional affine span needed to reproduce the five-point Gram table is greater than **36,000 squared clock units**. A common flat three-dimensional fit would give zero. We also tried every one of the 20 interior ports as the center of a sphere with the four marked rim ports on its surface. At even the best center the largest and smallest measured radius differ by **more than 71.39 clock units**. Both lower bounds include the native wave and return-solve uncertainties. Curved geodesic distances and centers between registered sites require additional measurements. These finite tests preclude certifying this resident-clock tetrahedron as the required independently measured proper ball of the entropy theorem. [Held geometry bounds](certificate.json)

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=code:.
pkg=model.experiments.interval_echo_einstein_station.equilibrium_refinement.area_survey
python -m "$pkg.measure"
python -m "$pkg.certificate"
python -m "$pkg.acquire"
python -m pytest -q --pyargs "$pkg"
```
