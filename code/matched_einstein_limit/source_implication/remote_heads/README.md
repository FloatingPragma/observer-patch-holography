# Reading both stored heads from one unmarked port

On the 80-cell board, launch one excess at unmarked port 5 and ask whether it is there at a chosen time. The population return curve on any open time interval identifies the entire pair of four-by-four reversal blocks, including their separate weights. The read requires neither a reversal label on individual preparations nor the tetrahedral area. The active-law generator and time calibration must be known. This is an operational way to measure the stored state and the expectation of any separately specified observable on it. Identifying a local null stress still requires a physical choice of energy, momentum flux, spatial stress and frame. [Exact ranks and count records](results.json) · [Independently prepared source](../full_algebra/PROOF.md)

## Exact rank

Let \(H_b\) be the integer matrix for twice the native flight generator in parity sector \(b\), up to its irrelevant scalar. Write \(A_{b,j}=(H_b^j)_{55}\), a four-by-four return amplitude block. The measured \(n\)-th derivative of the binary return effect is \((-i/2)^n\) times the direct sum of
\[
 B_{b,n}=\sum_{j=0}^n(-1)^j {n\choose j} A_{b,j}^{\mathsf T} A_{b,n-j}.
\]
In each block the even orders are real symmetric and the odd orders are real skew before the factor of \(i\). Record the ten upper-triangular and six strict upper-triangular entries in disjoint columns. Integer propagation through order 40 gives 32 independent rows. The same minor has determinant 184 modulo 1009 and 688 modulo 1013. A nonzero residue proves full real rank. Equality of population curves on any open interval gives equality of derivatives by analyticity. Every state on the specified block algebra is therefore distinguished. Coherences between parity blocks are invisible to this instrument. Turning off the squared-contrast term during the probe keeps the same ports and read but leaves rank exactly five. The even flight then acts as identity on both coins and the odd flight as identity on the right coin. Five independent derivative rows attain their combined dimension bound. [Integer flight and changed-law check](measure.py)

## A conserved parked collar charge and the cost of interrogation

Use the selected faithful state \(\chi\) as the frozen reference. The left head includes its reversal sector. Set \(\pi=\chi_L\otimes\chi_R\) and \(K_{\rm collar}=\log\pi-\log\chi\), an observable inferred from the reference probability curves. Erasing the pairing gives the tangent \(\pi-\chi\). Its two marginal variations vanish, so its entropy derivative is exactly the change of this one observable,
\[
 \delta S=\mathrm{Tr}[(\pi-\chi)K_{\rm collar}]
 =D(\pi\Vert\chi)+D(\chi\Vert\pi)>0.
\]
The two relative entropies at the numerical reference are 0.120213 and 0.110893 nats. The adjusting right rotation preserves entropy and both marginal states; its infinitesimal change of this frozen charge is zero. It follows that \(\delta\langle K_{\rm collar}\rangle=-\eta\,\delta A|_V\) for this physical pairing operation at the exact certified maximum. The charged observable commutes with the parked Hamiltonian \(16 I+P_{\rm odd}\) on the retained support. Thus its parked conservation and calibrated scalar area response follow on the same preparation. [Parked energy projection](../null_operator/README.md) · [Fixed-volume certificate](../full_algebra/PROOF.md)

There is a constraint when a probe flies. Let \(P\) be its sharp launch-port projector at port 5. Suppose a full conserved operator \(C\) is compatible with this event, \([C,P]=0\), and restricts on that port to \(K_{\rm collar}\), up to a state-independent scalar. Since \(C\) commutes with every free flight \(U(t)\), its restriction must commute with every return block \(PU(t)P\), its adjoint and hence every binary return effect. These effects span the full two-block algebra by the rank proof. On the parity blocks only a scalar in each block commutes with them all. The unpairing tangent preserves both block probabilities but changes \(K_{\rm collar}\) strictly, excluding this extension. The numerical Hilbert-Schmidt norm of its commutator with the time-two return effect is 0.0034434. With an excess released at port 5 for two units and removed with position unread, the native population replay changes mean head charge by \(+0.0000914690\) nats and preserves both parity probabilities. This last value uses the numerical source and flight evolution.

A total stress could exchange charge with a probe, include hopping terms or be insensitive to these sharp port events under physical coarse graining. A regional transport law would have to account for that exchange in independent records on the same observer balls. The proved parked observable leaves this flux and its contraction with a null direction undetermined. [Finite operator lemmas](../../../../Lean/ObserverPatchHolography/EinsteinBranch/StationarityReadout.lean)

## Finite counts on the selected equilibrium

Use 32 quarter-tick appointments and all 16 ordered launch/detector pairs from ports 5, 40, 0 and 38, including coincident pairs. Each setting is a binary position observation on a fresh preparation. Every miss counts. These are generated binomial records at the archived numerical equilibrium. The exact optimizer is within \(2\times10^{-7}\) in Hilbert-Schmidt norm of the certified center. Exact rank applies to either source. The finite-count intervals below concern the simulated preparation under the numerically evaluated native flight probabilities.

At \(2^{36}\) trials per setting, ordinary unweighted inversion gives parked energy \(16.202432\pm0.000093\) (99% sampling bound), but the smallest inferred state eigenvalue is \(-0.2564\). Full rank alone gives little precision in the least sensitive directions. The separate high-exposure experiment uses \(2^{52}\) trials per setting, a source-independent weighting, and distinct count banks for the modular reference, baseline and two operations. This is an aggregate simulation with 512 binary settings per bank. Individual terminal observations have been sampled from their known population probabilities.

| Change from the common baseline | Regional modular response, nats | Parked exchange energy |
|---|---:|---:|
| Unread left reversal sign | \(0.4094535\pm0.07935\) | \(0.29756325\pm0.00000072\) |
| Erase pairing, retain both one-head distributions | \(0.2311063\pm0.07935\) | \(-0.00000004\pm0.00000072\) |

The table uses simultaneous 99% sampling bounds. An independent held baseline and separate intervention banks give the two affine contrasts against the measured operator \(-\log\hat\chi\). At the baseline each contrast against the true \(-\log\chi\) equals the first entropy derivative. For the pairing operation the ideal energy difference is exactly zero. Entropy preservation by the volume-adjusting unitary makes these also the modular derivatives used in the signed fixed-volume test. The tomography uses single terminal clicks. The area uses a separate two-way first-detection instrument. Neither instrument measures the Ricci tensor on a null ray.

Let \(p_i^{\max}=\|E_i\|\) for terminal effect \(E_i\). The trace-constrained linear inverse is weighted by \(1/\sqrt{p_i^{\max}}\). For any one dual row \(a\), the variance of its linear estimate on \(m\) independent banks of \(N\) trials per setting is at most \(m\sum_i a_i^2p_i^{\max}/N\); each centered summand is bounded by \(\|a\|_\infty/N\). Bernstein's inequality bounds the two-sided sampling error by the square root of twice this variance bound times \(L\), plus \(2L\|a\|_\infty/(3N)\). Use \(L=\log(2\cdot31/0.005)\) for the 31 reference coordinates and \(L=\log(2\cdot4/0.005)\) for the four held quantities. The orthonormal reference coordinates give Hilbert-Schmidt radius 0.000866; the smallest reconstructed eigenvalue is 0.022252. The logarithm resolvent then bounds the modular-operator error by \(-\log(1-r/\lambda_{\min}(\hat\chi))\). A contrast of two densities costs at most twice this operator error. The archived count banks include the unsuccessful low-exposure inversion.

From the repository root, run the population and count replay with the scientific Python environment.

~~~sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.remote_heads.measure
~~~
