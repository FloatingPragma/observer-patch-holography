# Clocked receiver on the asymmetric repeated fans

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

On the three- and four-fan routes of the [asymmetric terminal control](../fan_path/TERMINAL_CONTROL.md), the hub and doublet receiver's records have positive partial transpose at both archived appointments. If observers at all other fans also verify their own local ring vacua and compare those outcomes, the accepted receiver records have negative partial transpose at both appointments. The probabilities and partial transpose entries below retain the original-launch denominator for positive signed-load preparations. These are simulated records on the supplied written-field guide. [Receipt](asymmetric_clock.json)

| Fans | Exposure | Receiver-only success range | Receiver-only Choi minimum | Success range with all remote vacuum reads | Choi minimum with all remote vacuum reads |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 27.0625 | 1.4723%–3.0518% | $+7.7130\times10^{-4}$ | 0.031753%–2.202311% | $-1.3222\times10^{-3}$ |
| 3 | 30.8125 | 1.4400%–2.5668% | $+6.1436\times10^{-4}$ | 0.028449%–1.165814% | $-9.1058\times10^{-4}$ |
| 4 | 27.0625 | 0.67343%–2.24261% | $+7.2113\times10^{-4}$ | 0.073796%–0.327778% | $-7.7764\times10^{-4}$ |
| 4 | 30.8125 | 0.88889%–2.03900% | $+7.8664\times10^{-4}$ | 0.022932%–0.519607% | $-5.4579\times10^{-4}$ |

The published receiver analyzes the two ring modes at the last occupied hub while holding a reference untouched at the launch. On separate trials, each other fan's observer measures whether its local ring degree of freedom has zero energy, and the receiver requires agreement from all of them. The appropriate local words are the spectator sectors behind the crossing. They can differ from the launch sectors. The audit checks that each such ring Laplacian is connected with uniform zero mode, so the local outcome is a rank-one projector. Joint acceptance gives a single operator on the two source amplitudes. The interval for its worst-input success excludes zero, and the partial transpose interval lies below zero. An observer needs all the remote vacuum verdicts to recover this heralded channel. Preparing a mobile doublet below the parked independent-ring threshold is a separate test.

## Independent error bound

Keep the original times $433/16$ and $493/16$ with the written bridge rate 4. The check reads the same source columns and local receiver calibrations and propagates under $H$ using Taylor polynomials on steps $\delta=1/16$. It checks that $2H$ is an exact symmetric integer sparse matrix, and its maximum absolute row sum determines an exact bound $M$ for $\|H\|_2$. The initial columns have exact dyadic coefficients.

Let $x=\delta M\le 5/2$. Each step keeps powers through degree $m=26$. The integral Taylor remainder for a self-adjoint generator, whose real-time evolution is unitary, has operator norm at most $x^{m+1}/(m+1)!$. This also applies to the Frobenius norm of the block of two inputs.

Write $b$ for the largest number of nonzero entries per row and $\gamma_k=ku/(1-ku)$ with binary64 unit roundoff $u=2^{-53}$. In either real component, a sparse row dot product, scalar division and multiplication for one term have forward error bounded by $\gamma_{2b+4}$ multiplying the corresponding absolute terms. The matrix of coefficient magnitudes has norm at most the same row-sum bound. Let $w_j=x^j/j!$. Relative to the Frobenius norm of the input columns, the total roundoff in the step is bounded by

$$
(1+\gamma_{m+1})\sum_{j=0}^{m} w_j(1+\gamma_{2b+4})^j-\sum_{j=0}^{m}w_j.
$$

The rational computation adds a $2^{-920}$ absolute underflow allowance. At most 4,096 rows, at most 64 terms per sparse product and $x\le 5/2$ keep the accumulated contribution of underflow below this allowance even if subnormals are flushed. If $\eta$ includes truncation, relative roundoff and this allowance, the exact solution has Frobenius norm $\sqrt2<2$ and after $s$ steps the error is at most $4s\eta$ for the verified $s\eta<1/2$. This follows by summing errors with growth no greater than $1+\eta$ per step.

Both receiver operations are contractions. If $\varepsilon$ bounds the propagated block error, the receiver's success effect changes by at most $3\varepsilon$. The Choi matrix is a quadratic product of that block with factor $1/2$. Its partial transpose preserves the Frobenius norm, so its least eigenvalue changes by at most $2\varepsilon$. The implementation adds an absolute $10^{-11}$ pad to both bounds for receiver contractions, construction of the matrices and the final dense eigensolves. Each small eigensolve is checked by its Gram defect and reconstruction residual before the computed values are used. The pad covers both residuals under polar correction with room for the binary64 products, projections on at most 64 spectator words, scalar normalizations, and reporting.

From the repository root, replay the probability enclosures with

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.asymmetric_clock
~~~
