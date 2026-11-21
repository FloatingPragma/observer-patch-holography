# Preserve the complete finite Gibbs response

This contribution to [the standing evidence audit #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
repairs evaluation of an existing finite Gibbs family. It follows the
[Hamiltonian assembly audit](HAMILTONIAN_ASSEMBLY.md): retaining an interaction
in the Hamiltonian does not ensure that a subsequent response calculation
retains its effect.

## Reproduced failures

For `A=diag(1,-1,1e-20)`, `B=diag(1,1,-2)` and zero multipliers, the state is
exactly `I/3`. Trace arithmetic gives `Cov(A,B)=-2e-20/3`. The old public
Duhamel evaluator returned zero: binary64 observable centering erased the
small component. Controls at `1e-100` and `1e-300` reproduced the failure.
A zero coupling in the Hessian cannot be inferred from this result.

For `T=diag(0,1e150)` and multiplier `7.45e-148`, the old evaluator returned
variance `4.940656458412465e-24`. Independent 500-digit evaluation of the
two-level partition function gives `2.822350730472012e-24`, a relative error
of approximately 75.05%. A subnormal excited population with few significant
bits was multiplied by the squared observable scale, producing an inaccurate,
normal-range answer. Positivity or an analytic derivative ratio misses this
error. The [maintainer's #1048 review](https://github.com/FloatingPragma/observer-patch-holography/pull/1048#pullrequestreview-5423083959)
identified the same mechanism in the entropy response API.

The state constructors also accepted such populations in individual states,
conditional sector states, sector weights and joint spectral weights. Two
normal factors can have a subnormal product, so local checks do not suffice.
The initial regression commit retains 31 failing cases before the repair.

## One response form, several existing uses

Fix finite Hermitian observables `T_a`, real multipliers `lambda_a`,
`H=sum_a lambda_a T_a`, `F(lambda)=log Tr exp(-H)` and `rho=exp(-H)/Tr exp(-H)`.
These are supplied data, not a selection of physical temperature.
Differentiating the exponential under the trace gives

```text
partial_a F = -Tr(rho T_a)
C_ab = partial_a partial_b F
     = integral_0^1 Tr(rho^s T_a rho^(1-s) T_b) ds
       - Tr(rho T_a) Tr(rho T_b).
```

This is the standard Bogoliubov/Kubo–Mori response; see
[Petz and Toth, equations (1)–(2)](https://dept.camden.rutgers.edu/math/files/Toth_036.pdf).
Let `p_i` be the thermal eigenvalues and `A_a=V* T_a V`. With logarithmic
mean `L(x,y)=(x-y)/(log x-log y)` and `L(x,x)=x`, expand the integral and pair
conjugate off-diagonal entries:

```text
C_ab = sum_{i<j} p_i p_j (A_a,ii-A_a,jj)(A_b,ii-A_b,jj)
       + 2 sum_{i<j} L(p_i,p_j) Re(A_a,ij conjugate(A_b,ij)).
```

The diagonal identity follows by expanding the pairwise sum and using
`sum_i p_i=1`. Both terms are real Gram forms with positive weights. For real
`v`, let `T(v)=sum_a v_a T_a`. Then `v^T C v=0` exactly when all diagonal
entries of `V* T(v) V` coincide and all off-diagonal entries vanish, equivalently
when `T(v)` is scalar. Thus the exact finite Hessian is strictly positive on
the observable span modulo identity. Duplicated constraints and scalar origins
give precisely the null directions in the existing MaxEnt uniqueness argument.
This does not identify exact null directions by thresholding a numerical Hessian.

If `S_a=sum_b M_ab T_b+c_a I`, old multipliers are `lambda=M^T mu` and
`C_S(mu)=M C_T(lambda) M^T`. Energy origins vanish; units act by congruence.
Centering and rescaling are legitimate algebraically but can still fail if
their intermediate arithmetic loses data.

The connection to the relative-entropy and projection audits is also exact:

```text
D(rho_lambda || rho_(lambda+delta))
  = F(lambda+delta)-F(lambda)-grad F(lambda) dot delta
  = (1/2) delta^T C(lambda) delta + O(||delta||^3).
```

The first equality substitutes the Gibbs logarithms; the second applies
Taylor's theorem to the finite analytic family. The same form is therefore
the susceptibility, MaxEnt objective Hessian and local relative-entropy
curvature. These identities are not independent tests of an implementation.

For `H=beta Z`, the commuting response is `C_ZZ=sech(beta)^2`, whereas the
transverse responses are `C_XX=C_YY=tanh(beta)/beta`. The first becomes
exponentially small; the others decrease inversely with the gap. Replacing
all quantum responses by a classical variance gives the wrong Hessian.

## Implementation and numerical boundary

`gibbs_response.py` evaluates the pairwise Gram formula. The public
`duhamel_covariance` starts with supplied entries and multipliers, removes
scalar origins, assembles and diagonalizes the thermal operator, and forms
weighted moments in a private multiprecision context. It never constructs
a binary64 state as an intermediate. Only the final covariance is converted.

Working precision has a 550-digit floor and increases with input dynamic
range. Evaluation with 40 additional digits must agree to relative `1e-25`
entrywise. Nonzero output underflow, overflow, insufficient relative output
precision and disagreement raise. There is no absolute correlation floor,
artificial spectral mass or Hessian ridge. Long-path controls retain responses
below `1e-180` after cancellation of terms of order `1e300`; the longer
unresolvable control must raise.

These are numerical safeguards, not interval bounds or a complete exact-zero
decision. Neither agreement of two precisions nor a returned zero proves an
algebraic identity. Very ill-conditioned inputs may require a separate
certified calculation. The kernel theorem concerns the exact finite family;
the optimizer's numerical rank and Hessian-resolution gates remain necessary.

### Audit correction: algebraic zeros must survive precision refinement

The maintainer-style audit found that the first PR head rejected ordinary
orthogonal observables at `I/d`: cancellation residues at two precisions
disagreed even though the exact covariance was zero. The first retained audit
regressions reproduced 27 failures across dimensions 3, 5 and 6, including
nonzero multipliers whose Hamiltonian sum is scalar and the Newton helper.
The same defect appeared inside degenerate thermal sectors, including after
exact real and complex dyadic changes of basis.

For a scalar Hamiltonian, the evaluator now computes
`Tr(A B)/d-Tr(A)Tr(B)/d^2` in rational arithmetic on the Hermitian parts of
the supplied entries, rounding only the result. Scalar-Hamiltonian membership
is itself checked by exact assembly. The internal uniform-spectrum Newton
case uses the same trace formula. This preserves both algebraic zeros and
representable small nonzero values without a tolerance-based exception.

For a nonuniform family that fails numerical precision agreement, the
public evaluator can certify individual zeros when the supplied Hamiltonian
has a rational spectrum. Its exact characteristic polynomial supplies all
eigenvalues and multiplicities; an incomplete rational spectrum supplies no
certificate. The exact spectral projectors are

```text
P_r = product_(s != r) (H-h_s I)/(h_r-h_s).
```

They are computed in canonical Gaussian-rational arithmetic, including
degenerate eigenspaces. Write `w_r=exp(-h_r)`, `d_r=Tr(P_r)`,
`a_r=Tr(P_r A)`, `b_r=Tr(P_r B)` and `Z=sum_r d_r w_r`. Then

```text
Z^2 C(A,B)
 = Z sum_(r,s) L(w_r,w_s) Tr(P_r A P_s B)
   - (sum_r w_r a_r)(sum_s w_s b_s).
```

This is the preceding integral formula with the exact spectral resolutions
inserted. For `r=s`, the logarithmic mean is `w_r`; for `r!=s` it is
`(w_r-w_s)/(h_s-h_r)`. Pair conjugate terms. The numerator is therefore a
finite sum of rational coefficients times `exp(-h_r-h_s)`. Grouping equal
exponents and proving every coefficient zero is a sufficient exact zero
certificate. This proof is independent of a selected eigenbasis and keeps
both the mean-product subtraction and interenergy coherences.

Only entries with that certificate are replaced by zero during a retried
precision check; every other entry must still agree and obey output-range
checks. This is not a general zero-decision procedure for irrational spectra.
Independent conditional-product and block-trace constructions verify the
certificates; nearby nonzero coherences down to `1e-300`, transverse variances
and an unsupported irrational spectrum prevent permissive acceptance.

Newton iteration reuses this Gram evaluator on its supplied binary64 spectrum
and vectors. It does not recover eigensolver accuracy lost upstream; replay
and residual checks still govern acceptance. The public covariance instead
recomputes the family at higher precision. Both interfaces state this distinction.

`gibbs_state` and `gibbs_sectors` still return binary64 states and now
conservatively refuse all subnormal thermal populations, including joint sector
populations. This is a numerical domain restriction, not a theorem that all
subnormals have poor precision. General state-validation APIs remain unchanged.
Existing dense-support checks remain. A weighted derivative can be resolved
when the corresponding binary64 state is not; the response API handles that
case without claiming to supply the state.

## Independent controls and impact

Controls use exact uniform-state trace arithmetic, 500/700-digit two-level
partition functions, shortest-walk expansions on coherent paths, and a separate
matrix-exponential Fréchet derivative implementation
([SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm_frechet.html)).
They cover complex noncommuting observables, extreme units, affine/unitary
changes of coordinates, scalar/duplicate null directions and all four
state-population boundaries. Both multiplier signs and gaps through 1000
are exercised. Existing malformed/masked/Boolean and lossy-conversion tests
continue to apply to both public entry points.

Seventeen isolated mutations are rejected: the old response, constant zero,
insufficient precision, removal of either Gram term, an arithmetic-mean kernel,
missing complex conjugation, early population rounding, discarded original
entries, disabled output/refinement checks and accepted subnormal states;
the audit also removes either exact-zero certificate, certifies every entry,
and deletes the mean-product or interenergy terms from the certificate.

| Consumer or surface | Impact |
| --- | --- |
| MaxEnt inference and #539 | Shared Gram evaluator in Newton steps. Live receipt regenerated; all five scientific checks remain true. |
| Existing nonclosure witness | Generic defects remain about `0.000181893249344` and `0.000250490780185` nats. The transverse product subfamily remains closed within numerical error. |
| Live numerical receipt | Multipliers/Hessian/closure values change at roundoff scale. Tiny residuals and derived bounds can change relatively more because they are already about machine precision. |
| Direct-sum collars | Refuses unresolved conditional, sector and joint populations; ordinary partition-weight and interaction controls still pass. |
| Entropy-response repair #1048 (merged) | Its conservative population refusal remains compatible. Its returned-state contract is unchanged here. |
| Information, geometry, #1049 and #1052 | No information, exact-support, recovery or entropy-balance theorem is changed. Combined integration is checked separately. |
| Papers, book, registry and public claims | Unchanged: finite response identities retain their stated premises. No physical energy, clock or thermal law is selected. |
| Frozen registrations and pinned evidence | Unchanged. Only the explicitly live #539 receipt is regenerated. |

Reproduce from the repository root with `PYTHONPATH=code`:

```text
python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
python code/maxent/maxent_closure_acceptance.py
```

The final command regenerates the live receipt, not frozen registrations.
This repair does not establish either physics task #1025/#1026.
