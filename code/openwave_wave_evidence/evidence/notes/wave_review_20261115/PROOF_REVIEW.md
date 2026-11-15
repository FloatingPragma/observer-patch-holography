# Written-proof assessment: the common free-wave construction

Subject: the committed `shared/unified_wave_submission_20261115` bundle at
`60519f4`. The reproduction receipt records the full commit and manifests.
This is a written review with independently implemented finite checks. It is
not a Lean build, external-group G6 confirmation or physical validation.

## Conclusion

I found no blocking mathematical error in the checked **free Gaussian
reconstruction** chain under its stated architecture, preparation, observable
dictionary and ordered limits. The chain is substantive: the nonlocal
covariance and its poles depend on the native transition law, while the curl
alone supplies neither. Its platform admissibility and its correspondence to
physical observations remain decisions for the full outside review.

This conclusion concerns rows 29/30 and the `c=0` source correspondence.
The complete interacting theory, orbital spectral proof and moving charged
apparatus have not received a full mathematical review here.

## Source map

All paths below are relative to the bundle.

| Step | Bound source | Review |
| --- | --- | --- |
| One-unit transition and Gaussian moments | `evidence/free/evidence/openwave/shared/native_occupation_fluctuation_limit/PROOF.md`, §§1–3 | Finite connected weighted graph, fair endpoints, disclosed uniform reference. The transition is the convex mixture of seam averaging projections. Stationary occupation correlations yield the covariance below. The simple isolated eigenvalue argument gives a fixed-graph CLT and convergence of every fixed mixed moment. |
| Native masses and controlled slice limit | `evidence/free/evidence/openwave/row_30_klein_gordon/work/2026-11-14_codex_slice_refinement/PROOF.md`, §§3–7 | The internal gaps and multiplicities follow from the icosahedral Laplacian. The harmonic recurrence, `asinh` pole, quadratic error bound and group-speed bound are consistent. Their units are the declared common coordinate/slice units. |
| Physical Maxwell quotient | `evidence/free/evidence/openwave/row_29_em_waves/work/2026-11-14_codex_record_maxwell/PROOF.md`, §§2–5 | Four independent occupation channels, local midpoint curl, odd electric/even magnetic reflection, and the stated positive-support domain give a positive rank-two quotient. Positivity is on the field-strength algebra, not on all four vector-potential components. |
| Volume, refinement and locality | `evidence/free/evidence/openwave/row_30_klein_gordon/work/2026-11-14_codex_continuum_reconstruction/PROOF.md`, §§1–4 | The fixed-graph moment limit precedes volume and spacing limits. The integrable four-dimensional infrared singularity, explicit shell bound and smooth test functions control the covariance limit. Wick pairing transfers fixed-moment errors. The resulting standard free scalar/Maxwell mass-shell representations have positive energy and a common causal cone. |
| Same-record clock | Same continuum directory, `SOURCE_CLOCK_JOIN.md`; first-hit directory, `PROOF.md` | The killed-cycle boundary problem gives the stated finite-period endpoint transform and cost. The finite-cycle transform differs from the infinite-line pole; their limit is controlled. The clock is a Euclidean stopping-record clock, not a unitary primitive update or laboratory clock. |
| Unified candidate | `evidence/unified/PROOF.md`, §§2–5; `instrument.py`; `PIN_SPEC.md`; `READOUT_SPEC.md` | At `c=0` the outer event always executes, regardless of complete inner outcomes. Marginalizing those outcomes preserves the full outer-history law, hence all field and stopping observables jointly. The conserved coupling bit is preparation, not the changed-law control. |

## Checks of the delicate steps

### 1. Contact term and finite windows

For the normalized centered occupation sum, reversibility gives

```text
P = I - L/(2W)
C_n = [Pi + 2 sum_{k=1}^{n-1}(1-k/n) Pi P^k Pi] / (4W)
C_infinity = L^+ - Pi/(4W).
```

On a nonconstant transfer eigenmode `z`, the difference is
`z(1-z^n)/(2W n(1-z)^2)`. This is a population expectation, not a
finite-sample confidence interval. The occupation-only candidate retains
the contact subtraction. Its finite moment-doubling target must retain it
too; optional auxiliary-coin completion is a different protocol.

The bound common member and readout specification make that choice explicit.
The earlier theorem's phrase "finite Green matrix" should always be read
through that specification, rather than silently targeting bare `L^+`.

### 2. The limiting procedure is constructive but not efficient

The finite-state analytic perturbation argument applies with the graph fixed.
It does not establish a uniform CLT rate as volume increases. The separate
exact-moment recurrence and doubling search provide a terminating finite
procedure for each fixed menu, using the actual contact-retaining target.

The volume bound uses
`(2j+1)^4-(2j-1)^4 = 64j^3+16j`; summing the infrared shell contribution is
bounded by a constant times `rho^2`. Away from the origin, finite derivative
bounds control Riemann sums. These ingredients justify a computable
diagonal, but do not supply a practical large-volume runtime estimate.

### 3. Field normalization and common speed

With `x=az`, the discrete inverse kernel supplies a factor `a^4` in the
Fourier integral. Multiplying each projected field by `a^-2` and smearing
with `a^4` therefore gives the stated continuum covariance. Both sectors
use the same external weight `a^-2` and the same slice convention.

For `S=lambda+|p|^2`,

```text
0 <= S - E_a^2 <= a^2 [sum_j p_j^4 + S^2]/12
E_a = 2 asinh(a sqrt(lambda+q_a)/2)/a.
```

The internal gaps are `0`, `5-sqrt(5)`, `6`, `5+sqrt(5)` with
multiplicities `1,3,5,3`. They are not empirical particle masses.

### 4. Maxwell positivity requires the actual quotient

Writing `C=i[k cross]`, `Q=|k|` and `P_T=I-kk^T/Q^2`,

```text
G = [[Q^2 P_T, -Q C], [-Q C, Q^2 P_T]]
  = [Q P_T, -C]^* [Q P_T, -C].
```

Here the right side denotes the Gram matrix of that 3-by-6 map.
Since `C^2=Q^2 P_T`, its rank is two for `Q>0`.
The polynomial Gaussian lift uses symmetric tensor powers of this form.
The constant potential mode is killed by curl. Finite differentiated
contacts cannot be dropped at the reflection plane; the proof explicitly
uses separated supports and then the controlled limit.

The finite `N=3` electric covariance has a small negative contact
eigenvalue. It does not establish exact finite reflection positivity.
Conversely, that retained artifact does not contradict the positive
limiting theorem.

### 5. Clock transfer and physical interpretation

Solving the cycle cut at its absorbing plane gives

```text
h_H(A) = cosh(omega(H/2-1))/cosh(omega H/2)
omega = 2 asinh(sqrt(A/w)/2)
E[tick attempts] = W(H-1)/w.
```

Successive hits compose by strong Markov and cyclic translation. At finite
period the clock energy `-log(h_H)/a` is not the infinite-line `omega/a`.
The source documents preserve this distinction and the cost divergence.

Continuing a positive Euclidean reconstruction into a Lorentzian free field
is legitimate mathematics. Whether the declared dictionary and limiting
observation family meet OpenWave's own-process criterion is the explicit
admission question. Primitive event-time unitarity, a fitted laboratory `c`
or an electron identity are not added here as extra wave-row requirements.

## Scope that the submission must preserve

- Four-axis routing, four independent channels, uniform preparations,
  refinement weights, spectral projection, tensor/scalar reflection and
  the curl/readout conventions are declared architecture and apparatus.
- The `100,000,000`-event covariance is an analytic population calculation.
  The common free record has four actual `32,768`-event tapes.
- Finite records corroborate source/decoder/control identities. The full
  Gaussian continuum and Lorentzian laws depend on the written limit proofs.
- Nonconstant modes have a scalar KG interpretation; identifying them with
  measured species or implementing an interacting detector is a further claim.
- The biased-coin control changes the actual source law. Turning the coupling
  preparation bit on/off is not that control.
- The static signed-source submission is a different candidate. Its charge
  and calibrated rest-energy results do not transfer to this construction
  merely because both call the same endpoint function.
