# Simulator contract: OPH issue #330 radial lift

## Scope

This contract separates four executable objects:

1. the source screen covariance;
2. the physical radial uniqueness theorem;
3. the branch-typed forward projection;
4. downstream observable transfer.

The first three may support an E4 primordial packet. The fourth is E5. A Planck-shape overlay, CAMB/CLASS run, or likelihood is never an ancestor of an E4 radial receipt.

## Required artifact split

Every run writes separate immutable artifacts:

```text
source/
  geometric_q_receipt.json
  source_release_amplitude.json
  screen_tilt.json
  screen_covariance.json
  physical_mode_basis.json
  source_embedding.json
  radial_dilation_intertwiner.json      # SOURCE_DILATION branch
  radial_cross_covariances.json         # RADIAL_TOMOGRAPHY branch
  radial_null_report.json
  radial_forward_residual.json
  radial_promotion_receipt.json

diagnostic/
  planck_shape_comparison.json
  camb_or_class_outputs/
  plots/
```

No file under `diagnostic/` may occur in the ancestor DAG of a file under `source/`.

## Branch selector

Exactly one of the following is declared:

```json
{"radial_branch": "SOURCE_DILATION"}
```

```json
{"radial_branch": "RADIAL_TOMOGRAPHY"}
```

```json
{"radial_branch": "PRIOR_CONTINUATION"}
```

`PRIOR_CONTINUATION` is never promoted as source-derived E4.

## Common finite forward operator

For a discrete physical spectral measure,

\[
A_{\ell j}=4\pi Z_q^2w_j^{\log k}|\Psi_{\ell}(k_j)|^2,
\qquad C=Ap.
\]

The run records:

- branch type and background-curvature status;
- \(k_j\), weights, units, and density-of-states convention;
- raw normalized radial window and hash;
- \(Z_q\) and theorem parent;
- each \(\Psi_\ell(k_j)\);
- all retained and held-out multipoles;
- numerical precision and quadrature convergence.

## `SOURCE_DILATION` receipt

The theorem-level object is not a fitted slope. The producer first exports the
scale-labelled source-refinement orbit and the commutative square

\[
D_{s,r}U_r=U_rR_{s,r},\qquad
C_{\zeta,r}=U_rC_{{\rm src},r}U_r^*.
\]

This gives the exact residual-transfer identity

\[
D_{s,r}^{-1}C_{\zeta,r}D_{s,r}-u_r(s)C_{\zeta,r}
=U_r\left(R_{s,r}^{-1}C_{{\rm src},r}R_{s,r}
-u_r(s)C_{{\rm src},r}\right)U_r^*.
\]

The finite objects must converge strongly on one embedded safe-band Hilbert
space, with a uniform covariance norm bound and a vanishing operator residual.
Only then does the continuum relation

\[
D_s^{-1}C_\zeta D_s=e^{-\theta s}C_\zeta
\]

hold on the common \(d\ln k\) physical mode basis.

The receipt records:

```json
{
  "receipt": "SCR330_RADIAL_DILATION_INTERTWINER_RECEIPT",
  "physical_mode_basis_id": "...",
  "safe_logk_band": [0.0, 0.0],
  "dilation_maps": ["..."],
  "scale_ratios": [1.0],
  "source_embedding_commutator_norms": [0.0],
  "screen_covariance_naturality_residual_norms": [0.0],
  "physical_operator_residual_norms": [0.0],
  "strong_covariance_cauchy_residuals": [0.0],
  "strong_dilation_cauchy_residuals": [0.0],
  "uniform_covariance_norm_bound": 0.0,
  "finite_to_continuum_passed": false,
  "off_band_leakage": 0.0,
  "max_absolute_log_residual": 0.0,
  "rms_log_residual": 0.0,
  "tolerance": 0.0,
  "refinement_sequence": ["..."],
  "passed": false
}
```

Pass conditions:

- the source mode space is the cofinal scale-labelled refinement orbit, not one angular cut;
- all maps are sourced from the same physical scale bridge;
- the embedding square \(D_{s,r}U_r=U_rR_{s,r}\) passes and its residual converges;
- source and coarse mode projectors commute within the declared residual;
- source covariance naturality and physical covariance naturality have the same residual under \(U_r\);
- the strong finite-to-continuum hypotheses and uniform covariance norm bound pass;
- covariance survival equals the screen source cocycle, not a separately chosen exponent;
- safe-band leakage and operator residual converge under refinement;
- finite diagonal checks \(\Delta^2(bk)=b^{-\theta}\Delta^2(k)\) pass;
- source DAG has no observational ancestor.

Negative controls:

- shuffled physical mode labels;
- wrong scale-ratio sign;
- a spectrum with a planted log-periodic wiggle;
- a separately fitted \(n_s\);
- a mode basis from a different source embedding.

Each control must fail.

## Thin-shell Mellin receipt

For `FlatExact` or `FlatAssumed` and a thin shell, compute

\[
I_\ell(\theta)=\frac{\sqrt\pi}{4}
\frac{\Gamma(1+\theta/2)}{\Gamma(3/2+\theta/2)}
\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}.
\]

The receipt records the convergence strip, every computed value, precision, and the arithmetic identity

\[
C_\ell^q=A_q\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}.
\]

The amplitude is source-derived from

\[
A_\zeta
=\frac{A_q}{\pi^{3/2}Z_q^2(k_\star R_\star)^\theta}
\frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}.
\]

The code must not solve this equation for \(A_\zeta\) using measured CMB \(C_\ell\).

## Finite-window receipt

The run computes the exact window transform and the certified bound

\[
\eta_{\ell,W}
=\frac{2\sqrt{J_\ell(\theta)}}{\theta}
\int W(dr)|r^{\theta/2}-R_\star^{\theta/2}|,
\]

\[
J_\ell(\theta)=I_\ell(\theta-2)
-\left[\ell(\ell+1)-\frac{\theta(\theta+1)}2\right]I_\ell(\theta),
\]

\[
|C_{\ell,W}-C_{\ell,R_\star}|
\le4\pi Z_q^2A_\zeta k_\star^\theta\eta
\left(2R_\star^{\theta/2}\sqrt{I_\ell}+\eta\right).
\]

The receipt stores \(I_\ell,J_\ell,\eta\), exact projected values, the bound, and the ratio of actual quadrature difference to the bound.

## `RADIAL_TOMOGRAPHY` receipt

The input must contain cross-covariances. Auto-spectra alone fail the receipt:

\[
C_\ell(r_i,r_j).
\]

The receipt records:

- radial basis and measure \(r^2dr\);
- finite spherical-Hankel matrix and unitarity residual;
- conjugated covariance \(H Q_\ell H^{-1}\);
- off-diagonal leakage in \(k\);
- recovered nonnegative multiplication spectrum;
- refinement convergence;
- reconstruction on held-out radii/windows.

A finite set of auto-windows is not tomography and must fail this receipt.

## Null-space receipt

For every finite operator, publish:

```json
{
  "shape": [0, 0],
  "singular_values": [],
  "rank_threshold": 0.0,
  "rank": 0,
  "nullity": 0,
  "condition_number_nonzero": 0.0,
  "right_null_basis_hash": "sha256:...",
  "resolution_kernel_hash": "sha256:..."
}
```

The raw right-null basis and resolution kernels are part of the evidence bundle, even when the source branch is one-dimensional. Their role is to prove that uniqueness comes from the source theorem, not from an accidentally square discretization.

## Prior-selected continuation

For a declared \(p_0,Q\), the exact representative is

\[
p_*=p_0+Q^{-1}A^T(AQ^{-1}A^T)^+(C-Ap_0).
\]

The run publishes \(p_0,Q,R_Q,N_Q\), positivity active set, and sensitivity to all declared prior variants. Its output type is `ConditionalRadialContinuation`.

## Forward residual

Every branch computes a forward residual without re-optimizing source parameters:

\[
r_\ell=C_\ell^q-4\pi Z_q^2\int d\nu(k)\Delta_\zeta^2(k)|\Psi_\ell(k)|^2.
\]

The residual artifact contains raw signed residuals, absolute and relative norms, numerical error budget, and held-out modes/windows. A failed residual blocks promotion.

## Curved branches

- `FlatExact`: exact flat formulas may be promoted if all other receipts pass.
- `FlatAssumed`: at most conditional physical status.
- `OpenCurved` / `ClosedCurved`: use the declared hyperspherical eigenfunctions and spectral measure; do not use the flat gamma conversion.
- `Unresolved`: fail physical radial promotion.

## E4 promotion

`SCR330_RADIAL_PROMOTION_RECEIPT` passes only if:

- all source/stress/clock/freezeout/phase/scale receipts pass;
- \(A_q,\theta,Z_q,W\) are source outputs;
- the complete null report and forward residual pass;
- exactly one uniqueness branch passes;
- the source DAG is clean;
- `physical_tt_te_ee_claim` is false.

## E5 firewall

A transfer artifact may consume the E4 primordial packet, but the E4 packet may not consume transfer outputs. TT/TE/EE, lensing, recombination, foregrounds, nuisance parameters, covariance, and likelihood live in a separate DAG and separate claim tier.

## Reference implementation

`oph_radial_lift_330.py` implements:

- exact \(I_\ell(\theta)\) and \(J_\ell(\theta)\);
- general-pivot/general-\(Z_q\) amplitude conversion;
- exact thin-shell gamma spectrum;
- finite-window quadrature and certified bound;
- finite projection matrix and SVD/null report;
- minimum-prior continuation;
- source-family forward residual;
- finite dilation residual checks;
- approximate-dilation shape bound;
- fail-closed receipts.

`pytest -q` must pass before a receipt bundle is accepted.
