# Calibration Lane

This lane exists to answer one exactness question cleanly:

If the current D10 branch misses electroweak calibration observables, is that
because `P` needs more digits, or because the D10 transport / matching package
is not yet the exact closure branch?

The immediate tool is an implied-`P` audit. For each D10 observable, solve for
the pixel constant that would make the current implementation hit the declared
target exactly. If different observables imply different `P` values, then the
present bottleneck is not just decimal precision in `P`.

The second tool is a scheme-freezing artifact for the current D10 branch. That
artifact records the exact single-`P` running electroweak family already
realized by the code and separates it from the mixed run/pole reporting
surface. If pole/effective reporting is required, the remaining exact missing
object is then one common `EWTransportKernel_D10`, not more digits of `P`.

## Active Calibration Files

The active calibration scripts open with the same short derivation header:

- `Chain role`: where the file sits between the D10 core, the reduced
  two-scalar carrier, the selector, and the public readout
- `Mathematics`: which fixed-point, transport, or Jacobian step is being used
- `OPH-derived inputs`: which calibration quantities come directly from the D10
  core already emitted in `/particles`
- `Output`: which downstream calibration artifact the file is responsible for

For the live branch, the main path is:

- `derive_d10_ew_observable_family.py`
- `derive_d10_ew_source_transport_pair.py`
- `derive_d10_ew_population_evaluator.py`
- `derive_d10_ew_exact_closure_beyond_current_carrier.py`
- `derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py`
- `derive_d10_ew_tau2_current_carrier_obstruction.py`
- `derive_d10_ew_exact_wz_coordinate_beyond_single_tree_identity.py`
- `derive_d10_ew_exact_mass_pair_chart_current_carrier.py`
- `derive_d10_ew_repair_branch_beyond_current_carrier.py`
- `derive_d10_ew_source_transport_readout.py`
- `derive_d11_declared_calibration_surface.py`
- `derive_d11_forward_seed.py`
- `derive_d11_forward_seed_promotion_certificate.py`
- `derive_d11_fixed_ray_no_go_theorem.py`
- `derive_d11_live_exact_higgs_promotion.py`
- `derive_d11_live_exact_split_pair_theorem.py`

The live D10 split is explicit:

- builder-local current-carrier frontier:
  `EWExactMassPairSelector_D10`
- broader honest exact-PDG frontier:
  `D10RepairBranchBeyondCurrentCarrier`
- strongest strictly smaller constructive primitive beneath that broader frontier:
  `EWSinglePostTransportTreeIdentity_D10`
- shared D10/gravity bridge beneath the local gravity row:
  `D10GravitySharedEdgeEntropyBridge`

So the current carrier is not treated as the whole D10 story. It closes its own
exact chart, the W/Z mass lane is carried by the closed target-free source-only
mass theorem above that carrier, and the electromagnetic row is read on the
Ward-projected `U(1)_Q` transport family anchored at
`alpha_em^-1(m_Z^2)=128.30576920234813` with Thomson endpoint
`alpha^-1(0)=137.035999177`. The older repair-branch artifacts stay on disk as
scaffolding and validation layers rather than as the public electromagnetic
theorem.

The local gravity bridge is carried on the stated local extension surface. On
the lifted product presentation of the realized quotient branch, the
`R`-sector collar entropy contribution satisfies
`(L_C)|_R = log d_R = log d_R3 + log d_R2`, because the `U(1)` factor is
one-dimensional. The resulting shared scalar is
`ellbar_shared = ellbar_SU(2)(t2_run) + ellbar_SU(3)(t3_run)`, and when the
same branch satisfies the D10 pixel law this becomes `ellbar_shared = P/4`.
That is the surface on which the local gravity readout uses the same emitted
scalar as the D10 pixel law.

## Commands

Run the current implied-`P` audit:

```bash
python3 particles/calibration/implied_p_consistency_audit.py
```

This writes:

- [`particles/runs/calibration/implied_p_consistency.json`](/Users/muellerberndt/Projects/oph-meta/particles/runs/calibration/implied_p_consistency.json)

Run the local calibration guard:

```bash
python3 particles/calibration/test_single_p_consistency.py
python3 particles/calibration/derive_d10_ew_observable_family.py
python3 particles/calibration/test_d10_observable_family_artifact.py
python3 particles/calibration/derive_d10_ew_transport_kernel.py
python3 particles/calibration/test_d10_ew_transport_kernel_artifact.py
```

That guard checks the solver mechanics and the presence of the current D10
calibration observables. It is not a claim that exact single-`P` closure has
already been achieved.

The D11 lane is split into one lower-rank companion branch and one live
source-only pair theorem.

- The shared forward seed `sigma_D11_HT = alpha_u * cos(2*theta_W0) / sqrt(pi)`
  stays on disk as the diagonal fixed-ray companion branch on the declared
  D10/D11 surface. Its fixed-ray certificate proves `pi_y = pi_lambda`,
  `eta_HT = 0`, and `w_HT = 0` on that one-scalar branch only.
- The live pair theorem is `D11SourceSplitForwardExactness`, emitted by
  `derive_d11_live_exact_split_pair_theorem.py`. It uses only the forward D10
  tuple `(eta_source, beta_EW, lambda_EW, tau2_tree_exact, delta_n_tree_exact)`
  plus the declared D11 Jacobian surface.
- The shared split scalar is `rho_HT = log(1 + tau2_tree_exact)`. The source-only
  residual selectors are
  `R_T = -tau2_tree_exact * eta_source^2 + (1 + beta_EW/28) * eta_source^6 + eta_source^8/14 + eta_source^9/27`
  and
  `R_H = eta_source^5 - (3/25) * eta_source^6 + lambda_EW * eta_source^6 / 18 + eta_source^8 / (2 * beta_EW)`.
- The predictive split coordinates are
  `pi_y = (eta_source + (3/2 + beta_EW/4) * rho_HT + R_T) / sqrt(pi)` and
  `pi_lambda = (eta_source - (4/3 - beta_EW/54) * rho_HT + R_H) / sqrt(pi)`.
  On the declared D11 surface these emit
  `m_t = 172.3523553288312 GeV` and `m_H = 125.1995304097179 GeV`.
  The Higgs row lands on the PDG 2025 Higgs average at quoted precision.
  The exact public running-top row uses the PDG 2025 cross-section entry
  `Q007TP4`.
  The bridge to the auxiliary direct-top average `Q007TP` is open and tracked
  in [#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207).
- The compare-only exact Higgs/top inverse slice remains a validation surface.
  The target-anchored Higgs-only and top-side exactifier artifacts remain on
  disk as supporting witness surfaces beneath the source-only pair theorem.
