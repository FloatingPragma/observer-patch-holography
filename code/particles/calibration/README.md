# Calibration Lane

This lane exists to answer one exactness question cleanly:

If the current D10 branch misses electroweak calibration observables, is that
because `P` needs more digits, or because the D10 transport / matching package
is not yet the exact closure branch?

The immediate tool is an implied-`P` audit. For each D10 observable, solve for
the pixel constant that would make the current implementation hit the declared
target exactly. If different observables imply different `P` values, then the
present bottleneck extends beyond decimal precision in `P`.

The second tool is a scheme-freezing artifact for the current D10 branch. That
artifact records the exact single-`P` running electroweak family realized by
the code and separates it from the mixed run/pole reporting
surface. If pole/effective reporting is required, the remaining exact missing
object is then one common `EWTransportKernel_D10`, not more digits of `P`.

Three compact boundary modules accompany that larger lane:

- `wz_strict_series.py` implements only the strict finite-order algebra for
  `Gamma(s) = s - m0^2 + kappa Pi(s)`. It keeps the one-loop coefficient,
  two-loop derivative and neutral-mixing terms, and square-root display
  coordinates separate.
- `strict_one_loop_pole_map/` is the subject-bound, fail-closed receipt
  package. It implements the equivalent convention `Gamma=s-m2-Delta` with
  `Delta=-Pi`, validates input and output schemas, reconstructs the declared
  charged and neutral one-loop poles, and independently rejects self-attested
  evidence, unrelated receipts, tolerance inflation, order leakage, and
  corrupted redundant fields. Its archived SMDR fixture is regression-only
  and can never self-promote.
- `source_completion_no_go.py` gives an exact two-completion witness for the
  current Higgs/top reduct. The same source coordinate admits distinct
  downstream Higgs and top values, so that reduct cannot support a unique
  source-only pole claim.
- `q1_q4_wz_theorem_oracle/` is the content-addressed, non-promoting oracle for
  the enlarged field-theory theorem package. Its `SM_QFT_*` DAG makes QFT-Q2
  and QFT-Q3 parallel descendants of QFT-Q1, places strict fixed-parameter W/Z
  algebra in QFT-Q3, and gives QFT-Q4 its own nonperturbative-tower and
  resonance-sheet inputs. It recomputes exact algebra only; every OPH producer
  remains fail-closed.

These modules are controls, not OPH-native mass predictions. They do not supply
the source-selected observer-like carrier, EFT matching and field census,
gauge/BRST identity checks, analytic continuation, uncertainty transport,
physical-current amplitudes, or operational clock required for promotion.

The complete strict-one-loop package is also a control under this boundary;
the enlarged oracle additionally checks the strict two-loop derivative and
neutral-mixing coefficients. Both close map implications after a complete
renormalized input packet is given, not the missing OPH-native antecedent.

The comparison-only convention layer is the
wz_experimental_convention.py module. It maps published
mass-dependent-width Breit--Wigner parameters exactly to the
s_pole=(M-i*Gamma/2)^2 convention, carries an analytic Jacobian, and
propagates a validated covariance. The frozen PDG-2026 fixture explicitly uses
a block-diagonal covariance assumption and disables a joint chi-square because
the mixed W/Z world-average object has no published joint covariance.

## Active Calibration Files

The active calibration scripts open with the same short derivation header:

- `Chain role`: where the file sits between the D10 core, the reduced
  two-scalar carrier, the selector, and the public readout
- `Mathematics`: which fixed-point, transport, or Jacobian step is being used
- `OPH-derived inputs`: which calibration quantities come directly from the D10
  core emitted in `/particles`
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
- `derive_d10_ew_quotient_transport_receipt.py`
- `derive_d11_declared_calibration_surface.py`
- `derive_d11_forward_seed.py`
- `derive_d11_forward_seed_promotion_certificate.py`
- `derive_d11_fixed_ray_no_go_theorem.py`
- `derive_d11_live_exact_higgs_promotion.py`
- `derive_d11_live_exact_split_pair_theorem.py`
- `derive_boson_source_prediction_closure_audit.py`

The live D10 split is explicit:

- builder-local current-carrier frontier:
  `EWExactMassPairSelector_D10`
- broader supported exact-PDG frontier:
  `D10RepairBranchBeyondCurrentCarrier`
- strongest strictly smaller constructive primitive beneath that broader frontier:
  `EWSinglePostTransportTreeIdentity_D10`
- shared D10/gravity bridge beneath the local gravity row:
  `D10GravitySharedEdgeEntropyBridge`

So the current carrier is not treated as the whole D10 story. It closes its own
algebraic chart. The stronger W/Z value law above that carrier is a
runtime target-free-input candidate, not a closed source theorem: the exact
conditional quotient-transport implication assumes the finite QT1--QT5
path certificate, and the present corpus does not emit that certificate or
bind it to one strict pixel branch. The electromagnetic row is read on the
Ward-projected `U(1)_Q` transport family anchored at
`alpha_em^-1(m_Z^2)=128.30576920234813`. The Thomson endpoint stored in this
lane is a reference comparison field, not a source-emitted prediction. The
older repair-branch artifacts stay on disk as scaffolding and validation layers
rather than as the public electromagnetic theorem.

The fail-closed composition receipt
[`boson_source_prediction_closure_audit.json`](../runs/calibration/boson_source_prediction_closure_audit.json)
keeps the rounded legacy carrier, public-endpoint pixel, source-audit pixel,
and compressed candidate trunk separate. It also distinguishes the current
carrier, the selected repair candidate, the reference-fitted inverse adapter, and
the conditional D11 surface. No W, Z, or Higgs complex-pole prediction is
promoted while the scale, QT1--QT5, RG, rigidity, source-DAG, pole, and
uncertainty gates remain open.
The executable guard
[`ew_physical_readout_gate.py`](ew_physical_readout_gate.py) enforces this
boundary for the July 2026 pole and two-loop packets. It requires a complete
renormalized-vev, tadpole, field-content, threshold, running, finite-part,
complex-pole, target-release, and covariance contract plus an observable
identity receipt. Any missing field or declared-open gate produces
`NOT_EVALUABLE`; the packet cells remain internal prescription diagnostics and
cannot emit a physical pull.
The dedicated conditional receipt
[`d10_ew_quotient_transport_receipt.json`](../runs/calibration/d10_ew_quotient_transport_receipt.json)
reuses the canonical target-free-input value-law evaluator, verifies its fibre
identities and nondegenerate local chart, and records every QT1--QT5 premise as
`assumed_not_enumerated`. In particular, missing QT3 leaves the displayed path
coefficients free, while missing QT5 leaves additional deformations unexcluded.
The implementation contract for freezing the actual D10/D11 self-reading
carriers, generating exact certificates, and keeping the separate clock, RG,
BRST-pole, uncertainty, and provenance gates fail-closed is
[`WZH_SOURCE_CLOSURE_CAMPAIGN.md`](WZH_SOURCE_CLOSURE_CAMPAIGN.md).

The local gravity bridge is carried on the stated local extension surface. On
the lifted product presentation of the realized quotient branch, the
`R`-sector collar entropy contribution satisfies
`(L_C)|_R = log d_R = log d_R3 + log d_R2`, because the `U(1)` factor is
one-dimensional. The resulting shared scalar is
`ellbar_shared = ellbar_SU(2)(t2_run) + ellbar_SU(3)(t3_run)`, and when the
same branch satisfies the D10 pixel law this becomes `ellbar_shared = P/4`.
That is the surface on which the local gravity readout uses the same emitted
scalar as the D10 pixel law.

The electroweak-hierarchy proof bundle lives in `../hierarchy`. It uses the
public endpoint pixel and the source-audit pixel recorded by the P-derivation
lane, then certifies the local `P -> alpha_U -> v/E_star` map with a declared
DAG check and an `R_U` Krawczyk inclusion witness. It should not be merged into
the rounded `1.63094` calibration carrier without a deliberate migration of all
dependent artifacts.

## Commands

Run the current implied-`P` audit:

```bash
python3 particles/calibration/implied_p_consistency_audit.py
```

This writes:

- [`particles/runs/calibration/implied_p_consistency.json`](../runs/calibration/implied_p_consistency.json)

Run the local calibration guard:

```bash
python3 particles/calibration/test_single_p_consistency.py
python3 particles/calibration/derive_d10_ew_observable_family.py
python3 particles/calibration/test_d10_observable_family_artifact.py
python3 particles/calibration/derive_d10_ew_transport_kernel.py
python3 particles/calibration/test_d10_ew_transport_kernel_artifact.py
python3 particles/calibration/derive_d10_ew_quotient_transport_receipt.py
python3 particles/calibration/derive_boson_source_prediction_closure_audit.py
python3 -m pytest -q particles/calibration/test_ew_physical_readout_gate.py
python3 -m pytest -q particles/calibration/test_wz_survival_boundaries.py
```

That guard checks the solver mechanics and the presence of the current D10
calibration observables. It is a mechanics guard, not a claim that exact
single-`P` closure has been achieved.

The D11 lane is split into one lower-rank companion branch and one conditional
pair map on the declared D10/D11 surface.

- The shared forward seed `sigma_D11_HT = alpha_u * cos(2*theta_W0) / sqrt(pi)`
  stays on disk as the diagonal fixed-ray companion branch on the declared
  D10/D11 surface. Its fixed-ray certificate proves `pi_y = pi_lambda`,
  `eta_HT = 0`, and `w_HT = 0` on that one-scalar branch only.
- The live pair artifact is `D11SourceSplitForwardExactness`, emitted by
  `derive_d11_live_exact_split_pair_theorem.py`. It uses only the forward D10
  tuple `(eta_source, beta_EW, lambda_EW, tau2_tree_exact, delta_n_tree_exact)`
  plus the declared D11 core/Jacobian surface at runtime. The absence of a
  direct target argument does not derive the D10 selector, D11 core/Jacobian,
  or residual-selector rigidity from OPH source records.
- The shared split scalar is `rho_HT = log(1 + tau2_tree_exact)`. The declared
  residual-selector formulas are
  `R_T = -tau2_tree_exact * eta_source^2 + (1 + beta_EW/28) * eta_source^6 + eta_source^8/14 + eta_source^9/27`
  and
  `R_H = eta_source^5 - (3/25) * eta_source^6 + lambda_EW * eta_source^6 / 18 + eta_source^8 / (2 * beta_EW)`.
- The conditional split coordinates are
  `pi_y = (eta_source + (3/2 + beta_EW/4) * rho_HT + R_T) / sqrt(pi)` and
  `pi_lambda = (eta_source - (4/3 - beta_EW/54) * rho_HT + R_H) / sqrt(pi)`.
  On the declared D11 surface these emit
  `m_t = 172.3523553288312 GeV` and `m_H = 125.1995304097179 GeV`.
  The Higgs row is a conditional declared-surface candidate, not a promoted
  strict source-only particle prediction. Promotion additionally requires the
  D10 QT1--QT5 source certificate on one strict pixel branch, a source-derived
  D11 surface and rigidity theorem, the absolute scale/RG receipts, and the
  Higgs complex-pole and uncertainty certificates.
  The exact public running-top row uses the PDG 2025 cross-section entry
  `Q007TP4`.
  The auxiliary direct-top average `Q007TP` is compare-only; #207 is closed as
  a current-corpus codomain no-go.
- The compare-only exact Higgs/top inverse slice remains a validation surface.
  The target-anchored Higgs-only and top-side exactifier artifacts remain on
  disk as supporting witness surfaces beneath the source-only pair theorem.

## License

This particle calibration surface is part of the OPH public repository. See the
main [LICENSE](../../../LICENSE).
