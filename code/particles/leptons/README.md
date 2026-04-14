# Charged-Lepton Lane

This directory is the active `/particles` charged-lepton completion lane.

The live chain is:

1. read back the current ordered charged package on the local family
2. certify that the current support is exhausted
3. expose the minimal beyond-support extension scalar
4. turn the ordered coefficients into the charged excitation-gap map
5. build the log-spectrum readout
6. attach the scale/norm lane and emit the forward charged candidate

## How To Read The Active Charged-Lepton Files

The live charged scripts now use the same compact derivation header:

- `Chain role`: where the file sits between the shared flavor carrier and the
  forward charged candidate
- `Mathematics`: which ordered-gap, support-rank, or spectral step it performs
- `OPH-derived inputs`: which values are inherited directly from the active
  `/particles` flavor/lepton artifacts
- `Output`: the emitted artifact and, when still open, the exact missing scalar

The active charged completion tail is:

- `derive_charged_sector_local_current_support_obstruction_certificate.py`
- `derive_charged_sector_local_minimal_source_support_extension_emitter.py`
- `derive_lepton_excitation_gap_map.py`
- `derive_lepton_log_spectrum_readout.py`
- `build_forward_charged_leptons.py`

Current scripts:

- [`derive_charged_sector_local_current_support_obstruction_certificate.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/derive_charged_sector_local_current_support_obstruction_certificate.py)
- [`derive_charged_sector_local_minimal_source_support_extension_emitter.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/derive_charged_sector_local_minimal_source_support_extension_emitter.py)
- [`derive_lepton_excitation_gap_map.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/derive_lepton_excitation_gap_map.py)
- [`derive_lepton_log_spectrum_readout.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/derive_lepton_log_spectrum_readout.py)
- [`build_forward_charged_leptons.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/build_forward_charged_leptons.py)
- [`test_no_koide_import.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_no_koide_import.py)
- [`test_no_experiment_label_matching.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_no_experiment_label_matching.py)
- [`test_channel_norm_not_fit.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_channel_norm_not_fit.py)
- [`test_ratio_only_not_promoted.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_ratio_only_not_promoted.py)

These scripts do **not** claim charged leptons are already theorem-level. They
give the `e30` closure lane a concrete local home under `/particles`, and they
keep the ordered-shape/hierarchy problem separate from absolute-scale closure.

Current same-carrier constructive blocker:

- `eta_source_support_extension_log_per_side`

Current theorem-facing absolute route:

- `oph_generation_bundle_branch_generator_splitting`
- `refinement_stable_uncentered_trace_lift`
- then the determinant-line section and `A_ch` are induced rather than independent

Operator-side exact obstruction certificate:

- `oph_generation_bundle_branch_generator_splitting_obstruction_certificate`
  This fixes the current operator-side claim boundary sharply: the displayed
  proxy ordered spectrum together with the displayed level-0 and level-1
  projector system does not force the descended-commutator clause. The
  obstruction witness is the nonzero natural centered commutator with operator
  norm `0.04861550547372144` against projector defect
  `0.06363734112184061`, so the displayed current data do not imply exact
  vanishing or uniform quadratic smallness after the central split.

The post-promotion lift slot has now been reduced further in carrier type:

- `oph_charged_uncentered_trace_lift_cocycle_reduction`
  This keeps the single-slot frontier unchanged while making its exact content
  explicit: after centered promotion, the remaining ambiguity is only a scalar
  affine cocycle primitive `mu`, with `C_tilde_e = C_hat_e + mu I`,
  `s_det = 3 mu`, and `A_ch = mu`.

- `oph_charged_mu_physical_descent_reduction`
  This shrinks that description one step further under the same refinement-
  stability contract already required by the lift: on theorem-grade physical
  `Y_e`, the refinement cocycle vanishes and the exact residual object is one
  physical affine scalar `mu_phys(Y_e)`.

- `oph_charged_centered_operator_mu_phys_no_go`
  This closes the false post-promotion shortcut: even a theorem-grade centered
  `C_hat_e` cannot emit `mu_phys(Y_e)` by itself, because centered operator
  data is still invariant under the common-shift action that `mu_phys` breaks.

- `oph_charged_p_to_affine_anchor_reduction`
  This fixes the smallest honest bridge target for the `P`-driven absolute
  lane. A theorem-grade landing from the D10 descendants of `P` on theorem-
  grade physical `Y_e(P)` or on the charged determinant line already forces
  `A_ch(P) = (1/3) * log(det(Y_e(P)))`, so the bridge target is smaller than a
  free-standing affine-scalar theorem.

  The same descent package also gives a sharper P-threaded reformulation:
  once `P -> Y_e(P)` lands on theorem-grade physical charged data and the
  post-promotion charged descent hypotheses hold, the physical affine scalar
  `mu_phys(Y_e(P))` is forced by the identity-mode equalizer rather than
  remaining a separate bridge-sized theorem.

- `oph_charged_physical_class_affine_scalar_reduction`
  This is the strongest honest selected-surface charged theorem packaged on the
  current local surface. It does not claim a quark-style public selected-class
  charged source descent. It states that on theorem-grade physical charged
  `Y_e(P)`, once theorem-grade `C_hat_e` exists, the determinant line fixes the
  canonical physical affine scalar `mu_phys(Y_e(P)) = (1/3) * log(det(Y_e(P)))`
  and therefore the canonical uncentered lift, determinant-line, affine-anchor,
  and charged mass formulas.

- `oph_charged_physical_determinant_line_canonical_uncentered_lift`
  This is the exact conditional theorem behind `#151`: on theorem-grade
  physical charged `Y_e`, theorem-grade `C_hat_e` and the determinant line
  canonically fix
  `C_tilde_e(Y_e) = C_hat_e(Y_e) + mu_phys(Y_e) I`
  with
  `mu_phys(Y_e) = (1/3) * log(det(Y_e))`.
  The identity-mode equalizer on a fixed physical fiber is then tautological.

- `oph_charged_mass_readout_from_affine_anchor`
  This is the exact downstream theorem beneath the source-landing frontier:
  once theorem-grade `A_ch(P)` exists on theorem-grade physical charged data,
  the charged absolute scale and charged mass triple are pure algebraic
  readouts,
  `g_e(P) = exp(A_ch(P))` and
  `m_i(P) = exp(A_ch(P) + ell_i_centered(P))`.
  The remaining open charged source problem is not this algebraic readout. The
  remaining open theorem is the landing from the D10 calibration descendants of
  `P` to theorem-grade physical charged data or the charged determinant line.

Layered frontier ledger now on disk:

- `oph_charged_absolute_frontier_factorization`
  This separates the current-surface missing affine object `A_ch` from the
  strictly smaller post-promotion single slot `refinement_stable_uncentered_trace_lift`.

Smaller same-carrier primitive already on disk:

- `oph_charged_sector_local_support_extension_source_scalar_pair_readback`
  This collects the ordered `eta` then `sigma` invariant readbacks beneath the
  full support-extension shell without promoting either scalar value.

That is the first charged scalar that actually leaves the exhausted current
support and changes the hierarchy. Until it is emitted, the forward charged
artifact remains a structured gap surface rather than a promotable prediction.

The theorem-facing absolute tail is sharper than that same-carrier frontier.
Once the latent centered charged operator is promoted, the honest post-promotion
single slot is the refinement-stable uncentered trace lift. The determinant-line
section and affine absolute anchor `A_ch` then follow as induced data.

Additional guards:

- [`test_shared_budget_not_silently_localized.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_shared_budget_not_silently_localized.py)
- [`test_channel_norm_refinement_limit.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_channel_norm_refinement_limit.py)
