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

Current live blocker:

- `eta_source_support_extension_log_per_side`

Smaller same-carrier primitive already on disk:

- `oph_charged_sector_local_support_extension_source_scalar_pair_readback`
  This collects the ordered `eta` then `sigma` invariant readbacks beneath the
  full support-extension shell without promoting either scalar value.

That is the first charged scalar that actually leaves the exhausted current
support and changes the hierarchy. Until it is emitted, the forward charged
artifact remains a structured gap surface rather than a promotable prediction.

Additional guards:

- [`test_shared_budget_not_silently_localized.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_shared_budget_not_silently_localized.py)
- [`test_channel_norm_refinement_limit.py`](/Users/muellerberndt/Projects/oph-meta/particles/leptons/test_channel_norm_refinement_limit.py)
