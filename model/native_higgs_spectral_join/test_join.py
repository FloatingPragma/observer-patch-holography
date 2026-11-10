"""Actual current/EOM ancestry and strict complete-record margins."""
import importlib.util
from pathlib import Path
from fractions import Fraction as F
import pytest
import sympy as s

spec = importlib.util.spec_from_file_location('native_higgs_spectral_join_public_model', Path(__file__).with_name('model.py'))
join = importlib.util.module_from_spec(spec); spec.loader.exec_module(join)


@pytest.fixture(scope='module')
def result():
    return join.certificate()


def test_full_current_parity_raw_L_and_Higgs_compiler(result):
    c = result['exact_full_source_operator_controls']
    assert c['actual_gauge_component_checks'] == 48
    assert c['actual_spin_component_checks'] == 12
    assert c['full_Yukawa_raw_lepton_commutator_zero']
    assert c['source_creation_probe_quantum_numbers']['raw_lepton_number_change'] == '-1'
    assert c['source_creation_probe_quantum_numbers']['Y'] == '1'
    assert len(c['all_local_occupation_sectors']) == 31


def test_same_actual_H_positive_cross_and_composite_weight(result):
    b = result['same_H_cross_residue_and_records']
    L, _ = map(F, b['excitation_energy_enclosure_for_that_exact_bin'])
    w, whi = map(F, b['elementary_unnormalized_bin_weight_interval'])
    gamma = F(b['real_elementary_composite_cross_residue_lower_exact'])
    assert gamma == L*w-F(720, 2**34) > F(28988, 100000)
    assert F(b['composite_unnormalized_bin_weight_lower_exact']) == gamma**2/whi > F(19997, 100000)


def test_preparation_inclusive_complete_interference_margin(result):
    b = result['same_H_cross_residue_and_records']
    assert F(b['unconditional_full_preparation_setting_hit_contrast_lower']) > F(7247, 100000)
    law = b['complete_detector_law']
    assert law['full_setting_spin_occupation_energy_record_count'] == 16
    assert law['potentially_nonzero_records_from_declared_empty_apparatus'] == 6
    assert law['preparation_partners_spin_references_and_all_filter_labels_retained']
    # Exact interference normalization: the B/2 coefficient and outer1/2
    # cannot be omitted when converting cross residue into complete hits.
    a, b2, cross = s.symbols('a b2 cross', real=True)
    plus = (a+b2/4+cross)/4; minus = (a+b2/4-cross)/4
    assert s.expand(plus-minus) == cross/2


def test_selector_scope_is_not_promoted_to_frozen_mass_identity(result):
    selector = result['same_H_cross_residue_and_records']['positive_cross_selector']
    assert selector['same_energy_elementary_composite_G1_doublet']
    assert not selector['is_the_maximum_elementary_residue_E_star_selector']
    assert not selector['maximum_elementary_E_star_common_composite_overlap_proved']
    assert not selector['physical_SM_weak_chirality_interactions_or_GeV_identity_proved']


def test_cross_inequality_detects_loss_of_positive_margin():
    assert join.cross_bound(F(1, 2), F(1, 4), F(1, 100), 1) > 0
    assert join.cross_bound(F(1, 2), F(1, 4), F(1, 2), 1) < 0
    with pytest.raises(ValueError): join.cross_bound(0, F(1, 4), 0, 1)


def test_wrong_neutral_composite_does_not_match_charged_current():
    neutral = join.poly.neutral_composite(0)
    assert not join.poly.commutator(join.poly.gauge_current(11), neutral)
    assert neutral != {}
    charged = join.poly.screened_eom('charged_lepton', 9, 0)
    assert join.poly.commutator(join.poly.gauge_current(11), charged) == charged


def test_frozen_research_science_is_replayed_exactly(result):
    comparison = result['independent_research_comparison']
    assert comparison['all_scientific_fields_equal'] and comparison['fresh_calculation_was_replayed']
    assert comparison['reference_sha256'] == '59902a655b4d314a32681445559858e076df1be29b522eb0398339b77cf4ff43'
