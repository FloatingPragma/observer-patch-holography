from copy import deepcopy

import pytest

from . import packet


def test_packet_replays_exact_source_and_retains_complex_frequency_domain():
    value=packet._produce()
    packet.validate(value)
    assert value['source_algebra']['selected_ground_symmetry']['exact_translation_cases']==4096
    assert value['spatial_bounds']['eta_zero_uniformity_claimed'] is False
    assert value['connected_field']['non_Gaussian_multitime_ground_field_proved'] is True
    assert value['Euclidean_hierarchy']['positive_time_reflection_is_word_Gram'] is True
    for section,key in [
        ('spatial_bounds','eta_zero_uniformity_claimed'),
        ('spatial_bounds','common_Lorentz_metric_or_Einstein_premise_completed'),
        ('spatial_bounds','limiting_state_or_symbol_values_computed'),
        ('connected_field','continuum_limit_or_scattering_proved'),
        ('connected_field','nonlinear_observable_encodings_excluded'),
        ('connected_field','local_state_faithfulness_assumed'),
        ('Euclidean_hierarchy','Euclidean_O4_covariance_inferred'),
        ('Euclidean_hierarchy','continuum_OS_axioms_completed'),
        ('Euclidean_hierarchy','limiting_hierarchy_values_computed'),
    ]:
        changed=deepcopy(value);changed[section][key]=True
        with pytest.raises(ValueError):packet.validate(changed)


def test_ordered_reflection_requires_the_odd_word_sign():
    value=packet._produce()
    changed=deepcopy(value)
    changed['Euclidean_hierarchy']['odd_word_cross_entry']='1/64'
    with pytest.raises(ValueError):packet.validate(changed)


def test_connected_time_bound_and_anomalous_pairings_cannot_be_promoted():
    value=packet._produce()
    for section,key,replacement in [
        ('uniform_bounds','time_radius','1'),
        ('uniform_bounds','variance_t_squared_lower','1/4'),
        ('local_identity','anomalous_pairings_retained',False),
    ]:
        changed=deepcopy(value)
        changed['connected_field'][section][key]=replacement
        with pytest.raises(ValueError):packet.validate(changed)


def test_fresh_worker_and_in_process_scientific_payload_agree():
    assert packet.report()==packet._produce()
