"""Adversarial joins against the literal common source/state/reference/clock."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

import pytest

from . import contract as c
from model.particle_field_state_response import state_response as warm
from model.particle_spatial_ground_state import ground as cold
from model.particle_gauss_ground_state import gauss
from model.particle_charged_residue import residue


@pytest.fixture(scope='module')
def primary():
    return c.member(0)


def forged_member_hashes(value):
    """An adversary can recompute every public checksum after substitution."""
    value = deepcopy(value)
    source_id = c.digest(value['source'])
    value['state']['source_id'] = source_id
    state_id = c.digest(value['state'])
    value['reference']['source_id'] = source_id
    value['reference']['state_id'] = state_id
    value['clock']['source_id'] = source_id
    value['ids'] = {'source_id': source_id, 'state_id': state_id,
                    'reference_id': c.digest(value['reference']),
                    'clock_id': c.digest(value['clock'])}
    value.pop('member_id', None)
    value['member_id'] = c.digest(value)
    return value


def test_member_is_an_actual_Gauss_plan_with_matching_action_characteristics(primary):
    c.verify_member(primary)
    plan = gauss.spatial_plan(side=primary['source']['side'],
        spacing=primary['source']['spacing'], link_degree=primary['source']['per_link_word_degree'],
        electric=primary['source']['electric_coefficient'], magnetic=primary['source']['magnetic_coefficient'],
        target_excess=primary['state']['target_ground_excess'])
    assert c.serial(plan) == primary['preparation_plan']
    assert plan['state_supported_in_source_Gauss_kernel']
    assert not plan['postselection']
    assert primary['matched_characteristics']['same_clock_gauge_Dirac_characteristic_match']


@pytest.mark.parametrize('substitute', ['warm', 'cold', 'native_Gauss', 'ground_trial'])
def test_native208_preparations_cannot_replace_the_spatial_state_even_with_forged_hashes(primary, substitute):
    value = deepcopy(primary)
    if substitute == 'warm':
        bound = warm.preparation_bound()
        value['state'].update(beta=c.serial(bound['beta']), source_only_Gauss_support=False,
                              amplitude='native208 warm positive Taylor polynomial')
    elif substitute == 'cold':
        bound = cold.polynomial_bound('4', -4, 22)
        value['state'].update(beta=c.serial(bound['beta']), Taylor_order=bound['Taylor_order'],
                              dyadic_scale=bound['dyadic_scale'], source_only_Gauss_support=False)
    elif substitute == 'native_Gauss':
        projector = gauss.native_projector()
        value['state'].update(native_physical_dimension=int(projector.trace()),
                              amplitude='P_N2 Q_native208', source_only_Gauss_support=True)
    else:
        trial = residue.ground_trial()
        value['state'].update(amplitude='native208 rational ground trial',
                              source_stationary=trial['source_density_exactly_stationary'],
                              ground_trial_error=c.serial(trial['half_trace_distance_upper']))
    with pytest.raises(ValueError, match='declared program'):
        c.verify_member(forged_member_hashes(value))


@pytest.mark.parametrize('substitute', ['clock_scale', 'reference_energy', 'reference_counterflow',
                                       'periodic_boundary', 'electric_coefficient', 'global_word_cutoff',
                                       'source_energy_shift'])
def test_clock_reference_and_source_changes_cannot_hide_behind_recomputed_ids(primary, substitute):
    value = deepcopy(primary)
    if substitute == 'clock_scale':
        value['clock']['Schrodinger'] = 'exp(-2i t H)'
    elif substitute == 'reference_energy':
        value['reference']['idle_H'] = 1
    elif substitute == 'reference_counterflow':
        value['reference']['idle_H'] = '-transpose(H_source)'
        value['reference']['joint_purification_stationary_under_source_only_clock'] = True
    elif substitute == 'periodic_boundary':
        value['source']['boundary'] = 'periodic_cube'
    elif substitute == 'electric_coefficient':
        alternative = cold.spatial_source_bounds(side=value['source']['side'], electric='1')
        value['source']['electric_coefficient'] = c.serial(alternative['electric_coefficient'])
    elif substitute == 'global_word_cutoff':
        value['source']['word_projection'] = 'global_total_word_degree'
    else:
        value['source']['energy_zero'] = 'subtract lowest source eigenvalue'
    with pytest.raises(ValueError, match='declared program'):
        c.verify_member(forged_member_hashes(value))


def test_a_receipt_from_another_constructed_family_member_is_not_admitted(primary):
    other = c.member(1)
    receipt = c.bind_result(other, 'specified source-only Gauss polynomial preparation',
                            other['preparation_plan'])
    with pytest.raises(ValueError, match='incompatible'):
        c.verify_result(primary, receipt)


def test_correctly_bound_actual_producer_result_is_admitted(primary):
    receipt = c.bind_result(primary, 'specified source-only Gauss polynomial preparation',
                            primary['preparation_plan'])
    c.verify_result(primary, receipt)


def test_direct_exact_producer_payload_roundtrips_without_fraction_string_ambiguity(primary):
    plan = gauss.spatial_plan(side=primary['source']['side'],
        spacing=primary['source']['spacing'], link_degree=primary['source']['per_link_word_degree'],
        electric=primary['source']['electric_coefficient'], magnetic=primary['source']['magnetic_coefficient'],
        target_excess=primary['state']['target_ground_excess'])
    receipt = c.bind_result(primary, 'specified source-only Gauss polynomial preparation', plan)
    c.verify_result(primary, receipt)


def test_implementation_bytes_change_source_identity_and_reject_stale_receipts(primary, monkeypatch):
    receipt = c.bind_result(primary, 'specified source-only Gauss polynomial preparation',
                            primary['preparation_plan'])
    implementation = c.ROOT/'model/particle_common_contract/contract.py'
    read_bytes = Path.read_bytes
    def substituted_bytes(path):
        value = read_bytes(path)
        return value+b'\n# isolated audit substitution\n' if path == implementation else value
    # No real file or inherited pin changes. This substitutes what the
    # source-definition hasher sees while retaining the same scalar bounds.
    monkeypatch.setattr(Path, 'read_bytes', substituted_bytes)
    changed = c.member(0)
    assert changed['preparation_plan'] == primary['preparation_plan']
    assert changed['ids']['source_id'] != primary['ids']['source_id']
    with pytest.raises(ValueError, match='incompatible'):
        c.verify_result(changed, receipt)


def test_recomputed_payload_hash_does_not_admit_a_different_preparation(primary):
    receipt = c.bind_result(primary, 'specified source-only Gauss polynomial preparation',
                            primary['preparation_plan'])
    wrong = cold.spatial_plan(side=primary['source']['side'],
        spacing=primary['source']['spacing'], link_degree=primary['source']['per_link_word_degree'],
        electric=primary['source']['electric_coefficient'], magnetic=primary['source']['magnetic_coefficient'],
        target_excess=primary['state']['target_ground_excess'])
    # Same configured action and thermal budget, but the actual producer
    # supplies the ambient full carrier rather than the physical Gauss state.
    receipt['result'] = c.serial(wrong)
    receipt['result_sha256'] = c.digest(receipt['result'])
    with pytest.raises(ValueError):
        c.verify_result(primary, receipt)


def test_recomputed_payload_hash_does_not_admit_unmatched_characteristics(primary):
    receipt = c.bind_result(primary, 'gauge and matter principal characteristic match',
                            primary['matched_characteristics'])
    alternative = c.characteristics.clock_family('1/7', '1')
    assert not alternative['same_clock_gauge_Dirac_characteristic_match']
    receipt['result'] = c.serial(alternative)
    receipt['result_sha256'] = c.digest(receipt['result'])
    with pytest.raises(ValueError):
        c.verify_result(primary, receipt)


def test_stationary_preparation_does_not_supply_an_exact_Gibbs_modular_clock(primary):
    assert primary['state']['source_stationary']
    assert not primary['state']['exact_Gibbs_state']
    assert not primary['clock']['modular_equals_minus_beta_times_source_time']
    forged = deepcopy(primary)
    forged['state']['exact_Gibbs_state'] = True
    forged['clock']['modular_equals_minus_beta_times_source_time'] = True
    with pytest.raises(ValueError):
        c.verify_member(forged_member_hashes(forged))
    # Provenance-only hashes must not turn this extra, unregistered state
    # theorem into an admitted consequence of the fixed program.
    with pytest.raises(ValueError):
        c.bind_result(primary, 'exact Gibbs modular clock from source stationarity',
                      {'source_stationary': True, 'modular_time': '-beta*s'})


def test_unknown_or_omitted_aggregate_certificate_domain_is_rejected():
    report = {name: {} for name in c.CONTROL_DOMAINS}
    c.validate_inventory(report)
    unknown = deepcopy(report)
    unknown['unclassified_relativistic_vacuum'] = {'proof': 'a matching hash'}
    with pytest.raises(ValueError, match='unclassified or omitted'):
        c.validate_inventory(unknown)
    missing = deepcopy(report)
    missing.pop('interacting_field_thermal_records')
    with pytest.raises(ValueError, match='unclassified or omitted'):
        c.validate_inventory(missing)


@pytest.mark.parametrize('field,value', [('physical_particle_assignments', ['photon']),
                                       ('source_predicted_GeV_masses', [{'electron': '.000511'}])])
def test_comparison_domains_cannot_promote_physical_particle_or_mass_claims(field, value):
    report = {name: {} for name in c.CONTROL_DOMAINS}
    report[field] = value
    with pytest.raises(ValueError, match='separately proved dictionary'):
        c.validate_inventory(report)
