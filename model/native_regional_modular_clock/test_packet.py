"""Owner replay and adversarial controls for the regional modular receipt."""
from dataclasses import replace
from fractions import Fraction as F
import json

import pytest

from . import certificate, collar, packet


@pytest.fixture(scope='module')
def fresh():
    return packet.produce()


def test_original_source_and_complete_seam_restrictions(fresh):
    native = fresh['same_reference_regional_certificate']
    assert native['native_source']['all_transfer_entries_authenticated'] == 4356
    assert native['native_source']['P_trace'] == '56'
    assert len(native['regional_seams']) == 30
    for seam in native['regional_seams']:
        assert set(seam['regional_density_blocks']) == {'0', '1', '2'}
        assert seam['number_commutator_HS_squared'] == '1/70560'
    assert native['conditional_sampler']['accepted_probability'] == '28/33'
    assert native['conditional_sampler']['rejected_probability'] == '5/33'


def test_actual_port_collar_and_number_preserving_scope(fresh):
    row = fresh['authenticated_finite_source_collars'][1]
    assert (row['source_level'], row['authenticated_cell_radius']) == (14, 4)
    assert (row['ports'], row['internal_seams'], row['crossing_seams'],
            row['actual_port_collar_distance']) == (348, 905, 17, 9)
    assert F(row['bounds']['regional_modular_to_larger_exchange_flow_upper']) < F(1068, 10**9)
    assert fresh['limiting_contract']['observable_algebra'] == (
        'regional observables preserving local occupation number')
    assert fresh['scope']['fixed_region_modular_overlap_compatibility'] is False
    assert fresh['scope']['physical_G1_to_G6_discharged'] is False
    assert fresh['custody']['all_inherited_producers_replayed'] is False


@pytest.mark.parametrize('mutation', (
    'omit_local_sector', 'falsify_coherence', 'reverse_clock',
    'remove_density_rate', 'claim_fixed_overlap', 'claim_preparation',
    'claim_all_level_routing', 'change_bound', 'change_source', 'boolean_to_integer'))
def test_scientific_and_typed_mutations_fail(fresh, mutation):
    retained = json.loads(packet.canonical(fresh))
    native = retained['same_reference_regional_certificate']
    if mutation == 'omit_local_sector':
        del native['regional_seams'][0]['regional_density_blocks']['0']
    elif mutation == 'falsify_coherence':
        native['regional_seams'][0]['number_commutator_HS_squared'] = '0'
    elif mutation == 'reverse_clock':
        retained['limiting_contract']['clock'] = 'Heisenberg exchange time +t'
    elif mutation == 'remove_density_rate':
        retained['limiting_contract']['density_condition'] = 'q/V -> p is sufficient'
    elif mutation == 'claim_fixed_overlap':
        retained['scope']['fixed_region_modular_overlap_compatibility'] = True
    elif mutation == 'claim_preparation':
        retained['scope']['global_transfer_reference_scalable_preparation_supplied'] = True
    elif mutation == 'claim_all_level_routing':
        retained['scope']['all_level_numerical_port_routing_proved'] = True
    elif mutation == 'change_bound':
        retained['authenticated_finite_source_collars'][1]['bounds'][
            'regional_modular_to_larger_exchange_flow_upper'] = '0'
    elif mutation == 'change_source':
        first = next(iter(retained['source_bindings']))
        retained['source_bindings'][first] = '0' * 64
    else:
        retained['scope']['physical_G1_to_G6_discharged'] = 0
    with pytest.raises(ArithmeticError, match='fresh owning replay'):
        packet.check_payload(retained, fresh)


def test_binding_inventory_excludes_self_receipt_and_research_notes(fresh):
    paths = set(fresh['source_bindings'])
    assert set(packet.OWNED_FILES) <= paths
    assert set(certificate.transfer.SOURCE_PATHS) <= paths
    assert not any(path.startswith('docs/') for path in paths)
    assert 'model/native_regional_modular_clock/results.json' not in paths
    assert packet.ROOT == certificate.transfer.ROOT


def test_mid_replay_source_change_is_rejected(fresh, monkeypatch):
    before = fresh['source_bindings']
    after = dict(before)
    after[next(iter(after))] = 'f' * 64
    calls = iter((before, after))
    monkeypatch.setattr(packet, 'source_bindings', lambda: next(calls))
    monkeypatch.setattr(certificate, 'certificate', lambda: fresh['same_reference_regional_certificate'])
    rows = iter(fresh['authenticated_finite_source_collars'])
    monkeypatch.setattr(packet.collar_certificate, 'row', lambda *_: next(rows))
    with pytest.raises(ArithmeticError, match='changed during'):
        packet.produce()


def test_original_kernel_mutation_fails_before_regional_claim(monkeypatch):
    native = certificate.transfer.native_counts()
    rows = [list(row) for row in native.counts]
    rows[0][0] += 1
    malformed = replace(native, counts=tuple(map(tuple, rows)))
    monkeypatch.setattr(certificate.transfer, 'native_counts', lambda: malformed)
    with pytest.raises(ArithmeticError, match='symmetric and stochastic'):
        certificate.certificate()


def test_density_convergence_without_rate_does_not_close_growing_bound():
    capacity_terms, errors = [], []
    for n in (10, 20, 40):
        V, v = n**6, n**2
        q, M = V // 2 + V // n, 3 * V
        capacity_terms.append(F(v*v, M))
        assert F(q, V) - F(1, 2) == F(1, n)
        errors.append(collar.uniform_modular_error(V, M, q, v, v, 1, F(1, 2)))
    assert all(a > b for a, b in zip(capacity_terms, capacity_terms[1:]))
    assert all(a < b for a, b in zip(errors, errors[1:]))
    assert errors[-1] > 100


def test_untyped_or_nonfinite_receipt_values_fail():
    for value in (float('nan'), float('inf'), 1.0, {1: 'wrong key'}):
        with pytest.raises(TypeError, match='exact typed'):
            packet.canonical(value)


def test_retained_receipt_matches_fresh_owner(fresh):
    assert packet.check_payload(packet.strict_json.load(packet.OUTPUT), fresh)


def test_duplicate_retained_keys_fail_closed(fresh, monkeypatch, tmp_path):
    path = tmp_path / 'ambiguous.json'
    path.write_text('{"schema":"first","schema":"second"}')
    monkeypatch.setattr(packet, 'produce', lambda: fresh)
    with pytest.raises(packet.strict_json.DuplicateKeyError):
        packet.verify(path)
