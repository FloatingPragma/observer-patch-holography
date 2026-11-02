from copy import deepcopy
import pytest

from . import packet


@pytest.fixture(scope='module')
def receipt():
    return packet._produce()


def test_cold_catalogue_replays_the_same_actual_reference(receipt):
    assert packet.report() == receipt


@pytest.mark.parametrize('key', ['physical_species_identified',
    'native_repair_fixed_algebra_extended', 'spontaneous_order_proved',
    'interacting_mass_shell_proved'])
def test_order_parameter_topology_cannot_be_promoted_to_particles(receipt, key):
    wrong = deepcopy(receipt)
    wrong['catalogue'][key] = True
    with pytest.raises(ValueError):
        packet.validate(wrong)


def test_all_valleys_and_actual_zero_order_are_retained(receipt):
    for section, key, changed in (
            ('clifford', 'weyl_valleys', 1),
            ('hedgehog', 'degree', 0)):
        wrong = deepcopy(receipt)
        wrong['catalogue'][section][key] = changed
        with pytest.raises(ValueError):
            packet.validate(wrong)
    wrong = deepcopy(receipt)
    wrong['actual_centered_ground']['onsite_mass_order_expectations'] = ['1']*4
    with pytest.raises(ValueError):
        packet.validate(wrong)


@pytest.mark.parametrize('path,changed', [
    (('clifford', 'bulk_gap'), '999'),
    (('exhaustion', 'constraint_rank'), 27),
    (('hedgehog', 'fermionic_zero_mode_proved'), True),
    (('antihedgehog', 'repair_protection_proved'), True),
    (('pure_spin_dirac_decomposition', 'physical_species_identified'), True),
    (('gapped_components', 'spin_texture', 'pi2'), '0'),
    (('source_anchor', 'all_four_anticommutators_zero'), 1),
])
def test_every_nested_theorem_and_scope_claim_is_replayed(receipt, path, changed):
    wrong = deepcopy(receipt)
    target = wrong['catalogue']
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = changed
    with pytest.raises(ValueError, match='classification'):
        packet.validate(wrong)


@pytest.mark.parametrize('corruption', ['empty', 'omitted', 'wrong_digest'])
def test_mass_manifest_must_be_complete_and_current(receipt, corruption):
    wrong = deepcopy(receipt)
    bindings = wrong['source_bindings']
    if corruption == 'empty':
        bindings.clear()
    elif corruption == 'omitted':
        bindings.pop(next(iter(bindings)))
    else:
        bindings[next(iter(bindings))] = '0'*64
    with pytest.raises(ValueError, match='manifest'):
        packet.validate(wrong)


@pytest.mark.parametrize('corruption', ['false_count', 'extra_claim', 'missing_field'])
def test_receipt_shape_and_zero_count_are_exact(receipt, corruption):
    wrong = deepcopy(receipt)
    if corruption == 'false_count':
        wrong['physical_comparison_count_added'] = False
    elif corruption == 'extra_claim':
        wrong['fermionic_zero_mode_proved'] = True
    else:
        wrong.pop('catalogue')
    with pytest.raises(ValueError):
        packet.validate(wrong)


def fake_parent():
    return {'schema': 'not-a-centered-ground-proof', 'all_pass': True,
            'source_bindings': packet.centered_packet.source_bindings(),
            'finite_gap': {'full_centered_gap_strict_lower': '1/128'},
            'charged_particles': {'onsite_mass_order': {
                'actual_unique_ground_order_expectations': ['0']*4}}}


def test_minimal_fabricated_parent_is_not_a_ground_proof(monkeypatch):
    fabricated = fake_parent()
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: fabricated)
    with pytest.raises(ValueError, match='invalid centered-source'):
        packet._produce()


@pytest.mark.parametrize('corruption', ['empty', 'omitted', 'wrong_digest'])
def test_parent_manifest_cannot_choose_which_sources_to_bind(monkeypatch, corruption):
    fabricated = fake_parent()
    bindings = fabricated['source_bindings']
    if corruption == 'empty':
        bindings.clear()
    elif corruption == 'omitted':
        bindings.pop(next(iter(bindings)))
    else:
        bindings[next(iter(bindings))] = '0'*64
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: fabricated)
    with pytest.raises(ValueError, match='centered reference source manifest'):
        packet._produce()
