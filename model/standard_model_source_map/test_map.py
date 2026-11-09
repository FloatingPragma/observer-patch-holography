from copy import deepcopy
from functools import lru_cache
import subprocess
import sys

import pytest

from . import algebra, packet


def test_source_weights_have_independent_occupation_count_formula():
    rows = algebra.source_components()
    assert len(rows) == 15
    for row in rows:
        mask = row['occupation_mask']
        colors = (mask & 7).bit_count()
        weak = (mask >> 3).bit_count()
        assert row['q6Y'] == str(-2*colors+3*weak)
        assert colors+weak in (2, 4)
    match = algebra.match_components()
    assert sorted(m for r in match.values() for m in r['source_masks']) == sorted(r['occupation_mask'] for r in rows)


@pytest.mark.parametrize('kind', ['missing', 'duplicate', 'charge', 'color_conjugation', 'multiplicity'])
def test_incomplete_or_incorrect_weight_matches_fail(kind):
    rows = deepcopy(list(algebra.source_components()))
    if kind == 'missing': rows.pop()
    elif kind == 'duplicate': rows[-1] = deepcopy(rows[0])
    elif kind == 'charge': rows[0]['Q'] = '0'
    elif kind == 'color_conjugation': rows[0]['color_weight'] = ['1', '0']
    else: rows[0]['weak_C2'] = '3/4'
    with pytest.raises(ValueError): algebra.match_components(rows)


def test_same_current_bosonic_directions():
    value = algebra.gauge_and_higgs()
    assert value['weak_raising_and_lowering_ranks_on_matter'] == [4, 4]
    assert value['Higgs_component_Q'] == ['1', '0']
    assert value['independent_color_current_directions'] == 8


@lru_cache(maxsize=1)
def built():
    return packet.build()


def test_export_separates_family_and_physical_identification():
    value = built()
    rows = value['rows']
    assert len(rows) == len({r['target'] for r in rows}) == 18
    assert value['physical_particle_map'] == []
    assert all(r['predicted_mass'] is None and r['identified_topological_defect'] is None for r in rows)
    assert all(r['source_family_index'] is None and r['physical_particle_identified'] is False for r in rows)
    leptons = [next(r for r in rows if r['target'] == n) for n in ('electron', 'muon', 'tau')]
    assert all(r['internal_components'] == leptons[0]['internal_components'] for r in leptons)
    assert len([r for r in rows if r['target_kind'] == 'neutrino_flavor_channel']) == 3
    assert next(r for r in rows if r['target'] == 'gluon')['color_components'] == 8


@pytest.mark.parametrize('kind', ['mass', 'physical', 'defect', 'family', 'source_pin',
                                 'mirror', 'Koide', 'coherent_phase', 'clock'])
def test_unsupported_export_promotions_rejected(kind):
    value = deepcopy(built())
    row = value['rows'][0]
    if kind == 'mass': row['predicted_mass'] = 0
    elif kind == 'physical': row['physical_particle_identified'] = True
    elif kind == 'defect': row['identified_topological_defect'] = 'hedgehog'
    elif kind == 'family': row['source_family_index'] = 1
    elif kind == 'source_pin': value['source_bindings'].pop(next(iter(value['source_bindings'])))
    elif kind == 'mirror':
        value['finite_source_attachments']['kinetic_and_mirror_attachment']['local_sine_regulator']['mirror_removal_or_weight_suppression_proved'] = True
    elif kind == 'Koide':
        value['finite_source_attachments']['source_action_masses']['policies'][0]['channels'][2]['Koide_diagnostic']['balanced_value_two_thirds_excluded'] = False
    elif kind == 'coherent_phase':
        value['finite_source_attachments']['coherent_Koide_source']['scope']['phase_or_mass_ratio_selected'] = True
    else:
        value['finite_source_attachments']['complete_feedback_limit']['physical_clock_or_vacuum_selected'] = True
    with pytest.raises(ValueError): packet.validate(value)


def test_generated_artifacts_match_computation():
    assert packet.OUTPUT.read_text() == packet.canonical(built())
    assert packet.MARKDOWN.read_text() == packet.render(built())
    assert packet.MASS_MARKDOWN.read_text() == packet.render_masses(built())


def test_export_preserves_all_source_duties_and_distinct_koide_program():
    value = built()
    assert value['configured_family_copies'] == 3
    source = value['finite_source_attachments']
    assert [p['policy'] for p in source['source_action_masses']['policies']] == ['always', 'C', 'E', 'CE']
    assert source['coherent_Koide_source']['scope']['different_configured_source_from_current_gated_Z_D'] is True
    assert source['source_action_masses']['physical_pole_mass_table'] == []
    assert source['fermion_operator_attachment']['physical_CAR_modes'] == 90


def test_scientific_cache_does_not_share_mutable_export_state():
    value = packet.build()
    value['finite_source_attachments']['source_action_masses']['policies'].clear()
    assert len(packet.build()['finite_source_attachments']['source_action_masses']['policies']) == 4


def test_fresh_complete_export_binds_its_runtime_source_imports():
    script = '''
from pathlib import Path
import sys
from model.standard_model_source_map import packet
value = packet.build()
missing = []
for module in tuple(sys.modules.values()):
    filename = getattr(module, '__file__', None)
    if not filename:
        continue
    path = Path(filename).resolve()
    if path.is_relative_to(packet.ROOT) and path.suffix == '.py':
        name = path.relative_to(packet.ROOT).as_posix()
        if name not in value['source_bindings']:
            missing.append(name)
assert not missing, sorted(set(missing))
'''
    subprocess.run([sys.executable, '-c', script], cwd=packet.ROOT, check=True,
                   capture_output=True, text=True)
