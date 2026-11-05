from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..centered_particle_spectrum import packet as ground_packet
from . import certificate, packet, source


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_complete_spatial_coverage_and_exact_source_band_projectors(result):
    s = result['spectrum']
    assert len(s['rows']) == s['positive_lines'] == 138
    assert s['visible_spatial_dimension_fixed_highest_spin'] == 552
    assert s['visible_N7_dimension_with_spin_ladder'] == 2208
    assert s['visible_N9_dimension_with_spin_and_eta_ladders'] == 2208
    assert s['orthogonal_union_N7_dimension'] == 4032
    assert s['unseen_N7_dimension_after_union'] == s['unseen_N9_dimension_after_ladder_union'] == 7408
    assert sum(s['actual_source_kinetic_band_counts'].values()) == 138
    assert min(s['actual_source_kinetic_band_counts'].values()) > 0
    adjacency = source.parent.ground.band_data()[0]
    for row in s['rows']:
        assert row['visible_spatial_rank'] == 4
        assert row['visible_states_per_charge_with_spin_ladder'] == 16
        assert len(row['residue_matrix_intervals']) == 8
        for x in range(8):
            for y in range(8):
                lo, hi = map(F, row['residue_matrix_intervals'][x][y])
                assert lo <= hi
                if x != y and adjacency[x, y] == 0:
                    assert lo == hi == 0
    assert result['source_spatial_commutant']['exact_commutator_rank'] == 62


def test_complete_covariance_and_exact_spin_eta_maps(result):
    assert result['source_covariance']['complete_operator_covariance_checks'] == 384
    ladder = result['spin_eta_ladders']
    assert [r['step_norm_squared'] for r in ladder['spin_lowering']] == [3, 4, 3]
    assert ladder['eta_raise_N7_to_N9_norm_squared'] == 1
    assert ladder['eta_raised_probe_is_bare_C_adjoint'] is False


def test_all_spatial_sum_rules_match_independent_source_and_retain_contact(result):
    for row in result['spectrum']['source_current_moment_matrices']:
        assert row['all_64_entries_checked'] is True
        for x in range(8):
            for y in range(8):
                a = tuple(map(F, row['spectral_matrix_intervals'][x][y]))
                b = tuple(map(F, row['independent_source_matrix_intervals'][x][y]))
                assert max(a[0], b[0]) <= min(a[1], b[1])
    d = result['spectrum']['origin_kinetic_contact_decomposition']
    assert F(d['centered_contact_first_moment'][0]) > F(35, 100)
    assert F(d['centered_contact_rate_norm_squared'][0]) > F(18, 100)
    assert F(d['kinetic_contact_interference'][0]) > F(2, 100)
    assert result['source_rate_algebra']['centered_contact_rate_norm_squared_is_probe_norm_squared_over_four'] is True


def test_spatial_propagation_matches_original_full_occupation_evolution(result):
    old = ground_packet.load_witness()
    g = source.parent.reduced('ground', tuple(map(int, old['ground_columns'])))
    vector = g['full_basis'] @ (old['ground_vectors'][:, 0].astype(float)/2.0**60)
    vector /= np.linalg.norm(vector)
    probes = np.column_stack([source.probe(x) @ vector for x in range(8)])
    e0 = float(old['ground_energies'][0])/2.0**56
    h = source.parent.centered(5, 2).astype(float)/16-e0*sparse.eye(1568)
    time = .5
    evolved = expm_multiply((-1j*time)*h, probes)
    direct = probes.conj().T @ evolved
    spectral = np.zeros((8, 8), complex)
    for row in result['spectrum']['rows']:
        energy = sum(map(float, row['energy_interval']))/2
        residue = np.array([[sum(map(float, entry))/2 for entry in line]
                            for line in row['residue_matrix_intervals']])
        spectral += np.exp(-1j*time*energy)*residue
    np.testing.assert_allclose(direct, spectral, rtol=0, atol=2e-10)
    assert abs(direct[0, 1]) > 1e-2
    assert abs(direct[0, 0].imag) > .1


def test_kinetic_source_change_is_rejected_even_if_symmetry_is_preserved(monkeypatch):
    previous = source.parent.ground.band_data()
    monkeypatch.setattr(source.parent.ground, 'band_data', lambda: (2*previous[0], *previous[1:]))
    with pytest.raises(ArithmeticError, match='actual complete sector source'):
        source.kinetic_source_anchor()


def test_missing_spatial_generators_does_not_fake_complete_commutant(monkeypatch):
    all_rows = source.parent.automorphisms()
    monkeypatch.setattr(source.parent, 'automorphisms', lambda: all_rows[:1])
    with pytest.raises(ArithmeticError, match='completely exhausted'):
        source.spatial_commutant_certificate()


def test_omitting_contact_changes_full_source_rate():
    maps = source.direct_rate_maps()
    assert (maps['full_rate16']-maps['kinetic_rate16']).nnz > 0
    assert not (maps['full_rate16']-maps['kinetic_rate16']-maps['centered_contact_rate16']).nnz


@pytest.mark.parametrize('bad', [True, -1, 8, 0.0])
def test_probe_domain_is_checked_before_cache_lookup(bad):
    source.probe(0)
    with pytest.raises(ValueError):
        source.probe(bad)


@pytest.mark.parametrize('change', ['rank', 'band', 'entry', 'moment', 'species', 'momentum', 'detector'])
def test_scientific_and_scope_mutations_fail_closed(result, change):
    damaged = deepcopy(result)
    s = damaged['spectrum']
    if change == 'rank':
        s['rows'][0]['visible_spatial_rank'] = 8
    elif change == 'band':
        s['rows'][0]['actual_source_kinetic_band_sign'] *= -1
    elif change == 'entry':
        s['rows'][0]['residue_matrix_intervals'][0][1] = ['0', '0']
    elif change == 'moment':
        s['source_current_moment_matrices'][2]['all_64_entries_checked'] = False
    else:
        s[{'species': 'physical_particle_species_identified',
           'momentum': 'spatial_labels_are_continuum_momentum',
           'detector': 'native_detector_execution_claimed'}[change]] = True
    with pytest.raises(ValueError, match='exact reconstruction'):
        packet.validate(damaged, result)


def test_retained_receipt_and_complete_source_manifest(result):
    retained = packet.strict_json.load(packet.OUTPUT)
    packet.validate(retained, result)
    for path in ('model/composite_spatial_response/PROOF.md',
                 'model/charged_composite_search/witness.npz',
                 'model/centered_particle_spectrum/witness.npz'):
        assert path in retained['source_bindings']


def test_string_path_cold_replay():
    assert packet.verify(str(packet.OUTPUT)) is True


def test_duplicate_json_rejected_before_proof_replay(tmp_path):
    wrong = tmp_path/'bad.json'
    wrong.write_text('{"all_pass":true,"all_pass":false}')
    with pytest.raises(ValueError):
        packet.verify(wrong)
