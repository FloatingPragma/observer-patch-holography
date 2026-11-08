"""Prevent kinetic labels or finite source energies becoming physical masses."""
import copy
from fractions import Fraction as F
import subprocess
import sys

import pytest
import sympy as s

from tools import strict_json
from . import algebra, bounds, charge_bridge as charge, packet


@pytest.fixture(scope='module')
def report():
    value = packet.build()
    assert packet.canonical(value) == packet.canonical(strict_json.load(packet.OUTPUT))
    assert packet.MARKDOWN.read_text() == packet.render(value)
    return value


def test_actual_source_census_and_empty_physical_table(report):
    assert len(report['free_branch_census']) == 8
    assert {tuple(r['kinetic_valley']) for r in report['free_branch_census']} == {
        (a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)}
    assert report['physical_particle_mass_table'] == []
    assert report['mass_ratios'] == []
    assert report['ordered_spin_branch']['invariant_rest_energy_squared'] is None
    assert not any(report['physical_admission'].values())


@pytest.mark.parametrize('mutation', [
    'physical_row', 'species', 'mass', 'missing_valley', 'missing_hole',
    'clock', 'spin_mass', 'ratio', 'native', 'topology', 'stationarity', 'type',
])
def test_promotion_or_dictionary_mutation_rejected(report, monkeypatch, mutation):
    value = copy.deepcopy(report)
    if mutation == 'physical_row': value['physical_particle_mass_table'] = [{'particle': 'electron', 'mass': '1'}]
    elif mutation == 'species': value['free_branch_census'][0]['empirical_species'] = 'neutrino'
    elif mutation == 'mass': value['free_branch_census'][0]['invariant_rest_energy_squared'] = '1.092'
    elif mutation == 'missing_valley': value['free_branch_census'].pop()
    elif mutation == 'missing_hole': value['free_branch_census'][0]['charge_branches'].pop()
    elif mutation == 'clock': value['free_branch_convention']['clock'] = 'uncentered original H'
    elif mutation == 'spin_mass': value['ordered_spin_branch']['invariant_rest_energy_squared'] = '0'
    elif mutation == 'ratio': value['zero_over_zero_ratio_defined'] = True
    elif mutation == 'native': value['finite_charge_operator_bridge']['native_donor_pulse_word_executed'] = True
    elif mutation == 'topology': value['finite_charge_operator_bridge']['real_space_topological_defect_identified'] = True
    elif mutation == 'stationarity': value['finite_charge_operator_bridge']['charge_conservation_implies_reference_stationarity'] = True
    else: value['verification_checks_pass'] = 1
    monkeypatch.setattr(packet, 'build', lambda: report)
    with pytest.raises(ValueError, match='admission changed'):
        packet.validate(value)


@pytest.mark.parametrize('edges', [2, 3, 4, 8, 31, 100])
def test_common_clock_obstruction_is_exact(edges):
    row = bounds.charged_clock_bound(edges, F(7, 3))
    coefficient = F(row['charged_energy_lower_per_exchange'])
    exact = F(row['rescaled_exchange_upper_if_a_charged_excitation_meets_cap'])
    simple = F(row['simpler_rescaled_exchange_upper'])
    assert coefficient == F(edges**8, 2)-edges**5
    assert coefficient >= F(edges**8, 4)
    assert exact*coefficient == F(7, 3) and exact <= simple


@pytest.mark.parametrize('edges', [True, 1, 2.0, -1])
def test_invalid_clock_domain_rejected(edges):
    with pytest.raises(ValueError): bounds.charged_clock_bound(edges)


def test_zero_control_cannot_define_a_mass_ratio():
    with pytest.raises(ValueError, match='zero invariant'):
        bounds.ratio(0, 0)
    assert bounds.ratio(F(1, 3), F(2, 3)) == F(1, 2)


@pytest.mark.parametrize('edges', [2, 4, 9, 40])
def test_finite_energy_charged_weight_bound(edges):
    exchange, cap = F(2, 3), F(3, 7)
    row = bounds.charged_weight_bound(edges, exchange, cap)
    gap = exchange*(F(edges**8, 2)-edges**5)
    assert F(row['charged_weight_upper']) == cap/gap
    assert F(row['charged_weight_upper']) <= F(row['simpler_charged_weight_upper'])
    # A two-level positive-clock control saturates C/Delta exactly.
    charged_probability = cap/gap
    assert charged_probability*gap == cap
    assert row['positive_ground_subtracted_clock'] == 'L=K-E_ground >= 0'


def test_charged_weight_zero_energy_and_clipping():
    assert bounds.charged_weight_bound(2, 1, 0)['charged_weight_upper'] == '0'
    assert bounds.charged_weight_bound(2, 1, 10**6)['charged_weight_upper'] == '1'
    with pytest.raises(ValueError): bounds.charged_weight_bound(2, 1, -1)
    with pytest.raises(ValueError): bounds.charged_weight_bound(2, 0, 1)


@pytest.mark.parametrize('modes', [1, 2, 3, 4])
def test_all_donor_inputs_and_charge_sectors(modes):
    dim = 1 << modes
    total = s.diag(*(i.bit_count() for i in range(2*dim)))
    for mode in range(modes):
        unitary = charge.donor_unitary(modes, mode)
        creator = charge.creator(modes, mode)
        assert unitary.T*unitary == s.eye(2*dim)
        assert total*unitary == unitary*total
        assert unitary[:dim, dim:] == creator
        assert unitary[dim:, dim:] == creator*creator.T


def test_real_source_charge_bridge_has_weight_without_particle_promotion(report):
    bridge = report['finite_charge_operator_bridge']
    assert bridge['actual_finite_reference_has_sharp_Q0']
    assert not bridge['sharp_charge_inferred_from_one_point_mean']
    assert all(row['addition_weight'] == row['removal_weight'] == '1/2'
               for row in bridge['actual_centered_reference_weights'])
    assert s.simplify(1-s.sympify(bridge['centered_ground_return_after_unread_complete_C_repair'])).is_positive
    assert not bridge['relativistic_particle_or_mass_identified']


def test_dropped_source_branch_cannot_pass_as_a_complete_charge_readout(monkeypatch):
    original = charge.full_domain.finite_branches
    def omit(**kwargs):
        rows = original(**kwargs)
        rows.pop((0, 0))
        return rows
    monkeypatch.setattr(charge.full_domain, 'finite_branches', omit)
    with pytest.raises(ArithmeticError, match='record slots'):
        charge.certificate()


def test_explicit_check_file_read_before_expensive_replays(tmp_path):
    path = tmp_path/'missing.json'
    run = subprocess.run([sys.executable, '-m', 'model.particle_mass_readout_gate', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0 and str(path) in run.stderr


def test_duplicate_receipt_keys_rejected(tmp_path):
    path = tmp_path/'bad.json'
    path.write_text('{"schema": "oph.particle_mass_readout_gate.v1", "schema": "other"}')
    run = subprocess.run([sys.executable, '-m', 'model.particle_mass_readout_gate', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0 and 'duplicate JSON key' in run.stderr
