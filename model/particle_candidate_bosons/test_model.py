"""Exact source covariance, full-menu normalization and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_literal_source_face_irreps_and_selection():
    row = m.symmetry_control()
    assert row['actual_binary_cube_rotations_checked'] == 48
    assert row['multiplet_dimensions'] == {'A1g': 1, 'Eg': 2, 'T1u': 3,
                                         'A1u': 1, 'Eu': 2, 'T1g': 3}
    assert len(row['loop_channels']) == 12
    assert [x['gauge_coordinate_conjugation'] for x in row['loop_channels']] == [1]*6+[-1]*6
    assert all(x['CAR_parity'] == 0 for x in row['loop_channels'])
    assert row['bare_scalar_even_to_A1u_overlap'] == row['bare_vector_to_A1u_overlap'] == '0'


def test_exact_laurent_identity_keeps_the_twelve_rows_in_one_menu():
    row = m.completeness_control()
    assert row['parent_sum_F_dagger_F'] == 'I'
    assert row['blank_anchor_Gram'] == 'I_12/12'
    assert row['success_rows']+row['failure_rows'] == 13
    assert row['extra_uniform_selector_divisor'] == 1
    assert not row['compressed_menu_assumed_exactly_saturated']


def _gauss_flows():
    """Independent bounded integer-flow basis for the charge-six control."""
    edges = m.reduced.cube_support()['interior_edges']
    incidence = s.zeros(7, 12)
    for i, (x, y) in enumerate(edges):
        if x < 7:
            incidence[x, i] = -1
        if y < 7:
            incidence[y, i] = 1
    pivots = incidence.rref()[1]
    free = [i for i in range(12) if i not in pivots]
    lift = -incidence[:, list(pivots)].inv()*incidence[:, free]
    states = []
    for values in product((-1, 0, 1), repeat=len(free)):
        assigned = lift*s.Matrix(values)
        if any(x not in (-1, 0, 1) for x in assigned):
            continue
        state = [0]*12
        for i, x in zip(pivots, assigned):
            state[i] = int(x)
        for i, x in zip(free, values):
            state[i] = x
        assert incidence*s.Matrix(state) == s.zeros(7, 1)
        states.append(tuple(state))
    return states


def test_discarded_word_gram_is_the_exact_completeness_deficit():
    states = _gauss_flows()
    assert len(states) == 69
    where = {q: i for i, q in enumerate(states)}
    projected = s.zeros(69)
    discarded = s.zeros(69)
    operators = []
    for terms in m.loop_polynomials():
        full = {}
        for j, q in enumerate(states):
            for word, a in terms.items():
                target = tuple(x+y for x, y in zip(q, word))
                full.setdefault(target, {})[j] = a
        C = s.zeros(69)
        for q, values in full.items():
            if q in where:
                for j, a in values.items():
                    C[where[q], j] = a
            else:
                for j, a in values.items():
                    for k, b in values.items():
                        discarded[j, k] += s.conjugate(a)*b
        projected += C.H*C
        operators.append(C)
    assert (projected+discarded).applyfunc(m.clean) == s.eye(69)
    assert discarded.applyfunc(m.clean) != s.zeros(69)
    # Numerical eigenvalues are only a regression on this omitted-native control.
    assert np.linalg.eigvalsh(np.array(projected.evalf(), complex)).max() < 1+1e-12
    assert all(C == C.H for C in operators)


def test_source_basis_normalization_error_is_rejected(monkeypatch):
    original = m.detector_basis()
    wrong = original.copy()
    wrong[:, 0] *= 2
    monkeypatch.setattr(m, 'detector_basis', lambda: wrong)
    m.loop_polynomials.cache_clear()
    m.completeness_control.cache_clear()
    try:
        with pytest.raises(ArithmeticError, match='sum to identity'):
            m.completeness_control()
    finally:
        m.loop_polynomials.cache_clear()
        m.completeness_control.cache_clear()


def test_native_derived_current_twirl_and_central_exception():
    row = m.native_current_control()
    assert row['actual_native_generator_rank'] == 12
    assert row['derived_su2_plus_su3_rank'] == 11
    assert not row['abelian_current_killed_by_adjoint_twirl_claimed']
    native = m.reduced.channels.native.native_currents()
    # These currents use increasing CAR masks; the packet helper uses a
    # different source chart. Lift each finite group element in this chart.
    subsets = [tuple(i for i in range(5) if mask >> i & 1)
               for mask in range(32) if mask.bit_count() in (2, 4)]
    def representation(g):
        return s.Matrix([[g.extract(a, b).det() if len(a) == len(b) else 0
                          for b in subsets] for a in subsets])
    rotations = {tuple(int(m.reduced.channels.group.rotation(q)[i, j])
                       for i in range(3) for j in range(3))
                 for q in m.reduced.channels.group.elements()}
    color = [representation(s.diag(s.Matrix(3, 3, r), s.eye(2))) for r in rotations]
    pauli = (s.eye(2), s.I*s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, 1], [-1, 0]]), s.I*s.diag(1, -1))
    weak = [representation(s.diag(s.eye(3), u)) for u in pauli]
    for j, current in enumerate(native):
        a = sum((u*current*u.H for u in color), s.zeros(15))/24
        b = (sum((u*a*u.H for u in weak), s.zeros(15))/4).applyfunc(m.clean)
        assert b == (current if j == 11 else s.zeros(15))


def test_uniform_full_complement_residue_bounds_are_outward():
    row = m.spectral_control()
    assert row['electric_anchor_outside_cluster_norm_upper'] == F(4, 143)
    assert row['each_actual_ground_to_multiplet_amplitude_strict_lower'] > F(2, 7)
    assert row['each_actual_ground_to_multiplet_residue_strict_lower'] == F(4, 49)
    assert row['multiplet_transition_energy_interval'] == [F(435, 7), F(447, 7)+F(1, 45675)]
    assert row['all_six_source_transition_ratios_upper'] == F(149, 145)
    for ground in (-F(1, 45675), F(0)):
        assert (F(447, 7)-ground)/(F(435, 7)-ground) <= F(149, 145)
    assert not row['six_multiplet_energies_proved_pairwise_distinct']
    assert not row['exact_H_replaced_by_a_twelve_state_compression']


@pytest.fixture(scope='module')
def bound():
    member = m.contract.member(0)
    return member, m.bind_member(member)


def test_original_complete_weights_do_not_condition_away_the_rare_sector(bound):
    member, value = bound
    r = value['result']
    assert r['member_id'] == member['member_id']
    assert all(r[k] == v for k, v in member['ids'].items())
    assert r['each_row_original_reported_weight_strict_lower_factored'] == {
        'base': 2, 'negative_exponent': 32406376, 'rational_multiplier': '1/13'}
    assert F(r['each_row_clock_weight_before_tiny_floor']) > F(1, 13)
    assert r['combined_twelve_row_original_reported_weight_strict_lower_multiplier'] == '12/13'
    assert sum(F(x) for x in r['multiplet_original_reported_weight_strict_lower_multipliers'].values()) == F(12, 13)
    assert not r['positive_window_normalized_by_success']


def test_zero_records_require_initial_symmetry_and_do_not_overstate_reference(bound):
    r = bound[1]['result']
    assert set(r['forbidden_original_record_categories'].values()) == {'0'}
    assert r['zero_records_require_the_stated_initial_symmetry_record']
    assert not r['zeroes_asserted_for_arbitrary_initial_native_sectors']
    assert not r['clock_leakage_breaks_these_exact_symmetry_zeroes']
    assert r['full_native_charge_conjugation']['CAR_modes'] == 480
    assert r['full_native_charge_conjugation']['empty_sector_maps_to_full_matter_sector']
    assert not r['full_native_charge_conjugation']['source_only_operation_fixes_retained_reference_state']
    assert not r['gauge_coordinate_conjugation_within_empty_sector']['is_full_native_charge_conjugation']


def test_no_physical_spin_species_or_mass_promotion(bound):
    r = bound[1]['result']
    assert r['physical_particle_map'] == r['physical_mass_table_GeV'] == []
    assert not r['intrinsic_Lorentz_spin_or_physical_CP_assigned']
    assert not r['uniform_refinement_particle_residue_proved']
    assert set(r['candidate_readouts']) == {'photon', 'gluon', 'W_plus_minus', 'Z', 'Higgs', 'A1u'}


def test_live_verifier_rejects_promoted_labels(bound):
    member, value = bound
    changed = deepcopy(value)
    changed['result']['physical_particle_map'] = ['Higgs']
    changed['result_sha256'] = m.contract.digest(changed['result'])
    with pytest.raises(ValueError, match='boson detector'):
        m.verify_bound(member, changed)
