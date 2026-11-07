from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import reference as r


@pytest.mark.parametrize('weight', [F(), F(1, 7), F(1, 2), F(5, 6), F(1)])
def test_full_quantum_stationary_state_including_memory_coherences(weight):
    state = r.stationary_entries(weight)
    assert sum(state[i, i] for i in range(132)) == 1
    assert r.apply_exact(state) == {key: value for key, value in state.items() if value}
    # An independent Q expectation reads both off-diagonal memory blocks.
    expected_charge = sum(state[i, 66+i]+state[66+i, i] for i in range(66))
    assert expected_charge == 2*weight-1
    pi = r.local.reference(r.parameters())
    assert tuple(state[i, i]+state[66+i, 66+i] for i in range(66)) == pi


def test_actual_complete_maps_conserve_charge_on_every_column_and_reference():
    total = [F()]*132
    slots = 0
    for row in r.branches():
        mapping = dict(row.cover_pairs())
        for i, j in mapping.items():
            total[i] += row.probability
            assert mapping[(i+66) % 132] == (j+66) % 132
        if row.abort:
            assert not mapping
        slots += 1
    assert slots == 720
    assert total == [F(1)]*132
    # Preservation on a state carrying an actual memory coherence.
    initial = {(0, 0): F(1, 2), (66, 66): F(1, 2),
               (0, 66): F(1, 2), (66, 0): F(1, 2)}
    out = r.apply_exact(initial)
    assert sum(out.get((i, 66+i), F())+out.get((66+i, i), F()) for i in range(66)) == 1


def test_reading_or_resetting_orientation_is_a_real_protection_failure():
    # Computational dephasing sends either sharp charge to the same state.
    plus, minus = r.stationary_entries(F(1)), r.stationary_entries(F())
    dephase = lambda state: {key: value for key, value in state.items() if key[0]//66 == key[1]//66 and value}
    assert dephase(plus) == dephase(minus)
    assert sum(dephase(plus).get((i, 66+i), F()) for i in range(66)) == 0


@pytest.mark.parametrize('weight,ratio', [(F(1, 4), F(1, 3)), (F(1, 2), F(1)), (F(3, 4), F(3))])
def test_charge_flip_energy_sign_and_exact_reverse(weight, ratio):
    assert r.modular_gap_ratio(7, 7, 1, -1, weight) == ratio
    assert r.modular_gap_ratio(7, 7, -1, 1, weight) == 1/ratio
    plus = r.stationary_entries(F(1))
    z_conjugate = {(i, j): (-1)**(i//66+j//66)*value for (i, j), value in plus.items()}
    assert z_conjugate == r.stationary_entries(F())


@pytest.mark.parametrize('weight,charge', [(F(), -1), (F(1), 1)])
def test_corner_clock_has_no_energy_for_absent_charge(weight, charge):
    assert r.modular_gap_ratio(0, 0, charge, charge, weight) == 1
    with pytest.raises(ValueError, match='support'):
        r.modular_gap_ratio(0, 0, charge, -charge, weight)


def test_every_raw_branch_is_modular_eigenoperator_on_complete_cover():
    # Diagonalize only the two-dimensional X, independently of the branch
    # producer; retain all native input/output columns and both charge blocks.
    pi = r.local.reference(r.parameters())
    p = r.parameters()
    for row in r.branches():
        ratios = set()
        for i, j in zip(row.inputs, row.outputs, strict=True):
            before, after = r.local.basis()[i], r.local.basis()[j]
            current = before[row.seam[0]]-after[row.seam[0]]
            ratios.add(pi[i]/pi[j])
            assert pi[i]/pi[j] == F(p.port_weights[row.seam[0]], p.port_weights[row.seam[1]])**current
        assert len(ratios) <= 1


def test_covariance_on_noncommuting_density_and_both_complete_charge_blocks():
    rng = np.random.default_rng(913)
    initial = rng.normal(size=132)+1j*rng.normal(size=132)
    initial /= np.linalg.norm(initial)
    x_eigenvectors = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    chart = np.kron(x_eigenvectors, np.eye(66))
    # A finite faithful reference with distinct charge and occupation gaps.
    probabilities = np.array([float(x) for x in r.local.reference(r.parameters())])
    phases = np.exp(0.37j*np.log(np.r_[probabilities/3, 2*probabilities/3]))
    unitary = (chart*phases)@chart.T
    for row in r.branches():
        if row.seam != (0, 5) or not row.inputs:
            continue
        operator = np.zeros((132, 132))
        for i, j in row.cover_pairs():
            operator[j, i] = np.sqrt(float(row.probability))
        first = operator@unitary@initial
        second = unitary@operator@initial
        # Compare complete subnormalized density outputs, so branch phases
        # cancel but probability loss cannot be hidden by normalization.
        assert np.linalg.norm(np.outer(first, first.conj())-np.outer(second, second.conj())) < 1e-14


def test_raw_current_and_car_join_is_sensitive_to_sign_corruption(monkeypatch):
    rows = list(r.branches())
    index = next(i for i, row in enumerate(rows) if -1 in row.signs)
    row = rows[index]
    rows[index] = replace(row, signs=tuple(1 for _ in row.signs))
    monkeypatch.setattr(r, 'branches', lambda: tuple(rows))
    with pytest.raises(ArithmeticError, match='CAR'):
        r.current_certificate()


def test_odd_exchange_has_holonomy_and_zero_telescope_energy():
    for law in ('C', 'E'):
        loop = r.protection.exchange_loop(law)
        initial = tuple(loop['initial_pair'])
        state = r.exchange.INDEX[initial]
        ratio, sign = F(1), 1
        for name, seam, coin, raw in loop['history']:
            row = next(row for row in r.branches() if (row.law, row.seam, row.coin, row.raw)
                       == (name, tuple(seam), coin, raw))
            i = row.inputs.index(state)
            final = row.outputs[i]
            ratio *= r.modular_gap_ratio(state, final, -1, -1, F(2, 3))
            sign *= row.signs[i]
            state = final
        assert r.exchange.PAIRS[state] == initial
        assert ratio == 1 and sign == -1


def test_complete_tape_cannot_distinguish_sharp_references():
    pi = r.local.reference(r.parameters())
    # Compare every two-event history on two actual seams, including zeros.
    first = [row for row in r.branches() if row.law == 'E' and row.seam == (0, 5)]
    second = [row for row in r.branches() if row.law == 'C' and row.seam == (0, 1)]
    for start in first:
        for end in second:
            next_map = {i: (j, sign) for i, j, sign in zip(end.inputs, end.outputs, end.signs, strict=True)}
            plus, minus = {}, {}
            for i, middle, sign in zip(start.inputs, start.outputs, start.signs, strict=True):
                if middle not in next_map:
                    continue
                final, sign2 = next_map[middle]
                probability = start.probability*end.probability*pi[i]
                plus[final] = plus.get(final, F())+probability
                minus[final] = minus.get(final, F())+probability*(sign*sign2)**2
            assert plus == minus


def test_sharp_sector_preparation_reuses_full_quantum_bound():
    value = r.preparation_certificate()
    assert value['attempts'] == 1600
    assert F(value['source_and_untouched_reference_half_trace_distance_upper']) < F(1, 10**6)
    assert value['preparation_archive_traced_only_for_this_comparison']


@pytest.mark.parametrize('bad', [0.5, 1, True, '1/2', F(-1), F(2)])
def test_original_exact_probability_required(bad):
    with pytest.raises(ValueError):
        r.stationary_entries(bad)


def test_certificate_and_scientific_boundaries():
    value = r.certificate()
    assert json.loads(json.dumps(value)) == value
    assert value['complete_stationary_space_complex_dimension'] == 2
    assert value['current_and_energy']['complete_record_instrument_covariant_under_same_modular_clock']
    assert not value['scope']['sharp_charge_sector_is_a_spatially_localized_source_defect']
    assert not value['scope']['interacting_Higgs_action_or_state_selected']
    assert not value['scope']['modular_clock_equals_repair_attempt_time_proved']
    thermo = value['same_corner_finite_thermodynamics']
    assert thermo['kernel'] == [['25/33', '8/33'], ['25/33', '8/33']]
    assert not thermo['physical_energy_calibration_proved']
    assert not thermo['refinement_uniform_low_temperature_family_produced']
    reverse = value['weighted_quantum_record_reversal']
    assert reverse['branches']['raw_slots'] == 720
    assert reverse['branches']['complete_reverse_channel_equals_original']
    assert not reverse['scope']['petz_reversal_is_physical_time_reversal']


def test_thermodynamic_adapter_rejects_wrong_coin_or_corner(monkeypatch):
    from . import gauge_join
    rows = [{'primitive_charge': charge, 'target_probability': '25/33',
             'source_probability': '8/33',
             'two_orthonormal_singlets_form_reducing_corner': True,
             'one_unread_attempt_exact_preparation_in_each_sharp_charge_corner': True}
            for charge in r.gauge_transport.NAMES]
    joined = {'Gauss_singlet_stationary_corners': rows}
    assert r.finite_thermodynamic_certificate(joined)['finite_four_law_core_inputs_constructed']
    rows[0]['source_probability'] = '1/2'
    with pytest.raises(ArithmeticError, match='differs'):
        r.finite_thermodynamic_certificate(joined)
    rows[0]['source_probability'] = '8/33'
    native = gauge_join.native_edge()
    monkeypatch.setattr(gauge_join, 'native_edge', lambda: dict(native, forward_probability=F(1, 2)))
    with pytest.raises(ArithmeticError, match='weights'):
        r.finite_thermodynamic_certificate(joined)
