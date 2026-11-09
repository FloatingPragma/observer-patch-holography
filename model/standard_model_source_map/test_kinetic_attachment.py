"""Source-current frames, physical charge conjugation and retained mirror weight."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import kinetic_attachment as ka
from model.common_source.chiral_candidate import kinetic


def test_hypercharge_fifth_moment_from_independent_five_mode_occupations():
    fundamental = (-2, -2, -2, 3, 3)
    charges = [sum(fundamental[j] for j in range(5) if mask >> j & 1)
               for mask in range(32) if mask.bit_count() in (2, 4)]
    actual = ka.mirror_invariant()
    expected = {}
    for q in charges:
        y = str(s.Rational(q, 6))
        expected[y] = expected.get(y, 0)+1
    assert actual['weights'] == expected
    assert sum(F(q, 6)**5 for q in charges) == F(5, 9)
    assert sum(F(q, 6)**3 for q in charges) == 0
    assert sorted(charges) != sorted(-q for q in charges)
    assert not actual['fifth_trace_is_a_four_dimensional_gauge_anomaly']


def test_every_source_valley_has_the_same_current_intertwiner():
    row = ka.current_frames()
    assert row['same_current_intertwining_identities'] == 96
    assert row['principal_symbol_axis_identities'] == 24
    assert len(row['native_matter_masks']) == 15
    assert sum(v['chirality'] for v in row['rows']) == 0
    assert row['inherited_configured_coordinate_speed'] == 'log(5)/2'
    assert not row['interacting_vacuum_or_continuum_attachment_inferred']


@pytest.mark.parametrize('valley', [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)])
def test_frame_against_independent_supplied_numeric_weyl_operator(valley):
    row = ka.valley_frame(valley)
    p = (F(1, 3), F(-2, 5), F(1, 7))
    shifted = tuple((-1)**b*float(x) for b, x in zip(valley, p, strict=True))
    source = kinetic.weyl_kinetic(shifted, speed=1)
    # The configured menu is left-handed; multiply by -chi for the canonical sign.
    target = -row['chirality']*kinetic.weyl_kinetic(tuple(map(float, p)), speed=1)
    change = np.asarray(row['map'], dtype=complex)
    assert np.max(np.abs(target@change-change@source)) < 1e-14
    assert np.max(np.abs(change.conj().T@change-np.eye(90))) == 0


def test_lorentz_current_transform_preserves_metric_not_pointwise_CAR_norm():
    row = ka.lorentz_currents()
    assert row['finite_spinor_current_identities'] == 24
    spinor = s.Matrix([1, 0])
    boost = s.diag(3, s.Rational(1, 3))
    assert boost.det() == 1
    assert (spinor.H*boost.H*boost*spinor)[0] == 9
    assert not row['pointwise_boost_is_unitary_on_equal_time_CAR']
    assert not row['finite_momentum_menu_closed_under_all_Lorentz_boosts']


def _cell(q, side):
    # Half-open cell about zero: -pi/2 <= k < pi/2, modulo the torus.
    return tuple(int(not (x < side//4 or x >= 3*side//4)) for x in q)


@pytest.mark.parametrize('side', [4, 8])
def test_exact_original_valley_probabilities_include_holes_at_opposite_momentum(side):
    counts = {b: 0 for b in product((0, 1), repeat=3)}
    nonzero = {b: 0 for b in counts}
    records = 0
    for q in product(range(side), repeat=3):
        zero = all(x in (0, side//2) for x in q)
        # Trace over both spin rows, separately for addition and removal.
        for addition in (False, True):
            momentum = q if addition else tuple((-x) % side for x in q)
            cell = _cell(momentum, side)
            counts[cell] += 1
            nonzero[cell] += int(not zero)
            records += 2
    result = ka.local_valley_weights(side)
    assert result['complete_charge_spin_momentum_records'] == records
    for row in result['rows']:
        b = tuple(row['node_pi_units'])
        assert F(counts[b], 2*side**3) == row['original_combined_weight']
        assert F(nonzero[b], 2*side**3) == row['original_strictly_positive_energy_weight']
    assert result['one_chirality_strictly_positive_energy_weight'] > 0


def test_full_four_row_car_instrument_is_complete_without_spin_postselection():
    # Two actual local CAR modes, with the second mode's occupied-prefix sign.
    ladders = []
    for j in range(2):
        a = s.zeros(4)
        for mask in range(4):
            if mask >> j & 1:
                a[mask ^ (1 << j), mask] = (-1)**((mask & ((1 << j)-1)).bit_count())
        ladders.append(a)
    kraus = [a/s.sqrt(2) for a in ladders]+[a.H/s.sqrt(2) for a in ladders]
    assert sum((a.H*a for a in kraus), s.zeros(4)) == s.eye(4)
    assert ladders[0]*ladders[1]+ladders[1]*ladders[0] == s.zeros(4)


def test_free_sine_symbol_projectors_give_exact_uniform_total_measure():
    pauli = ka.geometry.pauli()
    for q in product(range(4), repeat=3):
        h = -sum((s.sin(s.pi*x/2)*sigma for x, sigma in zip(q, pauli, strict=True)), s.zeros(2))
        energy2 = s.trace(h*h)/2
        if energy2 == 0:
            occupied = s.diag(1, 0)
        else:
            occupied = (s.eye(2)-h/s.sqrt(energy2))/2
        empty = s.eye(2)-occupied
        assert (occupied*occupied-occupied).applyfunc(s.simplify) == s.zeros(2)
        assert s.trace(occupied) == s.trace(empty) == 1
        assert (s.trace(occupied)+s.trace(empty))/2 == 1


def test_actual_valley_projector_fails_to_reduce_a_local_hypercharge_current():
    points = list(product(range(4), repeat=3))
    modes = [q for q in points if _cell(q, 4) == (0, 0, 0)]
    columns = [s.Matrix([s.I**sum(a*b for a, b in zip(q, x, strict=True))/8
                         for x in points]) for q in modes]
    fourier = s.Matrix.hstack(*columns)
    assert fourier.H*fourier == s.eye(8)
    projection = fourier*fourier.H
    site = s.zeros(64); site[0, 0] = 1
    commutator = projection*site-site*projection
    hs2 = s.expand(s.trace(commutator.H*commutator))
    assert hs2 == s.Rational(7, 32)
    source_y2 = sum(F(q, 6)**2 for q in (-4,)*3+(-3,)*2+(1,)*6+(2,)*3+(6,))
    row = ka.local_gauge_obstruction()
    assert hs2*s.Rational(source_y2) == s.Rational(row['one_spin_one_generation_commutator_HS_squared'])
    assert row['one_spin_one_generation_commutator_HS_squared'] == '35/48'


@pytest.mark.parametrize('bad', [True, 1, 2, 5, -4, F(4)])
def test_invalid_spatial_grids_are_not_silently_coerced(bad):
    with pytest.raises(ValueError):
        ka.local_valley_weights(bad)


def test_current_dictionary_mutation_is_rejected(monkeypatch):
    original = ka.matter.action_intertwiner
    def changed():
        matrix, rows = original()
        matrix = s.MutableDenseMatrix(matrix)
        matrix.row_swap(0, 1)
        return matrix, rows
    monkeypatch.setattr(ka.matter, 'action_intertwiner', changed)
    with pytest.raises(ArithmeticError, match='current permutation'):
        ka.current_frames()


def test_public_probability_values_do_not_alias_future_calls():
    row = ka.local_valley_weights()
    row['rows'][0]['original_combined_weight']._numerator = 99
    assert ka.local_valley_weights()['rows'][0]['original_combined_weight'] == F(1, 8)


def test_full_certificate_is_JSON_safe_without_changing_exact_helpers():
    import json
    value = ka.certificate()
    encoded = json.dumps(value, sort_keys=True, allow_nan=False)
    decoded = json.loads(encoded)
    weights = decoded['complete_free_local_valley_weights']
    assert all(row['original_combined_weight'] == '1/8' for row in weights['rows'])
    assert all(row['original_strictly_positive_energy_weight'] == '7/64' for row in weights['rows'])
    assert weights['one_chirality_combined_weight'] == '1/2'
    assert weights['one_chirality_strictly_positive_energy_weight'] == '7/16'
    assert ka.local_valley_weights()['one_chirality_combined_weight'] == F(1, 2)
