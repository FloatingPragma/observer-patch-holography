"""Same-source Ward coefficients, CAR signs and complete dressed records."""
from itertools import product

import pytest
import sympy as s

from . import gauge_attachment as ga


def _clean(state):
    return {key: s.expand(value) for key, value in state.items() if s.expand(value) != 0}


def _add(out, key, value):
    out[key] = out.get(key, 0)+value


def _bilinear(mask, p, q):
    if not mask >> q & 1:
        return None
    sign = (-1)**((mask & ((1 << q)-1)).bit_count())
    mask ^= 1 << q
    if mask >> p & 1:
        return None
    sign *= (-1)**((mask & ((1 << p)-1)).bit_count())
    return mask | (1 << p), sign


def _hop(state, axis):
    out = {}
    for (mask, label), amplitude in state.items():
        for (p, q, row, col), coefficient in ga.hopping_coefficients(axis).items():
            if label == col and (target := _bilinear(mask, p, q)) is not None:
                _add(out, (target[0], row), amplitude*coefficient*target[1])
    return _clean(out)


def _gauss(state, generator, endpoint):
    out = {}; labels = ga.link_labels(); lookup = {v: i for i, v in enumerate(labels)}
    for (mask, label), amplitude in state.items():
        for (a, b), value in generator.todok().items():
            for spin in range(2):
                target = _bilinear(mask, 30*endpoint+15*spin+a, 30*endpoint+15*spin+b)
                if target is not None:
                    _add(out, (target[0], label), amplitude*value*target[1])
            if label:
                x, y = labels[label]
                if endpoint == 0 and a == x:
                    _add(out, (mask, lookup[b, y]), -amplitude*value)
                if endpoint == 1 and b == y:
                    _add(out, (mask, lookup[x, a]), amplitude*value)
    return _clean(out)


def test_all_actual_currents_and_full_cubic_tensor():
    row = ga.anomaly_certificate()
    assert row['symmetric_cubic_tensor_components'] == 364
    assert row['linear_trace_components'] == 12
    assert row['weak_doublets_counting_color'] == 4
    assert not row['regulated_chiral_quantum_measure_constructed']
    # An individual source block has a nonzero abelian anomaly. Cancellation
    # uses every M15 block, not a traceless single-generator assumption.
    center = ga.currents()[-1]
    parts = [s.trace(center.extract(list(block), list(block))**3)
             for _, block in ga.native.native_blocks()]
    assert any(parts) and s.expand(sum(parts)) == 0


def test_full_coefficients_exhaust_all_axes_currents_and_link_rows():
    row = ga.ward_certificate()
    assert row['full_Fock_Ward_identities'] == 72
    assert row['electric_projector_Ward_identities'] == 24
    assert row['edge_Fock_dimension'] == str(2**60)
    assert row['coefficient_array_dimension'] == 3600
    assert row['wrong_orientation_endpoint_failures'] == 2
    assert not row['physical_Fock_space_truncated_to_one_particle']
    for axis in range(3):
        h = ga.coefficient_matrix(axis)
        assert h.H == h
        assert len(ga.hopping_coefficients(axis)) == 236


@pytest.mark.parametrize('axis', range(3))
def test_actual_sixty_mode_car_with_spectators_and_every_native_current(axis):
    # Sparse exact vectors in the full sixty-mode Fock space, not matrices
    # in a one-particle truncation. Occupied spectators exercise CAR signs.
    masks = [0, (1 << 60)-1, (1 << 31) | (1 << 2),
             (1 << 45) | (1 << 12) | (1 << 19),
             sum(1 << j for j in range(0, 60, 3))]
    for mask, label in product(masks, (0, 1, 40, 59)):
        state = {(mask, label): s.Integer(1)}
        for generator in ga.currents():
            for endpoint in (0, 1):
                assert _gauss(_hop(state, axis), generator, endpoint) == _hop(
                    _gauss(state, generator, endpoint), axis)


def test_literal_dressed_reduction_checks_the_owning_source():
    row = ga.dressed_reduction_certificate()
    assert row['literal_native_Gauss_state_checks'] == 480
    assert row['literal_Weyl_link_Hamiltonian_columns'] == 60
    assert row['restricted_polynomial_at_unit_parameters'] == 'H^2-H=d I/4'
    assert not row['finite_transport_frequency_identified_with_particle_mass']


@pytest.mark.parametrize('axis', range(3))
def test_all_block_position_spin_records_keep_original_weights(axis):
    row = ga.dressed_history(axis)
    assert len(row['records']) == 21
    assert row['measurement_complete_on_full_Fock']
    assert row['records'][-1]['outcome'] == 'other_occupation'
    assert sum(s.Rational(r['original_probability']) for r in row['records']) == 1
    assert any(s.Rational(r['original_probability']) == 0 for r in row['records'])
    for _, block in ga.native.native_blocks():
        selected = [r for r in row['records'] if r['block_indices'] == list(block)]
        assert len(selected) == 4
        assert sum(s.Rational(r['original_probability']) for r in selected) == s.Rational(len(block), 15)
    assert not row['initial_block_measurement_renormalized']
    assert not row['native_pulse_execution_claimed']
    # A sigma_z edge preserves spatial spin; sigma_x/y flip it on transport.
    for r in row['records']:
        if r['position'] == 0 and r['spin'] == (0 if axis != 2 else 1):
            assert r['original_probability'] == '0'


def test_zero_clock_keeps_zero_transport_and_all_twenty_one_slots():
    row = ga.dressed_history(0, duration=0)
    assert len(row['records']) == 21
    assert sum(s.Rational(r['original_probability']) for r in row['records'] if r['position'] == 0) == 0


def test_full_fock_record_projectors_include_all_other_occupations():
    blocks = [set(indices) for _, indices in ga.native.native_blocks()]
    for mask in [0, (1 << 60)-1, 3]+[1 << mode for mode in range(60)]:
        one = mask.bit_count() == 1
        mode = mask.bit_length()-1
        effects = [int(one and mode//30 == position and mode//15 % 2 == spin
                       and mode % 15 in block)
                   for block, position, spin in product(blocks, range(2), range(2))]
        effects.append(int(not one))
        assert len(effects) == 21 and sum(effects) == 1


def test_quantum_link_vacuum_cannot_be_replaced_by_identity_link():
    # Every off-diagonal link is nonzero nilpotent; none is a unitary Wilson
    # link or has a nonzero eigenvalue that could freeze its action to I.
    for index in range(1, 60):
        u = s.SparseMatrix(60, 60, {(index, 0): 1})
        assert u*u == s.zeros(60)
        assert u.H*u != s.eye(60)


@pytest.mark.parametrize('axis', [True, -1, 3, 1.0])
def test_axis_admission(axis):
    with pytest.raises(ValueError):
        ga.hopping_coefficients(axis)


@pytest.mark.parametrize('duration', [-1, 0.125, s.Symbol('t')])
def test_clock_admission(duration):
    with pytest.raises(ValueError):
        ga.dressed_history(0, duration=duration)


def test_independent_anomaly_component_detects_a_changed_current(monkeypatch):
    values = list(ga.currents()); values[-1] = values[-1]+s.I*s.eye(15)
    monkeypatch.setattr(ga, 'currents', lambda: tuple(values))
    with pytest.raises(ArithmeticError, match='anomaly'):
        ga.anomaly_certificate()


def test_bad_link_orientation_is_an_actual_nonzero_ward_operator():
    generator = ga.currents()[-1]
    h = ga.coefficient_matrix(2)
    broken = ga.gauss_coefficients(generator, 0, orientation=-1)
    assert (broken*h-h*broken).is_zero_matrix is False


def test_result_objects_do_not_share_mutable_record_state():
    first = ga.dressed_history(0)
    expected = repr(ga.dressed_history(0))
    first['records'][0]['block_indices'].append(999)
    assert repr(ga.dressed_history(0)) == expected
