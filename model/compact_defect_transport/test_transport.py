from fractions import Fraction as F
from itertools import product
from math import atan2, pi

import numpy as np
import pytest
from scipy import linalg as la

from . import bounds, source, transport


def test_all_continuing_records_preserve_two_degrees_and_complete_probability():
    states, rows, transition = source.continuing_protocol()
    assert len(states) == 3 and len(rows) == 12
    assert transition == [[F(1, 2), F(1, 4), F(1, 4)],
                          [0, F(3, 4), F(1, 4)], [0, F(1, 4), F(3, 4)]]
    stationary = [F(0), F(1, 2), F(1, 2)]
    assert [sum(stationary[i]*transition[i][j] for i in range(3))
            for j in range(3)] == stationary
    for state in states:
        for triangles, expected in zip(source.surfaces(), (-1, 1)):
            area = 0.
            for triangle in triangles:
                u, v, w = [np.asarray(state[x], dtype=float) for x in triangle]
                u, v, w = (x/la.norm(x) for x in (u, v, w))
                area += 2*atan2(np.linalg.det([u, v, w]), 1+u@v+v@w+w@u)
            assert abs(area/(4*pi)-expected) < 1e-12
    assert all(row['local_shell_degrees'] == [-1, 1] for row in rows)
    assert all(row['incident_triangles_checked'] == 12 for row in rows
               if row['point'] == [0, 0, 0])


def test_guards_fail_closed_and_rejections_preserve_the_literal_state():
    state = source.initial()
    for point, proposal, reason in [((0, 0, 0), (0, 0, 0), 'zero_direction'),
                                    ((0, 0, 1), (0, 0, 1), 'fixed_core'),
                                    ((0, 0, 0), (1, 0, -1), 'reference_budget'),
                                    ((3, 0, 0), (1, 0, 0), 'fixed_exterior')]:
        target, row = source.update(state, point, proposal)
        assert target == state and not row['accepted'] and row['reason'] == reason
    bad, row = source.update(state, (-1, -1, -1), (-1, -1, -F(1003, 3000)), eta=F(6, 25))
    assert bad == state and row['reason'] == 'shell_homotopy'
    for invalid in (0., True, F(-1), F(1, 4)):
        with pytest.raises(ValueError):
            source.validate(state, eta=invalid)
    malformed = dict(state)
    del malformed[0, 0, 0]
    with pytest.raises(ValueError):
        source.validate(malformed)


def test_erasing_texture_is_allowed_by_spectral_count_but_changes_winding():
    value = source.exact_control()
    assert value['unguarded_erasure_degrees'] == [0, 0]
    assert value['erasure_retains_four_zero_onsite_coordinates']
    state = source.initial()
    erased = {p: ((0, 0, 1) if any(v) else v) for p, v in state.items()}
    assert sum(not any(v) for v in erased.values()) == 2
    with pytest.raises(ValueError, match='budget'):
        source.validate(erased)


def test_rational_resolvent_bounds_and_round_domains():
    first = bounds.certificate()
    beta = F(first['spectral_projection_difference_norm_upper'])
    assert beta == F(3572880, 221206129) < 1
    assert F(first['single_quench_leakage_probability_upper']) == beta**2
    assert F(bounds.certificate(rounds=4)['uncontrolled_history_leakage_probability_upper']) == 16*beta**2
    assert bounds.certificate(rounds=0)['uncontrolled_history_leakage_probability_upper'] == '0'
    assert first['controlled_record_endpoint_band_leakage_probability'] == '0'
    for kwargs in ({'mass': 7}, {'eta': F(1, 5)}, {'rounds': True}, {'rounds': -1}):
        with pytest.raises(ValueError):
            bounds.certificate(**kwargs)


def test_direct_rotation_against_full_spectral_polar_decomposition():
    rng = np.random.default_rng(123)
    generator = rng.normal(size=(10, 10))+1j*rng.normal(size=(10, 10))
    generator = (generator+generator.conj().T)/2
    first = np.eye(10, 3, dtype=complex)
    second = la.expm(-.01j*generator)@first
    actual = transport.polar_transport(first, second)
    p, q = first@first.conj().T, second@second.conj().T
    s = q@p+(np.eye(10)-q)@(np.eye(10)-p)
    expected = s@la.inv(la.sqrtm(s.conj().T@s))
    assert la.norm(actual-expected) < 1e-12
    assert la.norm(actual@p@actual.conj().T-q) < 1e-12
    assert la.norm(transport.polar_transport(first, first)-np.eye(10)) < 1e-12
    with pytest.raises(ValueError):
        transport.polar_transport(first*2, second)


def test_adjacent_gate_replay_and_physical_path():
    points = tuple(product(range(-1, 2), repeat=3))
    path = transport.snake_path(points)
    rng = np.random.default_rng(456)
    unitary, _ = la.qr(rng.normal(size=(54, 54))+1j*rng.normal(size=(54, 54)))
    result = transport.local_compile(unitary, path)
    assert result['nearest_neighbor_two_mode_gates'] <= 54*53//2
    assert result['complete_replay_Frobenius_error_below_1e_minus_8']
    assert transport.local_compile(np.eye(54, dtype=complex), path)['nearest_neighbor_two_mode_gates'] == 0
    with pytest.raises(ValueError):
        transport.local_compile(unitary*2, path)
    with pytest.raises(ValueError):
        transport.local_compile(unitary, [0]*54)
