"""Native incidence, global topology, work accounting and complete outputs."""
from fractions import Fraction as F

import pytest
import sympy as s

from model.particle_observer_outputs import holonomy_transport as h


def test_actual_geometry_and_order_two_upstream_holonomies():
    data = h.complex_data(); _, flux = h._upstream()
    assert (len(data['edges']), len(data['faces'])) == (30, 20)
    assert data['punctures'] == (0, 13)
    for e in range(30):
        actual = h.curvature(1 << e)
        upstream = flux.face_holonomies({e: 3}, data['support']['boundary'], 20)
        assert upstream == [3*((actual >> f) & 1) for f in range(20)]
        assert actual.bit_count() == 2
    assert data['ranks']['full_H1_dimension'] == 0


def test_every_vertex_gauge_and_complete_flat_space():
    gauges = {h.gauge_update(v) for v in range(1 << 12)}
    assert len(gauges) == 1 << 11
    for links in gauges:
        assert h.curvature(links) == 0
        assert h.gauge_update(h.flat_gauge_potential(links)) == links
    with pytest.raises(ValueError, match='nonzero curvature'):
        h.flat_gauge_potential(1)


def test_every_pair_of_faces_reuses_a_complete_transport_path():
    data = h.complex_data(); gate, flux = h._upstream()
    for first in range(20):
        for last in range(first+1, 20):
            path = flux.dual_path(data['support'], first, last)
            upstream = gate.verify_local_link_path(data['support'], (first, last), path, 3)
            links = 0
            for face, edge in path:
                links ^= 1 << edge
                assert h.curvature(links) == (1 << first) ^ (1 << face)
            assert upstream == {e: 3 for e in h._indices(links)}


def test_regional_stokes_and_exact_protecting_operation_condition():
    for region in range(1 << 9):
        boundary = h.region_boundary(region)
        for e in range(30):
            assert h.region_parity(1 << e, region) == (boundary >> e) & 1
            assert h.preserves_region(1 << e, region) == (not ((boundary >> e) & 1))
        for v in range(12):
            assert h.preserves_region(h.gauge_update(1 << v), region)
    assert h.region_boundary((1 << 20)-1) == 0
    rows = h.graphical_path()['stages']
    assert rows[1]['boundary_holonomy_sign'] == -1
    assert all(r['boundary_holonomy_sign'] == -1 and r['update_preserves_region']
               for r in rows[2:-1])
    assert not rows[1]['update_preserves_region']
    assert not rows[-1]['update_preserves_region']


def _local_states(edge):
    f, g = h.complex_data()['edge_faces'][edge]
    spectator = next(k for k in range(20) if k not in (f, g))
    # Padding one untouched face realizes every local input within global
    # even parity; it is identical on both sides of each local update.
    return tuple(((a << f) | (b << g) | ((a ^ b) << spectator), work)
                 for a in range(2) for b in range(2) for work in range(2))


def test_all_local_gate_columns_are_involutions_and_conserve_total_energy():
    for edge in range(30):
        states = _local_states(edge)
        for gate in (h.work_exchange, h.move_or_hold):
            outputs = {gate(*state, edge) for state in states}
            assert outputs == set(states)
            for state in states:
                final = gate(*state, edge)
                assert gate(*final, edge) == state
                assert h.energy_units(*final) == h.energy_units(*state)


def test_complete_coherent_creation_matrix_and_hold_columns():
    edge = 0; states = _local_states(edge); index = {v: i for i, v in enumerate(states)}
    matrix = s.zeros(8)
    for j, state in enumerate(states):
        for out, (real, imag) in h.creation_pulse(*state, edge).items():
            matrix[index[out], j] = s.Rational(real.numerator, real.denominator)+s.I*s.Rational(imag.numerator, imag.denominator)
    assert matrix.H*matrix == s.eye(8)
    energy = s.diag(*(h.energy_units(*state) for state in states))
    assert matrix*energy == energy*matrix
    assert any(matrix[j, j] == 1 for j in range(8))
    with pytest.raises(ValueError, match='unit circle'):
        h.creation_pulse(0, 1, edge, F(1, 2), F(1, 2))


def test_complete_history_retains_no_birth_and_zero_labels():
    report = h.complete_history()
    assert report['face_and_work_readout_alphabet_size'] == 1 << 21
    assert report['odd_face_parity_pointer_labels_with_zero_projector'] == 1 << 20
    assert report['additional_outside_Gauss_outcome']
    rows = report['nonzero_complete_histories']
    assert {r['birth_record']: F(r['probability']) for r in rows} == {False: F(9, 25), True: F(16, 25)}
    assert all(r['records'][-1] == {'defect_faces': [], 'work_register': 1} for r in rows)
    assert report['all_other_histories_have_exact_zero_weight']
    assert not report['postselection_used']


def test_certificate_separates_native_memory_and_added_curvature():
    result = h.certificate()
    assert result['native_boundary']['actual_current_branches_checked'] == 4680
    assert result['native_boundary']['actual_native_current_face_holonomy'] == 0
    assert not result['native_boundary']['unextended_C_E_creates_these_spatial_defects']
    assert not result['complete_spaces']['nontrivial_flat_sphere_holonomy']
    assert not result['protection']['gauge_invariance_alone_protects_defects']
    assert result['physical_particle_assignments'] == []
    assert all(r['total_energy_in_J'] == 2 for r in result['graphical_path']['stages'])


@pytest.mark.parametrize('value', [1, -1, 1 << 20, True])
def test_invalid_or_odd_curvature_is_rejected(value):
    with pytest.raises(ValueError):
        h.work_exchange(value, 1, 0)
