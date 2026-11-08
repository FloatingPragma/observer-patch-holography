"""Independent finite-source Bloch transport and exact response-channel checks."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as sy

from model.source_family_compatibility import composite
from model.source_family_compatibility import momentum_response as source
from model.source_family_compatibility import nodal_channels as n


def test_actual_periodic_source_lifts_match_the_nodal_bloch_frame():
    points = tuple(product(range(4), repeat=3))
    index = {x: i for i, x in enumerate(points)}
    hopping = np.zeros((64, 64), dtype=np.int64)
    frame = np.zeros((64, 8), dtype=complex)
    for x in points:
        frame[index[x], n.INDEX[tuple(t % 2 for t in x)]] = 1j ** sum(x)
        for j in range(3):
            y = list(x)
            y[j] = (y[j] + 1) % 4
            i, k = index[x], index[tuple(y)]
            hopping[i, k] += source.signed_link(x, j)
            hopping[k, i] += source.signed_link(x, j)
    np.testing.assert_array_equal(frame.conj().T @ frame, 8 * np.eye(8))
    np.testing.assert_array_equal(hopping @ frame, np.zeros_like(frame))
    for j, expected in enumerate(n.magnetic_generators()):
        lift = np.zeros((64, 64), dtype=np.int64)
        for x in points:
            y = list(x)
            y[j] = (y[j] - 1) % 4
            lift[index[x], index[tuple(y)]] = source.magnetic_phase(x, j)
        np.testing.assert_array_equal(lift @ hopping, hopping @ lift)
        np.testing.assert_array_equal(lift @ frame, frame @ (-1j * expected))
    for (perm, signs), expected in zip(source.cubic_rotations(), n.cubic_lifts(), strict=True):
        lift = np.zeros((64, 64), dtype=np.int64)
        for x in points:
            y = tuple(t % 4 for t in source.rotate(x, perm, signs))
            lift[index[y], index[x]] = source.rotation_phase(x, perm)
        np.testing.assert_array_equal(lift @ hopping, hopping @ lift)
        np.testing.assert_array_equal(lift @ frame, frame @ expected)


def test_commutant_dimension_by_independent_exact_linear_elimination():
    identity = sy.eye(8)
    constraints = sy.Matrix.vstack(*[
        sy.kronecker_product(identity, sy.Matrix(g))
        - sy.kronecker_product(sy.Matrix(g).T, identity)
        for g in n.magnetic_generators() + n.cubic_lifts()
    ])
    kernel = constraints.nullspace()
    assert len(kernel) == 2
    assert len(n.commutant_basis(n.magnetic_generators())) == 8
    for b in n.commutant_basis(n.magnetic_generators() + n.cubic_lifts()):
        assert constraints * sy.Matrix(b.T.reshape(-1)) == sy.zeros(1728, 1)


def test_reciprocal_cell_transport_is_necessary_at_boundary_node():
    correct = n.cubic_lifts()
    wrong = []
    for perm, signs in source.cubic_rotations():
        matrix = np.zeros((8, 8), dtype=np.int64)
        for x in n.PARITIES:
            y = tuple(t % 2 for t in source.rotate(x, perm, signs))
            matrix[n.INDEX[y], n.INDEX[x]] = source.rotation_phase(x, perm)
        wrong.append(matrix)
    assert any(not np.array_equal(a, b) for a, b in zip(correct, wrong, strict=True))
    # A dimension count alone misses this error: both commutants have dimension two.
    wrong_basis = n.commutant_basis(n.magnetic_generators() + tuple(wrong))
    assert len(wrong_basis) == 2
    gamma = tuple(np.asarray(g) for g in n.gamma_matrices())
    volume = gamma[0] @ gamma[1] @ gamma[2]
    assert any(not np.array_equal(8 * b, np.trace(b) * np.eye(8)
                                 - np.trace(volume @ b) * volume) for b in wrong_basis)


@pytest.mark.parametrize('bad', [(), (np.eye(2) * .5,), (np.ones((2, 2)),)])
def test_exact_signed_orbit_input_contract(bad):
    with pytest.raises(ValueError):
        n.commutant_basis(bad)


def test_actual_full_ground_retarded_moments_obey_energy_reflection():
    data = composite.source((2, 1, 1), F(3, 2), F(4))
    numerator, denominator, _ = composite.ground(data, F(-9, 2))
    h = data['H'].astype(object)
    operators = [a.astype(object) for a in data['a']]
    original = [a.copy() for a in operators]
    grading = np.diag([1, 1, -1, -1]).astype(object)
    for order in range(7):
        moment = np.empty((4, 4), dtype=object)
        for i, a in enumerate(operators):
            for j, b in enumerate(original):
                value = a @ b.T + b.T @ a
                moment[i, j] = sum(int(x) * int(y) for x, y in
                                   zip(numerator.T.flat, value.flat, strict=True))
        np.testing.assert_array_equal(grading @ moment @ grading, (-1) ** order * moment)
        if order == 0:
            np.testing.assert_array_equal(moment, denominator * np.eye(4, dtype=int))
        operators = [a @ h - h @ a for a in operators]


def test_original_valley_transport_and_source_scope():
    row = n.certificate()
    assert row['nodal_little_group_commutant_dimension'] == 2
    assert row['including_hidden_spin_projector_ranks'] == [8, 8]
    assert [r['kinetic_degree'] for r in row['source_valley_channels']].count(1) == 4
    assert row['independent_after_source_energy_reflection'] == 1
    assert row['original_nodal_lorentzian_scores_equal']
    assert not row['uniform_growing_packet_scaling_supplied']
    assert not row['physical_spin_or_species_identification']
    assert not row['physical_mass_prediction']
    row['source_valley_channels'].clear()
    assert len(n.certificate()['source_valley_channels']) == 8


def test_scalar_memory_score_against_full_block_resolvent_with_real_memory():
    gamma = tuple(np.asarray(g, complex) for g in n.gamma_matrices())
    channel = 1j * gamma[0] @ gamma[1] @ gamma[2]
    identity = np.eye(8)
    zero = np.zeros((8, 8))
    alpha, coupling, hidden = F(4), F(3, 2), F(5, 4)
    clock = np.block([[zero, float(alpha) * identity, zero],
                      [float(alpha) * identity, zero, float(coupling) * identity],
                      [zero, float(coupling) * identity, float(hidden) * channel]])
    for eta in (F(1, 8), F(1, 2), F(2)):
        a = -coupling**2 * hidden / (eta**2 + hidden**2)
        b = coupling**2 * eta / (eta**2 + hidden**2)
        expected = n.lorentzian_score(alpha, eta, a, b)
        g = np.linalg.inv(1j * float(eta) * np.eye(24) - clock)[:8, :8]
        q = -float(eta) * (g - g.conj().T) / (2j)
        np.testing.assert_allclose(q, float(expected) * identity, atol=2e-15, rtol=0)
        assert expected == n.lorentzian_score(alpha, eta, -a, b)
        assert eta**2 / (alpha**2 + eta**2) <= expected <= 1
    with pytest.raises(ValueError):
        n.lorentzian_score(4, F(1, 2), 0, -1)
