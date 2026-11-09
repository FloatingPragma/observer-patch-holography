"""Independent finite-source controls for the perturbative coefficient audit."""
from fractions import Fraction as F
from itertools import combinations, product
import json

import numpy as np
import pytest

from model.high_field_particle import algebra
from model.native_ordered_ground import source


def tensor_label(label):
    pauli = {'I': np.eye(2, dtype=complex),
             'X': np.array([[0, 1], [1, 0]], dtype=complex),
             'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
             'Z': np.diag([1, -1]).astype(complex)}
    out = np.ones((1, 1), dtype=complex)
    for axis in reversed(label):  # The owning source uses least-significant-site bits.
        out = np.kron(out, pauli[axis])
    return out


@pytest.mark.parametrize('sites,edges', [
    (1, ()), (2, ((0, 1),)), (3, ((0, 1), (1, 2))),
    (4, ((0, 1), (0, 3), (1, 2), (2, 3))),
    (4, ((0, 1), (0, 1), (1, 2)))])
def test_full_source_equals_independent_pauli_tensor(sites, edges):
    j, b = F(2, 7), F(3, 5)
    actual, scale = algebra.transformed_source(sites, edges, j, b)
    target = np.zeros_like(actual, dtype=complex)
    for x in range(sites):
        label = ['I'] * sites
        label[x] = 'Z'
        target -= float(b * scale / 2) * tensor_label(label)
    for x, y in edges:
        target -= float(j * scale / 4) * np.eye(1 << sites)
        for axis, sign in [('X', 1), ('Y', -1), ('Z', -1)]:
            label = ['I'] * sites
            label[x] = label[y] = axis
            target += sign * float(j * scale / 4) * tensor_label(label)
    # All tensor entries in this test are exactly represented binary integers/halves.
    assert np.array_equal(actual, target)
    assert algebra.frame_audit(sites, edges)['full_source_staggered_charge_conserved']


def test_source_change_is_rejected(monkeypatch):
    original = source.matrix

    def wrong_source(*args, **kwargs):
        matrix, scale = original(*args, **kwargs)
        changed = matrix.copy()
        changed[0, 0] += 1
        return changed, scale

    monkeypatch.setattr(source, 'matrix', wrong_source)
    with pytest.raises(ArithmeticError):
        algebra.frame_audit(2, ((0, 1),))


@pytest.mark.parametrize('sites,edges', [
    (2, ((0, 1),)), (3, ((0, 1), (1, 2))),
    (4, ((0, 1), (0, 3), (1, 2), (2, 3))),
    (5, ((0, 1), (0, 2), (0, 3), (0, 4))),
    (6, algebra.periodic_edges((6,))), (8, source.cube()[0])])
def test_all_virtual_triples_from_independent_edge_incidence(sites, edges):
    result = algebra.expansion(sites, edges)
    multiplicities = {edge: edges.count(edge) for edge in set(edges)}
    expected = [[F(0) for _ in range(sites)] for _ in range(sites)]
    triples = set()
    for subset in combinations(range(sites), 3):
        amplitude = [F(multiplicities.get(tuple(y for y in subset if y != x), 0), 2)
                     if x in subset else F(0) for x in range(sites)]
        if any(amplitude):
            triples.add(sum(1 << x for x in subset))
        for x in range(sites):
            for y in range(sites):
                expected[x][y] -= amplitude[x] * amplitude[y] / 2
    assert result['absolute_second'] == tuple(map(tuple, expected))
    assert result['three_flip_intermediate_masks'] == tuple(sorted(triples))
    assert result['all_complementary_rows_checked'] == (1 << sites) - sites
    assert result['full_order_two_intertwiner_exact']


def test_dimer_exact_spectrum_expansion_and_ground_subtraction():
    result = algebra.expansion(2, ((0, 1),))
    # The exact even block has energies 1-e/2 +/- sqrt(1+e^2/4);
    # both odd states have absolute energy 1 after the field scalar shift.
    assert result['ground_first'] == -F(1, 2)
    assert result['ground_second'] == -F(1, 8)
    assert result['absolute_first'] == ((F(0), F(0)), (F(0), F(0)))
    assert result['absolute_second'] == ((F(0), F(0)), (F(0), F(0)))
    assert result['relative_second'] == ((F(1, 8), F(0)), (F(0), F(1, 8)))
    assert result['three_flip_intermediate_masks'] == ()
    # An exact rational Taylor control avoids a fitted or floating eigenvalue.
    e = F(1, 20)
    x = e * e / 4
    lower, upper = 1 + x / 2 - x * x / 8, 1 + x / 2
    assert lower > 0 and lower * lower <= 1 + x <= upper * upper
    # The positive correction would be lost by omitting ground subtraction.
    assert result['relative_second'] != result['absolute_second']


def test_irregular_graph_and_parallel_amplitudes():
    edges = ((0, 1), (0, 1), (1, 2), (1, 2), (1, 2))
    result = algebra.expansion(3, edges)
    assert result['degree'] == (2, 5, 3)
    assert result['squared_multiplicity_degree'] == (4, 13, 9)
    assert result['ground_second'] == -F(13, 8)
    assert result['relative_first'] == ((F(1), F(0), F(0)),
                                        (F(0), F(5, 2), F(0)),
                                        (F(0), F(0), F(3, 2)))
    assert result['relative_second'][0][2] == -F(3, 4)
    assert result['relative_second'][1][1] == F(13, 8)
    assert result['ground_second'] != -F(len(edges), 8)


def test_charge_parity_and_number_have_distinct_source_meanings():
    edges = ((0, 1), (1, 2))
    h, scale = algebra.transformed_source(3, edges)
    number = np.array([mask.bit_count() for mask in range(8)])
    assert np.any(h * (number[:, None] - number[None, :]))
    assert h[3, 0] == scale // 2  # Pair creation is present in the full source.
    result = algebra.frame_audit(3, edges)
    assert result['one_flip_charge_relative_to_product_vacuum'] == [-1, 1, -1]
    assert result['one_flip_parity_relative_to_product_vacuum'] == -1
    second = algebra.expansion(3, edges)['relative_second']
    assert second[0][2] == -F(1, 8)
    assert second[0][1] == second[1][2] == 0


def test_second_wave_contains_five_flip_rows():
    result = algebra.expansion(6, algebra.periodic_edges((6,)))
    assert result['second_wave_maximum_flip_count'] == 5


def test_cycle_fourier_coefficients_exact_without_eigenvalue_fit():
    result = algebra.expansion(4, algebra.periodic_edges((4,)))
    matrix = result['relative_second']
    roots = [1, 1j, -1, -1j]
    for wave, eigenvalue in [(0, F(0)), (1, F(1, 2)), (2, F(0)), (3, F(1, 2))]:
        vector = [roots[(wave * x) % 4] for x in range(4)]
        assert [sum(matrix[x][y] * vector[y] for y in range(4)) for x in range(4)] == [eigenvalue * value for value in vector]


def test_declared_periodic_cube_is_not_bulk_unit_edge_cube():
    edges = algebra.periodic_edges((2, 2, 2))
    assert len(edges) == 24 and len(set(edges)) == 12
    result = algebra.expansion(8, edges)
    assert result['degree'] == (6,) * 8
    assert result['squared_multiplicity_degree'] == (12,) * 8
    coordinates = tuple(product(range(2), repeat=3))
    for wave in coordinates:
        vector = [(-1) ** sum(k * x for k, x in zip(wave, point)) for point in coordinates]
        adjacency_value = 2 * sum((-1) ** k for k in wave)
        second_value = F(24 - adjacency_value ** 2, 8)
        matrix = result['relative_second']
        assert [sum(matrix[x][y] * vector[y] for y in range(8)) for x in range(8)] == [second_value * value for value in vector]
        assert second_value != F(12 - adjacency_value ** 2, 8)
    assert algebra.expansion(8, source.cube()[0])['degree'] == (3,) * 8


def test_bulk_path_support_hessian_and_nonconstant_symbol():
    result = algebra.cubic_paths()
    coefficients = result['second_order_coefficients']
    assert len(coefficients) == 19
    assert coefficients[(0, 0, 0)] == F(3, 4)
    assert coefficients[(2, 0, 0)] == -F(1, 8)
    assert coefficients[(1, 1, 0)] == -F(1, 4)
    assert result['at_zero_momentum'] == -3
    assert result['hessian_at_zero_momentum'] == ((F(3), F(0), F(0)),
                                                  (F(0), F(3), F(0)),
                                                  (F(0), F(0), F(3)))
    # k=(pi/2,pi/2,pi/2) gives phases in {+1,-1} on this even support.
    opposite = sum(value * (-1) ** (sum(offset) // 2)
                   for offset, value in coefficients.items())
    assert opposite == F(3, 2)
    assert opposite - result['at_zero_momentum'] == F(9, 2)


@pytest.mark.parametrize('sites,edges', [(0, ()), (9, ()),
                                        (3, ((0, 1), (1, 2), (0, 2))),
                                        (2, ((1, 0),)), (2, [(0, 1)])])
def test_source_domain_rejections(sites, edges):
    with pytest.raises(ValueError):
        algebra.expansion(sites, edges)


@pytest.mark.parametrize('shape', [(3,), (1, 2), (4, 4), (), [2, 2]])
def test_periodic_domain_rejections(shape):
    with pytest.raises(ValueError):
        algebra.periodic_edges(shape)


def test_actual_cube_spectrum_diagnostic_is_separate_from_exact_proof():
    edges = source.cube()[0]
    epsilon = F(1, 1000)
    h, scale = source.matrix(8, edges, epsilon, F(1))
    spectrum = np.linalg.eigvalsh(h.astype(float) / scale)
    actual_gaps = spectrum[1:9] - spectrum[0]
    expansion = algebra.expansion(8, edges)
    effective = (np.eye(8) + float(epsilon) * np.array(expansion['relative_first'], dtype=float)
                 + float(epsilon ** 2) * np.array(expansion['relative_second'], dtype=float))
    # A finite-source numerical control, not a certified remainder or band proof.
    assert np.max(np.abs(actual_gaps - np.linalg.eigvalsh(effective))) < 1e-8


def test_certificate_exact_json_scope_and_independent_returns():
    first = algebra.certificate()
    assert json.loads(json.dumps(first)) == first
    assert first['first_order_hopping_zero']
    assert first['infinite_lattice_band_uses_external_Yarotsky_smallness_theorem']
    for key in ('explicit_Yarotsky_threshold_supplied',
                'finite_coupling_remainder_computed_by_this_certificate',
                'joint_n_inverse_fifth_field_trajectory_used',
                'finite_periodic_double_edges_identified_with_bulk_unit_edges',
                'native_ground_preparation_executed_here',
                'measured_particle_inputs_used', 'physical_particle_assignment',
                'physical_mass_prediction', 'relativistic_dispersion_claimed'):
        assert first[key] is False
    first['finite_source_controls'][0]['frame']['exact_global_Hadamard_identity'] = False
    assert algebra.certificate()['finite_source_controls'][0]['frame']['exact_global_Hadamard_identity'] is True
