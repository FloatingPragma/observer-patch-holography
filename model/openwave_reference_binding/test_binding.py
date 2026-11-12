"""Independent string contraction, spectral and domain controls."""
from itertools import product
import numpy as np
import pytest
import sympy as s
from model.coherent_defect_wave.chiral_joint.oriented_links.model import endpoint_current, labels
from model.coherent_defect_wave.physical_attachment.chiral import representation
from . import binding, strings


def test_dense_three_edge_string_gauss_and_endpoint_orientation():
    path = strings.shortest_paths(7)[-1]
    state = strings.string_state('L', path)
    active = [sorted({key[1][j] for key in state}) for j in range(3)]
    dims = [2, *map(len, active), 2]
    rows = list(product(*(range(d) for d in dims))); lookup = {key:i for i,key in enumerate(rows)}
    vector = s.zeros(len(rows), 1)
    for (reference, links, matter), value in state.items():
        key = (reference, *(active[j].index(x) for j,x in enumerate(links)), matter)
        vector[lookup[key]] = value
    assert (vector.H*vector)[0] == 1 and len(rows) == 256
    identities = [s.eye(d) for d in dims]
    nodes = [0]+[edge[2] for edge in path]
    for current in (8, 9, 10, 11):
        t = representation()[2][current].extract((10, 11), (10, 11))
        for node in nodes:
            operator = s.zeros(256)
            factors = []
            if node == 0: factors.append((0, t))
            if node == 7: factors.append((4, -t.T))
            for j, (_, start, finish, forward) in enumerate(path):
                if node in (start, finish):
                    endpoint = (0 if forward else 1) if node == start else (1 if forward else 0)
                    factors.append((j+1, endpoint_current(current, endpoint).extract(active[j], active[j])))
            for axis, matrix in factors:
                matrices = list(identities); matrices[axis] = matrix
                operator += s.kronecker_product(*matrices)
            assert operator*vector == s.zeros(256, 1)
    assert any(not edge[3] for edge in path)
    assert strings.gauss_image('L', path, state, 7, 11, wrong_matter_sign=True)


def test_exact_doublet_minimum_and_singlet_intervals():
    doublet = binding.family_certificate('L'); singlet = binding.family_certificate('ec')
    assert doublet['shortest_string_saturates_global_constrained_minimum']
    assert not singlet['shortest_string_saturates_global_constrained_minimum']
    for row in doublet['rows']:
        assert list(map(s.Rational, row['binding_cost_relative_to_origin_interval'])) == [row['distance']]*2
    for row in singlet['rows'][1:]:
        low, high = map(s.Rational, row['binding_cost_relative_to_origin_interval'])
        assert 0 < low < row['distance'] < high
    assert len(doublet['rows']) == len(singlet['rows']) == 8


def test_exact_irrep_uses_actual_local_H_and_positive_metric():
    for family in ('ec', 'L'):
        row = binding.matter_irrep_certificate(family)
        assert row['actual_local_H_intertwining_exact'] and row['exact_Gram_positive']
        assert row['isolated_selected_root_count'] == 1
        assert row['native_irrep_dimension'] == (1 if family == 'ec' else 2)


def test_compressed_energy_is_not_a_stationary_full_H_eigenvalue():
    # A full-domain counterexample to replacing a compression by an invariant band.
    h = s.Matrix([[0, s.Rational(1, 7)], [s.Rational(1, 7), 1]])
    p = s.diag(1, 0); hopping = h-s.diag(0, 1)
    assert p*hopping*p == s.zeros(2) and h*p != p*h
    assert h.det() != 0 and min(np.linalg.eigvalsh(np.asarray(h, float))) < 0


def test_reference_cost_and_ground_origin_cannot_remove_binding():
    row = binding.family_certificate('L')['rows'][-1]
    distance = s.Integer(row['distance'])
    local, ground, cost = s.symbols('local ground cost', real=True)
    at_origin = local+cost-ground
    at_remote = local+cost+distance-ground
    assert s.simplify(at_remote-at_origin) == distance
    assert s.diff(at_remote-at_origin, cost) == s.diff(at_remote-at_origin, ground) == 0


def test_cloud_is_retained_and_path_domain_rejects_bad_cells():
    value = binding.certificate()
    for row in value['screened_cloud_counterexample']['records']:
        assert s.Rational(row['combined_origin_probe_and_remote_hop_squared_norm']) > 0
        assert row['single_unscreened_charge_projectors_all_annihilate_state']
        assert not row['flux_path_from_origin_to_remote_hole_exists']
    assert not value['scope']['screened_cloud_minimum_or_string_breaking_computed']
    with pytest.raises(ValueError): strings.shortest_paths(8)
    with pytest.raises(ValueError): strings.string_state('Q', ())
