"""Cut covariance, exact source representations and retained-record tests."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import model as m


def test_induced_graph_cut_and_complement_count_every_interior_edge_once():
    for side in (2, 3):
        cells = list(product(range(side), repeat=3))
        regions = [cells[:1], cells[:3], cells, [x for x in cells if sum(x) % 2]]
        total = 3*side**2*(side-1)
        for region in regions:
            a = m.cut(side, region)
            other = [x for x in cells if x not in region]
            complement_internal = len(m.cut(side, other)['interior_links']) if other else 0
            assert len(a['interior_links'])+len(a['crossing_links'])+complement_internal == total
            assert sum(d for _, d in a['inside_boundary']) == len(a['crossing_links'])
            assert a['CAR_modes'] == 60*len(region)


def test_all_cube_rotations_include_link_reversal_and_cut_factor_covariance():
    for side in (2, 3):
        r = m.cube_control(side)
        assert r['proper_rotations_checked'] == 24
        assert r['induced_link_cut_covariant']
    r = m.cut(2, ((0, 0, 0), (1, 0, 0)))
    assert len(r['interior_links']) == 1
    assert len(r['crossing_links']) == 4
    assert r['inside_boundary'] == (((0, 0, 0), 2), ((1, 0, 0), 2))


def test_bulk_cube_edge_entropy_bound_scales_with_boundary_not_volume():
    member = m.contract.member(1)
    region = tuple(product(range(2, 5), repeat=3))
    b = m.boundary_entropy_bound(member, region)
    assert len(b['crossing_links']) == 6*3**2
    assert len(b['interior_links']) == 3*3**2*2
    assert b['edge_entropy_upper_coefficient_of_log_dmax'] == 108
    assert b['lattice_cut_face_area'] == F(27, 2)
    assert b['edge_entropy_upper_per_lattice_area_coefficient_of_log_dmax'] == 8
    assert b['largest_native_irrep_dimension'] == 6


def test_no_cut_edge_term_for_entire_physical_box():
    member = m.contract.member(0)
    b = m.boundary_entropy_bound(member, tuple(product(range(2), repeat=3)))
    assert b['edge_entropy_upper_coefficient_of_log_dmax'] == 0
    assert b['boundary_irrep_dimension_upper_factored'] == {'base': 6, 'exponent': 0}
    assert b['inside_boundary_charge_bounds'] == []


def test_regional_energy_difference_contains_only_crossing_terms():
    member = m.contract.member(0)
    cells = list(product(range(2), repeat=3))
    for mask in range(1, 256):
        region = tuple(x for i,x in enumerate(cells) if mask & (1 << i))
        r = m.regional_energy_difference(member, region)
        assert r['boundary_difference_norm_upper'] <= r['uniform_cut_boundary_difference_norm_upper']
        assert all(0 < value < 1 for value in r['crossing_term_coefficients'].values())
        if mask == 255:
            assert r['boundary_difference_norm_upper'] == 0
            assert r['crossing_term_coefficients'] == {}


def test_actual_native_character_schmidt_values_for_every_representation():
    blocks = m.stress.channels.native_number_control()['blocks']
    assert [(x['Q6'], x['internal_dimension']) for x in blocks] == [(-4,3),(-3,2),(1,6),(2,3),(6,1)]
    for block in blocks:
        d = block['internal_dimension']
        # Orthonormal Haar coefficient basis: trace(A B) reverses the two indices.
        coefficients = np.zeros((d*d, d*d))
        for i, j in product(range(d), repeat=2):
            coefficients[i*d+j, j*d+i] = 1/d
        reduced = coefficients @ coefficients.T
        assert np.allclose(reduced, np.eye(d*d)/(d*d), atol=1e-15)
        assert abs(np.trace(reduced)-1) < 1e-14


def _entropy(rho):
    values = np.linalg.eigvalsh(rho)
    return -sum(x*np.log(x) for x in values if x > 1e-15)


def test_center_record_entropy_and_edge_term_retain_all_probabilities():
    weights = [.2, .3, .5]; dims = [1, 2, 6]
    multiplicity_states = [np.diag([.4,.6]), np.array([[.6,.1j],[-.1j,.4]]), np.ones((1,1))]
    blocks = [p*np.kron(t, np.eye(d)/d) for p,t,d in zip(weights,multiplicity_states,dims)]
    size = sum(len(b) for b in blocks); rho = np.zeros((size,size), complex)
    record = np.zeros((3*size,3*size), complex); offset = 0
    for k, block in enumerate(blocks):
        n = len(block); rho[offset:offset+n, offset:offset+n] = block
        start = k*size+offset; record[start:start+n, start:start+n] = block; offset += n
    alg = -sum(p*np.log(p) for p in weights)+sum(p*_entropy(t) for p,t in zip(weights,multiplicity_states))
    edge = sum(p*np.log(d) for p,d in zip(weights,dims))
    assert abs(_entropy(rho)-alg-edge) < 1e-13
    assert abs(_entropy(record)-_entropy(rho)) < 1e-13
    assert edge > 0


def test_retained_center_measurement_can_disturb_global_reference_coherence():
    # Two matching sectors have a coherent physical singlet superposition.
    psi = np.array([np.sqrt(.3), 0., 0., np.sqrt(.7)])
    full = np.outer(psi, psi)
    projectors = [np.diag([1.,1.,0.,0.]), np.diag([0.,0.,1.,1.])]
    post = sum(P@full@P for P in projectors)
    assert np.linalg.norm(post-full) > .6
    assert abs(_entropy(full)) < 1e-13
    assert abs(_entropy(post)+.3*np.log(.3)+.7*np.log(.7)) < 1e-13
    # Source marginal is unchanged while the full state loses coherence.
    assert np.allclose(np.trace(full.reshape(2,2,2,2), axis1=1,axis2=3),
                       np.trace(post.reshape(2,2,2,2), axis1=1,axis2=3))


def test_region_validation_and_member_custody():
    with pytest.raises(ValueError): m.cut(2, ())
    with pytest.raises(ValueError): m.cut(2, ((0,0,0), (0,0,0)))
    with pytest.raises(ValueError): m.cut(2, ((2,0,0),))
    member = m.contract.member(0); receipt = m.bind_member(member)
    m.verify_bound(member, receipt)
    floor = receipt['result']['physical_density_floor']
    assert floor['operator'] == 'P_G' and floor['base'] == 2 and floor['negative_exponent'] > 0
    assert receipt['result']['edge_entropy_positive_lower_factored'] == {
        'multiplier_of_physical_density_floor': 4, 'natural_log_argument': 108}
    bad = deepcopy(receipt); bad['result']['Einstein_equation_derived'] = True
    with pytest.raises(ValueError): m.verify_bound(member, bad)
    badmember = deepcopy(member); badmember['source']['spacing'] = '2'
    with pytest.raises(ValueError): m.bind_member(badmember)
