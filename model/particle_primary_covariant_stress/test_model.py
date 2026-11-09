"""Independent allocation, operator continuity and symmetry controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_physical_action_and_all_cube_rotations():
    for side in (2, 3):
        result = m.allocation_control(side)
        assert result['proper_rotations_checked'] == 24
        assert result['term_counts'] == {'onsite_matter': side**3,
            'electric': 3*side**2*(side-1), 'bond': 3*side**2*(side-1),
            'plaquette': 3*side*(side-1)**2}


def test_bulk_and_corner_supports_use_four_vertices_and_actual_three_owners():
    member = m.contract.member(1)
    bulk = m.cell_energy(member, (2, 2, 2))
    assert sum(row['id'][0] == 'plaquette' for row in bulk['terms']) == 12
    assert all(len(row['owner_support']) == 3 for row in bulk['terms']
               if row['id'][0] == 'plaquette')
    assert bulk['norm_upper'] == (180+F(1827, 2)*4+F(3, 7))*2
    corners = [m.cell_energy(member, x)['norm_upper']
               for x in product((0, 7), repeat=3)]
    assert len(set(corners)) == 1
    with pytest.raises(ValueError): m.cell_energy(member, (8, 0, 0))


def test_every_small_region_difference_is_exactly_boundary_improvement():
    cells = list(product(range(2), repeat=3))
    for mask in range(1, 256):
        region = tuple(x for i, x in enumerate(cells) if mask & (1 << i))
        r = m.regional_improvement_control(2, region)
        assert r['energy_difference_coefficients'] == r['inward_boundary_flow_coefficients']
        if mask == 255:
            assert r['crossing_routes'] == 0 and r['energy_difference_coefficients'] == {}
    with pytest.raises(ValueError): m.regional_improvement_control(2, (cells[0], cells[0]))


def test_owner_allocation_really_is_not_rotation_covariant():
    old = m._old_allocations(2)
    R = s.diag(-1, -1, 1)
    x = (0, 0, 0); y = m.channels.rotate_cell(x, 2, R)
    rotated = {m._term(k, [m.channels.rotate_cell(v, 2, R) for v in vs]): w
               for (k, vs), w in old[x].items()}
    assert rotated != old[y]


def test_noncommuting_operator_improvement_and_Ward_identity():
    old = m._old_allocations(2); new = m._allocations(2)
    terms = sorted(set().union(*(set(v) for v in new.values())))
    X = s.Matrix([[0, 1], [1, 0]]); Y = s.Matrix([[0, -s.I], [s.I, 0]])
    Z = s.diag(1, -1)
    ops = {t: (i+1)*[X, Y, Z][i % 3] for i, t in enumerate(terms)}
    def assigned(table):
        return {x: sum((w*ops[t] for t, w in weights.items()), s.zeros(2))
                for x, weights in table.items()}
    a = assigned(old); b = assigned(new)
    H = sum(a.values(), s.zeros(2))
    assert H == sum(b.values(), s.zeros(2))
    for x in a:
        delta = b[x]-a[x]
        Jdiff = sum((s.I*(b[x]*b[y]-b[y]*b[x]-a[x]*a[y]+a[y]*a[x])
                     for y in a), s.zeros(2))
        assert s.simplify(Jdiff-s.I*(delta*H-H*delta)) == s.zeros(2)


def test_C_odd_tangent_has_zero_linear_entropy_and_energy_but_not_zero_cost():
    rho = np.diag([.6, .4]); C = np.diag([1, -1])
    X = np.array([[0., 1.], [1., 0.]])
    K = -np.diag(np.log(np.diag(rho))); energy = np.diag([2., 5.])
    assert np.array_equal(C @ X @ C, -X)
    assert np.trace(X @ K) == np.trace(X @ energy) == 0
    perturbed = rho+.1*X
    values, vectors = np.linalg.eigh(perturbed)
    log = (vectors*np.log(values)) @ vectors.T
    assert np.trace(perturbed @ (log+K)) > .01


def test_complete_covariant_instrument_retains_entangled_reference():
    C = np.diag([1, -1]); H = np.diag([1., 3.])
    rho = np.array([[.6, .1j], [-.1j, .4]])
    ev, vec = np.linalg.eigh(rho)
    root = (vec*np.sqrt(ev))@vec.conj().T
    psi = root.reshape(-1); paired = np.kron(C, C)@psi
    effect = (np.eye(2)+H/4)/2
    effects = [effect, np.eye(2)-effect]
    sums = [np.zeros((4, 4), complex), np.zeros((4, 4), complex)]
    for effect in effects:
        K = np.diag(np.sqrt(np.diag(effect))); gate = np.kron(K, np.eye(2))
        u = gate@psi; v = gate@paired
        assert abs(np.vdot(u, u)-np.vdot(v, v)) < 1e-14
        assert np.linalg.norm(v-np.kron(C, C)@u) < 1e-14
        sums[0] += np.outer(u, u.conj()); sums[1] += np.outer(v, v.conj())
    assert abs(np.trace(sums[0])-1) < 1e-14
    assert np.linalg.norm(sums[0]-sums[1]) > .01


def test_member_binding_rejects_result_and_clock_mutation():
    member = m.contract.member(0); receipt = m.bind_member(member)
    m.verify_bound(member, receipt)
    bad = deepcopy(receipt); bad['result']['proper_metric_or_Lorentz_stress_derived'] = True
    with pytest.raises(ValueError): m.verify_bound(member, bad)
    badmember = deepcopy(member); badmember['source']['spacing'] = '2'
    with pytest.raises(ValueError): m.bind_member(badmember)
