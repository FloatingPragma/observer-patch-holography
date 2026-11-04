"""Independent quantum, geometric and long-wavelength texture controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from . import energy, packet


def rotation(q):
    w, x, y, z = np.asarray(q)/np.linalg.norm(q)
    return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])


def test_exact_matrix_identity_matches_a_nonproduct_quantum_ground():
    sigma = (np.array([[0., 1], [1, 0]]), np.array([[0., -1j], [1j, 0]]), np.diag([1., -1]))
    H = (sum(np.kron(S, S) for S in sigma)-np.eye(4))/4
    ground = np.array([0, 1., -1., 0])/np.sqrt(2)
    qs = [np.array([1., 2., 3., 4.]), np.array([2., -1., 4., 3.])]
    unitaries = []
    for q in qs:
        q /= np.linalg.norm(q)
        unitaries.append(q[0]*np.eye(2)-1j*sum(q[j+1]*sigma[j] for j in range(3)))
    rotated = np.kron(*unitaries)@ground
    actual = (rotated.conj()@H@rotated-ground.conj()@H@ground).real
    relative = rotation(qs[0]).T@rotation(qs[1])
    expected = np.trace(-np.eye(3)/4@(relative-np.eye(3)))
    assert abs(actual-expected) < 1e-12
    assert actual > .1


def test_arbitrary_rotation_lower_bound_and_local_frame_residual():
    rng = np.random.default_rng(121)
    for _ in range(30):
        R = rotation(rng.normal(size=4))
        a, b = .1, -.2
        kappa = -(a+b)
        actual = np.trace(np.diag([a, b, b])@(R-np.eye(3)))
        direction = kappa*(1-R[0, 0])
        residual = -b*(1+R[0, 0]-R[1, 1]-R[2, 2])
        assert abs(actual-direction-residual) < 1e-12
        assert residual >= -1e-12
        assert actual >= direction-1e-12
    # A nonuniform rotation of the local frame around e1 has no direction
    # change but does have energy. Local gauge quotienting would erase it.
    R = rotation([1, 1, 0, 0])
    assert R[0, 0] == 1
    assert np.trace(np.diag([a, b, b])@(R-np.eye(3))) > .3


def coplanar_energy(scale):
    coordinates = np.arange(-scale-1, scale+2)/scale
    p = np.maximum(0, 1-coordinates**2)**4
    theta = np.einsum('i,j,k->ijk', p, p, p)
    return sum(np.sum(1-np.cos(np.diff(theta, axis=j))) for j in range(3))/scale


def test_finite_amplitude_texture_has_the_predicted_static_scaling():
    from sympy import Rational
    gradient = float(Rational(energy.exact_identities()['coplanar_profile_gradient_integral']))
    limit = gradient/2
    values = [coplanar_energy(L) for L in (4, 8, 16, 32)]
    errors = [abs(value-limit) for value in values]
    assert errors[-1] < errors[0]/30
    assert errors[-1] < .01
    assert limit > .1


def radial_direction_energy(radius):
    axis = np.arange(-radius, radius+1, dtype=float)
    points = np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), axis=-1)
    length = np.linalg.norm(points, axis=-1)
    direction = np.divide(points, length[..., None], out=np.zeros_like(points), where=length[..., None] > 0)
    value = 0.
    for j in range(3):
        left, right = [slice(None)]*3, [slice(None)]*3
        left[j], right[j] = slice(None, -1), slice(1, None)
        left, right = tuple(left), tuple(right)
        inside = (length[left] >= 1) & (length[right] >= 1) & (length[left] <= radius) & (length[right] <= radius)
        increment = 1-np.sum(direction[left]*direction[right], axis=-1)
        value += np.sum(increment[inside])
    return value/radius


def test_radial_hedgehog_direction_energy_grows_linearly_with_radius():
    values = [radial_direction_energy(L) for L in (8, 16, 32)]
    assert all(x < y for x, y in zip(values, values[1:]))
    assert 11 < values[-1] < 4*np.pi
    assert abs(values[-1]-4*np.pi) < abs(values[0]-4*np.pi)/3


def test_symbolic_identity_topological_degree_and_rational_error():
    value = energy.certificate()
    assert value['identities']['compact_Hopf_lift_absolute_degree'] == 1
    assert energy.static_error(1, 128, 1, 2**48) < F(1, 10**6)
    assert value['local_unitary_orbit_only']
    assert not value['energy_minimized_over_all_states']
    assert not value['finite_mass_topological_particle_proved']


@pytest.mark.parametrize('args', [(1., 1, 1, 1), (1, 0, 1, 1), (1, 1, F(1, 2), 1), (True, 1, 1, 1)])
def test_static_bound_rejects_invalid_domain(args):
    with pytest.raises(ValueError):
        energy.static_error(*args)


def test_scope_promotion_and_stale_import_are_rejected(monkeypatch):
    value = deepcopy(packet.report())
    value['energy']['energy_minimized_over_all_states'] = True
    with pytest.raises(ValueError, match='proof, scope or source'):
        packet.validate(value)
    pins = packet.bindings()
    pins['model/spin_texture_energy/PROOF.md'] = '0'*64
    monkeypatch.setattr(packet, 'bindings', lambda: pins)
    with pytest.raises(ValueError, match='after import'):
        packet.report()
