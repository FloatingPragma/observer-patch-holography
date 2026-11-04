"""Independent same-source radial matrices and complete-Fock controls."""
from fractions import Fraction as F
from itertools import combinations

import mpmath
import numpy as np
import pytest
from scipy.linalg import eigh

from ..bounded_defect_index import finite
from ..full_fock_field.field import FieldGeometry, FullFock
from . import algebra, core


def chain_oracle(coupling=F(1, 256), mass=4):
    geometry = FieldGeometry((3, 1, 1), speed=1., interaction=float(coupling))
    carrier = FullFock(6)
    c = [carrier.annihilator(j).toarray() for j in range(6)]
    sigma = np.array([[0., 1.], [1., 0.]])
    h = geometry.single_h().toarray()+np.kron(np.diag([mass, 0, -mass]), sigma)
    H0 = sum(h[i, j]*c[i].conj().T@c[j] for i in range(6) for j in range(6))
    numbers = np.array([mask.bit_count() for mask in range(64)])
    H = carrier.hamiltonian(geometry).toarray()+np.diag(float(coupling)*(-numbers/2+3/4))
    H += sum((h[i, j]-geometry.single_h()[i, j])*c[i].conj().T@c[j]
             for i in range(6) for j in range(6))
    values, vectors = eigh(h)
    negative = values < -1e-9
    sea = vectors[:, negative]
    zero = np.kron(np.array([[-.5j/mass], [1.], [-.5j/mass]])/
                   np.sqrt(1+.5/mass**2), np.eye(2))
    columns = []
    for occupation in ((), (0,), (1,), (0, 1)):
        orbitals = np.column_stack([sea, zero[:, list(occupation)]]) if occupation else sea
        vector = np.zeros(64, complex)
        for occupied in combinations(range(6), orbitals.shape[1]):
            vector[sum(1 << j for j in occupied)] = np.linalg.det(orbitals[list(occupied), :])
        columns.append(vector)
    P = np.column_stack(columns)
    D0 = np.diag([(((mask >> 2) & 1)+((mask >> 3) & 1)-1)**2 for mask in range(64)])
    sectors = {}
    for number in range(7):
        rows = [mask for mask in range(64) if mask.bit_count() == number]
        energies, basis = eigh(H[np.ix_(rows, rows)])
        full = np.zeros((64, len(rows)), complex)
        full[rows] = basis
        sectors[number] = energies, full
    return h, H0, H, c, P, zero, D0, sectors


@pytest.mark.parametrize('mass', [4, 64])
def test_complete_radial_cube_zero_frame_is_quaternion_and_has_the_paired_sea(mass):
    h = finite.source_matrix(1, mass=mass).toarray()
    origin = [26, 27]
    outside = [j for j in range(54) if j not in origin]
    F0 = np.zeros((54, 2), complex)
    F0[origin] = np.eye(2)
    F0[outside] = -np.linalg.solve(h[np.ix_(outside, outside)], h[np.ix_(outside, origin)])
    gram = F0.conj().T@F0
    assert np.linalg.norm(gram-gram[0, 0]*np.eye(2)) < 1e-12
    zero = F0/np.sqrt(gram[0, 0].real)
    assert np.linalg.norm(h@zero) < 1e-12
    values, vectors = eigh(h)
    assert sum(abs(values) < 1e-10) == 2
    assert np.min(abs(values[abs(values) > 1e-10])) >= mass-2
    R = vectors[:, values < -1e-10]@vectors[:, values < -1e-10].conj().T
    weights, weighted = [], np.zeros((2, 2), complex)
    for site in range(27):
        W = zero[2*site:2*site+2]
        rho = np.trace(W.conj().T@W).real/2
        weights.append(rho)
        assert np.linalg.norm(W.conj().T@W-rho*np.eye(2)) < 1e-12
        weighted += W.conj().T@R[2*site:2*site+2, 2*site:2*site+2]@W
    assert abs(sum(weights)-1) < 1e-12
    ipr = np.dot(weights, weights)
    assert np.linalg.norm(weighted-(1-ipr)*np.eye(2)/2) < 1e-12
    assert np.linalg.norm(R[np.ix_(origin, origin)]-(1-weights[13])*np.eye(2)/2) < 1e-12
    bounds = core.core_bounds(mass=mass)
    assert weights[13] >= float(bounds['origin_weight_lower'])
    assert ipr >= float(bounds['ipr_lower'])
    # The calculation retains each complete spin matrix without replacing
    # it by a scalar spatial ansatz.
    assert zero.shape == (54, 2)


def test_actual_cube_source_symmetries_and_invalid_time_reversal_controls():
    value = algebra.certificate()
    assert value['source_symmetries']['one_body_Kramers_square'] == -1
    assert not value['source_symmetries']['plain_spin_time_reversal_commutes']
    assert not value['source_symmetries']['inversion_times_plain_spin_time_reversal_commutes']
    exact, _ = algebra.source_matrix()
    assert np.max(abs(np.array(exact, complex)-finite.source_matrix(1, mass=64).toarray())) < 1e-13


def test_scalar_mass_mutation_does_not_inherit_radial_particle_hole(monkeypatch):
    original = algebra.onsite
    def changed(point, mass):
        value = original(point, mass)
        return value+np.eye(2) if point == (0, 0, 0) else value
    monkeypatch.setattr(algebra, 'onsite', changed)
    with pytest.raises(ArithmeticError, match='symmetry'):
        algebra.symmetries()


@pytest.mark.parametrize('mass', [4, 64])
def test_full_filled_sea_contact_and_actual_origin_parity(mass):
    h, H0, H, _, P, zero, D0, _ = chain_oracle(F(1), mass)
    assert np.linalg.norm(P.conj().T@P-np.eye(4)) < 1e-12
    assert np.linalg.norm(h@zero) < 1e-12
    compressed = P.conj().T@(H-H0)@P
    weights = [np.linalg.norm(zero[2*x:2*x+2, 0])**2 for x in range(3)]
    ipr = np.dot(weights, weights)
    constant = compressed[0, 0]
    expected = np.diag([constant, constant-ipr/2, constant-ipr/2, constant])
    assert np.linalg.norm(compressed-expected) < 1e-11
    alpha = weights[1]
    parity = np.diag(P.conj().T@D0@P).real
    assert np.max(abs(parity-np.array([1+alpha**2, 1-alpha**2, 1-alpha**2, 1+alpha**2])/2)) < 1e-12
    # The local physical projector differs from an exact zero-orbital parity.
    assert 0 < parity[1] < .1 and parity[0] < 1


def test_exact_full_Fock_radial_Kramers_square_and_particle_hole():
    _, _, H, c, _, _, _, _ = chain_oracle(F(1, 64), 4)
    time = np.zeros((64, 64), complex)
    ph = np.zeros_like(time)
    for mask in range(64):
        occupied = [j for j in range(6) if mask & (1 << j)]
        images = [2*(2-j//2)+(j % 2 ^ 1) for j in occupied]
        sign = (-1)**sum(images[i] > images[j] for i in range(len(images)) for j in range(i+1, len(images)))
        for j in occupied:
            sign *= (-1)**((j//2-1) % 2)*(-1 if j % 2 == 0 else 1)
        time[sum(1 << j for j in images), mask] = sign
        vector = np.zeros(64, complex)
        vector[63] = 1
        for mode in reversed(occupied):
            coefficient = (-1)**((mode//2-1) % 2+mode % 2)
            vector = coefficient*c[mode ^ 1]@vector
        ph[:, mask] = vector
    expected_square = np.diag([(-1)**mask.bit_count() for mask in range(64)])
    assert np.array_equal(time@time.conjugate(), expected_square)
    assert np.array_equal(time@H.conjugate()@time.conj().T, H)
    assert np.array_equal(ph@H@ph.conj().T, H)


@pytest.mark.parametrize('mass,coupling', [(4, F(1, 64)), (4, F(1, 256)), (64, F(1, 65536))])
def test_all_number_ground_charged_poles_and_center_effect_bounds(mass, coupling):
    _, H0, H, c, P, zero, D0, sectors = chain_oracle(coupling, mass)
    bounds = core.core_bounds(coupling, mass=mass, dimensions=1)
    energies, ground = sectors[3][0], sectors[3][1][:, :2]
    eg = energies[0]
    assert abs(energies[1]-eg) < 1e-11
    assert all(v[0][0] >= eg-1e-11 for v in sectors.values())
    gap = sectors[4][0][0]-eg
    assert abs(sectors[2][0][0]-eg-gap) < 1e-11
    assert float(bounds['actual_charged_gap_lower']) <= gap <= float(bounds['actual_charged_gap_upper'])
    free_ground = P[:, 1:3]
    G = ground@ground.conj().T
    G0 = free_ground@free_ground.conj().T
    delta = float(bounds['projection_angle_bound'])
    assert np.linalg.norm(G-G0, 2) <= delta+1e-12
    ground_effect = np.trace(G@D0).real/2
    assert ground_effect <= float(bounds['center_parity_actual_ground_upper'])+1e-12
    assert abs(ground_effect-np.trace(G0@D0).real/2) <= delta+1e-12
    for number, plus in [(4, True), (2, False)]:
        vector = sectors[number][1][:, 0]
        C = np.outer(vector, vector.conjugate())
        charged_effect = np.trace(C@D0).real
        assert charged_effect >= float(bounds['center_parity_actual_charged_lower'])-1e-12
        assert charged_effect-ground_effect >= float(bounds['center_parity_actual_charged_minus_ground_lower'])-1e-12
        rows = np.array([(ground.conj().T@(a if plus else a.conj().T)@vector).ravel() for a in c])
        residue = rows@rows.conj().T/2 if plus else rows.conj()@rows.T/2
        assert np.linalg.norm(residue-zero@zero.conj().T/2, 2) <= 3*delta+1e-12
        low = zero.conj().T@residue@zero
        assert np.linalg.norm(low-np.trace(low)*np.eye(2)/2) < 1e-12
        assert min(np.diag(residue)[2:4].real) >= float(bounds['each_origin_spin_low_pole_weight_lower'])-1e-12
        assert min(sectors[number][0][1:])-eg >= float(bounds['other_charged_excitations_threshold_lower'])-1e-11
    full_free = np.linalg.eigvalsh(H0)
    shifted = np.linalg.eigvalsh(H)-full_free[0]
    compressed = np.linalg.eigvalsh(P.conj().T@(H-H0)@P)
    assert max(abs(shifted[:4]-compressed)) <= float(bounds['ordered_low_energy_error'])+1e-11


@pytest.mark.parametrize('bits', [2, 8])
def test_faithful_actual_source_Gibbs_bound_without_underflow(bits):
    _, _, H, _, _, _, _, _ = chain_oracle(F(1, 64), 4)
    values = np.linalg.eigvalsh(H)
    row = core.faithful_reference(bits, coupling=F(1, 64), mass=4, dimensions=1)
    ctx = mpmath.mp.clone()
    ctx.dps = 60
    beta = ctx.mpf(row['inverse_temperature'].numerator)/row['inverse_temperature'].denominator
    weights = [ctx.exp(-beta*ctx.mpf(float(value-values[0]))) for value in values]
    normalization = sum(weights)
    probabilities = [x/normalization for x in weights]
    assert min(probabilities) > 0
    trace_distance = sum(abs(p-(ctx.mpf('0.5') if j < 2 else 0)) for j, p in enumerate(probabilities))
    assert trace_distance <= ctx.mpf(1)/2**bits
    assert row['modular_to_Hamiltonian_time_factor'] == -row['inverse_temperature']
    assert not row['stationary_under_native_repair_inferred']


def test_exact_shared_mass64_witness_and_parent_domain():
    value = core.certificate()
    retained = core.core_bounds()
    assert retained['origin_weight_lower'] == F(7688, 7691)
    assert retained['ipr_lower'] == F(7688, 7691)**2
    assert retained['projection_angle_bound'] == F(27, 16252874)
    assert retained['exact_low_Kramers_ground_multiplicity'] == 2
    assert retained['center_parity_actual_ground_upper'] < F(1, 1000)
    assert retained['center_parity_actual_charged_lower'] > F(999, 1000)
    assert value['interaction_background_same_as_exact_bounded_radial_zero_modes']
    assert not value['finite_linear_texture_spectrum_substituted']


@pytest.mark.parametrize('arguments', [dict(coupling=0), dict(coupling=True), dict(coupling=.1),
                                      dict(coupling=F(1), mass=4), dict(coupling=F(1, 100), mass=4),
                                      dict(radius=0), dict(radius=True), dict(dimensions=2),
                                      dict(mass=2), dict(hopping_bound=1), dict(velocity=0)])
def test_undefined_or_uncontrolled_interacting_domains_rejected(arguments):
    with pytest.raises(ValueError):
        core.core_bounds(**arguments)
