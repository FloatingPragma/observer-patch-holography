"""Independent full-Fock controls for source state selection and entropy."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import eigh, expm

from tools import strict_json
from . import algebra, bounds, packet


def independent_source(sites, edges, c, U):
    """Tensor Jordan-Wigner construction, independent of the bit-hop audit."""
    z = np.diag([1., -1.])
    lower = np.array([[0., 1.], [0., 0.]])
    modes = []
    # Tensor mode order is reversed to match the source integer basis.
    for mode in range(2*sites):
        item = np.array([[1.]])
        for k in reversed(range(2*sites)):
            item = np.kron(item, z if k < mode else lower if k == mode else np.eye(2))
        modes.append(item)
    number = [a.T@a for a in modes]
    identity = np.eye(4**sites)
    defect = [(number[2*x]+number[2*x+1]-identity)@(
        number[2*x]+number[2*x+1]-identity) for x in range(sites)]
    N = sum(defect)
    T = np.zeros_like(N)
    for x, y, sign in edges:
        for spin in (0, 1):
            hop = modes[2*x+spin].T@modes[2*y+spin]
            T += sign*c/2*(hop+hop.T)
    return U/2*N+T, N, T


def exact_frame(G, N, nu):
    diagonal = np.diag(N)
    freq = diagonal[:, None]-diagonal[None, :]
    equal = freq == 0
    frame = np.eye(len(G))
    rotated = G.copy()
    for _ in range(5):
        off = np.where(equal, 0., rotated)
        S = np.divide(off, nu*freq, out=np.zeros_like(off), where=~equal)
        step = expm(S)
        frame, rotated = step@frame, step@rotated@step.T
    return frame, rotated


def spin_marginal(state, dictionary):
    count = int(round(np.sqrt(len(state))))
    tensor = (dictionary.T@state@dictionary).reshape(count, count, count, count)
    return np.einsum('abad->bd', tensor)


def thermal(H, beta):
    energy, vectors = eigh(H)
    weights = np.exp(-beta*(energy-energy[0]))
    weights /= weights.sum()
    return (vectors*weights)@vectors.T, energy[0]-np.log(np.exp(-beta*(energy-energy[0])).sum())/beta


def entropy(state):
    values = eigh(state, eigvals_only=True)
    positive = values[values > 1e-14]
    return float(-(positive*np.log(positive)).sum())


@pytest.mark.parametrize('sites,edges', [(2, ((0, 1, 1),)), (3, ((0, 1, -1), (1, 2, 1)))])
def test_full_source_dictionary_and_charge_extension_independently(sites, edges):
    T2, defect, dictionary = algebra.carrier(sites, edges)
    G, N, T = independent_source(sites, edges, 2, 6)
    assert np.array_equal(T, T2)
    assert np.array_equal(N, np.diag(defect.sum(axis=0)))
    assert np.allclose(G, 3*N+T2)
    assert np.array_equal(dictionary.T@dictionary, np.eye(4**sites))
    P = dictionary[:, :2**sites]
    assert not np.any(P.T@T2@P)
    assert set(np.diag(N)) == set(range(sites+1))
    assert all(algebra.audit(sites, edges)['identities'].values())


def test_actual_full_Fock_ground_and_thermal_references_converge():
    sites, edges, beta = 3, ((0, 1, -1), (1, 2, 1)), .8
    _, _, dictionary = algebra.carrier(sites, edges)
    P = dictionary[:, :2**sites]
    _, _, hopping = independent_source(sites, edges, 1, 1)
    spin = -(hopping@P).T@(hopping@P)
    target, target_F = thermal(spin, beta)
    lowest_spin = eigh(spin, eigvals_only=True)[0]
    ground_errors, thermal_errors = [], []
    for c in (8, 16, 32):
        G, N, _ = independent_source(sites, edges, c, c*c)
        Y, rotated = exact_frame(G, N, c*c/2)
        values, vectors = eigh(G)
        ground_vectors = vectors[:, values < values[0]+1e-9]
        ground = ground_vectors@ground_vectors.T/ground_vectors.shape[1]
        assert values[0] < 0
        assert np.trace(ground@N)/sites < (36*c/(c*c/2))**2
        ground_spin = spin_marginal(Y@ground@Y.T, dictionary)
        ground_errors.append(abs(np.trace(ground_spin@spin)-lowest_spin))
        source, source_F = thermal(G, beta)
        assert np.trace(source@G)/sites < np.log(4)/beta
        omega = spin_marginal(Y@source@Y.T, dictionary)
        thermal_errors.append(np.linalg.norm(omega-target, ord='nuc'))
        free_defect = np.trace(omega@spin)-entropy(omega)/beta-target_F
        relative = beta*(np.trace(omega@spin)-target_F)-entropy(omega)
        assert free_defect >= -1e-11
        assert abs(relative-beta*free_defect) < 1e-11
        assert thermal_errors[-1]**2 <= 2*relative+1e-10
        assert abs(source_F-target_F) < .08
        assert np.linalg.norm(rotated-Y@G@Y.T, 2) < 1e-9
    assert all(b < a*.4 for a, b in zip(thermal_errors, thermal_errors[1:]))
    # Ground-state selection is tested without choosing a vector in its doublet.
    assert max(ground_errors) < 1e-5


def test_nontrivial_modular_clock_has_source_sign_and_temperature():
    G, _, _ = independent_source(2, ((0, 1, 1),), 1, 3)
    beta, s = .7, .23
    rho, _ = thermal(G, beta)
    values, basis = eigh(rho)
    modular = (basis*np.exp(1j*s*np.log(values)))@basis.T
    A = np.diag([(-1.)**(mask & 1) for mask in range(16)])
    actual = modular@A@modular.conj().T
    expected = expm(-1j*beta*s*G)@A@expm(1j*beta*s*G)
    wrong_sign = expm(1j*beta*s*G)@A@expm(-1j*beta*s*G)
    assert np.linalg.norm(actual-expected, 2) < 1e-12
    assert np.linalg.norm(actual-wrong_sign, 2) > .1
    assert values[0] >= 4**(-2)*np.exp(-beta*2*(3/2+6))


def test_entropy_bound_retains_entangled_charge_spin_states():
    rng = np.random.default_rng(419)
    raw = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    state = raw@raw.conj().T
    state /= np.trace(state)
    tensor = state.reshape(4, 4, 4, 4)
    spin = np.einsum('abad->bd', tensor)
    charge = np.einsum('abcb->ac', tensor)
    assert abs(entropy(state)-entropy(spin)) <= entropy(charge)+1e-12
    # Discarding charge entropy would fail for this complete-carrier state.
    assert abs(entropy(state)-entropy(spin)) > .1


def test_rational_bounds_are_volume_uniform_where_claimed():
    first = bounds.reference_bounds(2**128, 2**256, support=8, volume=64)
    second = bounds.reference_bounds(2**129, 2**258, support=8, volume=64)
    large = bounds.reference_bounds(2**128, 2**256, support=8, volume=64**3)
    for key in ('ground_quadratic_form_error', 'spin_free_energy_defect_per_site'):
        assert 0 < first[key] < F(1, 10**12)
        assert F(49, 100) < second[key]/first[key] < F(51, 100)
        assert first[key] == large[key]
    assert first['finite_volume_trace_norm_error'] < large['finite_volume_trace_norm_error']
    assert first['thermal_dressed_onsite_defect_upper'] < F(1, 4)
    for value in (F(0), F(1, 9), F(2, 3), F(7, 2**400)):
        assert bounds.sqrt_upper(value)**2 >= value
    assert bounds.sqrt_upper(F(1, 9)) == F(1, 3)


@pytest.mark.parametrize('args', [(1, 1), (True, 2**80), (1., 2**80),
                                  (2**40, 2**80, 0), (2**40, 2**80, 1, True),
                                  (2**40, 2**80, 1, 2, 1)])
def test_invalid_domain_rejected(args):
    with pytest.raises(ValueError):
        bounds.reference_bounds(*args)


def test_receipt_replay_and_unfounded_promotions_rejected():
    retained = strict_json.load(packet.OUTPUT)
    packet.validate(retained)
    for change in ('source_bindings', 'all_pass', 'bounds'):
        altered = deepcopy(retained)
        altered[change] = False if change == 'all_pass' else True
        with pytest.raises(ValueError):
            packet.validate(altered)
    for flag in ('unique_thermodynamic_phase_identified',
                 'infinite_volume_KMS_identification_proved',
                 'source_state_selection_derived_from_repair'):
        altered = deepcopy(retained)
        altered['bounds'][flag] = True
        with pytest.raises(ValueError):
            packet.validate(altered)
