"""Independent complete-Fock spectral measures, currents and archive controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import eigh, block_diag

from . import algebra, bounds, response
from ..defect_field_bridge import full_domain


def oracle(profile, coupling):
    """All 64 occupation masks; no effective core Hamiltonian is substituted."""
    sigma = np.array([[0., 1.], [1., 0.]])
    if profile == 'texture':
        spatial = np.array([[1., -.5, 0], [-.5, 0, .5], [0, .5, -1.]])
        phi = np.array([1., 2., 1.])/np.sqrt(6)
        gap, ipr = F(6, 5), F(1, 2)
    elif profile == 'original':
        spatial = np.array([[0, .5j, 0], [-.5j, 0, .5j], [0, -.5j, 0]])
        phi = np.array([1., 0., 1.])/np.sqrt(2)
        gap, ipr = F(2, 3), F(1, 2)
    else:
        raise ValueError('unknown independent finite source')
    h = np.kron(spatial, sigma)
    annihilators = []
    for mode in range(6):
        c = np.zeros((64, 64), complex)
        for state in range(64):
            occupied = [j for j in range(6) if state & (1 << j)]
            if mode in occupied:
                c[state ^ (1 << mode), state] = (-1)**occupied.index(mode)
        annihilators.append(c)
    n = [c.conj().T@c for c in annihilators]
    H0 = sum(h[i, j]*(annihilators[i].conj().T@annihilators[j])
             for i in range(6) for j in range(6))
    contacts = [(n[2*x]-.5*np.eye(64))@(n[2*x+1]-.5*np.eye(64)) for x in range(3)]
    H = H0+float(coupling)*sum(contacts)
    sectors = {}
    for number in range(7):
        masks = [mask for mask in range(64) if mask.bit_count() == number]
        energies, vectors = eigh(H[np.ix_(masks, masks)])
        embedded = np.zeros((64, len(masks)), complex)
        embedded[masks, :] = vectors
        sectors[number] = energies, embedded
    zero = np.kron(phi[:, None], np.eye(2))
    estimate = bounds.residue_bounds(
        sites=3, coupling=coupling, gap_lower=gap, inverse_participation=ipr,
        detector_probability=F(2, 3) if profile == 'texture' else F(1, 2))
    return h, H0, H, annihilators, n, contacts, sectors, zero, estimate


def residues(annihilators, ground, plus, minus):
    add = np.array([(ground.conj().T@c@plus).ravel() for c in annihilators])
    remove = np.array([(ground.conj().T@c.conj().T@minus).ravel() for c in annihilators])
    return add@add.conj().T/2, remove.conj()@remove.T/2


@pytest.mark.parametrize('profile', ['texture', 'original'])
@pytest.mark.parametrize('coupling', [F(1, 1024), F(1, 4096)])
def test_all_sector_native_Lehmann_poles_residues_and_high_energy_weight(profile, coupling):
    _, H0, H, annihilators, _, _, sectors, zero, estimate = oracle(profile, coupling)
    eg = sectors[3][0][0]
    ground = sectors[3][1][:, :2]
    assert abs(sectors[3][0][1]-eg) < 1e-12
    assert min(values[0] for values, _ in sectors.values()) >= eg-1e-12
    plus, minus = sectors[4][1][:, 0], sectors[2][1][:, 0]
    charged = sectors[4][0][0]-eg
    assert abs(sectors[2][0][0]-eg-charged) < 1e-12
    assert float(estimate['charged_gap_lower']) <= charged <= float(estimate['charged_gap_upper'])
    Mplus, Mminus = residues(annihilators, ground, plus, minus)
    projector = zero@zero.conj().T
    error = float(estimate['pole_matrix_operator_norm_error'])
    for M in (Mplus, Mminus):
        assert np.linalg.norm(M-projector/2, 2) < error
        values = np.linalg.eigvalsh(M)
        assert np.max(abs(values[:4])) < 1e-12
        assert values[4] > 0 and abs(values[4]-values[5]) < 1e-12
    compressed = zero.conj().T@Mplus@zero
    z = compressed.trace().real/2
    assert np.linalg.norm(compressed-z*np.eye(2)) < 1e-12
    assert np.linalg.norm(zero.conj().T@Mminus@zero-compressed) < 1e-12
    assert .5-error <= z <= .5+1e-12
    assert min(sectors[2][0][1]-eg, sectors[4][0][1]-eg) >= float(
        estimate['other_charged_excitations_threshold_lower'])-1e-12
    for spin in (0, 1):
        c = sum(zero[i, spin].conjugate()*annihilators[i] for i in range(6))
        add_weights = np.sum(abs(ground.conj().T@c@sectors[4][1])**2, axis=0)/2
        rem_weights = np.sum(abs(ground.conj().T@c.conj().T@sectors[2][1])**2, axis=0)/2
        assert abs(add_weights.sum()-.5) < 1e-12
        assert abs(rem_weights.sum()-.5) < 1e-12
        assert add_weights[1:].sum() <= error and rem_weights[1:].sum() <= error
        for time in (0., 1., 100., 1000.):
            actual = (add_weights@np.exp(-1j*(sectors[4][0]-eg)*time)
                      +rem_weights@np.exp(1j*(sectors[2][0]-eg)*time))
            assert abs(actual-2*z*np.cos(charged*time)) <= 2*error+1e-11
    # The unperturbed sector projector is compared independently of Slater conventions.
    masks = [mask for mask in range(64) if mask.bit_count() == 3]
    _, free_vectors = eigh(H0[np.ix_(masks, masks)])
    free = np.zeros((64, 2), complex)
    free[masks] = free_vectors[:, :2]
    G = ground@ground.conj().T
    G0 = free@free.conj().T
    assert np.linalg.norm(G-G0, 2) <= float(estimate['projection_angle_bound'])
    assert np.linalg.norm((G-G0)/2, 'nuc') <= float(
        estimate['normalized_ground_trace_norm_error'])
    assert np.linalg.norm(H@ground-ground*eg) < 1e-11


def test_centered_contact_energy_current_is_retained_in_the_continuity_law():
    h, _, H, c, numbers, contacts, _, _, _ = oracle('texture', F(1, 1024))
    local, free = [], []
    for x in range(3):
        term = sum(h[i, j]*(c[i].conj().T@c[j]) for i in (2*x, 2*x+1)
                   for j in (2*x, 2*x+1))
        for y in range(3):
            if y == x:
                continue
            for i in (2*x, 2*x+1):
                for j in (2*y, 2*y+1):
                    term += (h[i, j]*c[i].conj().T@c[j]
                             +h[j, i]*c[j].conj().T@c[i])/2
        free.append(term)
        local.append(term+float(F(1, 1024))*contacts[x])
    assert np.count_nonzero(sum(local)-H) == 0
    energy_current = [[1j*(a@b-b@a) for b in local] for a in local]
    interaction_current = energy_current[0][1]-1j*(free[0]@free[1]-free[1]@free[0])
    assert np.linalg.norm(interaction_current) > 0
    assert np.linalg.norm(energy_current[0][2]) > 0  # distance-two support retained
    for x in range(3):
        assert np.count_nonzero(1j*(H@local[x]-local[x]@H)+sum(energy_current[x])) == 0
        nx = numbers[2*x]+numbers[2*x+1]
        outgoing = np.zeros_like(H)
        for y in range(3):
            if y != x:
                for i in (2*x, 2*x+1):
                    for j in (2*y, 2*y+1):
                        outgoing += 1j*(h[i, j]*c[i].conj().T@c[j]
                                        -h[j, i]*c[j].conj().T@c[i])
        assert np.count_nonzero(1j*(H@nx-nx@H)+outgoing) == 0


def test_complete_archive_distance_retains_actual_native_event_outputs():
    vector = np.zeros(16, complex)
    vector[[1, 4]] = [1/np.sqrt(2), 1j/np.sqrt(2)]
    rho = np.outer(vector, vector.conjugate())
    branches = full_domain.finite_branches(sites=2, flavors=2, law='E', seam=(0, 1))
    blocks = []
    for row in branches.values():
        K = np.zeros((16, 16), complex)
        for origin, output, sign in row:
            K[output, origin] = sign/np.sqrt(2)
        blocks.append(K@rho@K.conj().T)
    assert abs(sum(np.trace(block).real for block in blocks)-1) < 1e-12
    p = np.exp(-.3)
    actual = block_diag(p*rho, *((1-p)*block for block in blocks))
    ideal = block_diag(rho, np.zeros((512, 512)))
    assert abs(np.linalg.norm(actual-ideal, 'nuc')/2-(1-p)) < 1e-12
    assert np.linalg.norm(p*rho/p-rho) < 1e-12  # conditioning would erase the positive distance
    assert 1-p > .25


def test_original_site_dark_control_and_exact_source_currents():
    h, _, _, _, _, _, _, zero, _ = oracle('original', F(1, 1024))
    assert np.count_nonzero(zero[2:4]) == 0
    assert np.linalg.norm(h@zero) < 1e-14
    certificate = algebra.certificate()
    assert len(certificate['complete_native_record_currents']) == 6
    assert all(row['every_branch_number_and_current_checked']
               for row in certificate['complete_native_record_currents'])


def test_two_proved_cube_sources_and_nonzero_probability_without_underflow():
    for profile in ('linear_texture_cube', 'unmodified_odd_cube'):
        value = response.source_result(profile)
        assert value['full_Fock_dimension'] == 2**54
        assert value['detector_pole_certified_visible']
        assert value['each_zero_mode_low_pole_weight_lower'] > F(49, 100)
        assert value['projection_angle_bound'] < F(1, 216)
        archive = bounds.poisson_archive(edges=54, rate=F(1),
                                         time=1/value['charged_gap_lower'])
        assert archive['negative_log_no_event_probability'] > 10**7
        assert archive['no_event_probability'].startswith('exp(-')
        assert not archive['autonomous_spectral_protection_inferred']
    packet = response.certificate()
    for window in packet['declared_positive_rate_continuing_windows'].values():
        assert F(window['negative_log_no_event_probability']) == F(1, 1000)
        assert F(window['rate_per_seam']) > 0
        assert window['every_event_and_abort_output_retained']


@pytest.mark.parametrize('field,value', [
    ('sites', True), ('sites', 2), ('coupling', 1.), ('coupling', F(1)),
    ('gap_lower', 0), ('inverse_participation', 0), ('detector_probability', -1)])
def test_invalid_or_uncontrolled_pole_parameters_fail(field, value):
    arguments = dict(sites=27, coupling=F(1, 65536), gap_lower=F(2, 3),
                     inverse_participation=F(1, 8), detector_probability=F(1, 8))
    arguments[field] = value
    with pytest.raises(ValueError):
        bounds.residue_bounds(**arguments)


def test_unknown_source_and_invalid_clock_are_not_silently_substituted():
    with pytest.raises(ValueError):
        response.source_result('selected_physical_particle')
    for arguments in (dict(edges=0, rate=1, time=1),
                      dict(edges=54, rate=-1, time=1),
                      dict(edges=54, rate=1, time=1.)):
        with pytest.raises(ValueError):
            bounds.poisson_archive(**arguments)
