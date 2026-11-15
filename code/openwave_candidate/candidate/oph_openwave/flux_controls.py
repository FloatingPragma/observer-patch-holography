"""Bare compact-state and static supplied-flux controls, without spin promotion.

The flux is a static Z6 cochain on one twelve-port carrier. Its one-axis
transport convention is supplied. Neither these doublets nor bare compact
states identify physical electron states in the joined integer-record model.
"""
import numpy as np
from scipy.linalg import expm
from sympy import Matrix

from .geometry import CARRIER_SEAMS, RING_PORTS, carrier, flat_ring


def static_flux_laplacian(flux):
    """Construct the carrier Laplacian with supplied oriented Z6 hop phases."""
    lap = np.zeros((12, 12), complex)
    eye = np.eye(12, dtype=complex)
    for a, b in CARRIER_SEAMS:
        column = eye[:, a] - np.exp(2j * np.pi * flux.get((a, b), 0) / 6) * eye[:, b]
        lap += np.outer(column, column.conj())
    return lap


def transport_rotation():
    """Return the chosen +120 degree class-3 transport lift about one axis."""
    sigma = {0: 11, 11: 5, 5: 0, 1: 10, 10: 4, 4: 1,
             2: 9, 9: 7, 7: 2, 3: 8, 8: 6, 6: 3}
    rotation = np.zeros((12, 12))
    for v in range(12):
        rotation[sigma[v], v] = -1 if v in {1, 5, 8, 9} else 1
    return rotation


def run():
    """Recompute scoped finite identities and controls without worker imports."""
    c = carrier()
    incidence = c.incidence
    lap = incidence @ incidence.T
    equations = np.vstack([lap - 6 * np.eye(12), np.eye(12)[list(RING_PORTS)]])
    compact_basis = Matrix(equations.astype(int).tolist()).nullspace()
    if len(compact_basis) != 2:
        raise ValueError('carrier compact eigenspace dimension changed')
    psi = np.array(compact_basis[0], dtype=float).ravel()
    psi /= np.linalg.norm(psi)
    screen = flat_ring(2, 2)
    extended = np.zeros(screen.ports)
    extended[:12] = psi
    compact_error = float(np.max(np.abs(screen.incidence @ screen.incidence.T @ extended - 6 * extended)))
    bare_h = c.hamiltonian
    time = .37
    initial = np.r_[psi, np.zeros(30)]
    measured = expm(-1j * time * bare_h) @ initial
    predicted = np.r_[np.cos(np.sqrt(6) * time) * psi,
                      -1j * np.sin(np.sqrt(6) * time) * (incidence.T @ psi) / np.sqrt(6)]
    breathing_error = float(np.max(np.abs(measured - predicted)))
    stationary = np.r_[psi, incidence.T @ psi / np.sqrt(6)] / np.sqrt(2)
    stationary_error = float(np.max(np.abs(bare_h @ stationary - np.sqrt(6) * stationary)))
    antipodal_flux = {(0, 1): 3, (0, 5): 3, (6, 8): 3, (7, 8): 3, (1, 7): 3}
    gauged = static_flux_laplacian(antipodal_flux)
    rotation = transport_rotation()
    commute_error = float(np.max(np.abs(rotation @ gauged - gauged @ rotation)))
    values = np.linalg.eigvalsh(gauged)
    doublet_error = float(np.max(np.abs(values[::2] - values[1::2])))
    expected_ground = (9 - np.sqrt(3) - np.sqrt(20 + 6 * np.sqrt(3))) / 2
    # One edge with class 3 has adjacent, nonantipodal flux faces.
    adjacent_values = np.linalg.eigvalsh(static_flux_laplacian({(0, 1): 3}))
    adjacent_split = float(adjacent_values[1] - adjacent_values[0])
    permutation = abs(rotation)
    ungauged_defect = float(np.max(np.abs(permutation @ gauged - gauged @ permutation)))
    fermions = [[a, k] for a in range(6) for k in range(6) if a * k % 6 == 3]
    return {
        'scope': 'Bare ring compact modes and static supplied one-axis Z6 flux controls; no joined-state stability or physical electron spin/statistics claim.',
        'checks': {
            'bare_compact_space_dimension_two': len(compact_basis) == 2,
            'bare_ring_compact_eigenstate': compact_error < 1e-12,
            'bare_port_transit_breathing': breathing_error < 1e-12,
            'bare_stationary_eigenvector': stationary_error < 1e-12,
            'static_flux_rotation_commutes': commute_error < 1e-12,
            'antipodal_flux_six_doublets': doublet_error < 1e-12,
            'class_three_ground': bool(abs(values[0] - expected_ground) < 1e-12),
            'transport_cube_minus_identity': bool(np.array_equal(np.linalg.matrix_power(rotation, 3), -np.eye(12))),
            'z6_five_fermion_sectors': len(fermions) == 5,
        },
        'mutations': {
            'ports_only_state_not_stationary': float(np.sum(abs(measured[12:]) ** 2)) > .1,
            'rotation_cube_not_rephase_invariant': bool(np.array_equal(np.linalg.matrix_power(-rotation, 3), np.eye(12))),
            'ungauged_rotation_rejected': ungauged_defect > .1,
            'nonantipodal_doublet_transfer_rejected': adjacent_split > .1,
            'even_flux_fermion_assignment_rejected': all(k % 2 for a, k in fermions),
        },
        'full_criterion_met': False,
        'supplied': {
            'compact_control': 'bare mass-zero incidence generator, ring torus 2x2',
            'time': time,
            'flux': 'class3 on antipodal faces; five class3 seams',
            'rotation': 'parallel-transport lift about one fixed face axis',
            'statistics': 'Z6 ribbon-algebra phase labels only; no exchange dynamics',
        },
        'measurements': {
            'compact_residual': compact_error, 'breathing_error': breathing_error,
            'stationary_error': stationary_error, 'commutator_error': commute_error,
            'doublet_error': doublet_error, 'ground': float(values[0]),
            'adjacent_flux_ground_split': adjacent_split,
            'ungauged_rotation_defect': ungauged_defect, 'fermion_sectors': fermions,
        },
        'physical_rows': {
            '4': 'Bare eigenstates do not transfer stability to joined integer-record particle states.',
            '6': 'One-axis kinetic labels supply no field-carried angular momentum on an electron.',
            '7': 'Rephasable one-axis transport sign and antipodal doublets supply no physical spin/statistics identification.',
        },
    }
