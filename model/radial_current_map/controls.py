"""Independent complete 64-state controls; floats serve only diagnostics."""
from fractions import Fraction as F

import numpy as np
from scipy.linalg import eigh

from ..full_fock_field.field import FieldGeometry
from ..radial_interacting_core.core import core_bounds


def chain_control(mass=4, coupling=F(1, 256)):
    modes, dimension = 6, 64
    c = []
    for mode in range(modes):
        a = np.zeros((dimension, dimension), complex)
        for mask in range(dimension):
            if mask & (1 << mode):
                a[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
        c.append(a)
    number = [a.conj().T@a for a in c]
    eye = np.eye(dimension)
    sigma = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    geometry = FieldGeometry((3, 1, 1), speed=1., interaction=float(coupling))
    h = geometry.single_h().toarray()+np.kron(np.diag([mass, 0, -mass]), sigma[0])
    H = sum((h[i, j]*c[i].conj().T@c[j] for i in range(modes) for j in range(modes)),
            np.zeros((dimension, dimension), complex))
    centered = sum(((number[2*x]-.5*eye)@(number[2*x+1]-.5*eye) for x in range(3)),
                   np.zeros_like(H))
    H += float(coupling)*centered
    pair = sum(((-1)**((x-1) % 2)*c[2*x].conj().T@c[2*x+1].conj().T for x in range(3)),
               np.zeros_like(H))
    uniform_pair = sum((c[2*x].conj().T@c[2*x+1].conj().T for x in range(3)), np.zeros_like(H))
    N = sum(number, np.zeros_like(H)); z = (N-3*eye)/2
    if (np.linalg.norm(H@pair-pair@H, 2) > 1e-11
            or np.linalg.norm(pair@pair.conj().T-pair.conj().T@pair-2*z, 2) > 1e-12):
        raise ArithmeticError('the complete chain eta algebra failed')
    sectors = {}
    for n in range(7):
        rows = [mask for mask in range(dimension) if mask.bit_count() == n]
        energies, vectors = eigh(H[np.ix_(rows, rows)])
        full = np.zeros((dimension, len(rows)), complex); full[rows] = vectors
        sectors[n] = energies, full
    ground = sectors[3][1][:, :2]
    minus = sectors[2][1][:, 0]
    plus = pair@minus
    eigplus = sectors[4][1][:, 0]
    low = np.column_stack([minus, *ground.T, plus])
    P = low@low.conj().T
    eta_residual = max(np.linalg.norm(pair@ground), np.linalg.norm(pair.conj().T@ground),
                       np.linalg.norm(pair@plus), np.linalg.norm(pair.conj().T@minus),
                       abs(np.vdot(plus, plus)-1),
                       np.linalg.norm(np.outer(plus, plus.conj())-np.outer(eigplus, eigplus.conj())))
    if eta_residual > 1e-9 or np.linalg.norm(low.conj().T@low-np.eye(4)) > 1e-9:
        raise ArithmeticError('the full interacting chain has the wrong eta multiplets')
    for site in range(3):
        added = np.array([plus.conj()@c[2*site+j].conj().T@ground for j in range(2)])
        removed = np.array([minus.conj()@c[2*site+j]@ground for j in range(2)])
        expected = (-1)**((site-1) % 2)*np.array([[0, 1], [-1, 0]])@removed
        if np.linalg.norm(added-expected) > 1e-9:
            raise ArithmeticError('the actual interacting pole amplitudes lost eta covariance')
    d0 = (number[2]+number[3]-eye)@(number[2]+number[3]-eye)
    spin = [sum((p[i, j]*c[2+i].conj().T@c[2+j] for i in range(2) for j in range(2)),
                np.zeros_like(H)) for p in sigma]
    for charge in (minus, plus):
        for S in spin:
            leakage = (eye-P)@S@charge
            expected = 1-np.vdot(charge, d0@charge).real
            if (np.linalg.norm(P@S@charge) > 1e-9
                    or abs(np.vdot(leakage, leakage).real-expected) > 1e-10 or expected <= 0):
                raise ArithmeticError('the complete charged-state spin leakage failed')
    row = core_bounds(coupling, mass=mass, dimensions=1)
    gap = sectors[2][0][0]-sectors[3][0][0]
    if not float(row['actual_charged_gap_lower']) <= gap <= float(row['actual_charged_gap_upper']):
        raise ArithmeticError('the complete chain failed its source gap bound')
    if np.linalg.norm(H@uniform_pair-uniform_pair@H, 2) < .1:
        raise ArithmeticError('the wrong unstaggered pairing countercontrol passed')
    uncentered = H+float(coupling)*N/2
    if np.linalg.norm(uncentered@pair-pair@uncentered-float(coupling)*pair, 2) > 1e-10:
        raise ArithmeticError('the uncentered-contact countercontrol failed')
    return {'sites': 3, 'all_Fock_states': dimension, 'all_number_sectors_checked': 7,
            'mass': str(mass), 'coupling': str(coupling),
            'same_literal_kinetic_and_centered_contact': True,
            'eta_commutators_and_actual_low_representation_pass': True,
            'all_sites_and_ground_columns_addition_removal_amplitudes_intertwine': True,
            'all_three_center_spin_currents_leak_outside_actual_quartet': True,
            'wrong_unstaggered_pairing_rejected': True,
            'uncentered_contact_has_nonzero_pair_commutator': True,
            'numeric_tolerance': '1e-9', 'numeric_eigenvalues_used_for_theorem_admission': False}


def certificate():
    return {'complete_Fock_controls': [chain_control(), chain_control(64, F(1, 65536))],
            'floating_residuals_are_diagnostics_only': True,
            'retained_raw_floats_in_cold_receipt': False}
