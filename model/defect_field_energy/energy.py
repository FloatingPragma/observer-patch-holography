"""Local energy operators on the specified spatial CAR action.

The four-mode oracle uses the exact corrected quasiparticle chart of the
joint source experiment. Its identity reference is stated explicitly.
"""
from __future__ import annotations

import math

import numpy as np

from ..defect_field_bridge import joint
from ..defect_field_interactions import interaction
from ..defect_field_propagation import propagation

OPERATOR_GUARD = 2e-9


def nambu(h, delta):
    """Centered quadratic F(B)=1/2(c†,c) B (c,c†)^T."""
    h, delta = np.asarray(h, complex), np.asarray(delta, complex)
    if (h.ndim != 2 or h.shape[0] != h.shape[1] or delta.shape != h.shape
            or np.linalg.norm(h-h.conj().T) > OPERATOR_GUARD
            or np.linalg.norm(delta+delta.T) > OPERATOR_GUARD):
        raise ValueError('Hermitian h and antisymmetric Delta on one finite CAR domain required')
    return np.block([[h, delta], [delta.conj().T, -h.T]])


def normal_order_data(B, scalar=0.):
    """Return h,Delta,constant with every centered-CAR scalar retained."""
    B = np.asarray(B, complex)
    if B.ndim != 2 or B.shape[0] != B.shape[1] or len(B) % 2:
        raise ValueError('an even-dimensional Nambu matrix required')
    modes = len(B)//2
    h, delta = B[:modes, :modes], B[:modes, modes:]
    if np.linalg.norm(B-nambu(h, delta)) > OPERATOR_GUARD:
        raise ValueError('matrix is not a Hermitian centered CAR quadratic')
    return h, delta, complex(scalar)-np.trace(h)/2


def component_energy(indices, phase=0., *, sites=12, spacing=1., speed=1.,
                     reference='bare_propagation'):
    """The same complete spatial block, with local energies and drive derivative.

    Local h_x contains all onsite pairing and half each incident hopping.
    The default scalar is zero, matching the literal propagation source
    factors. The optional common-static reference subtracts11 cell vacuum
    energies on12 sites and requires its separate scalar source phase.
    """
    if isinstance(phase, (bool, str)) or not math.isfinite(float(phase)):
        raise ValueError('a finite external Higgs phase required')
    if reference not in ('bare_propagation', 'common_static'):
        raise ValueError('an explicit supported scalar reference convention required')
    indices = tuple(indices)
    h, _ = propagation.block(indices, sites=sites, spacing=spacing, speed=speed)
    sub = np.ix_(indices, indices)
    pairing = interaction.higgs_pairing(phase)[sub]
    derivative = sum((1j*row['creation_phase_sign']
                      *np.exp(1j*row['creation_phase_sign']*phase)*row['pairing'][sub]
                      for row in interaction.channel_pairings().values()),
                     np.zeros_like(pairing))
    delta = np.kron(np.eye(sites), pairing)
    delta_derivative = np.kron(np.eye(sites), derivative)
    modes_per_site, modes = len(indices), len(h)
    local, local_derivative = [], []
    for site in range(sites):
        support = np.arange(site*modes_per_site, (site+1)*modes_per_site)
        hi, di, dpi = np.zeros_like(h), np.zeros_like(h), np.zeros_like(h)
        hi[support, :] += h[support, :]/2
        hi[:, support] += h[:, support]/2
        di[np.ix_(support, support)] = pairing
        dpi[np.ix_(support, support)] = derivative
        local.append(nambu(hi, di))
        local_derivative.append(nambu(np.zeros_like(h), dpi))
    base_pairing = interaction.higgs_pairing(0)[sub]
    cell_vacuum = -float(np.linalg.svd(base_pairing, compute_uv=False).sum())/2
    scalar = -(sites-1)*cell_vacuum if reference == 'common_static' else 0.
    B = nambu(h, delta)
    dB = nambu(np.zeros_like(h), delta_derivative)
    if (np.linalg.norm(sum(local)-B) > OPERATOR_GUARD
            or np.linalg.norm(sum(local_derivative)-dB) > OPERATOR_GUARD):
        raise ArithmeticError('local energies do not reconstruct the same spatial action')
    return {'indices': indices, 'sites': sites, 'modes': modes, 'B': B,
            'local_B': tuple(local), 'dB_dphi': dB,
            'local_dB_dphi': tuple(local_derivative),
            'reference_scalar': scalar, 'local_reference_scalar': scalar/sites,
            'single_cell_vacuum': cell_vacuum, 'reference_convention': reference}


def energy_current(data, source, destination):
    """Centered quadratic matrix of i[h_source,h_destination]."""
    if any(type(site) is not int or not 0 <= site < data['sites'] for site in (source, destination)):
        raise ValueError('two actual site indices required')
    left, right = data['local_B'][source], data['local_B'][destination]
    return 1j*(left@right-right@left)


def continuity_residual(data, site, phase_rate=0.):
    """Compare local Heisenberg derivative plus explicit work to incoming flow."""
    local = data['local_B'][site]
    power = phase_rate*data['local_dB_dphi'][site]
    direct = 1j*(data['B']@local-local@data['B'])+power
    resolved = sum((energy_current(data, other, site) for other in range(data['sites'])),
                   np.zeros_like(local))+power
    return float(np.linalg.norm(direct-resolved))


def quadratic_operator(B, scalar=0.):
    """Independent complete Jordan--Wigner realization for at most six modes."""
    h, delta, constant = normal_order_data(B, scalar)
    modes = len(h)
    if not 1 <= modes <= 6:
        raise ValueError('explicit complete-Fock oracle is bounded to six modes')
    dimension = 1 << modes
    operators = []
    for mode in range(modes):
        a = np.zeros((dimension, dimension), complex)
        for mask in range(dimension):
            if mask >> mode & 1:
                a[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
        operators.append(a)
    result = constant*np.eye(dimension, dtype=complex)
    for i, a in enumerate(operators):
        for j, b in enumerate(operators):
            result += h[i, j]*a.conj().T@b
            result += (delta[i, j]*a.conj().T@b.conj().T
                       +delta[i, j].conjugate()*b@a)/2
    return result


def ring_current_oracle():
    """Nonzero J_0->2 of the actual12-ring, on its exact four-mode support.

    The intermediate site1 has disappeared from the commutator operator;
    no state or particle sector is projected. Fermionic routing represents
    the same CAR subalgebra on the four endpoint modes.
    """
    data = component_energy((27, 78))
    current = energy_current(data, 0, 2)
    support = np.array([0, 1, 4, 5])
    full = np.r_[support, support+data['modes']]
    outside = np.setdiff1d(np.arange(2*data['modes']), full)
    if np.linalg.norm(current[outside]) > OPERATOR_GUARD or np.linalg.norm(current[:, outside]) > OPERATOR_GUARD:
        raise ArithmeticError('ring-current support is not the complete four endpoint modes')
    bare = quadratic_operator(current[np.ix_(full, full)])
    R = small_energy_oracle()['R_qp_to_bare']
    return {'operator': R.conj().T@bare@R, 'bare_operator': bare,
            'B': current[np.ix_(full, full)], 'ring_sites': (0, 2),
            'original_internal_component': (27, 78), 'support_is_operator_identity_not_state_projection': True}


def small_energy_oracle(phase=0., *, reference_scalar=0.):
    """Return the complete16D energy algebra used by the joint experiment.

    Mode order is (site0:+,-;site1:+,-), little-endian occupation bits.
    ``J[y,x]=i[h_y,h_x]`` is energy current from site y to site x.
    This is the sum of the selected source factors, not an assertion that
    their executed ordered product equals its simultaneous exponential.
    """
    if any(isinstance(value, (bool, str)) or not math.isfinite(float(value))
           for value in (phase, reference_scalar)):
        raise ValueError('finite Higgs phase and scalar energy reference required')
    phase, reference_scalar = float(phase), float(reference_scalar)
    ladders = joint.ladders()
    hopping_bare = sum((coefficient*ladders[a].conj().T@ladders[b]
                       +coefficient.conjugate()*ladders[b].conj().T@ladders[a]
                       for a, b, coefficient in ((0, 2, .5j), (1, 3, -.5j))),
                      np.zeros((16, 16), complex))
    local_chart = joint.ideal_chart()
    full_chart = np.kron(local_chart, local_chart)
    hopping = full_chart.conj().T@hopping_bare@full_chart
    mass = interaction.bare_component_chart((27, 78))['energy']
    local_pairing = interaction.pair_hamiltonian(mass, phase)
    derivative = mass*(-math.sin(phase)*np.diag([-1., 0., 0., 1.]))
    derivative = derivative.astype(complex)
    derivative[0, 3] = derivative[3, 0] = mass*math.cos(phase)
    onsite = (np.kron(np.eye(4), local_pairing), np.kron(local_pairing, np.eye(4)))
    differentiated = (np.kron(np.eye(4), derivative), np.kron(derivative, np.eye(4)))
    local = tuple(value+hopping/2+reference_scalar*np.eye(16)/2 for value in onsite)
    H = local[0]+local[1]
    currents = np.zeros((2, 2, 16, 16), complex)
    for source in range(2):
        for destination in range(2):
            currents[source, destination] = 1j*(local[source]@local[destination]
                                                -local[destination]@local[source])
    return {'H': H, 'h_sites': local, 'J': currents,
            'dH_dphi': differentiated[0]+differentiated[1],
            'dh_sites_dphi': differentiated, 'H_hopping': hopping,
            'H_onsite': onsite[0]+onsite[1], 'R_qp_to_bare': full_chart,
            'reference_scalar': reference_scalar, 'energy': mass,
            'basis': 'joint corrected static QP chart; site outer, local(+,-) inner',
            'component': (27, 78),
            'reference_scope': 'selected four-mode factor; onsite vacuum constants included, optional global identity shift stated separately'}
