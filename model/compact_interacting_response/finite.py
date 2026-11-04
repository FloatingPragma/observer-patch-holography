"""Finite CAR readouts retaining the compact band's nonzero kinetic energy."""
from functools import lru_cache

import numpy as np
from scipy import linalg as la, sparse as sp

from model.compact_defect_pair import spectral
from model.defect_current_bridge.construction import native_kernel
from model.interacting_defect_response.finite import numpy_car
from model.mcpr_native_memory_bridge.construction import design, direct_memory
from model.common_source import native_family_record_memory as memory_law


def numerical_block(mass, coupling, e1, e2):
    """The exact reducing four-mode slice of the static compact mean source."""
    if not np.isfinite([mass, coupling, e1, e2]).all() or mass <= 0 or min(coupling, e1, e2) < 0:
        raise ValueError('positive mass and nonnegative finite coupling and spatial energies required')
    c = numpy_car(4)
    occupations = [a.T@a for a in c]
    number = sum(occupations, np.zeros((16, 16)))
    pairing = sum((c[i].T@c[j].T+c[j]@c[i] for i, j in ((0, 1), (2, 3))), np.zeros((16, 16)))
    interaction = (number-2*np.eye(16))@(number-2*np.eye(16))
    kinetic = e1*(occupations[0]-occupations[1])+e2*(occupations[2]-occupations[3])
    return kinetic+mass*pairing+coupling*interaction, interaction, c


def expected_spectrum(mass, coupling, e1, e2):
    r = np.hypot(mass, coupling)
    values = [0., 4*coupling, 2*coupling-2*r, 2*coupling+2*r]
    values += [s*e1+t*e2 for s in (-1, 1) for t in (-1, 1)]
    values += [coupling+s*e+t*mass for e in (e1, e2) for s in (-1, 1) for t in (-1, 1)]
    return np.sort(values)


def block_readout(mass, coupling, e1, e2):
    h, interaction, c = numerical_block(mass, coupling, e1, e2)
    energies, vectors = la.eigh(h)
    expected = expected_spectrum(mass, coupling, e1, e2)
    if np.max(abs(energies-expected)) > 1e-11:
        raise ArithmeticError('complete compact four-mode spectrum failed')
    if energies[1]-energies[0] < 1e-10:
        raise ValueError('the readout control requires a nondegenerate slice ground')
    ground = vectors[:, 0]
    weights = abs(vectors.conj().T@c[0]@ground)**2+abs(vectors.conj().T@c[0].T@ground)**2
    poles = []
    for energy, weight in zip(energies-energies[0], weights, strict=True):
        if weight <= 1e-12:
            continue
        if poles and abs(energy-poles[-1][0]) < 1e-10:
            poles[-1][1] += float(weight)
        else:
            poles.append([float(energy), float(weight)])
    if abs(sum(weight for _, weight in poles)-1) > 1e-10:
        raise ArithmeticError('the complete addition/removal sum rule failed')
    paired = 2*coupling-2*np.hypot(mass, coupling)
    blocked = -e1-e2
    single = coupling-mass-max(e1, e2)
    regime = min((paired, 'paired'), (blocked, 'blocked'), (single, 'single_odd'))[1]
    if regime == 'blocked':
        analytic_poles = [[e1+coupling-mass, .5], [e1+coupling+mass, .5]]
    elif regime == 'paired':
        r = np.hypot(mass, coupling)
        analytic_poles = sorted([[2*r-coupling+t*mass+s*e1, (1-t*mass/r)/4]
                                 for t in (-1, 1) for s in (-1, 1)])
    else:
        analytic_poles = None
    if analytic_poles is not None:
        merged = []
        for energy, weight in analytic_poles:
            if weight <= 1e-12:
                continue
            if merged and abs(energy-merged[-1][0]) < 1e-10:
                merged[-1][1] += weight
            else:
                merged.append([energy, weight])
        if np.shape(poles) != np.shape(merged) or np.max(abs(np.array(poles)-merged)) > 1e-10:
            raise ArithmeticError('independent pole and residue formula failed')
    n0, n2 = c[0].T@c[0], c[2].T@c[2]
    connected = float(ground.conj()@(n0@n2)@ground-(ground.conj()@n0@ground)*(ground.conj()@n2@ground))
    return {'mass_coefficient': format(mass, '.15g'), 'coupling': format(coupling, '.15g'),
            'spatial_energies': [format(e1, '.15g'), format(e2, '.15g')],
            'complete_Fock_energies': [format(x, '.12g') for x in energies],
            'slice_ground_regime': regime, 'slice_ground_energy': format(energies[0], '.12g'),
            'addition_plus_removal_poles_and_weights': [[format(x, '.12g') for x in pair] for pair in poles],
            'connected_cross_pair_density': format(connected, '.12g'),
            'pairing_and_quartic_commutator_nonzero': bool(la.norm(h@interaction-interaction@h) > 1e-12),
            'complete_spectrum_and_readout_formula_error_below': '1e-10',
            'slice_is_full_core_ground_sector': False}


@lru_cache(maxsize=1)
def compact_frame():
    source = spectral.source()
    _, energies, w = spectral.eigensystem(source)
    return source, energies, w


def source_controls():
    source, energies, w = compact_frame()
    e1, e2 = sorted(float(x) for x in energies if x > 0)
    memory = design('forward', 0)
    k = native_kernel(memory['mass_matrix_over_v'])
    h, c = source['h'], source['C']
    a = w.conj().T@(h@w)
    s = w.conj().T@(c@w.conjugate())
    # The tensor identities evaluate every spatial row and all 45 internal
    # coordinates while avoiding the redundant full 61740-square matrix.
    kinetic_error = la.norm(h@w-w@a)
    pair_error = la.norm(c@w.conjugate()-w@s)*sp.linalg.norm(k)
    core_delta = np.kron(s, k.conjugate().toarray())
    ac = np.kron(a, np.eye(45))
    core = np.block([[ac, core_delta], [core_delta.conj().T, -ac.conjugate()]])
    expected = sorted(x+y for x in energies for y in np.r_[-la.svdvals(k.toarray()), la.svdvals(k.toarray())])
    spectral_error = la.norm(la.eigvalsh(core)-expected)
    if max(kinetic_error, pair_error, spectral_error) > 1e-9:
        raise ArithmeticError('compact source/memory quadratic inclusion failed')
    full_blocks = [block_readout(m, 1/128, e1, e2) for m in memory['masses_over_v']]
    dropped = [block_readout(m, 1/128, 0., 0.) for m in memory['masses_over_v']]
    if any(row['slice_ground_regime'] != 'blocked' for row in full_blocks):
        raise ArithmeticError('retained source witness changed its ground sector')
    if any(row['slice_ground_regime'] != 'paired' for row in dropped):
        raise ArithmeticError('dropping compact kinetics failed its negative control')
    # All 21 charged singular channels and the twelve neutral modes enter this
    # full-core variational bracket. It does not diagonalize the 2**180 Fock space.
    charged = list(memory['masses_over_v'])
    charged += [mass for mass in (6.**-6/np.sqrt(2), 6.**-3/np.sqrt(2), 1/np.sqrt(2)) for _ in range(3)]
    charged += [mass for mass in (6.**-6/np.sqrt(2), 6.**-4/np.sqrt(2), 6.**-2/np.sqrt(2)) for _ in range(3)]
    free_ground = -2*sum(max(m, e1)+max(m, e2) for m in charged)-3*(e1+e2)
    even_blocks = 2*sum(int(m > e1)+int(m > e2) for m in charged)
    variational_upper = free_ground+even_blocks/128
    blocked_upper = -45*(e1+e2)
    return {'source': 'compact opposite-degree pair, radius=3, M=64, half-separation=1',
            'physical_annihilators': h.shape[0]*45, 'localized_annihilators': 180,
            'full_Nambu_rows': 2*h.shape[0]*45, 'localized_Nambu_rows': 360,
            'same_native_internal_dimension': 45, 'same_native_currents': 12,
            'memory_source_branch': 'forward', 'memory_inverse_branch': 0,
            'memory_target_ancestry': 'retained target-informed MCPR design; no runtime observed lepton masses',
            'finite_spatial_energies': [format(e1, '.15g'), format(e2, '.15g')],
            'actual_native_memory_kernel_inclusion_checked': True,
            'all_row_kinetic_pairing_and_spectral_errors_below': '1e-9',
            'interaction': 'lambda*(N_W-90)^2', 'declared_lambda_over_v': '1/128',
            'charged_lepton_slices': full_blocks,
            'full_180_mode_vacuum_bounds': {
                'charged_singular_channels_with_color_multiplicity': len(charged),
                'pairing_blocks_in_even_quadratic_ground': even_blocks,
                'quadratic_ground_energy_lower': format(free_ground, '.12g'),
                'quadratic_ground_trial_energy_upper': format(variational_upper, '.12g'),
                'fully_blocked_trial_energy_upper': format(blocked_upper, '.12g'),
                'best_displayed_ground_upper': format(min(variational_upper, blocked_upper), '.12g'),
                'exact_interacting_ground_energy_computed': False},
            'drop_spatial_kinetics_changes_every_slice_ground': True,
            'zero_spatial_control_poles': [row['addition_plus_removal_poles_and_weights'] for row in dropped]}


def complete_memory_control():
    """Complete native records on a nonlinear probe retaining compact kinetics."""
    _, energies, _ = compact_frame()
    e1, e2 = sorted(float(x) for x in energies if x > 0)
    design_row = design('forward', 0)
    law = direct_memory(design_row)
    transition, enabled = law['transition'], law['enabled_transition']
    fixed, _, _ = numerical_block(.25, .5, e1, e2)
    pairing, _, _ = numerical_block(1., 0., 0., 0.)
    gain = design_row['gain_over_v']
    projectors = [float(np.asarray(p, float)[0, 0]) for p in memory_law._projectors(False)]
    mean = fixed+gain*law['coefficient_matrix'][0, 0]*pairing
    initial = np.zeros((16, 16), complex)
    initial[0, 0] = 1
    unitary = la.expm(-1j*mean)
    target = unitary@initial@unitary.conj().T
    errors = []
    for steps in (32, 64, 128):
        off = la.expm(-1j*fixed/steps)
        on = np.array([la.expm(-1j*(fixed+gain*x*pairing)/steps) for x in projectors])
        blocks = law['stationary'][:, None, None]*initial
        for _ in range(steps):
            disabled = np.einsum('jk,kab->jab', transition-enabled, blocks)
            selected = np.einsum('jk,kab->jab', enabled, blocks)
            blocks = off@disabled@off.conj().T+on@selected@on.conj().transpose(0, 2, 1)
        errors.append(float(np.sum(abs(la.eigvalsh(blocks.sum(axis=0)-target)))))
    if not errors[0] > errors[1] > errors[2] > 0 or errors[0]/errors[-1] < 3.5:
        raise ArithmeticError('complete compact-memory probe failed its weak-step trend')
    return {'probe_Fock_rows': 16, 'memory_inputs': 12,
            'fine_records_per_memory_input': law['retained_records_per_input'],
            'spatial_kinetic_term_retained': True, 'noncommuting_quartic_term_retained': True,
            'steps': [32, 64, 128], 'trace_norm_errors': [format(x, '.10g') for x in errors],
            'all_disabled_and_enabled_outcomes_summed': True,
            'probe_is_full_compact_core': False,
            'probe_is_invariant_slice_of_each_family_record': False}
