"""Hopping derivatives with the complete compact spectral complement retained.

The inequalities in ``bounds`` are rational analytic bounds. Eigenvectors,
derivative residuals and thermal readouts are floating numerical controls.
"""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
from scipy import linalg as la, sparse as sp

from model.compact_defect_pair.spectral import source
from model.compact_interacting_response.finite import numerical_block
from model.interacting_defect_response.thermodynamics import gibbs, response


def bounds(mass=64):
    if isinstance(mass, bool) or not isinstance(mass, int) or mass < 8:
        raise ValueError('integer texture amplitude at least eight required')
    gap = F(mass-6)
    derivative = F(6)/gap
    return {
        'hopping_parameter_domain': ['1/2', '3/2'],
        'band_dimension': 4,
        'band_absolute_energy_upper': '3',
        'complement_absolute_energy_lower': str(mass-3),
        'band_complement_gap_lower': str(gap),
        'spectral_projector_derivative_norm_upper': str(derivative),
        'Kato_generator_norm_upper': str(derivative),
        'localized_number_derivative_norm_upper': str(180*derivative),
        'quartic_force_norm_upper_per_coupling': str(32400*derivative),
        'localized_hopping_force_norm_upper': '360',
    }


@lru_cache(maxsize=1)
def source_parts():
    value = source(radius=3, mass=64, half_separation=1)
    h = value['h']
    onsite = sp.block_diag([h[2*i:2*i+2, 2*i:2*i+2]
                           for i in range(len(value['points']))], format='csr')
    hopping = (h-onsite).toarray()
    return value, onsite.toarray(), hopping


@lru_cache(maxsize=5)
def moving_band(parameter=1.):
    if not np.isfinite(parameter) or not .5 <= parameter <= 1.5:
        raise ValueError('hopping parameter must lie in [1/2,3/2]')
    value, onsite, hopping = source_parts()
    h = onsite+parameter*hopping
    energies, vectors = la.eigh(h, driver='evr')
    low = np.flatnonzero(abs(energies) < 32)
    high = np.flatnonzero(abs(energies) >= 32)
    if len(low) != 4:
        raise ArithmeticError('compact spectral count failed')
    w, q = vectors[:, low], vectors[:, high]
    low_energies = energies[low]
    numerator = q.conj().T@hopping@w
    denominator = low_energies[None, :]-energies[high, None]
    x = q@(numerator/denominator)
    projector = w@w.conj().T
    derivative = x@w.conj().T+w@x.conj().T
    kato = x@w.conj().T-w@x.conj().T
    compressed_force = w.conj().T@hopping@w
    return {'h': h, 'T': hopping, 'C': value['C'], 'energies': low_energies,
            'W': w, 'P': projector, 'dP': derivative, 'K': kato,
            'A': np.diag(low_energies), 'J': compressed_force,
            'energy_derivatives': compressed_force.diagonal().real,
            'complement_energies': energies[high], 'X': x}


def source_control():
    row = moving_band()
    h, t, p, dp, k, w, c = (row[key] for key in ('h', 'T', 'P', 'dP', 'K', 'W', 'C'))
    norm = la.norm
    residuals = {
        'projector_tangent': norm(dp-dp@p-p@dp),
        'differentiated_spectral_commutator': norm(h@dp-dp@h+t@p-p@t),
        'Kato_transport_equation': norm(k@p-p@k-dp),
        'Kato_skew_adjoint': norm(k+k.conj().T),
        'pairing_compatibility': norm(c@dp.conjugate()-dp@c),
        'Kato_pairing_compatibility': norm(c@k.conjugate()-k@c),
        'moving_frame_commutator_cancels': norm(w.conj().T@(h@k-k@h)@w),
    }
    epsilon = 1/1024
    minus, plus = moving_band(1-epsilon), moving_band(1+epsilon)
    projector_difference = (plus['P']-minus['P'])/(2*epsilon)
    energy_difference = (plus['energies']-minus['energies'])/(2*epsilon)
    residuals['independent_projector_centered_difference'] = norm(projector_difference-dp)
    residuals['independent_energy_centered_difference'] = norm(energy_difference-row['energy_derivatives'])
    if max(residuals.values()) > 1e-7:
        raise ArithmeticError('complete moving-band source derivative failed')
    actual_norm = float(la.svdvals(row['X'])[0])
    if not 0 < actual_norm < float(F(bounds()['spectral_projector_derivative_norm_upper'])):
        raise ArithmeticError('moving-projector norm or its analytic majorant failed')
    return {
        'one_body_dimension': h.shape[0], 'all_complement_modes_retained': len(row['complement_energies']),
        'band_energies': [format(float(v), '.12g') for v in row['energies']],
        'hopping_energy_derivatives': [format(float(v), '.12g') for v in row['energy_derivatives']],
        'projector_derivative_norm': format(actual_norm, '.12g'),
        'maximum_numerical_identity_or_difference_error_below': '1e-7',
        'floating_diagnostics_are_exact_certificates': False,
        'rational_uniform_bounds': bounds(),
    }


def force_control():
    """Exact three-mode control detecting the omitted moving-projector force."""
    import sympy as sy
    from model.interacting_defect_response.algebra import annihilators

    modes = annihilators(3)
    occupations = [a.T*a for a in modes]
    number = occupations[0]+occupations[1]
    tangent = modes[2].T*modes[0]+modes[0].T*modes[2]
    centered = number-sy.eye(8)
    interaction_force = centered*tangent+tangent*centered
    # A one-body rotation of modes zero and two moves the rank-two band.
    generator = modes[2].T*modes[0]-modes[0].T*modes[2]
    interaction = centered*centered
    if interaction_force != generator*interaction-interaction*generator:
        raise ArithmeticError('moving quartic derivative is not its commutator')
    if interaction_force == sy.zeros(8):
        raise ArithmeticError('the frozen-projector negative control is trivial')
    # Algebraic expansion to first order supplies an independent derivative.
    symbol = sy.Symbol('s', real=True)
    polynomial = (centered+symbol*tangent)*(centered+symbol*tangent)
    direct = polynomial.diff(symbol).subs(symbol, 0)
    if direct != interaction_force:
        raise ArithmeticError('direct noncommuting square derivative failed')
    return {'complete_Fock_rows': 8,
            'moving_quartic_force_rank': int(interaction_force.rank()),
            'moving_quartic_force_squared_trace': str(sy.trace(interaction_force*interaction_force)),
            'direct_polynomial_and_CAR_commutator_derivatives_agree': True,
            'freezing_the_localization_projector_omits_nonzero_force': True}


def thermal_control():
    from model.mcpr_native_memory_bridge.construction import design
    row = moving_band()
    positive = row['energies'][2:]
    velocity = row['energy_derivatives'][2:]
    mass = design('forward', 0)['masses_over_v'][2]
    coupling, beta = 1/128, 64.
    h, _, operators = numerical_block(mass, coupling, *positive)
    occupation = [a.conj().T@a for a in operators]
    force = velocity[0]*(occupation[0]-occupation[1])+velocity[1]*(occupation[2]-occupation[3])
    value = response(h, force, beta)
    epsilon = 1/1024
    m = moving_band(1-epsilon)['energies'][2:]
    p = moving_band(1+epsilon)['energies'][2:]
    hm = numerical_block(mass, coupling, *m)[0]
    hp = numerical_block(mass, coupling, *p)[0]
    vm, vp = gibbs(hm, beta), gibbs(hp, beta)
    estimates = {
        'free_energy_derivative': (vp['free_energy']-vm['free_energy'])/(2*epsilon),
        'energy_derivative': (vp['mean_energy']-vm['mean_energy'])/(2*epsilon),
        'entropy_derivative': (vp['entropy']-vm['entropy'])/(2*epsilon),
    }
    errors = [abs(estimates['free_energy_derivative']-value['mean_perturbation']),
              abs(estimates['energy_derivative']-value['energy_derivative']),
              abs(estimates['entropy_derivative']-value['entropy_derivative']),
              abs(value['entropy_derivative']-beta*(value['energy_derivative']-value['mean_perturbation']))]
    if max(errors) > 1e-6 or abs(value['mean_perturbation']) < 1e-6:
        raise ArithmeticError('same-source thermal hopping work failed')
    return {'complete_Fock_rows': 16,
            'source_mass': format(mass, '.12g'), 'coupling': '1/128', 'beta': '64',
            'source_positive_spatial_energies': [format(float(x), '.12g') for x in positive],
            'mean_hopping_force': format(value['mean_perturbation'], '.12g'),
            'energy_derivative': format(value['energy_derivative'], '.12g'),
            'entropy_derivative': format(value['entropy_derivative'], '.12g'),
            'work_subtracted_first_law_error_below': '1e-6',
            'independent_source_parameter_centered_differences_checked': True,
            'declared_reducing_slice_is_full_180_mode_Gibbs_state': False,
            'supplied_hopping_drive_is_a_metric_variation': False}
