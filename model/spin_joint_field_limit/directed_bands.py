"""Directed original-probability bands and a complete periodic source control."""
from fractions import Fraction as F

import numpy as np

from model.native_ordered_ground import source
from model.centered_particle_spectrum.certificate import sqrt_bounds
from . import arctan_bands


def directed_band(relative_radius=F(1, 4), momentum_scale=F(1)):
    if type(relative_radius) not in (int, F) or not 0 < relative_radius <= F(1, 4):
        raise ValueError('exact relative cap radius in (0,1/4] required')
    if type(momentum_scale) not in (int, F) or momentum_scale <= 0:
        raise ValueError('positive exact momentum scale required')
    delta, rho = F(relative_radius), F(momentum_scale)
    parent = arctan_bands.continuum_band()
    r = 3*rho/2
    energy = [r*x for x in parent['scaled_energy_band']]
    return {
        'reference': 'joint positive-field ground and H_B/(Jp) clock',
        'cap_centres': ['k0', '-k0'], 'centre_norm': rho,
        'cap_radius': delta*rho, 'annular_parameter': r,
        'momentum_norm_interval': [(1-delta)*rho, (1+delta)*rho],
        'energy_interval': energy,
        'each_charge_and_direction_original_probability_lower':
            parent['each_original_charge_probability_lower']/2,
        'energy_over_momentum_interval': [energy[0]/((1+delta)*rho),
                                          energy[1]/((1-delta)*rho)],
        'complete_records': ['minus_negative_cap', 'minus_positive_cap',
                             'minus_outside', 'plus_negative_cap',
                             'plus_positive_cap', 'plus_outside', 'complement'],
        'profile_changes_when_cap_shrinks': True,
        'bounded_creator_has_no_assumed_momentum_support': True,
        'source_even_collective_internal_charge': True,
        'original_odd_CAR_identification': False,
        'energy_shell_or_plane_wave_atom_inferred': False,
        'effective_finite_source_to_continuum_error_claimed': False,
    }


def translation_score(cap_radius, translation_norm, original_probability):
    inputs = (cap_radius, translation_norm, original_probability)
    if any(type(x) not in (int, F) or x < 0 for x in inputs):
        raise ValueError('nonnegative exact spatial and original-probability inputs required')
    radius, length, weight = map(F, inputs)
    if weight > F(1, 4):
        raise ValueError('a bounded creator branch has original probability at most 1/4')
    phase = radius*length
    return {'original_probability': weight,
            'translation_phase_vector_squared_error_upper': min(F(4), phase**2)*weight,
            'original_interference_score_defect_upper': min(F(2), phase**2/2)*weight,
            'success_probability_divided_out': False,
            'time_phase_or_dispersion_inferred': False}


def _permutation(sites, transform):
    return np.array([sum(((word >> x) & 1) << transform(x)
                         for x in range(sites)) for word in range(1 << sites)])


def finite_operators():
    """Literal six-site source, in its native Marshall frame."""
    sites = 6
    edges = ((0, 1), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5))
    h, scale = source.matrix(sites, edges, F(1), F(1, 32))
    size = 1 << sites
    identity = np.eye(size, dtype=np.int64)
    words = np.arange(size)
    z = [np.diag(1-2*((words >> x) & 1)) for x in range(sites)]
    x = [identity[words ^ (1 << site)] for site in range(sites)]
    # -i Y = X Z in this real matrix convention.
    minus_i_y = [x[site] @ z[site] for site in range(sites)]
    twice_f = (2, -1, -1, 2, -1, -1)
    epsilon = (1, -1, 1, -1, 1, -1)
    charge = sum((epsilon[i]*x[i] for i in range(sites)), np.zeros_like(identity))
    az = sum((epsilon[i]*twice_f[i]*z[i] for i in range(sites)), np.zeros_like(identity))
    ay = sum((twice_f[i]*minus_i_y[i] for i in range(sites)), np.zeros_like(identity))
    creators = {q: az+q*ay for q in (-1, 1)}  # B_q = numerator/32.
    translate = identity[:, _permutation(sites, lambda i: (i+2) % sites)]
    invert = identity[_permutation(sites, lambda i: (-i) % sites)]
    odd = identity[_permutation(sites, lambda i: (i+3) % sites)]
    zero = np.zeros_like(identity)
    for unitary in (translate, invert, odd):
        if not np.array_equal(unitary @ h, h @ unitary):
            raise ArithmeticError('declared symmetry is not a source symmetry')
    if not np.array_equal(translate @ translate @ translate, identity):
        raise ArithmeticError('even translation must have order three')
    if not np.array_equal(invert @ translate @ invert, translate.T):
        raise ArithmeticError('site inversion does not reverse the actual translation')
    if not np.array_equal(charge @ h, h @ charge):
        raise ArithmeticError('declared charge does not commute with source clock')
    if not np.array_equal(translate @ charge, charge @ translate):
        raise ArithmeticError('spatial translation does not preserve charge')
    if not np.array_equal(odd @ charge @ odd, -charge):
        raise ArithmeticError('odd source symmetry must exchange the charges')
    for q, b in creators.items():
        if not np.array_equal(charge @ b-b @ charge, 2*q*b):
            raise ArithmeticError('creator has the wrong source charge')
        if not np.array_equal(invert @ b @ invert, b):
            raise ArithmeticError('even profile is not inversion invariant')
        if not np.array_equal(odd @ b @ odd, -creators[-q]):
            raise ArithmeticError('charge pairing changes the source preparation')
        if not np.array_equal(b+translate @ b @ translate.T
                              +translate.T @ b @ translate, zero):
            raise ArithmeticError('profile has an unwanted zero momentum component')
    if not np.array_equal(creators[1].T, creators[-1]):
        raise ArithmeticError('opposite creators are not adjoints')
    # Fourier projection of ||A||<=1/2 gives ||B_q||<=1/2. This independent
    # exact row bound also certifies that the complete complement is positive.
    effects = sum((b.T @ b for b in creators.values()), zero)
    if max(sum(abs(int(v)) for v in row) for row in effects) > 512:
        raise ArithmeticError('complete creation complement is not certified')
    return {'H_integer': h, 'H_scale': scale, 'twice_charge': charge,
            'creator_numerators': creators, 'creator_denominator': 32,
            'even_translation': translate, 'inversion': invert,
            'charge_pairing': odd, 'complement_numerator': 1024*identity-effects,
            'sites': sites, 'edges': edges, 'p': F(1, 4), 'B': F(1, 32)}


def ground_proposal(h, scale, bits=48):
    if bits != 48:
        raise ValueError('the literal source certificate uses 48-bit proposals')
    energies, vectors = np.linalg.eigh(h/scale)
    if vectors[:, 0].sum() < 0:
        vectors[:, 0] *= -1
    denominator = 1 << bits
    return (np.rint(energies*denominator).astype(np.int64),
            np.rint(vectors*denominator).astype(np.int64))


def verify_ground(h, scale, energies, vectors, bits=48):
    """Full-basis exact residual admission; no floating eigenpair is trusted."""
    size = len(h)
    if bits != 48:
        raise ValueError('the literal source certificate uses 48-bit proposals')
    if (energies.shape != (size,) or vectors.shape != (size, size)
            or energies.dtype != np.int64 or vectors.dtype != np.int64):
        raise ValueError('complete literal integer eigensystem required')
    e, v, d = energies.astype(object), vectors.astype(object), 1 << bits
    if any(a > b for a, b in zip(e, e[1:])):
        raise ValueError('ordered complete eigenvalue proposal required')
    bound = F(195, 32)
    if any(abs(F(int(a), d)) > bound for a in e):
        raise ValueError('source energy exceeds its local norm bound')
    gram = v.T @ v-d*d*np.eye(size, dtype=object)
    delta = F(max(sum(abs(int(x)) for x in row) for row in gram), d*d)
    if delta >= F(1, 10**6):
        raise ArithmeticError('full basis is not certified independent')
    residual = h.astype(object) @ v*d-scale*v*e[None, :]
    error = sqrt_bounds(F(sum(int(x)**2 for x in residual.flat), scale**2*d**4))[1]
    radius = error+2*bound*delta
    e0, e1 = F(int(e[0]), d), F(int(e[1]), d)
    norm = int(v[:, 0] @ v[:, 0])
    r0 = sqrt_bounds(F(sum(int(x)**2 for x in residual[:, 0]), scale**2*d*d*norm))[1]
    separation = e1-radius-e0
    if separation <= r0:
        raise ArithmeticError('source ground is not isolated')
    ground_error = 2*r0/separation
    if ground_error > F(1, 10**8):
        raise ArithmeticError('ground-vector error exceeds admitted precision')
    return {'ground_vector_error': ground_error, 'ground_energy': (e0-r0, e0+r0),
            'gap_lower': separation-r0, 'source_norm_bound': bound,
            'maximum_energy_upper': F(int(e[-1]), d)+radius,
            'full_basis_residual_radius': radius}


def _ratio(interval, positive):
    if positive[0] <= 0:
        raise ArithmeticError('original probability has no positive floor')
    values = [a/b for a in interval for b in positive]
    return min(values), max(values)


def finite_control():
    row = finite_operators()
    h, scale = row['H_integer'], row['H_scale']
    energies, vectors = ground_proposal(h, scale)
    checked = verify_ground(h, scale, energies, vectors)
    g = vectors[:, 0].astype(object)
    norm = int(g @ g)
    b = row['creator_numerators'][1].astype(object)
    t = row['even_translation'].astype(object)
    image = b @ g
    values = []
    hj = image.copy()
    for power in range(3):
        # For real images and H commuting with T, either nontrivial Z3
        # character has moment (v*H^jv-v*H^jTv)/3.
        value = F(int(image @ hj-image @ (t @ hj)), 3*32**2*norm*scale**power)
        error = checked['ground_vector_error']*checked['source_norm_bound']**power/2
        values.append((value-error, value+error))
        hj = h.astype(object) @ hj
    weight = values[0]
    raw_mean, raw_second = _ratio(values[1], weight), _ratio(values[2], weight)
    maximum_mean_square = max(x*x for x in raw_mean)
    minimum_mean_square = 0 if raw_mean[0] <= 0 <= raw_mean[1] else min(x*x for x in raw_mean)
    variance = ((raw_second[0]-maximum_mean_square)/row['p']**2,
                (raw_second[1]-minimum_mean_square)/row['p']**2)
    e0 = checked['ground_energy']
    mean = ((raw_mean[0]-e0[1])/row['p'], (raw_mean[1]-e0[0])/row['p'])
    return {'actual_intervals': {'each_original_charge_direction_probability': weight,
                                'scaled_energy_mean': mean,
                                'scaled_energy_variance': variance},
            'ground_checks': checked,
            'complete_record_count': 13,
            'nonzero_charge_direction_records': 4,
            'zero_probability_records_retained': 8,
            'source_instance': 'six-site periodic bipartite ring',
            'p': row['p'], 'B': row['B'], 'J': F(1),
            'actual_translation_step': 2, 'actual_translation_order': 3,
            'translation_eigenphases': ['0', '2*pi/3', '-2*pi/3'],
            'source_charge_values': [-1, 1],
            'all_four_success_probabilities_exactly_equal': True,
            'all_input_effect_sum_identity': True,
            'finite_control_proves_continuum_floor': False,
            'finite_translation_labels_identified_with_continuum_momentum': False,
            'physical_particle_or_mass_assignment': False}


def certificate():
    """Fixed admitted intervals, insensitive to the eigensolver's last digits."""
    result = finite_control()
    admitted = {
        'each_original_charge_direction_probability': [F(481981, 10**8), F(481982, 10**8)],
        'scaled_energy_mean': [F(6097426, 10**6), F(6097427, 10**6)],
        'scaled_energy_variance': [F(48956, 10**6), F(48957, 10**6)],
    }
    for key, interval in admitted.items():
        actual = result['actual_intervals'][key]
        if not interval[0] <= actual[0] <= actual[1] <= interval[1]:
            raise ArithmeticError('actual source result escapes admitted '+key)
    checked = result['ground_checks']
    ground_interval = [F(-4309116, 10**6), F(-4309115, 10**6)]
    if not ground_interval[0] <= checked['ground_energy'][0] <= checked['ground_energy'][1] <= ground_interval[1]:
        raise ArithmeticError('actual ground escapes the fixed energy interval')
    if checked['gap_lower'] < F(2, 3):
        raise ArithmeticError('actual finite gap lacks its declared floor')
    if 4*(checked['maximum_energy_upper']-checked['ground_energy'][0]) > 18:
        raise ArithmeticError('complete positive source energy interval is too small')
    if admitted['scaled_energy_variance'][0]/admitted['scaled_energy_mean'][1]**2 <= F(1, 800):
        raise ArithmeticError('directional finite spectral width is not certified positive')
    result.pop('actual_intervals')
    result.pop('ground_checks')
    result.update({
        'certified_intervals': admitted,
        'ground_energy_interval': ground_interval,
        'ground_vector_error_upper': F(1, 10**8),
        'complete_positive_scaled_energy_band': [F(8, 3), F(18)],
        'each_original_probability_lower': F(1, 208),
        'relative_energy_variance_lower': F(1, 800),
        'full_integer_residual_and_Gram_admission': True,
        'floating_solver_trusted': False,
        'complement_probability_interval':
            [1-4*admitted['each_original_charge_direction_probability'][1],
             1-4*admitted['each_original_charge_direction_probability'][0]],
        'native_compiler_templates_reused': True,
        'complete_native_pulse_tape_executed': False,
    })
    return {'directed_joint_source_band': directed_band(),
            'finite_periodic_source_control': result}
