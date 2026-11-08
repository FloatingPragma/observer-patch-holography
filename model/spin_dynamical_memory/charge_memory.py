"""Complete same-clock memory within each actual conserved spin-charge channel."""
from fractions import Fraction as F
import numpy as np
from model.protected_defect_catalog import charge
from model.centered_particle_spectrum.certificate import sqrt_bounds, decimal_bounds
from . import finite
from .algebra import memory_interval


def build():
    operator = charge.operators()
    ground, ground_error, ground_proof = charge.ground_input()
    witness, checked, h, scale, own_ground = finite.verified_source()
    if (not (ground == own_ground).all() or not (h == operator['H']).all()
            or scale != operator['scale'] or checked != ground_proof):
        raise ArithmeticError('charge and memory must use exactly the same source and state')
    attachment = charge.full_source_charge()
    unit = 1 << finite.certificate.BITS
    vectors = witness['vectors'].astype(object)
    energies = [F(int(x), unit) for x in witness['energies']]
    origin = energies[0]
    elo, ehi = map(F, checked['ground_energy_interval'])
    delta = F(checked['complete_polar_Weyl_radius'])+max(origin-elo, ehi-origin)
    gram = F(checked['Gram_defect'])
    norm = sum(int(x)**2 for x in ground)
    hg = h.astype(object)@ground
    ground_energy = F(int(ground@hg), scale*norm)
    rows = []
    total_weight = F()
    identity = np.eye(256, dtype=np.int64)
    local_x = identity[np.arange(256)^1]
    completeness = np.zeros((256, 256), dtype=np.int64)
    for q in (-1, 1):
        if not np.array_equal(operator['twice_creation'][q],
                operator['local_Z']@(identity-q*operator['signs'][0]*local_x)):
            raise ArithmeticError('local X pointer followed by Z must prepare this exact charged branch')
        completeness += operator['twice_creation'][q].T@operator['twice_creation'][q]
        twice_a = operator['twice_creation'][q].astype(object)
        image = twice_a@ground
        image_norm = sum(int(x)**2 for x in image)
        if image_norm == 0:
            raise ArithmeticError('actual charge channel has zero source weight')
        weight = F(image_norm, 4*norm)
        total_weight += weight
        wl, wh = max(F(), weight-2*ground_error), min(F(1), weight+2*ground_error)
        if wl <= F(4, 25):
            raise ArithmeticError('charge herald probability not certified')
        normalized_error = 4*ground_error*sqrt_bounds(F(norm, image_norm))[1]
        ha = h.astype(object)@image
        ha_mean = F(int(image@ha), scale*image_norm)
        mean = ha_mean-ground_energy
        variance = F(int(ha@ha), scale**2*image_norm)-ha_mean**2
        mean_error = 28*(normalized_error+ground_error)
        variance_error = 6*14**2*normalized_error
        if variance-variance_error <= 0:
            raise ArithmeticError('same-charge complement memory not certified')
        coefficients = vectors.T@image
        weights = [F(int(x)**2, unit**2*image_norm) for x in coefficients]
        eta = F(1, 4)
        response = sum(w/(e-origin+eta) for w, e in zip(weights, energies))
        error = delta/eta**2+(gram*(2+gram)+2*normalized_error)/eta
        gl, gh = response-error, response+error
        ml, mh = memory_interval((mean-mean_error, mean+mean_error), (gl, gh), eta)
        if ml <= 0:
            raise ArithmeticError('resolved memory not strictly positive')
        rows.append({'source_spin_charge': q,
            'original_channel_probability_interval': decimal_bounds(wl, wh, 12),
            'normalized_branch_vector_error': str(normalized_error),
            'conditional_mean_excitation_energy_interval': decimal_bounds(mean-mean_error, mean+mean_error, 12),
            'conditional_memory_strength_interval': decimal_bounds(variance-variance_error, variance+variance_error, 12),
            'eta': str(eta), 'conditional_resolvent_interval': decimal_bounds(gl, gh, 12),
            'conditional_memory_interval': decimal_bounds(ml, mh, 12),
            'original_complete_charge_and_energy_score_interval': decimal_bounds(wl*eta*gl, wh*eta*gh, 12),
            'all_other_charge_or_complement_score_interval': decimal_bounds(1-wh*eta*gh, 1-wl*eta*gl, 12),
            'all_spectral_levels_retained': len(weights),
            'ground_energy_origin_and_clock_identical_to_full_response': True,
            'physical_charge_or_particle_assignment': None})
    if total_weight != 1:
        raise ArithmeticError('charge preparation decomposition lost original probability')
    if not np.array_equal(completeness, 4*identity):
        raise ArithmeticError('both local creation outcomes must form a complete instrument')
    return {'full_source_charge_attachment_replayed': attachment,
        'creation_instrument': 'measure local X0, then apply Z0; Kq=Z0(I-q*epsilon0*X0)/2',
        'complete_creation_outcome_count': 2, 'ideal_nonzero_charges': [-1, 1],
        'both_local_creation_outcomes_retained': True,
        'independent_fine_charge_audit_record_count': 256,
        'fine_X_audit_claimed_to_prepare_coherent_charge_branches': False,
        'native_local_creation_instrument_tape_executed_here': False,
        'channels': rows,
        'full_response_memory_explained_only_by_different_charge_means': False,
        'physical_particle_count_assigned': 0}
