"""Full finite-source zero-field control with an exact spectral certificate."""
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from model.native_ordered_ground import certificate, source
from model.centered_particle_spectrum.certificate import decimal_bounds, sqrt_bounds
from .algebra import memory_interval

ROOT = Path(__file__).resolve().parents[2]
WITNESS = 'model/spin_dynamical_memory/zero_field_witness.npz'


def matrix():
    raw, scale = source.matrix()
    h = raw.copy()
    for mask in range(256):
        for site in range(8):
            h[mask^(1 << site), mask] += scale//8
    if scale % 8:
        raise ArithmeticError('exact source selection-field removal required')
    return h, scale


def propose(path=ROOT/WITNESS):
    """Untrusted finite proposal; all output is independently certified below."""
    h, scale = matrix()
    e, v = np.linalg.eigh(h/scale)
    if v[:, 0].sum() < 0:
        v[:, 0] *= -1
    unit = 1 << certificate.BITS
    np.savez_compressed(path, energies=np.rint(e*unit).astype(np.int64),
                        vectors=np.rint(v*unit).astype(np.int64))


def verify_witness(witness):
    h, scale = matrix()
    if set(witness) != {'energies', 'vectors'}:
        raise ValueError('complete zero-field witness required')
    if witness['energies'].shape != (256,) or witness['vectors'].shape != (256, 256):
        raise ValueError('all 256 zero-field levels and vectors required')
    if any(a.dtype != np.dtype('int64') for a in witness.values()):
        raise ValueError('literal integer zero-field witness required')
    unit = 1 << certificate.BITS
    energies = witness['energies'].astype(object)
    vectors = witness['vectors'].astype(object)
    if any(a > b for a, b in zip(energies, energies[1:])):
        raise ValueError('ordered zero-field spectrum required')
    if any(abs(F(int(x), unit)) > 12 for x in energies):
        raise ValueError('source zero-field norm bound exceeded')
    gram = vectors.T@vectors-unit**2*np.eye(256, dtype=object)
    defect = F(max(sum(abs(int(x)) for x in row) for row in gram), unit**2)
    if defect >= F(1, 10**6):
        raise ValueError('zero-field full Gram defect too large')
    residual = h.astype(object)@vectors*unit-scale*vectors*energies[None, :]
    epsilon = sqrt_bounds(F(sum(int(x)**2 for x in residual.flat), scale**2*unit**4))[1]
    radius = epsilon+24*defect
    e0, e1 = F(int(energies[0]), unit), F(int(energies[1]), unit)
    if e1-e0 <= 2*radius:
        raise ValueError('zero-field ground not uniquely isolated')
    norm = F(sum(int(x)**2 for x in vectors[:, 0]), unit**2)
    ground_residual = sqrt_bounds(F(sum(int(x)**2 for x in residual[:, 0]),
                                   scale**2*unit**4)/norm)[1]
    state_error = 2*ground_residual/(e1-radius-e0)
    if state_error >= F(1, 10**8):
        raise ValueError('zero-field ground proposal too loose')
    return {'dimension': 256, 'Gram_defect': str(defect), 'polar_radius': str(radius),
        'ground_energy_interval': decimal_bounds(e0-ground_residual, e0+ground_residual, 13),
        'gap_interval': decimal_bounds(e1-radius-e0-ground_residual,
                                      e1+radius-e0+ground_residual, 12),
        'ground_vector_error': str(state_error),
        'complete_spectrum_proves_unique_finite_ground': True,
        'floating_eigensolver_trusted': False}


def finite_field_continuity(gap_lower, field, eta):
    gap_lower, field, eta = map(F, (gap_lower, field, eta))
    delta = 4*field
    if eta <= 0 or field < 0 or gap_lower <= 2*delta:
        raise ValueError('positive resolution and isolated same-box field perturbation required')
    vector_error = 2*delta/(gap_lower-delta)
    response_error = 2*delta/eta**2+2*vector_error/eta
    mean_error = 2*delta+24*vector_error
    common_response_floor = 1/(24+2*delta+eta)
    memory_error = mean_error+response_error/common_response_floor**2
    return {'field': str(field), 'generator_perturbation_norm': str(delta),
        'ground_vector_difference_bound': str(vector_error),
        'resolvent_difference_bound': str(response_error),
        'mean_difference_bound': str(mean_error),
        'memory_difference_bound': str(memory_error)}


def build():
    with np.load(ROOT/WITNESS, allow_pickle=False) as archive:
        witness = {key: archive[key] for key in archive.files}
    proof = verify_witness(witness)
    h, scale = matrix()
    flip = np.arange(256)^255
    if not np.array_equal(h[np.ix_(flip, flip)], h):
        raise ArithmeticError('zero-field transverse parity changed')
    unit = 1 << certificate.BITS
    vectors = witness['vectors'].astype(object)
    ground = vectors[:, 0]
    norm = sum(int(x)**2 for x in ground)
    a = ground*np.array([1-2*(mask&1) for mask in range(256)], dtype=object)
    ha = h.astype(object)@a
    hg = h.astype(object)@ground
    state_error = F(proof['ground_vector_error'])
    ha_mean = F(int(a@ha), scale*norm)
    mean = ha_mean-F(int(ground@hg), scale*norm)
    variance = F(int(ha@ha), scale**2*norm)-ha_mean**2
    mean_error, variance_error = 48*state_error, 6*12**2*state_error
    energies = [F(int(x), unit) for x in witness['energies']]
    origin = energies[0]
    lo, hi = map(F, proof['ground_energy_interval'])
    delta = F(proof['polar_radius'])+max(origin-lo, hi-origin)
    gram = F(proof['Gram_defect'])
    shape = gram*(2+gram)+2*state_error
    coefficients = vectors.T@a
    weights = [F(int(c)**2, unit**2*norm) for c in coefficients]
    eta = F(1, 4)
    response = sum(w/(e-origin+eta) for w, e in zip(weights, energies))
    error = delta/eta**2+shape/eta
    ml, mh = memory_interval((mean-mean_error, mean+mean_error),
                             (response-error, response+error), eta)
    if ml <= 0 or variance-variance_error <= 0:
        raise ArithmeticError('zero-field dynamical memory not certified')
    continuation = finite_field_continuity(F(proof['gap_interval'][0]), F(1, 2**40), eta)
    if F(continuation['memory_difference_bound']) >= ml/2:
        raise ArithmeticError('same-box small-field continuation failed')
    return {'source': 'same 12-edge open cube, J=1, B=0, unique finite ground',
        'certificate': proof, 'all_spectral_levels_retained': len(weights),
        'mean_excitation_energy_interval': decimal_bounds(mean-mean_error, mean+mean_error, 12),
        'memory_strength_interval': decimal_bounds(variance-variance_error, variance+variance_error, 12),
        'eta': str(eta), 'memory_interval': decimal_bounds(ml, mh, 12),
        'quantified_small_field_continuation': continuation,
        'memory_explained_only_by_selection_field_torque': False,
        'finite_B_to_zero_limit_identified_as_infinite_ordered_ground': False}
