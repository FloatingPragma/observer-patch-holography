"""Certified filtered weights from the complete FA98 integer witness."""
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from model.native_ordered_ground import certificate, source
from model.centered_particle_spectrum.certificate import decimal_bounds
from .filters import energy_filter, fourier_bounds

ROOT = Path(__file__).resolve().parents[2]
WITNESS = 'model/native_ordered_ground/witness.npz'


def finite_readout():
    with np.load(ROOT/WITNESS, allow_pickle=False) as archive:
        witness = {k: archive[k] for k in archive.files}
    verified = certificate.verify(witness)
    vectors = witness['vectors'].astype(object)
    ground = vectors[:, 0]
    norm = sum(int(x)**2 for x in ground)
    scale = 1 << certificate.BITS
    energies = [F(int(x), scale) for x in witness['energies']]
    origin = energies[0]
    e0lo, e0hi = map(F, verified['ground_energy_interval'])
    delta = F(verified['complete_polar_Weyl_radius'])+max(origin-e0lo, e0hi-origin)
    defect = F(verified['Gram_defect'])
    state_error = F(verified['normalized_ground_vector_error'])
    a, b, d = F(1), F(2), F(1, 4)
    control = fourier_bounds(a, b, d)
    moment = control['kernel_first_absolute_moment_upper']
    normalization = control['kernel_L1_upper']+1
    probability_error = 2*moment*delta+defect*(2+defect)+2*state_error
    energy_error = delta*(1+56*moment)+28*defect*(2+defect)+56*state_error
    odd = set(source.cube()[1])
    probes = {
        'local_Z0': ([1-2*(mask & 1) for mask in range(256)], 1),
        'uniform_Z_over_8': ([sum(1-2*((mask >> x) & 1) for x in range(8)) for mask in range(256)], 8),
        'staggered_Z_over_8': ([sum((-1 if x in odd else 1)*(1-2*((mask >> x) & 1)) for x in range(8)) for mask in range(256)], 8),
    }
    rows = []
    for name, (diagonal, divisor) in probes.items():
        coefficients = vectors.T @ (ground*np.array(diagonal, dtype=object))
        weights = [F(int(x)**2, scale**2*norm*divisor**2) for x in coefficients]
        filtered = [w*energy_filter(e-origin, a, b, d)**2
                    for w, e in zip(weights, energies, strict=True)]
        s = sum(filtered, F())
        numerator = sum((w*(e-origin) for w, e in zip(filtered, energies, strict=True)), F())
        lower, upper = max(F(), s-probability_error), min(F(1), s+probability_error)
        if lower <= 0:
            raise ArithmeticError('no certified positive finite filter weight')
        nlo, nhi = max(F(), numerator-energy_error), numerator+energy_error
        rows.append({'probe': name,
                     'filtered_weight_interval': decimal_bounds(lower, upper, 11),
                     'ideal_LCU_success_probability_interval': decimal_bounds(lower/normalization**2, upper/normalization**2, 11),
                     'conditional_excitation_energy_interval': decimal_bounds(nlo/upper, nhi/lower, 10),
                     'all_rejected_outcomes_probability_interval': decimal_bounds(1-upper/normalization**2, 1-lower/normalization**2, 11)})
    return {'source': 'FA98 open eight-site cube, J=1, B=1/4, true unique H_B ground',
            'complete_witness_dimension': 256, 'filter': {'a': '1', 'b': '2', 'd': '1/4'},
            'time_LCU_normalization': str(normalization),
            'weight_operator_and_state_error': str(probability_error),
            'energy_numerator_operator_and_state_error': str(energy_error),
            'readout_unit': 'dimensionless excitation energy of the selected H_B clock',
            'rows': rows, 'new_native_tape_executed': False,
            'finite_reference_identified_with_infinite_ordered_state': False,
            'relativistic_rest_mass_or_physical_species_identified': False}
