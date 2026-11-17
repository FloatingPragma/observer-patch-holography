"""Exact obstruction to identifying binary parity recording with number noise.

The laboratory premise is the supplied Lindblad operator sqrt(gamma)*N,
gamma=18/s, in Appendix D of https://arxiv.org/html/2603.04744v1. Its data
deposit is https://doi.org/10.5281/zenodo.18752287. The paper attributes the
noise to trap-drive fluctuations; its two spin measurements select successful
preparation/evolution, rather than recording occupation parity after every
repair. These premises do not identify a native C/E instrument with that
apparatus. No detector records are fitted or statistically tested here.

For fixed occupation-diagonal binary projectors, unread recording multiplies
|m><n| by 1[b(m)=b(n)]. A nonnegative sum of generators r*(D_b-I) therefore
has decay rates that are nonnegative sums of cut metrics. Every such rate
matrix obeys the triangle inequality. Number diffusion has rates
gamma*(m-n)^2/2 and violates that inequality on levels 0,1,2 whenever gamma>0.
The exact separation and its sharp entrywise error floor below add a physical
mapping obstruction to the existing C/E parity-refinement identity.

The comparison excludes only this specified generator class, not a complete
native law. Rotated/noncommuting or adaptive operations, weak measurements,
extra states and effective reduced dynamics are outside its scope. Seconds
belong to the published apparatus noise model, not native repair steps.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import math
from numbers import Integral
from pathlib import Path

import numpy as np

from model.anharmonic_observations.source import phase_diffusion
from model.observer.source_discrimination import certificate
from model.observer.scalar_response.experiment import canonical, compare_receipt
from model.predictions import observer_instrument as native_instrument
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


PAPER_URL = 'https://arxiv.org/html/2603.04744v1'
DATA_URL = 'https://doi.org/10.5281/zenodo.18752287'
PUBLISHED_RATE = F(18)
PAIRS = ((0, 1), (1, 2), (0, 2))
ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    'model/identifiability_checks/dephasing_observability.py',
    'model/anharmonic_observations/source.py',
    'model/anharmonic_observations/protocol.json',
    'model/observer/native.py',
    'model/observer/source_discrimination/certificate.py',
    'model/observer/source_discrimination/circuit.py',
    'model/predictions/observer_instrument.py',
)


def _nonnegative_exact(value, name):
    if isinstance(value, bool) or not isinstance(value, (Integral, F)):
        raise ValueError(f'{name} must be an exact rational, not a float')
    value = F(value)
    if value < 0:
        raise ValueError(f'{name} must be nonnegative')
    return value


def cut_rates(bits, rate=F(1)):
    """Three off-diagonal decay rates of r*(D_b-I), all exactly rational."""
    bits = tuple(bits)
    if len(bits) != 3 or any(type(bit) is not int or bit not in (0, 1) for bit in bits):
        raise ValueError('three binary occupation labels are required')
    rate = _nonnegative_exact(rate, 'recording rate')
    return tuple(rate * int(bits[a] != bits[b]) for a, b in PAIRS)


def number_rates(rate=PUBLISHED_RATE):
    """Rates of D[sqrt(gamma)*N] on |0><1|, |1><2|, |0><2|."""
    rate = _nonnegative_exact(rate, 'number diffusion rate')
    return tuple(rate * F((a-b)**2, 2) for a, b in PAIRS)


def generator_separation(rate=PUBLISHED_RATE):
    """An exact separating functional and the sharp uniform rate error.

    For any number of nonnegative recording rates, linearity extends the
    eight verified inequalities to their entire conic hull. If all three
    target rates have error at most epsilon, the triangle inequality gives
    gamma <= 3*epsilon. Two cuts, each of rate 5*gamma/6, attain equality.
    Thus the minimax error is gamma/3; no numerical optimizer is involved.
    """
    rate = _nonnegative_exact(rate, 'number diffusion rate')
    target = number_rates(rate)
    cuts = []
    for bits in product((0, 1), repeat=3):
        values = cut_rates(bits)
        margin = values[0] + values[1] - values[2]
        if margin < 0:
            raise AssertionError('a binary cut violated its triangle inequality')
        cuts.append({'bits': list(bits), 'rates': [str(v) for v in values],
                     'triangle_margin': str(margin)})
    target_margin = target[0] + target[1] - target[2]
    best_rate = 5*rate/6
    first, second = cut_rates((0, 1, 1), best_rate), cut_rates((0, 0, 1), best_rate)
    closest = tuple(a+b for a, b in zip(first, second, strict=True))
    errors = tuple(a-b for a, b in zip(closest, target, strict=True))
    floor = rate/3
    if target_margin != -rate or max(map(abs, errors)) != floor:
        raise AssertionError('the exact separating functional or sharpness failed')
    return {
        'levels': [0, 1, 2], 'coherence_pairs': [list(pair) for pair in PAIRS],
        'number_decay_rates_per_apparatus_second': [str(v) for v in target],
        'all_binary_cuts': cuts,
        'separating_functional': 'Gamma01+Gamma12-Gamma02',
        'target_functional_exact': str(target_margin),
        'positive_recording_generators_can_equal_target': rate == 0,
        'sharp_maximum_entry_error_per_apparatus_second': str(floor),
        'attaining_cuts': [[0, 1, 1], [0, 0, 1]],
        'attaining_recording_rate_each': str(best_rate),
        'attaining_decay_rates': [str(v) for v in closest],
        'attaining_errors': [str(v) for v in errors],
        'class': 'sum r_k*(D_bk-I), r_k>=0, fixed occupation-diagonal binary projective recordings',
        'bound_scope': 'three generator entries, not a trace/diamond distance, fitted rate, or experimental confidence bound',
    }


def finite_time_witness(rate=PUBLISHED_RATE, duration=F(1, 1000)):
    """The same-parity 0/2 coherence separates noise from partial pinching.

    On (|0>+|2>)/sqrt(2), any mixture (1-q)*I+q*D_parity leaves the
    plus-projector probability at one. Number diffusion gives
    (1+exp(-2*gamma*t))/2. The rational gap lower bound follows from
    exp(x)>=1+x, x=2*gamma*t. The duration is an illustrative apparatus
    hold time, not a retained experiment or a conversion from repair ticks.
    """
    rate = _nonnegative_exact(rate, 'number diffusion rate')
    duration = _nonnegative_exact(duration, 'apparatus duration')
    exponent = 2*rate*duration
    difference_lower = exponent/(2*(1+exponent))
    # Use the actual supplied noise implementation with an untouched spin.
    dimension = 4
    rho = np.zeros((2*dimension, 2*dimension), dtype=complex)
    rho[np.ix_((0, 2), (0, 2))] = .5
    actual = phase_diffusion(rho, dimension, float(rate), float(duration))
    expected = (1+math.exp(-float(exponent)))/2
    actual_probability = float(np.trace(rho@actual).real)
    residual = abs(actual_probability-expected)
    if residual > 2e-15:
        raise AssertionError('the supplied noise channel differs from its exact formula')
    return {
        'state': 'spin_0 tensor (|0>+|2>)/sqrt(2)',
        'readout': 'projector onto the prepared state',
        'illustrative_apparatus_duration_seconds': str(duration),
        'number_dephasing_exponent_exact': str(exponent),
        'number_dephasing_plus_probability_formula': '(1+exp(-2*gamma*t))/2',
        'number_dephasing_plus_probability_numeric': expected,
        'partial_parity_plus_probability_exact': '1',
        'all_partial_parity_strengths_covered': '0<=q<=1',
        'probability_gap_lower_bound_exact': str(difference_lower),
        'strict_separation': exponent > 0,
        'supplied_channel_formula_residual': residual,
        'floating_residual_is_rigorous_error_bound': False,
        'prepared_state_is_retained_laboratory_trial': False,
    }


def native_parity_control():
    """Check actual C/E instruments, without identifying loads with phonons."""
    proof = certificate.parity_refinement_certificate()
    states = native_instrument.STATES

    def index(*ports):
        return states.index(tuple(ports.count(port) for port in range(12)))

    x, y, z = index(2, 3), index(0, 1), index(1, 2)
    witnesses = []
    for label, right in (('same_parity_totals_0_2', y), ('opposite_parity_totals_0_1', z)):
        unit = {(x, right): (F(1), F(0))}
        current = native_instrument.instrument_density('C', unit)
        contrast = native_instrument.instrument_density('E', unit)
        current_weight = sum((value[0] for block in current.values()
                              for value in block.values()), F(0))
        contrast_weight = sum((value[0] for block in contrast.values()
                               for value in block.values()), F(0))
        expected = (F(1), F(1)) if right == y else (F(1, 2), F(0))
        if (current_weight, contrast_weight) != expected:
            raise AssertionError('actual native matrix-unit parity control failed')
        witnesses.append({
            'label': label, 'source_ports': [[2, 3], [0, 1] if right == y else [1, 2]],
            'C_output_coefficient_sum_exact': str(current_weight),
            'E_output_coefficient_sum_exact': str(contrast_weight),
            'coefficient_sum_is_probability': False,
        })
    return {
        'full_encoded_parity_refinement': proof,
        'actual_native_instrument_matrix_unit_controls': witnesses,
        'seam_loads_identified_with_laboratory_phonon_numbers': False,
        'scope': 'E adds seam-parity recording to C within each old current/coin block; no apparatus identification',
    }


def run():
    protocol = strict_json.load(ROOT/'model/anharmonic_observations/protocol.json')
    apparatus = protocol['apparatus']
    if (F(str(apparatus['dephasing_per_second'])) != PUBLISHED_RATE
            or apparatus['dephasing_operator'] != 'sqrt(18) * N'
            or protocol['primary_paper'] != PAPER_URL
            or protocol['dataset_doi'] != '10.5281/zenodo.18752287'):
        raise AssertionError('the retained laboratory premise differs from this comparison')
    return {
        'schema': 'oph.dephasing_observability.v1',
        'scope': 'exact channel-mapping obstruction and native source controls, not an empirical C/E verdict',
        'primary_source': {
            'paper': PAPER_URL, 'data': DATA_URL,
            'locations': ['Methods II.1', 'Methods II.2', 'Appendix C', 'Appendix D'],
            'supplied_noise_operator': 'sqrt(18)*N',
            'rate_per_apparatus_second': str(PUBLISHED_RATE),
            'published_noise_mechanism': 'fluctuations of the radio-frequency trapping signal',
            'binary_spin_records': 'postselection after preparation and complete evolution, not per-step occupation-parity records',
            'final_readout': 'spin readout after a controlled displacement estimates the characteristic function',
        },
        'generator_obstruction': generator_separation(),
        'finite_time_obstruction': finite_time_witness(),
        'native_source_control': native_parity_control(),
        'not_excluded': ['rotated or noncommuting projectors', 'adaptive recording',
                         'weak measurements', 'coherent coupled evolution',
                         'extra states and effective reduced dynamics'],
        'physical_mapping_established': False,
        'laboratory_candidate_likelihoods_computed': False,
        'C_or_E_rejected_by_data': False,
        'additional_requirements': [
            'an operational encoding linking native load and accessible oscillator preparations',
            'a physical implementation/identification of the recorded channel and retained records',
            'an independently calibrated duration and noise model for that implementation',
            'candidate-specific predictions for the actual deposited detector settings',
        ],
        'source_bindings': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                            for path in SOURCE_PATHS}),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Binary-recording/number-dephasing obstruction replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
