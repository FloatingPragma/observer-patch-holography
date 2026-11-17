"""Complete declared-service perturbation bounds, including source admission."""
from fractions import Fraction as F

CAP = 2**22


def budget(k=6, turns=2, calibrations=2):
    n, b, c = 24, 2, 64
    admission_attempts = c*(n-1)+6*b*b*(c*(n-1)+1)+n*c
    q = 256*2**k
    terms = {'admission_native_and_controller_services': 9*admission_attempts,
        'admission_initialization_and_finalization': 4096,
        'calibration_services': 9*37568*calibrations,
        'preparation_services': 6000,
        'rotation_services': 2*(192*turns)*(9*21+8),
        'marker_services': 9*(96*turns)*16,
        'counters_and_largest_ramsey_probe': 10*(2*(q-1)+128*2**k),
        'counter_choice_qft_and_analyzers': 2048}
    if sum(terms.values()) > CAP:
        raise ValueError('Complete source service cap exceeded')
    return {'admission_attempts_upper': admission_attempts, 'terms': terms,
            'total_upper': sum(terms.values()), 'uniform_upper': CAP}


def envelope(bits, copies):
    delta = F(1, 2**bits)
    error = (2*copies*CAP+1)*delta
    floor = F(3, 100)**copies
    if error >= floor:
        raise ValueError('Conditional perturbation certificate unavailable')
    return {'per_service_operator_norm_error': delta,
            'joint_initial_trace_norm_error': delta,
            'full_output_trace_norm_error': error,
            'acceptance_floor': floor,
            'conditional_binary_error': error/floor,
            'conditional_visibility_error': 2*error/floor}


CONTRACT = ('Each service is the complete declared bounded source/inverse, '
    'memory, record or next-address controller map, with its full Stinespring '
    'operator-norm error bounded by delta. This is not an elementary hardware '
    'gate count. The initial trace-norm error is one bound on the FULL JOINT '
    'input, not inferred from separate marginals. The hybrid bound includes '
    'admission, calibration, every failed prefix and one chosen readout. '
    'The exact rational-phase circuit is the nominal reference. These '
    'perturbation hypotheses are separate from its certified arithmetic.')
