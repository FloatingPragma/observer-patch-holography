"""Independent finite preparation calculation for the joint defect-wave clock."""
from collections import deque
from fractions import Fraction as F

from model.coherent_seam_repair.protection import BoundaryFlux


def executed_history():
    """Retain both filter flags and both detector outcomes after one idle tick.

    This starts from one specified pure basis state. Its success probability
    is distinct from the uniform-ensemble reference preparation probability.
    """
    from .model import JointObserver
    observer = JointObserver()
    initial = observer.basis_state(observer.initial_word, 0)
    rows = []
    for preparation, state in observer.preparation_filter(initial):
        evolved = observer.step(state)
        for detector, branch in observer.wave_basis_measurement(evolved, 0):
            rows.append({'preparation': preparation['outcome'],
                         'autonomous_ticks_after_preparation': 1,
                         'detector_coordinate': detector['coordinate'],
                         'detector_outcome': detector['outcome'],
                         'joint_probability': str(branch.norm_squared())})
    if sum(F(row['joint_probability']) for row in rows) != 1:
        raise ArithmeticError('preparation and readout histories are incomplete')
    return {'initial_basis': [observer.initial_word, 0, 3, 0, 0],
            'input_is_uniform_ensemble': False,
            'outcomes': rows,
            'total_probability': '1'}


def fixed_word_counts():
    """Trace contribution from the actual five controlled transport events.

    The wave port is unchanged by these events. At each port, compose the
    incident link moves in sorted source-seam order and count fixed link words.
    """
    law = BoundaryFlux()
    initial = 1 << law.edges.index((0, 1))
    words, queue = {initial}, deque([initial])
    while queue:
        word = queue.popleft()
        for edge in law.interior:
            result = law.transport(word, edge)[0]
            if result not in words:
                words.add(result)
                queue.append(result)
    counts = []
    for port in range(12):
        fixed = 0
        for word in sorted(words):
            result = word
            for edge in law.interior:
                if port in law.edges[edge]:
                    result = law.transport(result, edge)[0]
            fixed += result == word
        counts.append(fixed)
    return tuple(sorted(words)), tuple(counts)


def certificate():
    """Exact trace and complete same-U filter, with no dense 72000 matrix."""
    words, counts = fixed_word_counts()
    if len(words) != 10 or counts != (6, 6, 10, 10, 6, 0, 10, 10, 10, 6, 10, 6):
        raise ArithmeticError('fixed-boundary transport domain changed')
    dimension = len(words) * 12 * 5 * 2 * 60
    port_trace = 5 * 2 * F(-29, 30) * (10 * F(1, 2) + 50)
    trace = port_trace * sum(counts)
    normalization = F(5, 4) * dimension + trace
    success = 4 * normalization / (9 * dimension)
    eigenvalue_lower = 1 / (4 * normalization)
    if (dimension, trace, normalization, success, eigenvalue_lower) != (
            72000, -47850, 42150, F(281, 1080), F(1, 168600)):
        raise ArithmeticError('joint reference calculation changed')
    # Coefficients of I, U and U* in K_a* K_a, using U*U=I.
    success_effect = (F(5, 9), F(2, 9), F(2, 9))
    failure_effect = (F(4, 9), F(-2, 9), F(-2, 9))
    if tuple(a + b for a, b in zip(success_effect, failure_effect)) != (1, 0, 0):
        raise ArithmeticError('complete filter normalization failed')
    return {
        'schema': 'oph.coherent_defect_wave.preparation.v1',
        'joint_dimension': dimension,
        'fixed_boundary_link_words': list(words),
        'control_fixed_words_by_port': list(counts),
        'fixed_word_sum': sum(counts),
        'wave_port_diagonal_trace': str(port_trace),
        'joint_unitary_trace': str(trace),
        'stationary_reference_normalization': str(normalization),
        'reference_eigenvalue_lower': str(eigenvalue_lower),
        'success_probability_from_uniform_basis_ensemble': str(success),
        'failure_probability_from_uniform_basis_ensemble': str(1-success),
        'all_input_success_probability_lower': '1/9',
        'complete_filter_effect_coefficients': {
            'basis': ['I', 'U', 'U*'],
            'success': list(map(str, success_effect)),
            'failure': list(map(str, failure_effect)),
            'sum': ['1', '0', '0'],
        },
        'fresh_preparation_failure_tail': '(799/1080)^k',
        'executed_pure_input_preparation_readout': executed_history(),
        'stationarity_domain': 'success-conditioned reduced joint source algebra, preparation label unread',
        'scope': {'finite_same_U_reference': True,
                  'complete_specified_filter': True,
                  'arbitrary_effect_tomography': False,
                  'joint_archive_stationarity': False,
                  'Born_from_classical_repair': False,
                  'physical_energy_or_GeV_mass': False},
    }
