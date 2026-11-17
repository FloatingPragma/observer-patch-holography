"""Rebuild source neighborhoods, complete channel enclosures and controls."""
from fractions import Fraction as F
import hashlib
from math import isfinite
from pathlib import Path

from ..native_action import incident_refinement as clocks
from ..observer.scalar_response.experiment import canonical, compare_receipt as _compare_receipt
from ..observer_time_fields import dual
from . import construction as c

ROOT = Path(__file__).resolve().parents[2]
EXPORT_BITS = 64


def _runtime_diagnostics(value):
    """Retained diagnostics must satisfy the producer's own acceptance gates."""
    from .oracle import OUTPUT_TOLERANCE
    try:
        rows = value['runtime_controls']['executions']
        if type(rows) is not list or len(rows) != 2:
            raise ValueError('both complete runtime controls are required')
        for row in rows:
            for key in ('full_output_distance_numeric', 'reduced_output_distance_numeric',
                        'zero_padding_nuclear_norm', 'checkpoint_distance', 'probability_mass'):
                number = row[key]
                if type(number) is not float or not isfinite(number):
                    raise ValueError('finite original floating runtime diagnostics required')
                if key == 'probability_mass':
                    valid = abs(number-1) <= OUTPUT_TOLERANCE
                elif key in ('zero_padding_nuclear_norm', 'checkpoint_distance'):
                    valid = number == 0
                else:
                    valid = 0 <= number <= OUTPUT_TOLERANCE
                if not valid:
                    raise ValueError(f'runtime diagnostic exceeds its producer gate: {key}')
    except (KeyError, TypeError) as error:
        raise ValueError('complete retained runtime diagnostics required') from error


def compare_receipt(expected, actual):
    """Check scientific diagnostic domains before tolerant numerical replay."""
    _runtime_diagnostics(expected)
    _runtime_diagnostics(actual)
    _compare_receipt(expected, actual)


def dyadic(value, upper=False):
    scale = 1 << EXPORT_BITS
    units = -((-value.numerator*scale)//value.denominator) if upper else (
        value.numerator*scale)//value.denominator
    return str(F(units, scale))


def interval(pair):
    return [dyadic(pair[0]), dyadic(pair[1], True)]


def small_controls():
    from . import oracle
    rows = []
    # Full native configuration recursion, with all four marks and remote waits.
    edges = ((0, 1), (0, 2), (1, 2), (2, 3), (3, 4))
    for law in ('C', 'E'):
        exact = oracle.exact_history(5, edges, root=(0, 1, 2), law=law, cap=2)
        bounded = c.output(tuple(range(5)), edges, edges[:4], law=law, cap=2, bits=40)
        entries = 0
        for record, columns, coherence in zip(bounded.records, bounded.columns,
                                               bounded.coherence, strict=True):
            matrix = exact[record]
            for i in range(8):
                for j in range(8):
                    if i == j:
                        lo, hi = columns[i % 2][i//2]
                        if not lo/2 <= matrix[i][j] <= hi/2:
                            raise ArithmeticError('native oracle escaped a population interval')
                    else:
                        target = coherence/2 if (i, j) in ((2, 5), (5, 2)) else F(0)
                        if matrix[i][j] != target:
                            raise ArithmeticError('native oracle disagrees on quantum disturbance')
                    entries += 1
        rows.append({'law': law, 'full_configuration_ports': 5, 'cap_per_interval': 2,
                     'complete_quantum_entries_compared': entries})
    return rows


def run():
    from . import oracle
    ports, edges, incident, observed, exits, mapping = clocks.geometry()
    counts = clocks.exit_counts(ports, edges, observed, exits, 48)
    spatial, remainder = clocks.exit_bound(counts, 2, len(incident))
    outputs = [c.output(ports, edges, incident, law=law) for law in ('C', 'E')]
    masks = tuple(sum(c.category(e) == m for e in incident) for m in range(4))
    pi = tuple(F(count, len(incident)) for count in masks)
    p = F(len(incident), len(edges))
    q = 1-(1-p)**outputs[0].cap
    lc, le = pi[0]+(pi[1]+pi[2])/2, pi[0]
    distance = (q*q*(lc*lc-le*le)+q*(1-q)*(lc-le))/2
    word_mass = q*q*pi[1]*pi[0]
    export_error = F(128, 1 << EXPORT_BITS)
    if outputs[0].columns != outputs[1].columns:
        raise ArithmeticError('C and E no longer share the native population law')
    if (any(spatial+out.timeout+out.arithmetic_error+export_error >= F(7779, 10**13)
            for out in outputs) or 2*spatial >= F(2731, 10**13)):
        raise ArithmeticError('the retained complete-output budget was not established')
    channels = []
    for out in outputs:
        channels.append({'law': out.law,
            'arithmetic_half_trace_allowance': dyadic(out.arithmetic_error, True),
            'full_screen_to_export_midpoint_allowance': dyadic(
                spatial+out.timeout+out.arithmetic_error+export_error, True),
            'coherences': [{'record': [status, list(word)],
                'coherence_multiplier': interval((coherence, coherence))}
                for (status, word), coherence in zip(out.records, out.coherence, strict=True)]})
    paths = set((*clocks.screen.SOURCE_PATHS, *clocks.stopped_control.SOURCE_PATHS,
        *clocks.clock_histories.SOURCE_PATHS, *dual.SOURCE_PATHS,
        'model/native_action/incident_refinement.py',
        'model/primitive_observer/runtime.py', 'model/observer/scalar_response/experiment.py',
        'model/common_source/protocol.py',
        'tools/strict_json.py'))
    paths.update(str(path.relative_to(ROOT)).replace('\\', '/')
                 for path in Path(__file__).parent.glob('*.py') if not path.name.startswith('test_'))
    return {'schema': 'oph.coherent_screen_observer.v1',
        'geometry': {'levels': [3, 5], 'full_ports': [15360, 245760],
            'full_seams': [40320, 645120], 'ball_ports': len(ports), 'ball_seams': len(edges),
            'incident_seams': len(incident), 'mask_counts': list(masks),
            'source_address_map_sha256': hashlib.sha256(canonical(mapping).encode()).hexdigest(),
            'first_exit_counts_through_48': list(counts)},
        'contract': {'logical_basis_pair_masks': [1, 2], 'environment_density': '1/3',
            'scheduled_ticks': 2, 'feedback_phase_after_first_mask': 1,
            'complete_capped_records': 21, 'bell_block_dimension': 8,
            'population_dual_dimension': 1+len(ports)+len(ports)*(len(ports)-1)//2,
            'arbitrary_logical_reference': True, 'arbitrary_environment_correlations': False,
            'native_cp_semantics_executed': True, 'coherent_pulse_word_emitted': False,
            'midpoint_table_claimed_cp': False},
        'budgets': {'cap_per_interval': outputs[0].cap, 'integer_grid_bits': outputs[0].bits,
            'outward_export_bits': EXPORT_BITS, 'export_midpoint_allowance': str(export_error),
            'one_screen_spatial_allowance': dyadic(spatial, True),
            'spatial_series_remainder': dyadic(remainder, True),
            'two_completed_screen_allowance': dyadic(2*spatial, True),
            'complete_timeout_probability': interval((outputs[0].timeout, outputs[0].timeout)),
            'ball_attempts_maximum': 2*outputs[0].cap,
            'ball_attempts_expected_interval': interval((q*(1+q)/p, q*(1+q)/p)),
            'completed_full_attempts_expected': [str(F(2*m, len(incident))) for m in (40320, 645120)],
            'phase_gates_maximum': 1, 'phase_gates_expected_interval': interval((q*pi[1], q*pi[1])),
            'preparation_service': 'one logical/reference input and V-2 independent Bernoulli sites',
            'storage_service': 'V source binary sites, untouched reference, two masks and stop status',
            'physical_service_errors_certified': False},
        'quantum_witnesses': {'C_unread_exchange_per_completed_tick': str(lc),
            'E_unread_exchange_per_completed_tick': str(le),
            'complete_bell_C_E_distance': interval((distance, distance)),
            'completed_bell_C_E_distance': str((lc*lc-le*le)/2),
            'C_word_1_0_negative_partial_transpose_expectation': interval((-word_mass/4, -word_mass/4)),
            'C_word_1_0_witness_violation_lower': dyadic(word_mass/4),
            'E_word_1_0_separable': True,
            'word_1_0_divided_by_its_probability': False},
        'channels': channels,
        'shared_population_blocks': [{'record': [status, list(word)],
            'columns': [[interval(pair) for pair in row] for row in columns]}
            for (status, word), columns in zip(outputs[0].records, outputs[0].columns, strict=True)],
        'independent_exact_controls': small_controls(),
        'runtime_controls': oracle.run(),
        'source_bindings': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                            for path in sorted(paths)}}
