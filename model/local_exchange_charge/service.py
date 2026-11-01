"""Complete local dual-rail source services for the ordered-cover charge.

The four data bits encode two sites by 00=empty, 01=A, 10=B.  The
collision value 11 is retained.  Nine-bit repair execution reuses the
existing source compiler and its emitted Gray macros.  Preparation is a
nonunital, recorded local channel; its archive is traced before the native
history comparison.  Numerical admissions are not hardware certificates.
"""
from functools import lru_cache
from itertools import product
from types import MappingProxyType
import math

import numpy as np

from ..common_source import native_routing as routing
from ..defect_field_bridge import full_domain as full
from ..exchange_sector import statistics_service as source
from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import CompleteRun, execute_complete
from ..reference_observation.observation import frozen


AB, BA = 9, 6
READOUT_RECORDS = tuple(product(range(6), repeat=4))
ADMISSION = full.SOURCE_BOUND
HISTORY_ADMISSION = 2e-6


def _admit(value, ceiling, label):
    result = float(value)
    if not math.isfinite(result) or result < 0 or result > ceiling:
        raise ArithmeticError(label+' exceeded its finite nonnegative admission')
    return result


def _complete_error(run, target):
    failed = _admit(run.failure_probability, ADMISSION**2, 'complete source failure weight')
    return _admit(np.linalg.norm(run.success-target)+math.sqrt(failed), ADMISSION,
                  'complete source output error')


def _law(law):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('C or E required')
    return law


def _native_record(record):
    if (type(record) is not tuple or len(record) != 2 or
            any(type(x) is not int for x in record) or
            record[0] not in (0, 1) or not 0 <= record[1] < 16):
        raise ValueError('exact binary coin and raw integer0..15 required')
    return record


def _array(value):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked source amplitudes are not admitted')
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all():
        raise ValueError('finite numeric source amplitudes required')
    with np.errstate(over='ignore', invalid='ignore'):
        result = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(result).all():
        raise ValueError('source amplitudes must stay finite after conversion')
    return result


def _public_run(run):
    # Byte-backed data are immutable, but ndarray.shape metadata are not.
    # Detach every public snapshot so metadata changes cannot poison caches.
    return CompleteRun(run.word, run.bits, frozen(run.initial),
                       MappingProxyType({key: frozen(value) for key, value in run.primitives.items()}),
                       frozen(run.success), frozen(run.failure_weights), frozen(run.failure_gram))


def independent_target(law):
    """All 512 output addresses from an independent two-site endpoint rule."""
    _law(law)
    result = np.zeros((512, 16), complex)
    for coin, mask in product(range(2), range(16)):
        left, right = mask & 3, mask >> 2
        if left == 3 or right == 3:
            raw, output = 8, mask
        else:
            nl, nr = int(left != 0), int(right != 0)
            current = nl-((nl+nr)//2+coin*((nl+nr) % 2))
            output = (right if current == -1 else 0) if current else left
            output += 4*((left if current == 1 else 0) if current else right)
            raw = {0: 0, 1: 1, -1: 2}[current]
            if law == 'E':
                raw += 4*((nl+nr) % 2)
        result[output+16*coin+32*raw, mask] = 1/math.sqrt(2)
    return result


def repair_word(law):
    """Expand the predecessor's existing emitted schedule, without a new compiler.

    This is the literal word underlying full.source_isometry(2,2,law).
    Its cached Gray matrix is expanded so numerical first-export failures
    can also be reconstructed by execute_complete.
    """
    _law(law)
    assignments, quarters = [], [0]*512
    for coin, mask in product(range(2), range(16)):
        raw, output, sign, _ = full.branch_action(mask, law, coin, sites=2,
                                                 flavors=2, seam=(0, 1))
        origin = mask+16*coin
        assignments.append((origin, output+16*coin+32*raw))
        quarters[origin] = 2 if sign == -1 else 0
    images = routing.complete_permutation(tuple(assignments))
    word = [compiler.pulse('Y', -math.pi/4, [4])]+source.diagonal_word(tuple(quarters))
    macro = routing._canonical_word(False)
    for target, pattern in routing.permutation_edges(images):
        before = routing._conjugation_word(pattern)
        wires = [bit for bit in range(9) if bit != target]+[target]
        word += before+compiler.remap(macro, wires)+compiler.inverse(before)
    return word


def repair_run(law):
    """Every successful column and every lazy first-failure Kraus output."""
    _law(law)
    return _public_run(_repair_run(law))


@lru_cache(maxsize=2)
def _repair_run(law):
    # Both executions consume the existing source words.  The compact replay
    # and the literal prefix replay are independently compared below.
    compact, _, counts = full.source_isometry(2, 2, law)
    word = repair_word(law)
    if compiler.resources(word)['native_pulses'] != counts['native_pulses']:
        raise ArithmeticError('expanded predecessor pulse count changed')
    run = execute_complete(word, 9, np.eye(512, dtype=complex)[:, :16])
    target = independent_target(law)
    _complete_error(run, target)
    _admit(np.linalg.norm(run.success-compact), ADMISSION, 'compact/literal source comparison')
    return run


def repair_kraus(law):
    """All 32 coin/raw slots, including every collision and padded zero slot."""
    rows = repair_run(law).success.reshape(16, 2, 16, 16)
    return {(coin, raw): rows[raw, coin] for coin, raw in product(range(2), range(16))}


def _repair_grid(law):
    rows = repair_kraus(law)
    if type(rows) is not dict or set(rows) != set(product(range(2), range(16))):
        raise ArithmeticError('complete coin/raw repair grid required, including zeros and collision8')
    actual = np.zeros((512, 16), complex)
    for key, matrix in rows.items():
        _native_record(key)
        value = _array(matrix)
        if value.shape != (16, 16):
            raise ArithmeticError('complete sixteen-column repair block required')
        coin, raw = key
        actual[32*raw+16*coin:32*raw+16*coin+16] = value
    _admit(np.linalg.norm(actual-independent_target(law)), ADMISSION, 'complete raw repair grid')
    return rows


def repair_record_kraus(law, record):
    """All 6^5 native measurement records of coin and four raw qubits."""
    _law(law)
    if (type(record) is not tuple or len(record) != 5 or
            any(type(x) is not int or not 0 <= x < 6 for x in record)):
        raise ValueError('five complete six-outcome native readout records required')
    if any(x > 1 for x in record):
        return np.zeros((16, 16), complex)
    return _repair_grid(law)[record[0], sum(record[j+1] << j for j in range(4))]


def analyzer_word():
    """Diagonalize AB/BA exchange, fixing every other four-bit basis vector."""
    chart = [event for target in (1, 2, 3) for event in compiler.cnot(0, target)]
    basis = chart+compiler._gray_rotation(6, 7, -math.pi/4, 4)+compiler.inverse(chart)
    return compiler.inverse(basis)


def ideal_basis():
    result = np.eye(16)
    result[:, BA] = 0; result[:, AB] = 0
    result[BA, BA] = result[AB, BA] = 1/math.sqrt(2)
    result[AB, AB] = 1/math.sqrt(2); result[BA, AB] = -1/math.sqrt(2)
    return result


def analyzer_run():
    return _public_run(_analyzer_run())


@lru_cache(maxsize=1)
def _analyzer_run():
    run = execute_complete(analyzer_word(), 4, np.eye(16))
    _complete_error(run, ideal_basis().T)
    return run


def _record_mask(record):
    if (type(record) is not tuple or len(record) != 4 or
            any(type(x) is not int or not 0 <= x < 6 for x in record)):
        raise ValueError('four complete six-outcome native readout records required')
    return None if any(x > 1 for x in record) else sum(x << j for j, x in enumerate(record))


def measurement_certificate(bits=4):
    if type(bits) is not int or not 1 <= bits <= 5:
        raise ValueError('one through five measured code qubits required')
    maps = compiler.transfer_kraus(2, coherent=False)
    if (len(maps) != 6 or not np.array_equal(sum(k.conj().T@k for k in maps), np.eye(6))):
        raise ArithmeticError('complete six-outcome source readout changed')
    if any(np.any(m[:, :2]) for m in maps[2:]):
        raise ArithmeticError('noncode readout outcome ceased to be a retained zero')
    return {'local_qubit_readouts': bits, 'native_readout_load_events': bits,
            'native_readout_export_events': bits, 'complete_raw_slots': 6**bits,
            'valid_code_slots': 2**bits, 'retained_noncode_zero_slots': 6**bits-2**bits}


def readout_kraus(record):
    mask = _record_mask(record)
    result = np.zeros((16, 16), complex)
    if mask is not None:
        result[mask] = analyzer_run().success[mask]
    return result


def readout_value(record):
    mask = _record_mask(record)
    return 1 if mask == BA else -1 if mask == AB else 0


def _selector_weights(reference):
    if type(reference) is not bool:
        raise ValueError('explicit reference/defect preparation switch required')
    return (0.5, 0.5) if reference else (0.75, 0.25)


def preparation_run(mask, selector):
    if type(mask) is not int or mask not in (AB, BA) or type(selector) is not int or selector not in (0, 1):
        raise ValueError('AB/BA input and an exact binary selector required')
    return _public_run(_preparation_run(mask, selector))


@lru_cache(maxsize=4)
def _preparation_run(mask, selector):
    target = BA if selector == 0 else AB
    reset = [compiler.pulse('X', math.pi/2, [j]) for j in range(4) if (mask ^ target) >> j & 1]
    word = reset+compiler.inverse(analyzer_word())
    run = execute_complete(word, 4, np.eye(16)[:, mask])
    expected = ideal_basis()[:, target]
    _complete_error(run, expected[:, None])
    return run


def preparation_kraus(record, selector, *, reference=False):
    """Local CPTP replacement channel, after tracing its retained archive.

    A complete computational measurement is followed by conditional AB/BA
    reset and the compiled exchange-basis word.  A supplied independent
    classical selector has weights3/4,1/4 (reference1/2,1/2).  All 2592
    record/selector slots exist, including noncode zeros.  No unitary
    preparation from the central uniform state is asserted.
    """
    weights = _selector_weights(reference)
    if type(selector) is not int or selector not in (0, 1):
        raise ValueError('exact binary selector required')
    mask = _record_mask(record)
    result = np.zeros((16, 16), complex)
    if mask is not None:
        vector = preparation_run(mask, selector).success[:, 0] if mask in (AB, BA) else np.eye(16)[:, mask]
        result[:, mask] = math.sqrt(weights[selector])*vector
    return result


def preparation_failure_effect(*, reference=False):
    result = np.zeros((16, 16), complex)
    for mask, selector in product((AB, BA), range(2)):
        result[mask, mask] += _selector_weights(reference)[selector]*preparation_run(mask, selector).failure_probability
    return result


def local_reference():
    """Actual two-port reduction of the uniform twelve-port oneA/oneB sector."""
    weights = np.zeros(16)
    weights[0] = 90/132
    weights[[1, 2, 4, 8]] = 10/132
    weights[[AB, BA]] = 1/132
    return np.diag(weights)


def prepared_density(*, reference=False):
    initial = local_reference()
    # Noncode records are explicitly addressable zeros; all valid records
    # and both supplied selector outcomes are included in this contraction.
    result = np.zeros((16, 16), complex)
    for record, selector in product(READOUT_RECORDS, range(2)):
        k = preparation_kraus(record, selector, reference=reference)
        result += k@initial@k.conj().T
    return result


def history_amplitude(preparation_record, selector, first, second, final, columns, *, reference=False):
    """An individually addressable complete successful history; never normalize."""
    value = _array(columns)
    if value.ndim != 2 or value.shape[0] != 16 or not np.isfinite(value).all():
        raise ValueError('complete four-bit arbitrary-reference columns required')
    _native_record(first); _native_record(second)
    return (readout_kraus(final)@_repair_grid('E')[second]@_repair_grid('C')[first]
            @preparation_kraus(preparation_record, selector, reference=reference)@value)


def end_to_end(*, reference=False):
    """Preparation archive traced, then all C/E records and charge readouts."""
    initial = local_reference(); prepared = prepared_density(reference=reference)
    abort = {'preparation': float(np.trace(preparation_failure_effect(reference=reference)@initial).real),
             'C': 0., 'E': 0., 'analyzer': 0.}
    rows = []
    c_run, e_run, analyzer = repair_run('C'), repair_run('E'), analyzer_run()
    abort['C'] = float(np.trace(c_run.failure_gram@prepared).real)
    c_grid, e_grid = _repair_grid('C'), _repair_grid('E')
    for first, kc in c_grid.items():
        after_c = kc@prepared@kc.conj().T
        abort['E'] += float(np.trace(e_run.failure_gram@after_c).real)
        for second, ke in e_grid.items():
            before = ke@after_c@ke.conj().T
            abort['analyzer'] += float(np.trace(analyzer.failure_gram@before).real)
            observed = analyzer.success@before@analyzer.success.conj().T
            diagonal = observed.diagonal().real.copy()
            rows.append({'C': first, 'E': second, 'all16_code_probabilities': diagonal,
                         'remaining1280_readout_probabilities': 'exact zero: initialized valid-code readout'})
    totals = np.sum([row['all16_code_probabilities'] for row in rows], axis=0)
    _admit(abs(float(np.sum(totals))+sum(abort.values())-1), HISTORY_ADMISSION,
           'complete end-to-end normalization')
    epsilon = 0. if reference else .5
    expected = np.array([(1+epsilon)/132, (1-epsilon)/132, 130/132])
    outcomes = np.array([totals[BA], totals[AB], np.sum(totals)-totals[BA]-totals[AB]])
    _admit(np.max(abs(outcomes-expected)), HISTORY_ADMISSION, 'source charge readout')
    for weight in abort.values():
        _admit(weight, 1e-20, 'complete history failure weight')
    return {'rows': rows, 'charge_outcomes_plus_minus_zero': outcomes,
            'exact_target_plus_minus_zero': expected, 'source_first_failure_mass': abort,
            'complete_mass': float(np.sum(totals))+sum(abort.values())}


def three_port_reference():
    """Ambient 64-state reduction of the same twelve-port uniform sector."""
    weights = np.zeros(64)
    weights[0] = 72/132
    for site, flavor in product(range(3), range(2)):
        weights[1 << (2*site+flavor)] = 9/132
    for first, second in product(range(3), repeat=2):
        if first != second:
            weights[(1 << (2*first))+(1 << (2*second+1))] = 1/132
    return np.diag(weights)


def on_three_ports(matrix, sites):
    """Ordinary tensor-local qutrit/dual-rail embedding, including all padding.

    This is a tensor chart, not a CAR hopping route through the spectator.
    No global orientation register or spectator-dependent sign is applied.
    """
    value = _array(matrix)
    if (value.shape != (16, 16) or not np.isfinite(value).all() or
            type(sites) is not tuple or len(sites) != 2 or
            any(type(x) is not int or not 0 <= x < 3 for x in sites) or len(set(sites)) != 2):
        raise ValueError('complete two-site matrix and distinct three-port sites required')
    result = np.zeros((64, 64), complex)
    occupied_mask = (3 << (2*sites[0])) | (3 << (2*sites[1]))
    for column in range(64):
        local = ((column >> (2*sites[0])) & 3)+4*((column >> (2*sites[1])) & 3)
        spectator = column & ~occupied_mask
        for output in range(16):
            row = spectator+((output & 3) << (2*sites[0]))+((output >> 2) << (2*sites[1]))
            result[row, column] = value[output, local]
    return result


def _embed_checked(matrix, sites):
    """Independent Kronecker oracle checks all64 input/output tensor columns."""
    result = on_three_ports(matrix, sites)
    spectator = next(site for site in range(3) if site not in sites)
    indices = [4*(((mask >> (2*sites[0])) & 3)+4*((mask >> (2*sites[1])) & 3))
               +((mask >> (2*spectator)) & 3) for mask in range(64)]
    independent = np.kron(np.asarray(matrix), np.eye(4))[np.ix_(indices, indices)]
    if not np.array_equal(result, independent):
        raise ArithmeticError('ambient embedding changed a spectator or a complete matrix column')
    return result


def moving_history(law='C', *, reference=False):
    """Actual source matrices on all 64 rows: ports(0,1,5), seam(0,5).

    The initial preparation core is(0,1), and charge is read on(1,5).
    The preparation archive is traced in both compared experiments.  All
    32 native coin/raw outputs and 1296 terminal readout slots remain.
    """
    _law(law)
    initial = three_port_reference(); prepared = np.zeros((64, 64), complex)
    for mask, selector in product(range(16), range(2)):
        record = tuple((mask >> j) & 1 for j in range(4))
        k = _embed_checked(preparation_kraus(record, selector, reference=reference), (0, 1))
        prepared += k@initial@k.conj().T
    prep_fail = _embed_checked(preparation_failure_effect(reference=reference), (0, 1))
    repair = repair_run(law); analyzer = analyzer_run()
    repair_fail = _embed_checked(repair.failure_gram, (0, 2))
    analysis = _embed_checked(analyzer.success, (1, 2))
    analyzer_fail = _embed_checked(analyzer.failure_gram, (1, 2))
    abort = {'preparation': float(np.trace(prep_fail@initial).real),
             'repair': float(np.trace(repair_fail@prepared).real), 'analyzer': 0.}
    rows = []
    for label, local in _repair_grid(law).items():
        k = _embed_checked(local, (0, 2))
        before = k@prepared@k.conj().T
        abort['analyzer'] += float(np.trace(analyzer_fail@before).real)
        observed = analysis@before@analysis.conj().T
        # All 64 ambient rows enter, including invalid qutrit or charge rows.
        probabilities = observed.diagonal().real.reshape(16, 4).sum(axis=1)
        rows.append({'native_record': label, 'all16_code_probabilities': probabilities,
                     'remaining1280_readout_probabilities': 'exact zero: initialized valid-code readout'})
    totals = np.sum([row['all16_code_probabilities'] for row in rows], axis=0)
    actual = np.array([totals[BA], totals[AB], np.sum(totals)-totals[BA]-totals[AB]])
    epsilon = 0. if reference else .5
    expected = np.array([(1+epsilon/2)/132, (1-epsilon/2)/132, 130/132])
    mass = float(np.sum(totals))+sum(abort.values())
    _admit(abs(mass-1), HISTORY_ADMISSION, 'complete moving normalization')
    _admit(np.max(abs(actual-expected)), HISTORY_ADMISSION, 'complete moving charge probabilities')
    for weight in abort.values():
        _admit(weight, 1e-20, 'complete moving failure weight')
    return {'rows': rows, 'charge_outcomes_plus_minus_zero': actual,
            'exact_target_plus_minus_zero': expected, 'source_first_failure_mass': abort,
            'complete_mass': mass, 'prepared_density': prepared}


def moving_amplitude(law, preparation_record, selector, native_record, final_record,
                     columns, *, reference=False):
    """Actual full 64-row post-measurement output for any retained history."""
    value = _array(columns)
    if value.ndim != 2 or value.shape[0] != 64 or not np.isfinite(value).all():
        raise ValueError('complete six-bit arbitrary-reference columns required')
    _law(law); _native_record(native_record)
    prep = _embed_checked(preparation_kraus(preparation_record, selector, reference=reference), (0, 1))
    move = _embed_checked(_repair_grid(law)[native_record], (0, 2))
    read = _embed_checked(readout_kraus(final_record), (1, 2))
    return read@move@prep@value


def _resources(run):
    word = [compiler.pulse(kind, angle, targets) for kind, angle, targets in run.word]
    actual = compiler.resources(word)
    return source._failure_account({key: actual[key] for key in source.COUNT_KEYS})


def certificate():
    measurement = measurement_certificate()
    repairs = {}
    for law in ('C', 'E'):
        run = repair_run(law)
        repairs[law] = {'data_input_columns': 16, 'output_addresses': 512,
                        'coin_raw_slots': 32, 'collision_raw': 8,
                        'literal_source_pulses': len(run.word),
                        'retained_first_failure_slots': 4*sum(len(t)==1 for _, _, t in run.word)+2*sum(len(t)==2 for _, _, t in run.word),
                        'all_prefix_failure_vectors_reconstructible': True,
                        'complete_Stinespring_Frobenius_admission': ADMISSION,
                        'numerical_failure_gram_trace_below': 1e-20,
                        'pulse_resources': _resources(run),
                        'record_readout': measurement_certificate(5)}
        _admit(run.failure_probability, 1e-20, 'numerical export leakage')
    defect, reference = end_to_end(), end_to_end(reference=True)
    moving = {law: moving_history(law) for law in ('C', 'E')}
    moving_reference = {law: moving_history(law, reference=True) for law in ('C', 'E')}
    for row in (defect, reference, *moving.values(), *moving_reference.values()):
        _admit(sum(row['source_first_failure_mass'].values()), 1e-20, 'complete history export leakage')
    max_prep = max(len(preparation_run(m, s).word) for m, s in product((AB, BA), range(2)))
    prep_resources = _resources(preparation_run(AB, 0))
    analysis_resources = _resources(analyzer_run())
    return {'schema': 'oph.local_exchange_charge.source.v1',
            'ambient_dual_rail_data_bits': 4, 'repair_work_bits': 9,
            'moving_witness_maximum_simultaneous_logical_bits': 11,
            'reused_six_level_native_processor_additional': True,
            'repairs': repairs,
            'preparation': {'kind': 'complete local measurement and conditional replacement',
                            'selector': 'supplied independent classical mixture',
                            'defect_selector_weights': ['3/4', '1/4'],
                            'reference_selector_weights': ['1/2', '1/2'],
                            'complete_record_selector_slots': 2*6**4,
                            'maximum_conditional_native_pulses': max_prep,
                            'longest_conditional_pulse_resources': prep_resources,
                            'independent_selector_draws': 1,
                            'selector_preparation_physical_cost_calibrated': False,
                            'preparation_archive_traced_before_comparison': True,
                            'preparation_archive_claimed_indistinguishable': False},
            'measurement': measurement,
            'measurement_is_terminal_demolition': True,
            'global_charge_protection_after_terminal_readout_claimed': False,
            'complete_terminal_poststates_individually_reconstructible': True,
            'analyzer_native_pulses': len(analyzer_run().word),
            'analyzer_pulse_resources': analysis_resources,
            'end_to_end_complete_C_E_readout_slots': 32*32*6**4,
            'end_to_end_physical_C_E_readout_slots': 6**14,
            'end_to_end_longest_native_pulse_count': max_prep+sum(repairs[x]['literal_source_pulses'] for x in ('C', 'E'))+len(analyzer_run().word),
            'end_to_end_readout_load_and_export_events_each': 18,
            'end_to_end_blank_coin_record_qubit_initializations': 10,
            'input_uniform_reference_and_selector_are_supplied_initializations': True,
            'record_archive_format': 'each six-outcome measurement stored as an integer0..5; three bits suffice per value',
            'end_to_end_record_archive_bits_including_preparation_selector': 55,
            'complete_history_probability_error_admission': HISTORY_ADMISSION,
            'complete_history_error_composition': 'four complete local dilation errors <=2e-7 give amplitude error <=(1+2e-7)^4-1; probability error <=e(2+e)<2e-6',
            'defect_exact_charge_probabilities': ['1/88', '1/264', '65/66'],
            'reference_exact_charge_probabilities': ['1/132', '1/132', '65/66'],
            'source_first_failure_mass_below': 1e-20,
            'complete_normalization_error_below': HISTORY_ADMISSION,
            'moving_witness': {'ports': [0, 1, 5], 'initial_core': [0, 1],
                               'native_seam': [0, 5], 'readout_region': [1, 5],
                               'laws_executed': ['C', 'E'], 'ambient_rows_retained': 64,
                               'complete_coin_raw_readout_slots_per_law': 32*6**4,
                               'complete_physical_native_readout_slots_per_law': 6**9,
                               'measurement_load_and_export_events_each': 13,
                               'blank_coin_record_qubit_initializations': 5,
                               'record_archive_bits_including_preparation_selector': 40,
                               'longest_pulse_resources_by_law': {
                                   law: source._failure_account(source._add(prep_resources, repairs[law]['pulse_resources'], analysis_resources))
                                   for law in ('C', 'E')},
                               'defect_exact_charge_probabilities': ['5/528', '1/176', '65/66'],
                               'reference_exact_charge_probabilities': ['1/132', '1/132', '65/66'],
                               'transported_charge_expectation': '1/264',
                               'initial_core_charge_expectation': '1/132',
                               'complete_probability_error_admission': HISTORY_ADMISSION,
                               'source_first_failure_mass_below': 1e-20,
                               'sector_projection_used': False},
            'all_source_failures_and_native_raw_outcomes_retained': True,
            'source_failures_are_terminal': True, 'success_conditioning_used': False,
            'full_twelve_port_source_array_executed': False,
            'local_operator_certificate_covers_arbitrary_untouched_references': True,
            'floating_admission_is_outward_or_hardware_proof': False}
