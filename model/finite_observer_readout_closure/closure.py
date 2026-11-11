"""Full-source witnesses, with no physical-time or sampled-rank promotion."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from model.relative_modular_dynamics import construction as modular
from model.joint_clock_process import construction as clocks
from model.joint_clock_process import clock_change
from model.relative_observer_geometry import temporal_records as temporal
from model.particle_common_contract.provenance import import_closure

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def matrix_record(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def state_digest(state):
    rows = [[i, j, str(value[0]), str(value[1])] for (i, j), value in sorted(state.items())]
    return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()


@lru_cache(maxsize=1)
def periodic_y():
    # Complete exact square-root effects on all M66, including coherences
    # outside the 86-dimensional protocol algebra.
    a, b = modular.PROBE
    identity = s.eye(66)
    support, y = s.zeros(66), s.zeros(66)
    support[a, a] = support[b, b] = 1
    y[a, b], y[b, a] = -s.I, s.I
    outside = identity-support
    kraus = tuple((support+sign*y)/2+outside/s.sqrt(2) for sign in (1, -1))
    assert y.H == y and y*y == support
    assert sum((k.H*k for k in kraus), s.zeros(66)) == identity
    for k, sign in zip(kraus, (1, -1), strict=True):
        assert k.H*k == (identity+sign*y)/2
    degree = modular.events.embedded()[0]
    assert degree[a] == degree[b]
    cases = []
    for law in ('C', 'E'):
        columns = modular.event_columns(law)
        assert len(columns) == 66**2
        for i in range(66):
            for j in range(66):
                assert sum(coefficient for u, v, coefficient in columns[66*i+j] if u == v) == int(i == j)
        marked = modular.marked(law)
        assert max(marked.no_current_counts[i][i] for i in range(66)) == 52
        harmonic = modular.modular_clock_audit(law)
        prepared = modular.prepared_states(law)
        assert list(map(modular.trace, prepared)) == [F(28, 33), F(5, 66), F(5, 66)]
        finite = []
        for quarter in (-1, 0, 1, 2):
            for cap in (4, 32, 128):
                rows = [{'record': [0, 0, -1], 'probability': '28/33',
                         'terminal_state_sha256': state_digest(prepared[0])}]
                plus = F(0)
                for flag in (1, 2):
                    initial = modular.modular_action(prepared[flag], quarter)
                    for n in range(1, cap+1):
                        arrived = modular.apply(law, initial, first_at=n)
                        for bit in (0, 1):
                            raw = bit if flag == 1 else 1-bit
                            output = modular.y_read(arrived, raw)
                            probability = modular.trace(output)
                            assert probability >= 0
                            plus += probability*int(bit == 0)
                            rows.append({'record': [flag, n, bit], 'probability': str(probability),
                                         'terminal_state_sha256': state_digest(output)})
                    tail = modular.apply(law, initial, no_current=cap)
                    rows.append({'record': [flag, cap+1, -1], 'probability': str(modular.trace(tail)),
                                 'terminal_state_sha256': state_digest(tail)})
                assert len(rows) == 4*cap+3 and sum(F(row['probability']) for row in rows) == 1
                infinite = F(5, 33)*F(harmonic['exact_quarter_turn_conditional_plus_probabilities'][str(quarter)])
                tail_probability = sum(F(row['probability']) for row in rows if row['record'][0] and row['record'][2] == -1)
                assert 0 <= infinite-plus <= tail_probability <= F(5, 33)*F(13, 15)**cap
                finite.append({'quarter_turn': quarter, 'attempt_cap_after_herald': cap,
                    'complete_record_count': len(rows), 'records': rows,
                    'complete_original_plus_probability': str(plus),
                    'original_infinite_event_plus_probability': str(infinite),
                    'original_timeout_probability': str(tail_probability),
                    'uniform_original_tail_upper': str(F(5, 33)*F(13, 15)**cap)})
        runtime = modular.covariance_experiment(law, 'phase_before_event', cap=4)
        reference = next(row for row in finite if row['quarter_turn'] == 1 and row['attempt_cap_after_herald'] == 4)
        # Compare by record identity; the producer places both timeouts last.
        assert {tuple(row['record']): row['probability'] for row in runtime['complete_records']} == {
            tuple(row['record']): row['probability'] for row in reference['records']}
        cases.append({'law': law, 'harmonic': harmonic, 'finite_complete_records': finite,
                      'independent_primitive_runtime_quarter_turn': runtime})
    return {'dashboard_row_index': 5, 'theorem': 'complete native first-event modular Y instrument',
        'carrier_dimension': 66, 'full_detector_Kraus_identity_exact': True,
        'all4356_matrix_units_of_stopped_map_trace_identity_exact': True,
        'all_untouched_references_admitted': True,
        'preparation_failure_and_both_herald_signs_retained': True, 'cases': cases,
        'finite_law_premises': ['declared native C/E instrument and uniform scheduler',
            'declared coherent preparation, modular phase and complete Y detector services'],
        'native_pulse_or_hardware_error_closed': False,
        'native_attempts_generate_modular_phase': False, 'physical_time_identified': False,
        'recommended_status': 'Theorem + executed',
        'recommended_scope': 'finite declared quantum instrument; ideal phase target and finite timeout envelope'}


@lru_cache(maxsize=1)
def joint_clocks():
    rows = []
    for law in ('C', 'E'):
        groups = clocks.source_groups(law)
        identity = sum((s.Rational(weight)*s.Matrix(k).H*s.Matrix(k)
                        for group in groups for weight, k in group), s.zeros(4))
        assert identity == s.eye(4)
        assert all(weight >= 0 for group in groups for weight, _ in group)
        rows.append({'law': law, 'all_input_native_Kraus_completeness_exact': True,
            'passive_clock_conversion': clock_change.certificate(law),
            'feedback_clock_conversion': clock_change.certificate(law, True),
            'complete_three_attempt_runtime': clocks.execute_experiment(law),
            'complete_reference_feedback_runtime': clocks.execute_experiment(law, feedback=True, reference=True, cap=2)})
    assert [r['complete_three_attempt_runtime']['complete_field_plus_probability'] for r in rows] == ['35/54', '14/27']
    counts = clock_change.relative_count_prediction()
    return {'dashboard_row_index': 7, 'theorem': 'joint marked-clock instrument and remote-count generating law',
        'full_work_operator_dimension': 16, 'cases': rows, 'complete_relative_count_law': counts,
        'finite_law_premises': ['actual four-port native support', 'independent three-edge scheduler',
            'declared coherent preparation, readout and record-conditioned feedback services'],
        'arbitrary_reference_extension_from_full_map_identity': True,
        'compiled_pulse_or_hardware_error_closed': False, 'physical_affine_clock_map_proved': False,
        'recommended_status': 'Theorem + executed',
        'recommended_scope': 'finite stochastic clock comparison with complete conditional quantum outputs'}


@lru_cache(maxsize=1)
def cyclic_transfer():
    source, kernel = temporal.core.source(), temporal.core.completed()
    temporal.core.dictionary_certificate()
    spectral = temporal.core.spectral_certificate()
    p = s.Rational(source.tick_probability)
    d = s.eye(66)-s.Matrix(source.remote_counts)/source.denominator
    b = s.Matrix(source.incident_counts)/source.denominator
    k = s.Matrix(kernel.counts)/kernel.denominator
    metric = p*d.inv(method='DM')/66
    assert k*d == b and metric*k == k.T*metric
    f = s.Matrix([s.Rational(state[0])-s.Rational(1, 6) for state in source.states])
    assert d*f == p*f
    moments = temporal.source_correlations(20)
    representation = temporal.recover(tuple(row[0][0]-F(1, 36) for row in moments[:15]))
    columns, value = [], f
    for _ in range(7):
        columns.append(value)
        value = k*value
    inclusion = s.Matrix.hstack(*columns)
    gram = s.Matrix(representation.gram)
    transfer = s.Matrix(representation.transfer)
    assert inclusion.T*metric*inclusion == gram
    assert k*inclusion == inclusion*transfer
    assert temporal._positive(gram)
    assert temporal._positive(gram*transfer) and temporal._positive(gram-gram*transfer)
    annihilator = sum((s.Rational(coefficient)*column for coefficient, column
                      in zip(representation.polynomial, columns+[value], strict=True)), s.zeros(66, 1))
    assert annihilator == s.zeros(66, 1)
    # Full-source powers beyond the inference window are independent of the
    # seven-dimensional recurrence. The inclusion proves equality at every lag.
    withheld = []
    value = f
    for lag in range(33):
        if lag >= 15:
            actual = (f.T*value)[0]/66
            predicted = s.Rational(temporal.predict(representation, lag))
            assert actual == predicted
            withheld.append({'lag': lag, 'full_source_centered_correlation': str(actual),
                             'record_only_prediction': str(predicted)})
        value = k*value
    return {'dashboard_row_index': 17, 'theorem': 'exact cyclic inclusion and all-lag recurrence of the native field records',
        'full_native_states': 66, 'cyclic_dimension': 7, 'inference_lags': [0, 14],
        'source_positive_transfer_bounds': {key: str(value) if isinstance(value, F) else value for key, value in spectral.items()},
        'source_metric': 'p (I-A)^(-1)/66, p=3/10',
        'full_source_inclusion': matrix_record(inclusion), 'record_Gram': matrix_record(gram),
        'record_transfer': matrix_record(transfer),
        'full_source_intertwiner_KC_equals_CR': True, 'full_source_annihilator_vector_zero': True,
        'positive_definite_Gram_and_strict_positive_contraction_exact': True,
        'recurrence_coefficients_ascending': list(map(str, representation.polynomial)),
        'independent_unused_lags': withheld,
        'finite_law_premises': ['uniform total-two preparation and independent native scheduler',
            'initial local bits and signed incident current records', 'exact model correlation laws'],
        'sampled_data_certifies_exact_rank': False, 'physical_Hamiltonian_or_time_identified': False,
        'recommended_status': 'Theorem + executed',
        'recommended_scope': 'finite native-law reconstruction with an exact all-lag inclusion witness'}


def bindings():
    seeds = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    seeds.update(modular.SOURCE_PATHS)
    seeds.update(clocks.SOURCE_PATHS)
    seeds.update(temporal.SOURCE_PATHS)
    seeds.update(('model/relative_modular_dynamics/README.md', 'model/joint_clock_process/README.md',
                  'model/relative_observer_geometry/README.md'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    pins = bindings()
    result = {'schema': 'oph.finite_observer_readout_closure.v1', 'source_bindings': pins,
        'periodic_Y': periodic_y(), 'joint_event_clocks': joint_clocks(), 'unused_lag_correlations': cyclic_transfer(),
        'scope': {'new_empirical_observations': 0, 'physical_time_or_Universe_attachment_closed': False,
                  'dashboard_status_changes_applied': False}}
    if bindings() != pins:
        raise ValueError('finite observable sources changed during replay')
    return result


def write():
    value = certificate()
    (HERE/'results.json').write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    return value
