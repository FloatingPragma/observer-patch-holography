"""Exact quantum-sector reduction and complete cross-carrier readouts.

Native attempts are discrete instrument operations, not modular or physical
durations. The graph, schedule, coherent preparation and analyzer are declared.
All probability oracles use rational arithmetic; runtime residuals are diagnostics.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
import hashlib
from itertools import product
import math
from pathlib import Path

import numpy as np

from ..federated_observer import source
from ..observer import native
from ..observer.scalar_response.experiment import canonical
from ..primitive_observer import runtime

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*source.SOURCE_PATHS,
    'model/federated_quantum_sector/__init__.py',
    'model/federated_quantum_sector/__main__.py',
    'model/federated_quantum_sector/construction.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/observer/quantum_mixing.py',
    'model/native_quantum_history/construction.py',
    'Lean/EventAlgebra/SlotLocalitySurface.lean',
    'Lean/ObserverPatchHolography/Locality/DependencyCone.lean',
    'Lean/SourceLaw/ReflectionPositivity.lean'))))
CURRENT = (-1, 0, 1)


def candidate(law):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('native candidate C or E required')
    return law


def _mask(value):
    if type(value) is not int or not 0 <= value < 1 << 4096:
        raise ValueError('a nonnegative binary configuration of at most 4096 ports required')
    return value


def _seam(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(v) is not int for v in value)
            or not 0 <= value[0] < value[1] < 4096):
        raise ValueError('two distinct ordered binary port addresses required')
    return value


def native_image(mask, seam, coin):
    """Call the original endpoint law and retain both candidate record labels."""
    mask, (a, b) = _mask(mask), _seam(seam)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError('an exact native placement coin required')
    values = ((mask >> a) & 1, (mask >> b) & 1)
    target, (_, current) = native.native(values, (0, 1), coin)
    result = mask & ~((1 << a) | (1 << b))
    result |= target[0] << a | target[1] << b
    return result, current, values[0]-values[1]


def native_unit(law, left, right, seam, current=None):
    """One outcome's exact action on |left><right|, with the coin forgotten.

    C groups by coin/current; E additionally resolves endpoint imbalance.
    Visible current is retained in both cases. Every hidden alternative is
    summed as a CP map, never added as a coherent amplitude.
    """
    candidate(law); _mask(left); _mask(right); _seam(seam)
    if current is not None and (type(current) is not int or current not in CURRENT):
        raise ValueError('a signed current or the complete unread instrument required')
    result = defaultdict(F)
    for coin in (0, 1):
        a = native_image(left, seam, coin)
        b = native_image(right, seam, coin)
        label = 1 if law == 'C' else 2
        if a[label] == b[label] and (current is None or a[1] == current):
            if a[1] != b[1]:
                raise ArithmeticError('one hidden Kraus label gives incompatible visible currents')
            result[a[0], b[0]] += F(1, 2)
    return dict(result)


def difference(left, right, ports):
    _mask(left); _mask(right)
    if type(ports) is not int or not 1 <= ports <= 4096 or max(left, right) >= 1 << ports:
        raise ValueError('a bounded binary carrier containing both configurations required')
    return tuple(((left >> p) & 1)-((right >> p) & 1) for p in range(ports))


def block_killing(graph, signed, law):
    """The scalar loss of a signed-difference block on the source federation."""
    candidate(law)
    if (type(graph) is not source.Federation or type(signed) is not tuple
            or len(signed) != graph.ports
            or any(type(v) is not int or v not in (-1, 0, 1) for v in signed)):
        raise ValueError('a source federation and its complete ternary signed difference required')
    if law == 'C':
        return F(sum(abs(signed[a]-signed[b]) for a, b in graph.seams), 2*len(graph.seams))
    return F(sum(signed[a] != signed[b] for a, b in graph.seams), len(graph.seams))


def sector_column(graph, left, right, law):
    """Independent killed-exclusion formula, without calling native repair."""
    signed = difference(left, right, graph.ports)
    result = defaultdict(F, {(left, right): 1-block_killing(graph, signed, law)})
    unit = F(1, 2*len(graph.seams))
    for a, b in graph.seams:
        if signed[a] == signed[b] == 0 and ((left >> a) & 1) != ((left >> b) & 1):
            toggle = (1 << a) | (1 << b)
            result[left, right] -= unit
            result[left ^ toggle, right ^ toggle] += unit
    return {key: value for key, value in result.items() if value}


def source_column(graph, left, right, law):
    difference(left, right, graph.ports)
    result = defaultdict(F)
    for seam in graph.seams:
        for key, value in native_unit(law, left, right, seam).items():
            result[key] += value/len(graph.seams)
    return dict(result)


def local_certificate():
    """All matrix units and all visible outcomes, with no state sampling."""
    rows = []
    for law, left, right in product(('C', 'E'), range(4), range(4)):
        signed = difference(left, right, 2)
        complete = defaultdict(F)
        for current in CURRENT:
            output = native_unit(law, left, right, (0, 1), current)
            for (a, b), value in output.items():
                assert difference(a, b, 2) == signed
                assert a.bit_count() == left.bit_count() and b.bit_count() == right.bit_count()
                complete[a, b] += value
            rows.append([law, left, right, current,
                [[a, b, str(v)] for (a, b), v in sorted(output.items())]])
        loss = F(abs(signed[0]-signed[1]), 2) if law == 'C' else F(signed[0] != signed[1])
        expected = defaultdict(F, {(left, right): 1-loss})
        if signed == (0, 0) and left in (1, 2):
            expected[left, right] -= F(1, 2)
            expected[3-left, 3-right] += F(1, 2)
        assert dict(complete) == {key: v for key, v in expected.items() if v}
        # Equality is valid after tensoring with every untouched reference.
    for law in ('C', 'E'):
        trace = defaultdict(F); unit = defaultdict(F)
        for i, j in product(range(4), repeat=2):
            for (a, b), value in native_unit(law, i, j, (0, 1)).items():
                assert native_unit(law, a, b, (0, 1)).get((i, j), F(0)) == value
                if a == b:
                    trace[i, j] += value
                if i == j:
                    unit[a, b] += value
        expected = {(i, i): F(1) for i in range(4)}
        assert dict(trace) == dict(unit) == expected
    return {'endpoint_matrix_units_per_candidate': 16, 'visible_current_outcomes': 3,
        'all_outcome_matrix_units_checked': len(rows),
        'complete_exact_local_map_sha256': hashlib.sha256(canonical(rows).encode()).hexdigest(),
        'trace_preserving_unital_and_number_preserving': True,
        'unread_channel_self_adjoint_on_operator_space': True,
        'each_visible_outcome_preserves_signed_difference': True,
        'every_untouched_reference_extension': True}


def exchange_rate(level, first, second, law):
    graph = source.federation(level)
    candidate(law)
    if (type(first) is not int or type(second) is not int
            or not 0 <= first < second < graph.ports):
        raise ValueError('two distinct ordered source port addresses required')
    degree_sum = len(graph.neighbours[first])+len(graph.neighbours[second])
    loss = (F(degree_sum, 2*len(graph.seams)) if law == 'C' else
            F(degree_sum-int(second in graph.neighbours[first]), len(graph.seams)))
    return 1-loss


def phase_interval(rate, count, digits=12):
    """Outward rational interval for (1+rate^count)/2 using integer powers."""
    if type(rate) is not F or not 0 <= rate <= 1 or type(count) is not int or not 0 <= count <= 20000:
        raise ValueError('a rational contraction and bounded integer attempt count required')
    if type(digits) is not int or not 0 <= digits <= 15:
        raise ValueError('a bounded exact decimal reporting grid required')
    scale = 10**digits
    a, b = rate.numerator**count, rate.denominator**count
    numerator, denominator = scale*(a+b), 2*b
    lower = numerator//denominator
    upper = -(-numerator//denominator)
    return {'lower': str(F(lower, scale)), 'upper': str(F(upper, scale))}


def graph_predictions():
    rows = []
    for level in (0, 1, 2):
        graph = source.federation(level)
        a, b = graph.inter_seams[0]
        signed = difference(1 << a, 1 << b, graph.ports)
        # The source and reduction are independently constructed on this
        # nonclassical matrix unit and on common-particle configurations.
        spectator = next(p for p in range(graph.ports) if p not in (a, b))
        for law in ('C', 'E'):
            rate = exchange_rate(level, a, b, law)
            assert rate == 1-block_killing(graph, signed, law)
            for common in (0, 1 << spectator):
                assert source_column(graph, (1 << a)|common, (1 << b)|common, law) == sector_column(
                    graph, (1 << a)|common, (1 << b)|common, law)
            count = len(graph.seams)
            rows.append({'level': level, 'ports': graph.ports, 'seams': count,
                'cross_carrier_field_ports': [a, b], 'candidate': law,
                'coherence_retention_per_uniform_attempt': str(rate),
                'declared_native_attempts': count,
                'plus_probability_exact_expression': f'(1+({rate})^{count})/2',
                'plus_probability_rational_interval': phase_interval(rate, count),
                'uniform_global_scheduler_not_four_seam_tape': True})
    return rows


def _zero(size):
    return [[F(0) for _ in range(size)] for _ in range(size)]


def _identity(size):
    return [[F(i == j) for j in range(size)] for i in range(size)]


def _exact_apply(groups, density, reference):
    """Rational full matrix contraction; work is the major tensor factor."""
    dim = len(density)
    output = []
    for group in groups:
        result = _zero(dim)
        for weight, matrix in group:
            entries = [(i, j, value) for i, row in enumerate(matrix)
                       for j, value in enumerate(row) if value]
            for i, a, left in entries:
                for j, b, right in entries:
                    for r, s in product(range(reference), repeat=2):
                        result[i*reference+r][j*reference+s] += weight*left*right*density[a*reference+r][b*reference+s]
        output.append(result)
    return output


def _runtime_instrument(groups, *, visible=False, resources=(), binding=None):
    return runtime.Instrument(kraus=tuple(tuple(np.asarray(matrix, float)*math.sqrt(float(weight))
        for weight, matrix in group) for group in groups), visible=visible,
        resources=resources, binding=binding)


def native_groups(law, vertices, seam):
    """Exact finite invariant support of the full tensor-local instrument."""
    candidate(law); _seam(seam)
    if (type(vertices) is not tuple or tuple(sorted(set(vertices))) != vertices
            or any(type(v) is not int or not 0 <= v < 4096 for v in vertices)
            or not 2 <= len(vertices) <= 8 or not set(seam) <= set(vertices)):
        raise ValueError('two through eight active source ports containing the selected seam required')
    maps = {}
    for column, vertex in enumerate(vertices):
        for coin in (0, 1):
            target, current, imbalance = native_image(1 << vertex, seam, coin)
            if target.bit_count() != 1 or target.bit_length()-1 not in vertices:
                raise ArithmeticError('the source tape leaks from its exact one-particle support')
            key = current, coin, current if law == 'C' else imbalance
            maps.setdefault(key, _zero(len(vertices)))[vertices.index(target.bit_length()-1)][column] = F(1)
    # Explicit zero Kraus maps retain impossible visible outcomes.
    return tuple(tuple((F(1, 2), matrix) for key, matrix in sorted(maps.items()) if key[0] == current)
                 or ((F(1), _zero(len(vertices))),) for current in CURRENT)


def apparatus(vertices, pair, mode):
    if type(mode) is not str or mode not in ('phase', 'dephase', 'identity'):
        raise ValueError('a declared coherent intervention required')
    a, b = (vertices.index(p) for p in pair)
    size = len(vertices)
    if mode == 'dephase':
        first = _zero(size); first[a][a] = F(1)
        second = _identity(size); second[a][a] = F(0)
        return (((F(1), first), (F(1), second)),)
    matrix = _identity(size)
    if mode == 'phase':
        matrix[a][a] = F(-1)
    return (((F(1), matrix),),)


def analyzer(vertices, pair):
    """Complete exchange analyzer, with a retained fair outcome off its pair.

    First resolve pair versus complement, then measure X on the pair and
    use a fair coin on the complement. This fixes the postmeasurement map.
    """
    a, b = (vertices.index(p) for p in pair)
    size = len(vertices)
    outside = _identity(size); outside[a][a] = outside[b][b] = F(0)
    result = []
    for sign in (1, -1):
        projector = _zero(size)
        projector[a][a] = projector[b][b] = F(1, 2)
        projector[a][b] = projector[b][a] = F(sign, 2)
        result.append(((F(1), projector), (F(1, 2), outside)))
    return tuple(result)


def experiment(law, sign=1, mode='identity', reference=False):
    candidate(law)
    if type(sign) is not int or sign not in (-1, 1) or type(reference) is not bool:
        raise ValueError('a phase sign and Boolean reference selection required')
    graph, vertices, original, _, _ = source.tape_data()
    # c,d are the two internal endpoints, not the selected crossing endpoints.
    a, b = original[0]
    pair = tuple(v for v in vertices if v not in (a, b))
    tape = (original[1], original[0], original[2], original[0])
    rd = 2 if reference else 1
    size = len(vertices)*rd
    initial = _zero(size)
    i, j = vertices.index(pair[0])*rd, vertices.index(pair[1])*rd+int(reference)
    initial[i][i] = initial[j][j] = F(1, 2)
    initial[i][j] = initial[j][i] = F(sign, 2)
    native_maps = [native_groups(law, vertices, seam) for seam in tape]
    intervention = apparatus(vertices, pair, mode)
    identity = apparatus(vertices, pair, 'identity')
    read = analyzer(vertices, pair)
    events = []
    for n, groups in enumerate(native_maps):
        events.append(runtime.Event(f'native_{n}', instruments=(_runtime_instrument(groups,
            visible=True, resources=(('native_attempts', 1), ('current_records', 1)),
            binding={'candidate': law, 'source_seam': list(tape[n]), 'visible_current': list(CURRENT)}),)))
        if n == 0:
            # Apply the declared intervention exactly when the first current is zero.
            events.append(runtime.Event('record_controlled_intervention', control=0,
                instruments=tuple(_runtime_instrument(g, resources=(('declared_interventions', 1),),
                    binding={'mode': mode if current == 0 else 'identity'})
                    for current, g in zip(CURRENT, (identity, intervention, identity), strict=True))))
    events.append(runtime.Event('exchange_read', instruments=(_runtime_instrument(read,
        visible=True, resources=(('exchange_readouts', 1),),
        binding={'outcomes': ['plus', 'minus'], 'pair': list(pair),
                 'instrument': 'pair projection followed by X read; fair complete complement'}),)))
    default = [[sum((initial[i*rd+r][j*rd+r] for r in range(rd)), F(0))
                for j in range(len(vertices))] for i in range(len(vertices))]
    program = runtime.Program(f'federated_coherent_{law}_{mode}_{sign}_{rd}', len(vertices),
        np.asarray(default, float), tuple(events),
        preparation_resources=(('declared_coherent_preparations', 1),),
        binding={'graph_level': 0, 'global_ports': graph.ports, 'active_one_particle_ports': list(vertices),
            'all_other_global_ports_vacuum': True, 'fixed_four_seam_tape': [list(e) for e in tape],
            'arbitrary_reference_support_embedding': True})
    exact = {(): initial}
    for n, groups in enumerate(native_maps):
        next_branches = {}
        for word, density in exact.items():
            for outcome, matrix in enumerate(_exact_apply(groups, density, rd)):
                if n == 0:
                    matrix = _exact_apply(intervention if CURRENT[outcome] == 0 else identity, matrix, rd)[0]
                next_branches[word+(outcome,)] = matrix
        exact = next_branches
    terminal = {word+(outcome,): matrix for word, density in exact.items()
                for outcome, matrix in enumerate(_exact_apply(read, density, rd))}
    return program, terminal, initial


def execute_experiment(law, sign=1, mode='identity', reference=False):
    program, expected, initial = experiment(law, sign, mode, reference)
    supplied = np.asarray(initial, float) if reference else None
    result = runtime.execute(program, initial=supplied)
    prefix = runtime.execute(program, initial=supplied, stop_after=3)
    restored = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(prefix))
    continued = runtime.execute(program, checkpoint=restored)
    residual = runtime.complete_distance(result, continued)
    exact_mass = F(0); plus = F(0); rows = []
    current_law = defaultdict(F)
    for branch in result.branches:
        exact = expected[branch.record]
        mass = sum((exact[i][i] for i in range(len(exact))), F(0))
        exact_mass += mass
        current_law[branch.record[:-1]] += mass
        if branch.record[-1] == 0:
            plus += mass
        residual = max(residual, float(np.linalg.norm(branch.density-np.asarray(exact, float), 'fro')))
        expected_resources = {'native_attempts': 4, 'current_records': 4,
                              'declared_interventions': 1, 'exchange_readouts': 1}
        if not reference:
            expected_resources['declared_coherent_preparations'] = 1
        assert dict(branch.resources) == expected_resources
        rows.append({'word': list(branch.record), 'probability': str(mass)})
    assert exact_mass == 1 and len(result.branches) == 162 and residual < 2e-11
    encoded = [[list(word), [[str(v) for v in row] for row in matrix]]
               for word, matrix in sorted(expected.items())]
    return {'candidate': law, 'preparation_sign': sign, 'intervention': mode,
        'untouched_reference_dimension': result.reference_dimension,
        'preparation_resource_scope': ('externally supplied joint reference density; preparation cost unpriced'
            if reference else 'one declared coherent preparation'),
        'complete_record_words': len(rows), 'positive_words': sum(row['probability'] != '0' for row in rows),
        'exact_plus_probability': str(plus), 'exact_total_probability': '1',
        'complete_record_law': rows,
        'all_terminal_matrix_entries_sha256': hashlib.sha256(canonical(encoded).encode()).hexdigest(),
        'current_only_law_sha256': hashlib.sha256(canonical(
            [[list(word), str(mass)] for word, mass in sorted(current_law.items())]).encode()).hexdigest(),
        'complete_state_and_checkpoint_residual_diagnostic': residual,
        'source_pulse_words_or_physical_service_error_certified': False}


def run():
    certificate = local_certificate()
    predictions = graph_predictions()
    executed = [execute_experiment(law, sign, mode, reference)
                for law in ('C', 'E')
                for sign, mode, reference in ((1, 'identity', False), (-1, 'identity', False),
                    (1, 'phase', False), (1, 'dephase', False), (1, 'identity', True))]
    for law in ('C', 'E'):
        rows = [row for row in executed if row['candidate'] == law]
        assert len({row['current_only_law_sha256'] for row in rows}) == 1
        assert rows[3]['exact_plus_probability'] == rows[4]['exact_plus_probability'] == '1/2'
    return {'schema': 'oph.federated_quantum_sector.v1',
        'local_all_input_certificate': certificate,
        'uniform_scheduler_exchange_field_predictions': predictions,
        'complete_cross_carrier_experiments': executed,
        'scope': {'all_binary_number_sectors_and_arbitrary_references': True,
            'native_records_alone_detect_initial_coherence': False,
            'declared_coherent_interventions_and_analyzer': True,
            'new_relativistic_action_or_physical_clock_inserted': False,
            'continuum_relativistic_field_or_empirical_match_established': False},
        'source_bindings': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}}
