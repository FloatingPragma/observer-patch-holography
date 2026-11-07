"""Exact low-degree effects and integer native warmup budgets.

The bound concerns a signed record score after ignoring preparation records.
It does not assert trace-distance preparation of the entire source density.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, isqrt

from ..native_correlation_sectors import sectors as s, records as r
from ..native_equilibration import construction as mixing


def polynomial(values, ports):
    """Validate an exact squarefree polynomial of degree at most four."""
    s.integer(ports, 4)
    if type(values) is not dict or len(values) > 100000:
        raise ValueError('a bounded exact coefficient dictionary required')
    result = {}
    for subset, value in values.items():
        if (type(subset) is not tuple or len(subset) > 4
                or any(type(i) is not int or not 0 <= i < ports for i in subset)
                or tuple(sorted(set(subset))) != subset
                or type(value) not in (int, F)):
            raise ValueError('canonical subsets and exact rational coefficients required')
        if value:
            result[subset] = F(value)
    return result


def evaluate(values, occupied, ports):
    values = polynomial(values, ports)
    if (type(occupied) is not tuple or any(type(i) is not int or not 0 <= i < ports for i in occupied)
            or tuple(sorted(set(occupied))) != occupied):
        raise ValueError('an increasing original occupation tuple required')
    selected = set(occupied)
    return sum((value for subset, value in values.items() if set(subset) <= selected), F())


def multiply(left, right, ports):
    left = polynomial(left, ports); right = polynomial(right, ports)
    if len(left)*len(right) > 1000000:
        raise ValueError('bounded exact polynomial multiplication required')
    result = defaultdict(F)
    for a, x in left.items():
        for b, y in right.items():
            merged = tuple(sorted(set(a) | set(b)))
            if len(merged) > 4:
                raise ValueError('product exceeds the proved fourth-degree domain')
            result[merged] += x*y
    return polynomial(dict(result), ports)


def seam_average(values, edge, ports):
    values = polynomial(values, ports)
    if (type(edge) is not tuple or len(edge) != 2
            or any(type(i) is not int for i in edge) or not 0 <= edge[0] < edge[1] < ports):
        raise ValueError('an original ordered seam required')
    a, b = edge; result = defaultdict(F)
    for subset, value in values.items():
        target = tuple(sorted(b if i == a else a if i == b else i for i in subset))
        result[subset] += value/2; result[target] += value/2
    return polynomial(dict(result), ports)


def evolve(family, values, steps):
    """Bounded exact sparse dual evolution; no expanded-history claim."""
    s.fair_graph(family); s.integer(steps)
    if steps > 12:
        raise ValueError('at most twelve exact sparse evolution steps required')
    ports = len(family.weights); values = polynomial(values, ports)
    for _ in range(steps):
        result = {}
        for degree in range(5):
            part = {a: v for a, v in values.items() if len(a) == degree}
            if part:
                result.update({a: v/F(2*len(family.seams))
                               for a, v in s.counts_step(family, part, degree).items()})
        values = polynomial(result, ports)
    return values


def score_effect(family, first, second, lag):
    """Pull back the complete reflected four-E score to its input cut."""
    return cross_effect(family, (first, second), (first, second), lag)


def cross_effect(family, past, future, lag):
    """The same complete pullback for distinct past and future seam pairs."""
    s.fair_graph(family); ports = len(family.weights)
    if type(past) is not tuple or type(future) is not tuple or len(past) != 2 or len(future) != 2:
        raise ValueError('two original pairs of probe seams required')
    phi = s.contrast(*past, ports)
    psi = s.contrast(*future, ports)
    if any(edge not in family.seams for edge in (*past, *future)):
        raise ValueError('all probe seams must belong to the native graph')
    evolved = evolve(family, psi, lag)
    value = multiply(phi, evolved, ports)
    # The disjoint fixed-seam averages commute; both measurement backactions
    # are included before the initial input cut.
    return seam_average(seam_average(value, past[1], ports), past[0], ports)


def uniform_mean(values, ports, occupied):
    s.degree_data(ports, occupied); values = polynomial(values, ports)
    result = F()
    for subset, value in values.items():
        degree = len(subset)
        if degree <= occupied:
            probability = F(1)
            for i in range(degree):
                probability *= F(occupied-i, ports-i)
            result += value*probability
    return result


def polynomial_moments(values, ports, occupied):
    """Exact uniform-slice mean and centered L2 norm from sparse coefficients."""
    s.degree_data(ports, occupied); values = polynomial(values, ports)
    if len(values)**2 > 1000000:
        raise ValueError('at most one million exact Gram entries evaluated')
    mean = uniform_mean(values, ports, occupied); second = F()
    for a, left in values.items():
        for b, right in values.items():
            degree = len(set(a) | set(b)); probability = F(1)
            if degree > occupied:
                continue
            for i in range(degree):
                probability *= F(occupied-i, ports-i)
            second += left*right*probability
    centered = second-mean*mean
    if centered < 0:
        raise ArithmeticError('exact polynomial Gram norm is negative')
    return {'mean': mean, 'second_moment': second, 'centered_l2_squared': centered}


def snapshot_effect(family, first, second, lag):
    """Retrospective snapshot effect; its native readout is a separate premise."""
    s.fair_graph(family); ports = len(family.weights)
    first = polynomial(first, ports); second = polynomial(second, ports)
    if max((len(a) for a in (*first, *second)), default=0) > 2:
        raise ValueError('two fields of degree at most two required')
    return multiply(first, evolve(family, second, lag), ports)


def variance(ports, occupied):
    s.integer(ports, 4); s.degree_data(ports, occupied)
    if not 2 <= occupied <= ports-2:
        raise ValueError('nonzero contrast variance requires 2 <= q <= V-2')
    return F(4*occupied*(occupied-1)*(ports-occupied)*(ports-occupied-1),
             ports*(ports-1)*(ports-2)*(ports-3))


def ceiling_sqrt(value):
    if type(value) is not F or value < 0:
        raise ValueError('an exact nonnegative rational required')
    answer = isqrt(value.numerator//value.denominator)
    return answer + int(answer*answer*value.denominator < value.numerator)


def polynomial_budget(ports, occupied, gap_lower, centered_l2_squared_upper,
                      normalization, precision_bits=10):
    """Bound a specified degree-four effect with proved exact norm inputs.

    Supplying a bound is not a proof of that bound. polynomial_moments derives
    it for bounded sparse inputs; separate analytic bounds cover larger ones.
    """
    s.integer(ports, 4); s.degree_data(ports, occupied)
    if (type(centered_l2_squared_upper) is not F or centered_l2_squared_upper < 0
            or type(normalization) is not F or normalization <= 0):
        raise ValueError('exact nonnegative norm bound and positive normalization required')
    if type(gap_lower) is not F or not 0 < gap_lower <= 1:
        raise ValueError('a proved exact positive contraction gap at most one required')
    if type(precision_bits) is not int or not 1 <= precision_bits <= 1000:
        raise ValueError('one through 1000 original precision bits required')
    degree = min(4, occupied, ports-occupied)
    rank = comb(ports, degree)-1
    prefactor_squared = rank*centered_l2_squared_upper/(normalization*normalization)
    prefactor = ceiling_sqrt(prefactor_squared)
    blocks = (prefactor-1).bit_length()+precision_bits if prefactor else 0
    block_length = (gap_lower.denominator+gap_lower.numerator-1)//gap_lower.numerator
    error = F(prefactor, 1 << blocks)
    if (prefactor*prefactor < prefactor_squared or block_length*gap_lower < 1
            or error > F(1, 1 << precision_bits)):
        raise ArithmeticError('integer observable preparation certificate failed')
    return {'ports': ports, 'occupation': occupied, 'correlation_degree_cutoff': degree,
            'mean_zero_polynomial_rank': rank,
            'centered_l2_squared_upper': str(centered_l2_squared_upper),
            'normalization': str(normalization),
            'gap_lower_per_attempt': str(gap_lower),
            'normalized_prefactor_squared': str(prefactor_squared),
            'normalized_prefactor_integer_upper': prefactor,
            'half_contraction_block_attempts': block_length, 'halving_blocks': blocks,
            'warmup_attempts': block_length*blocks,
            'normalized_score_error_upper': str(error),
            'unnormalized_reference_score_trace_norm_upper': str(error*normalization)}


def preparation_budget(ports, occupied, gap_lower, precision_bits=10):
    """Sufficient integer count for the normalized four-E score error."""
    v = variance(ports, occupied)
    return {**polynomial_budget(ports, occupied, gap_lower, v, v, precision_bits),
            'stationary_variance': str(v)}


def prepared_program(level, occupied, warmup, lag):
    s.integer(level); s.integer(warmup); s.integer(lag)
    if level > 2:
        raise ValueError('this attachment evaluates native federation levels zero through two')
    if warmup+lag+4 > s.source.interpreter.MAX_ATTEMPTS:
        raise ValueError('warmup plus complete probe exceeds the native attempt ceiling')
    base = r.probe_program(level, occupied, lag)
    prefix = ([{'op': 'repair', 'selector': {'kind': 'uniform'},
                'law': 'mixture', 'repeat': warmup}] if warmup else [])
    return s.source.interpreter.program(federation={'level': level, 'total': occupied},
                                        commands=prefix+base['commands'], overflow='export_oldest')


def source_budgets(level, precision_bits=10, lag=3):
    s.integer(level); s.integer(lag)
    if level > 2 or lag > s.source.interpreter.MAX_ATTEMPTS-4:
        raise ValueError('evaluated federation level and bounded complete probe required')
    family, first, second, _ = s.probe_data(level)
    graph = mixing.graph_certificate(len(family.weights), family.seams)
    gap = F(graph['quantum_gap_lower_per_attempt'])
    rows = []
    for occupied in (2, len(family.weights)//3, len(family.weights)//2):
        budget = preparation_budget(len(family.weights), occupied, gap, precision_bits)
        attempts = budget['warmup_attempts']+lag+4
        admissible = attempts <= s.source.interpreter.MAX_ATTEMPTS
        program = prepared_program(level, occupied, budget['warmup_attempts'], lag) if admissible else None
        # The comparison is with the separately proved global-state bound
        # d*r^n/v for this same scalar score. Only logarithmic counts are kept.
        dimension = comb(len(family.weights), occupied)
        global_prefactor = (dimension*F(1)/variance(len(family.weights), occupied))
        global_integer = (global_prefactor.numerator+global_prefactor.denominator-1)//global_prefactor.denominator
        global_blocks = (global_integer-1).bit_length()+precision_bits
        rows.append({**budget, 'lag_attempts': lag, 'fixed_E_read_attempts': 4,
                     'complete_attempts': attempts, 'single_program_admitted': admissible,
                     'program_digest': s.source.interpreter.digest(program) if program else None,
                     'global_state_comparison_halving_blocks': global_blocks,
                     'global_state_comparison_warmup_attempts': global_blocks*budget['half_contraction_block_attempts'],
                     'full_expanded_history_executed': False})
    return {'level': level, 'graph_certificate': graph,
            'probe_seams': [list(first), list(second)], 'budgets': rows}


def source_effect_control(warmup=2, lag=2, law='E'):
    """Direct complete H66 source adjoints, including every raw zero slot.

    Full matrix-unit effects are diagonal because every raw map is a partial
    injection. The separate explicit small-matrix test checks this assertion.
    """
    s.integer(warmup); s.integer(lag)
    if warmup > 4 or lag > 4 or law not in ('C', 'E', 'mixture'):
        raise ValueError('bounded native H66 oracle parameters required')
    interpreter = s.source.interpreter
    family = s.source.local_family(s.source.local.Parameters(port_weights=(1,)*12))
    first, second = next((a, b) for a in family.seams for b in family.seams if not set(a) & set(b))
    pair = [{'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(e)}, 'law': 'E', 'repeat': 1}
            for e in (first, second)]
    uniform = {'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': law, 'repeat': 1}
    commands = [uniform]*warmup+pair+[uniform]*lag+pair
    program = interpreter.program(biased=False, commands=commands, overflow='export_oldest')
    states = s.source.local.basis()
    marker = {(i, i, i, i): F(1) for i in range(66)}
    effect = [F(1)]*66; branch_count = 0; zero_count = 0
    for position in reversed(range(len(commands))):
        updated = [F()]*66
        outcomes = interpreter._outcomes(program, commands[position], marker)
        for outcome in outcomes:
            branch_count += 1; zero_count += int(not outcome.density)
            score = 1
            if warmup <= position < warmup+2:
                score = (outcome.record.label-2*outcome.record.current
                         if outcome.record.label is not None else 0)
            elif position >= warmup+2+lag:
                score = outcome.record.label if outcome.record.label is not None else 0
            for (target, origin, target2, origin2), value in outcome.density.items():
                if target != target2 or origin != origin2:
                    raise ArithmeticError('native diagonal marker changed its reference code')
                updated[origin] += value*score*effect[target]
        effect = updated
    expected = evolve(family, score_effect(family, first, second, lag), warmup)
    configurations = tuple(tuple(i for i, bit in enumerate(state) if bit) for state in states)
    if effect != [evaluate(expected, x, 12) for x in configurations]:
        raise ArithmeticError('complete native H66 score differs from its degree-four pullback')
    return {'law': law, 'warmup_attempts': warmup, 'lag_attempts': lag,
            'source_columns': 66, 'checked_raw_slots': branch_count,
            'zero_operator_slots': zero_count, 'program_digest': interpreter.digest(program),
            'stationary_score': str(sum(effect, F())/66),
            'maximum_initial_basis_score': str(max(effect)),
            'minimum_initial_basis_score': str(min(effect)),
            'effect_coefficients': len(expected), 'degree_upper': max(map(len, expected), default=0)}


def prepared_h66_control(precision_bits=10, lag=2):
    """All-input exact adjoint propagation for the admitted fair H66 source."""
    from math import lcm
    s.integer(lag)
    if lag > 4:
        raise ValueError('at most four evaluated lag attempts required')
    budget = preparation_budget(12, 2, F(11, 240), precision_bits)
    if budget['warmup_attempts'] > 1000:
        raise ValueError('at most 1000 exact H66 marginal warmup attempts evaluated')
    interpreter = s.source.interpreter
    family = s.source.local_family(s.source.local.Parameters(port_weights=(1,)*12))
    first, second = next((a, b) for a in family.seams for b in family.seams if not set(a) & set(b))
    states = s.source.local.basis()
    configurations = tuple(tuple(i for i, bit in enumerate(state) if bit) for state in states)
    effect = score_effect(family, first, second, lag)
    vector = tuple(evaluate(effect, x, 12) for x in configurations)
    denominator = lcm(*(x.denominator for x in vector))
    values = tuple(int(x*denominator) for x in vector)
    center = sum(vector, F())/66
    command = {'op': 'repair', 'selector': {'kind': 'uniform'}, 'law': 'mixture', 'repeat': 1}
    program = interpreter.program(biased=False, commands=[command])
    marker = {(i, i, i, i): F(1) for i in range(66)}
    counts = [defaultdict(int) for _ in range(66)]
    for outcome in interpreter._outcomes(program, command, marker):
        for (target, origin, target2, origin2), weight in outcome.density.items():
            if target != target2 or origin != origin2:
                raise ArithmeticError('source adjoint count lost its diagonal marker')
            counts[origin][target] += 60*weight
    if any(sum(row.values()) != 60 or any(v.denominator != 1 for v in row.values()) for row in counts):
        raise ArithmeticError('actual source does not supply sixty integer population counts')
    counts = tuple(tuple((i, int(v)) for i, v in sorted(row.items())) for row in counts)
    for _ in range(budget['warmup_attempts']):
        values = tuple(sum(weight*values[j] for j, weight in row) for row in counts)
        denominator *= 60
    error = max(abs(F(value, denominator)-center) for value in values)/variance(12, 2)
    grid = 1 << 80
    directed_upper = F(-((-error.numerator*grid)//error.denominator), grid)
    if directed_upper > F(budget['normalized_score_error_upper']):
        raise ArithmeticError('all-input native propagation violates its analytic preparation budget')
    pair = [{'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(e)}, 'law': 'E', 'repeat': 1}
            for e in (first, second)]
    commands = [{**command, 'repeat': budget['warmup_attempts']}]+pair
    if lag:
        commands.append({**command, 'repeat': lag})
    commands += pair
    full_program = interpreter.program(biased=False, commands=commands, overflow='export_oldest')
    return {'budget': budget, 'lag_attempts': lag, 'source_columns': 66,
            'stationary_score': str(center), 'normalized_stationary_score': str(center/variance(12, 2)),
            'all_input_normalized_error_upper_80bit': str(directed_upper),
            'complete_program_digest': interpreter.digest(full_program),
            'complete_attempts': full_program['declared_resources']['attempts_upper'],
            'propagation': 'exact integer source adjoint marginal, with preparation records summed',
            'expanded_history_executed': False}
