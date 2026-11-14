"""Native feedback response from passive clock/field records, without a fit.

A finite amplitude, occupation-controlled repair has an exact locally balanced
likelihood identity. Native preparation establishes that balance; finite
warm-up bounds equilibrium forecasts. Complete records retain clock timeouts.
Parameters count native events.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as sp

from ..laboratory_clock import event_geometry as inherited
from ..observer import native
from ..record_laws import burn_in
from ..observer.scalar_response.experiment import compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*inherited.SOURCE_PATHS, *burn_in.SOURCE_PATHS,
    'model/predictions/repair_response.py', 'tools/strict_json.py'))))
PULSE = (0, 1)
# The impossible negative-current and coin-one positive-current labels are
# present in the readout alphabet and receive exact zero probability.
MARKS = ((0, -1, 0), *((1, coin, current) for coin in (0, 1) for current in (-1, 0, 1)))
READ_PORTS = (0, 1, 2)


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original integer inside the declared range is required')
    return value


@lru_cache(maxsize=1)
def source():
    kernel = inherited.transfer.native_counts()
    if PULSE not in kernel.seams:
        raise ArithmeticError('the pulse seam is absent from the native carrier')
    return kernel


def pulse_weights(state):
    """Postpulse density under endpoint balance, for each complete pulse mark.

    These formulas use only the disclosed endpoint bits. No counterfactual
    future, hidden state search or fitted response parameter is consulted.
    """
    if (type(state) is not tuple or len(state) != 12 or
            any(type(x) is not int or x not in (0, 1) for x in state)):
        raise ValueError('twelve original binary occupations are required')
    a, b = (state[j] for j in PULSE)
    result = {mark: F(0) for mark in MARKS}
    result[0, -1, 0] = F(1-a)
    result[1, 0, 1] = F((1-a)*b, 2)
    result[1, 1, 0] = F(a, 2)
    result[1, 0, 0] = F(a*b, 2)
    if sum(result.values()) != 1+F(b-a, 2):
        raise ArithmeticError('the complete pulse density lost its exact identity')
    return result


def _prior(prior):
    states = source().states
    if prior is None:
        prior = (F(1, len(states)),)*len(states)
    if (type(prior) is not tuple or len(prior) != len(states)
            or any(type(p) is not F or p < 0 for p in prior) or sum(prior) != 1):
        raise ValueError('a complete normalized exact native prior is required')
    return tuple(inherited.exact_core._fraction(p) for p in prior)


def pulse_distribution(*, reverse_forecast=False, prior=None):
    """Actual probability law, or an unnormalized passive weighted measure."""
    states, prior = source().states, _prior(prior)
    if type(reverse_forecast) is not bool:
        raise ValueError('a Boolean calculation mode is required')
    result = defaultdict(F)
    for state, weight in zip(states, prior, strict=True):
        if reverse_forecast:
            for mark, density in pulse_weights(state).items():
                result[mark, state] += weight*density
        elif state[PULSE[0]] == 0:
            result[(0, -1, 0), state] += weight
        else:
            for coin in (0, 1):
                after, (_, current) = native.native(state, PULSE, coin)
                mark = (1, coin, current)
                before_bits = tuple(state[j] for j in READ_PORTS)
                if pulse_readback(before_bits, mark) != tuple(after[j] for j in READ_PORTS):
                    raise ArithmeticError('the pulse record lost the disclosed initial fields')
                result[mark, after] += weight/2
    return {key: value for key, value in result.items() if value}


def pulse_readback(before_bits, mark):
    """Postpulse bits from the pre-pulse disclosure and retained pulse current."""
    if (type(before_bits) is not tuple or len(before_bits) != 3
            or any(type(x) is not int or x not in (0, 1) for x in before_bits)
            or type(mark) is not tuple or len(mark) != 3
            or any(type(x) is not int for x in mark) or mark not in MARKS):
        raise ValueError('three disclosed bits and a complete pulse mark are required')
    trigger, coin, current = mark
    if trigger != before_bits[0]:
        raise ValueError('the pulse trigger must read the disclosed port-zero bit')
    if not trigger:
        return before_bits
    # Checking the deterministic native endpoint prevents a fabricated mark
    # from acquiring observer-access status by a formal current subtraction.
    pair = before_bits[:2]+(0,)*10
    after, expected = native.native(pair, PULSE, coin)
    if expected[1] != current:
        raise ValueError('the retained pulse mark is incompatible with the disclosed bits')
    return after[0], after[1], before_bits[2]


def _reading(mark, initial, tape, terminal, status, clock_b):
    """Every field is reconstructed from disclosed bits and retained currents."""
    fields = list(initial)
    for edge, current in tape:
        for i, port in enumerate(READ_PORTS):
            fields[i] += current*(int(port == edge[1])-int(port == edge[0]))
    if tuple(fields) != tuple(terminal[j] for j in READ_PORTS):
        raise ArithmeticError('the visible current tape did not reconstruct the field')
    a_ticks = mark[0]+int(status == 'tick')
    b_ticks = mark[0]*int(clock_b in PULSE)+sum(clock_b in edge for edge, _ in tape)
    return mark, initial, tape, status, (a_ticks, b_ticks), tuple(fields)


def stopped_law(cap=2, *, clock_b=1, reverse_forecast=False, intervention=True, prior=None):
    """Complete union tape until the first postpulse port-zero tick or timeout.

    The pulse counts as a tick iff the trigger is one. The independent
    attempt watchdog limits execution but its counter is not a clock readout.
    """
    cap = _integer(cap, 1, 2)
    clock_b = _integer(clock_b, 1, 2)
    if type(intervention) is not bool or type(reverse_forecast) is not bool:
        raise ValueError('Boolean intervention and forecast choices are required')
    if reverse_forecast and not intervention:
        raise ValueError('the forecast concerns the declared intervention')
    prior = _prior(prior)
    initial_law = (pulse_distribution(reverse_forecast=reverse_forecast, prior=prior) if intervention else
                   {((0, -1, 0), state): mass for state, mass in zip(source().states, prior) if mass})
    initial_mass = sum(initial_law.values())
    populations = {}
    for (mark, state), mass in initial_law.items():
        endpoint = tuple(state[j] for j in READ_PORTS)
        before = (endpoint[0]+mark[2], endpoint[1]-mark[2], endpoint[2])
        # In the intervention arm before is the originally disclosed field
        # record, recoverable from endpoint+mark. No extra field read occurs.
        initial = pulse_readback(before, mark) if intervention else endpoint
        key = mark, initial, (), state
        populations[key] = mass
    result = defaultdict(F)
    for step in range(cap):
        future = defaultdict(F)
        for (mark, initial, tape, state), mass in populations.items():
            for edge in source().seams:
                for coin in (0, 1):
                    terminal, (_, current) = native.native(state, edge, coin)
                    read = tape+((edge, current),) if any(j in edge for j in READ_PORTS) else tape
                    if 0 in edge:
                        result[_reading(mark, initial, read, terminal, 'tick', clock_b)] += mass/60
                    elif step+1 == cap:
                        result[_reading(mark, initial, read, terminal, 'timeout', clock_b)] += mass/60
                    else:
                        future[mark, initial, read, terminal] += mass/60
        populations = future
    if sum(result.values()) != initial_mass:
        raise ArithmeticError('the complete stopped law lost initial measure mass')
    return dict(result)


def preparation_defect(prior):
    """Exact marked-measure L1 defect; no normalization or half-L1 convention.

    Only the moving mark differs: A_mu(m,s)-W_mu(m,s) =
    [mu(swap s)-mu(s)]/2 when (s_0,s_1)=(0,1).
    Consequently L1(A_mu-W_mu)=TV(mu,swap mu)/2 <= TV(mu,uniform).
    A common complete continuation contracts this L1 norm at any finite cap.
    """
    prior = _prior(prior)
    states = source().states
    mass = dict(zip(states, prior))
    swapped = tuple(mass[(s[1], s[0])+s[2:]] for s in states)
    tv = sum((abs(p-F(1, 66)) for p in prior), F())/2
    swap_tv = sum((abs(p-q) for p, q in zip(prior, swapped)), F())/2
    actual = pulse_distribution(prior=prior)
    forecast = pulse_distribution(prior=prior, reverse_forecast=True)
    l1 = sum((abs(actual.get(k, F())-forecast.get(k, F()))
              for k in actual.keys() | forecast.keys()), F())
    for mark, state in actual.keys() | forecast.keys():
        expected = ((mass[(state[1], state[0])+state[2:]]-mass[state])/2
                    if mark == (1, 0, 1) and state[:2] == (0, 1) else F())
        if actual.get((mark, state), F())-forecast.get((mark, state), F()) != expected:
            raise ArithmeticError('the exact preparation defect lost its moving-branch identity')
    if l1 != swap_tv/2 or l1 > tv:
        raise ArithmeticError('the sharp preparation bound failed')
    return {'source_probability_TV_to_uniform': tv, 'source_probability_TV_to_swap': swap_tv,
            'marked_forecast_L1_error': l1, 'forecast_total_mass': sum(forecast.values()),
            'forecast_is_probability': sum(forecast.values()) == 1,
            'exact_marked_identity': swap_tv == 0,
            'any_complete_continuation_L1_error_upper': l1,
            'any_readout_event_absolute_error_upper': l1}


def warmup_prior(initial, attempts):
    """Exact diagonal warm-up, with its entire preparation tape discarded."""
    initial = _source_state(initial)
    attempts = _integer(attempts, 0, 1000)
    kernel = source()
    row = [int(state == initial) for state in kernel.states]
    transitions = [tuple((j, n) for j, n in enumerate(values) if n)
                   for values in kernel.counts]
    for _ in range(attempts):
        after = [0]*66
        for value, outgoing in zip(row, transitions):
            for j, multiplicity in outgoing:
                after[j] += value*multiplicity
        row = after
    return tuple(F(n, 60**attempts) for n in row)


def local_balancing_preparation(prior):
    """One actual fair repair gives endpoint-swap balance on the binary sector.

    Preparation coin/current are not made available to the subsequent menu.
    Its deterministic incident tick is included as a known clock offset.
    """
    prior = _prior(prior)
    output = defaultdict(F)
    for state, mass in zip(source().states, prior):
        for coin in (0, 1):
            after, _ = native.native(state, PULSE, coin)
            output[after] += mass/2
    balanced = tuple(output[s] for s in source().states)
    if not preparation_defect(balanced)['exact_marked_identity']:
        raise ArithmeticError('one native local preparation failed endpoint balance')
    return {'prior': balanced, 'native_preparation_attempts': 1,
            'incident_tick_offsets_on_read_ports': tuple(int(p in PULSE) for p in READ_PORTS),
            'preparation_coin_current_unavailable_to_future_protocol': True}


def preparation_budget(attempts, coverage=F(1)):
    """Inherited binary-sector warm-up bound in repair attempts, not seconds.

    The exact pulse bound improves the triangle bound between two normalized
    arms: weighting passive records from the same warmed preparation costs
    at most delta in L1, even though that weighted measure may have mass !=1.
    """
    attempts = _integer(attempts, 0, 1000)
    coverage = inherited.exact_core._fraction(coverage)
    delta = burn_in.binary_bound(attempts, coverage)
    return {'warmup_attempts': attempts, 'exogenous_seam_coverage': coverage,
            'preparation_half_trace_distance_upper': delta,
            'complete_actual_readout_probability_TV_to_equilibrium_upper': delta,
            'same_preparation_weighted_forecast_L1_error_upper': delta,
            'same_preparation_any_event_forecast_error_upper': delta,
            'warmup_tape_and_controller_discarded': True,
            'binary_sector_support_required': True,
            'forecast_renormalized': False,
            'stationary_initial_preparation_required': False}


def density_response(attempts):
    """All twelve mean-field responses, from endpoint sums and passive covariance."""
    attempts = _integer(attempts, 0, 8)
    kernel = source()
    distribution = defaultdict(F)
    for (_, state), mass in pulse_distribution().items():
        distribution[state] += mass
    row = tuple(distribution[state] for state in kernel.states)
    for _ in range(attempts):
        row = inherited.transfer._row_times(row, kernel.counts)
    actual = tuple(sum((mass*(state[c]-F(1, 6)) for state, mass in zip(kernel.states, row)), F())
                   for c in range(12))
    laplacian = sp.zeros(12)
    for a, b in kernel.seams:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    transfer = (sp.eye(12)-laplacian/60)**attempts
    predicted = tuple(F(5, 66)*F(transfer[1, c]-transfer[0, c]) for c in range(12))
    if actual != predicted:
        raise ArithmeticError('the finite-amplitude response differs from passive covariance')
    return actual


def discounted_response(z=1):
    """Geometrically discounted response is the same native spatial Green dipole."""
    if type(z) not in (int, F) or z <= 0:
        raise ValueError('a positive exact discount parameter is required')
    z = inherited.exact_core._fraction(z)
    kernel = source()
    c = sp.Matrix(kernel.counts)
    mass = defaultdict(F)
    for (_, state), value in pulse_distribution().items():
        mass[state] += value
    delta = sp.Matrix(1, 66, [sp.Rational(mass[s]-F(1, 66)) for s in kernel.states])
    field = sp.Matrix(66, 12, lambda i, j: sp.Rational(kernel.states[i][j])-sp.Rational(1, 6))
    oracle = delta*((60+sp.Rational(z))*sp.eye(66)-c).inv(method='DM')*field
    laplacian = sp.zeros(12)
    for a, b in kernel.seams:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    dipole = sp.zeros(1, 12); dipole[0, 1], dipole[0, 0] = 1, -1
    predicted = sp.Rational(5, 66)*dipole*(laplacian+sp.Rational(z)*sp.eye(12)).inv()
    if oracle != predicted:
        raise ArithmeticError('the full configuration resolvent lost the spatial response')
    return tuple(F(value) for value in predicted)


def _profile(values):
    if type(values) is not tuple or len(values) != 12:
        raise ValueError('a fixed twelve-port source profile is required')
    return tuple(inherited.exact_core._fraction(value) for value in values)


def _source_state(state):
    if (type(state) is not tuple or len(state) != 12
            or any(type(x) is not int or x not in (0, 1) for x in state)
            or sum(state) != 2):
        raise ValueError('an original twelve-port binary source state with total two is required')
    return state


def activation(state, edge, profile, epsilon):
    """Local controller probability; the subsequent native endpoint is unchanged."""
    profile = _profile(profile)
    epsilon = inherited.exact_core._fraction(epsilon)
    state = _source_state(state)
    if (type(edge) is not tuple or len(edge) != 2
            or any(type(x) is not int for x in edge) or edge not in source().seams):
        raise ValueError('an original native oriented seam is required')
    a, b = edge
    bias = epsilon*(profile[b]-profile[a])
    if abs(bias) > F(1, 2):
        raise ValueError('the certified activation probability must lie in [1/4,3/4]')
    return (1+bias*(state[a]-state[b]))/2


def controlled_matrices(profile, epsilon=F(0)):
    """Full native slot kernel and actual executed-repair stopping kernel.

    Idle/control slots do not tick this particular aggregate dynamical clock.
    Successful native attempts tick even at zero current. Controller service
    coin writes, if implemented natively, have their own retained cost count.
    """
    profile, epsilon = _profile(profile), inherited.exact_core._fraction(epsilon)
    kernel = source(); index = {state: i for i, state in enumerate(kernel.states)}
    success, success_derivative = sp.zeros(66), sp.zeros(66)
    activity_derivative = []
    for i, state in enumerate(kernel.states):
        derivative = F(0)
        for a, b in kernel.seams:
            probability = activation(state, (a, b), profile, epsilon)
            slope = (profile[b]-profile[a])*(state[a]-state[b])/2
            derivative += slope/30
            for coin in (0, 1):
                after, _ = native.native(state, (a, b), coin)
                j = index[after]
                success[i, j] += sp.Rational(probability)/60
                success_derivative[i, j] += sp.Rational(slope)/60
        activity_derivative.append(sp.Rational(derivative))
    activity = success*sp.ones(66, 1)
    diag = sp.diag(*activity)
    first = sp.diag(*activity_derivative)
    k = sp.Matrix(kernel.counts)/60
    embedded = sp.diag(*(1/a for a in activity))*success
    slot = sp.eye(66)-diag+success
    derivative_slot = success_derivative-first
    derivative_embedded = 2*success_derivative-2*first*k
    uniform = sp.ones(1, 66)/66
    if (uniform*success_derivative != sp.zeros(1, 66)
            or slot*sp.ones(66, 1) != sp.ones(66, 1)
            or embedded*sp.ones(66, 1) != sp.ones(66, 1)):
        raise ArithmeticError('the native marked control or Palm normalization changed')
    return {'slot': sp.ImmutableMatrix(slot), 'embedded': sp.ImmutableMatrix(embedded),
            'success': sp.ImmutableMatrix(success), 'activity': tuple(F(a) for a in activity),
            'slot_derivative_at_zero': sp.ImmutableMatrix(derivative_slot),
            'embedded_derivative_at_zero': sp.ImmutableMatrix(derivative_embedded),
            'activity_derivative': tuple(F(a) for a in activity_derivative)}


def controlled_response(profile, count, *, clock='executed'):
    """Exact derivative of every terminal configuration weight, at epsilon=0."""
    profile = _profile(profile)
    count = _integer(count, 0, 8)
    if clock not in ('executed', 'slot'):
        raise ValueError('the clock is either executed dynamical repairs or control slots')
    matrices = controlled_matrices(profile)
    key = 'embedded' if clock == 'executed' else 'slot'
    p, dp = matrices[key], matrices[key+'_derivative_at_zero']
    uniform, derivative = sp.ones(1, 66)/66, sp.zeros(1, 66)
    for _ in range(count):
        derivative = derivative*p+uniform*dp
    laplacian = sp.zeros(12)
    for a, b in source().seams:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    k = sp.eye(12)-laplacian/60
    h = sp.Matrix(profile)
    vector = (2*k*(sp.eye(12)-k**count)*h if clock == 'executed'
              else 2*(sp.eye(12)-((sp.eye(12)+k)/2)**count)*h)
    field = sp.Matrix(66, 12, lambda i, j: source().states[i][j]-sp.Rational(1, 6))
    predicted = (field*vector).T/66
    if derivative != predicted or sum(derivative) != 0:
        raise ArithmeticError('the full native response lost its stopping-clock correction')
    means = derivative*field
    if means != sp.Rational(5, 33)*vector.T:
        raise ArithmeticError('the source susceptibility normalization changed')
    return {'terminal_probability_derivative': tuple(F(v) for v in derivative),
            'occupation_mean_derivative': tuple(F(v) for v in means)}


def controlled_history(profile, epsilon, initial, *, ticks=1, slot_cap=2):
    """Complete marked finite execution, including idle slots and timeout fields.

    Initial occupations are disclosed. Every selected seam, execution flag,
    coin and current is retained. The terminal field is reconstructed by those
    currents, including when the cap fires; no timeout is discarded or zeroed.
    """
    profile, epsilon = _profile(profile), inherited.exact_core._fraction(epsilon)
    ticks, slot_cap = _integer(ticks, 1, 2), _integer(slot_cap, 1, 2)
    initial = _source_state(initial)
    populations = {((), initial, 0): F(1)}; result = defaultdict(F)
    for step in range(slot_cap):
        future = defaultdict(F)
        for (tape, state, completed), mass in populations.items():
            for edge in source().seams:
                probability = activation(state, edge, profile, epsilon)
                branches = [(0, -1, 0, state, 1-probability)]
                for coin in (0, 1):
                    after, (_, current) = native.native(state, edge, coin)
                    branches.append((1, coin, current, after, probability/2))
                for active, coin, current, after, weight in branches:
                    record = tape+((edge, active, coin, current),)
                    count = completed+active
                    if count == ticks or step+1 == slot_cap:
                        fields = list(initial)
                        for (a, b), _, _, j in record:
                            fields[a] -= j; fields[b] += j
                        if tuple(fields) != after:
                            raise ArithmeticError('the full retained tape lost field readback')
                        key = record, 'complete' if count == ticks else 'timeout', count, tuple(fields)
                        result[key] += mass*weight/30
                    else:
                        future[record, after, count] += mass*weight/30
        populations = future
    if sum(result.values()) != 1:
        raise ArithmeticError('the complete controlled execution lost probability')
    return dict(result)


def native_service_coin_word(state, coins):
    """Fair-bit source using a recorded equal-occupation seam, with charged cost.

    The field/current interface discloses the equal pair. Each service write
    is an actual native repair, preserves the load and is separately counted.
    The word parameter specifies a branch, not chosen random coin outcomes.
    """
    state = _source_state(state)
    if type(coins) is not tuple or any(
            type(coin) is not int or coin not in (0, 1) for coin in coins):
        raise ValueError('a native binary state and an original binary coin word are required')
    edge = next((edge for edge in source().seams if state[edge[0]] == state[edge[1]]), None)
    if edge is None:
        raise ArithmeticError('the source triangle supplied no equal binary pair')
    records = []
    for coin in coins:
        after, record = native.native(state, edge, coin)
        if after != state or record != (coin, 0):
            raise ArithmeticError('a supposedly neutral coin service changed the carrier')
        records.append((edge, *record))
    return {'records': tuple(records), 'branch_probability': F(1, 2**len(coins)),
            'native_service_events': len(coins), 'dynamical_clock_ticks': 0,
            'terminal_load': state}


def dyadic_activation(probability, bits):
    """An executable finite fair-bit controller, with explicit approximation."""
    probability = inherited.exact_core._fraction(probability)
    bits = _integer(bits, 1, 10000)
    if not 0 <= probability <= 1:
        raise ValueError('a probability is required')
    denominator = 2**bits
    numerator = (probability*denominator).numerator//(probability*denominator).denominator
    return {'probability': F(numerator, denominator), 'absolute_error_upper': F(1, denominator),
            'fair_bits_per_slot': bits, 'native_service_events_if_equal_seam_realized': bits,
            'uniform_seam_selection_service_is_additional': True}


def likelihood_control(profile, epsilon, *, ticks=1, slot_cap=4):
    """Exact stopped likelihood/score moments from complete native branches.

    Six propagated moments retain correlations with the configuration and
    completed clock count. Timeout uses its actual terminal field. The coarse
    marks specify seam, activation, native coin/current; additional controller
    fair-bit words are separately recorded service resources.
    """
    profile, epsilon = _profile(profile), inherited.exact_core._fraction(epsilon)
    ticks, slot_cap = _integer(ticks, 1, 2), _integer(slot_cap, 1, 8)
    states, seams = source().states, source().seams
    # (mass, score, score^2, likelihood, likelihood^2, likelihood*score).
    moments = {(state, 0): (F(1, 66), F(0), F(0), F(1, 66), F(1, 66), F(0))
               for state in states}
    energy = sum(((profile[b]-profile[a])**2 for a, b in seams), F())
    variance = energy/30
    x = epsilon**2*slot_cap*variance
    if x >= 1:
        raise ValueError('the rational uniform likelihood allowance needs epsilon^2*q*v<1')
    for _ in range(slot_cap):
        future = defaultdict(lambda: [F(0)]*6)
        for (state, completed), values in moments.items():
            if completed == ticks:
                future[state, completed] = [a+b for a, b in zip(future[state, completed], values)]
                continue
            p, ps, pss, pl, pll, pls = values
            for edge in seams:
                a, b = edge
                probability = activation(state, edge, profile, epsilon)
                score = (profile[b]-profile[a])*(state[a]-state[b])
                branches = [(state, completed, -score, F(1, 60))]
                for coin in (0, 1):
                    after, _ = native.native(state, edge, coin)
                    branches.append((after, completed+1, score, F(1, 120)))
                if probability != (1+epsilon*score)/2:
                    raise ArithmeticError('the retained activation score changed')
                for after, count, step_score, weight in branches:
                    factor = 1+epsilon*step_score
                    row = (p, ps+step_score*p, pss+2*step_score*ps+step_score**2*p,
                           factor*pl, factor**2*pll, factor*(pls+step_score*pl))
                    target = future[after, count]
                    for j, value in enumerate(row):
                        target[j] += weight*value
        moments = dict(future)
    total = tuple(sum((row[j] for row in moments.values()), F()) for j in range(6))
    p, ps, pss, pl, pll, pls = total
    remainder_squared = pll+1+epsilon**2*pss-2*pl-2*epsilon*pls+2*epsilon*ps
    polynomial_bound = (1+epsilon**2*variance)**slot_cap-1-slot_cap*epsilon**2*variance
    rational_bound = x**2/(2*(1-x))
    timeout = sum((row[0] for (_, count), row in moments.items() if count < ticks), F())
    biased_timeout = sum((row[3] for (_, count), row in moments.items() if count < ticks), F())
    field = {state: sum(((F(n)-F(1, 6))*h for n, h in zip(state, profile)), F())
             for state in states}
    baseline_variance = sum((row[0]*field[state]**2 for (state, _), row in moments.items()), F())
    derivative = sum((row[1]*field[state] for (state, _), row in moments.items()), F())
    mean = sum((row[3]*field[state] for (state, _), row in moments.items()), F())
    source_variance = sum((field[state]**2/F(66) for state in states), F())
    centered = tuple(h-sum(profile)/12 for h in profile)
    if (p != 1 or pl != 1 or ps != 0 or baseline_variance != source_variance
            or source_variance != F(5, 33)*sum(h*h for h in centered)
            or not 0 <= remainder_squared <= polynomial_bound <= rational_bound
            or (mean-epsilon*derivative)**2 > baseline_variance*remainder_squared
            or biased_timeout**2 > pll*timeout):
        raise ArithmeticError('the exact stopped likelihood or correlated timeout bound failed')
    unlimited = controlled_response(profile, ticks)['occupation_mean_derivative']
    unlimited_derivative = sum((h*d for h, d in zip(profile, unlimited)), F())
    maximum_field = max(abs(v) for v in field.values())
    derivative_tail_squared = 8*maximum_field**2*ticks*variance*timeout
    if (derivative-unlimited_derivative)**2 > derivative_tail_squared:
        raise ArithmeticError('the stopped derivative comparison lost its clock correlation')
    return {'ticks': ticks, 'slot_cap': slot_cap, 'epsilon': epsilon,
            'baseline_score_mean': ps, 'baseline_score_second_moment': pss,
            'likelihood_mean': pl, 'likelihood_second_moment': pll,
            'linear_likelihood_remainder_second_moment': remainder_squared,
            'linear_likelihood_remainder_polynomial_upper': polynomial_bound,
            'linear_likelihood_remainder_rational_upper': rational_bound,
            'conditional_score_variance_upper': variance,
            'baseline_terminal_field_variance': baseline_variance,
            'capped_field_derivative': derivative, 'unbounded_clock_field_derivative': unlimited_derivative,
            'capped_derivative_error_squared_upper': derivative_tail_squared,
            'biased_terminal_field_mean': mean, 'field_nonlinear_error': mean-epsilon*derivative,
            'baseline_timeout': timeout, 'biased_timeout': biased_timeout,
            'timeout_fields_are_terminal_fields': True}


def continuum_response_budget(level, epsilon=F(1, 10**6)):
    """All-level bounds for a fixed bounded flat coordinate at density one-third.

    These are analytic resource and error bounds, not allocated large screens.
    The positive limiting response is a controlled mean, not a physical force
    identification or a proof of a unique nonlinear finite-epsilon law.
    """
    level = _integer(level, 0, 100)
    epsilon = inherited.exact_core._fraction(epsilon)
    if not 0 < epsilon <= F(1, 100):
        raise ValueError('the retained finite-signal regime is 0<epsilon<=1/100')
    scale, cells = 4**level, 20*4**level
    ports, seams = 12*cells, 630*scale
    number = ports//3
    chi = F(number*(ports-number), ports*(ports-1))
    ticks, slots = 2*seams*scale, 8*seams*scale
    variance = F(slots, chi*seams*scale)  # normalized scaled full-lift energy <=1
    if not variance <= 36 or epsilon**2*variance > F(1, 2):
        raise ArithmeticError('the uniform nonlinear response allowance failed')
    # A>=3/4096, ||coordinate||^2>=3/16, 1-exp(-x)>=x/(1+x).
    linear_lower = F(9, 32792)
    nonlinear_allowance = epsilon**2*36
    limiting_lower = epsilon*linear_lower-nonlinear_allowance
    required_denominator = 16*ports*slots**2
    bits = (required_denominator-1).bit_length()
    return {'level': level, 'ports': ports, 'seams': seams, 'particles': number,
            'susceptibility': chi, 'normalized_profile_maximum_upper': F(1),
            'uniform_positivity_radius': F(1, 4),
            'executed_dynamical_repairs_at_tau_one': ticks,
            'selection_slot_cap': slots, 'baseline_expected_slots_without_cap': 2*ticks,
            'biased_expected_slots_without_cap_upper': 4*ticks,
            'score_variance_upper': variance,
            'baseline_timeout_dyadic_exponent': ticks//2,
            'baseline_timeout_upper': '2^(-baseline_timeout_dyadic_exponent)',
            'biased_timeout_probability_squared_upper_factor': F(2),
            'capped_linear_derivative_error_squared_upper_factor': F(8*ports*ticks, chi**2*seams*scale),
            'capped_linear_derivative_error_squared_upper': 'factor * 2^(-baseline_timeout_dyadic_exponent)',
            'epsilon': epsilon, 'limiting_linear_response_lower': linear_lower,
            'limiting_nonlinear_mean_error_upper': nonlinear_allowance,
            'liminf_mean_response_lower': limiting_lower,
            'liminf_lower_is_positive': limiting_lower > 0,
            'dyadic_activation_bits_per_slot': bits,
            'maximum_native_coin_service_events': bits*slots,
            'additional_uniform_seam_selection_calls': slots,
            'dyadic_mean_error_squared_upper': F(1, 64*slots**2*ports*chi),
            'all_bounds_include_terminal_fields_on_timeout': True,
            'complete_nonlinear_continuum_law_identified': False,
            'clock_is_a_primitive_port_incident_clock': False}


def controlled_response_run():
    profile = (F(1), F(-1))+(F(0),)*10
    rows = []
    for count in (1, 2, 4):
        tick = controlled_response(profile, count)
        slot = controlled_response(profile, 2*count, clock='slot')
        defect = tick['occupation_mean_derivative'][0]-slot['occupation_mean_derivative'][0]
        if defect == 0:
            raise ArithmeticError('naive slot-to-clock substitution lost its finite control')
        rows.append({'executed_repairs': count, 'correct_response': tick,
                     'double_slot_substitution_response': slot,
                     'naive_clock_substitution_port_zero_defect': defect})
    initial = next(state for state in source().states if state[0] and not state[1])
    law = controlled_history(profile, F(1, 8), initial, ticks=1, slot_cap=2)
    zero = sum((p for key, p in law.items() if any(r[1] and r[3] == 0 for r in key[0])), F())
    idle = sum((p for key, p in law.items() if any(not r[1] for r in key[0])), F())
    timeout = sum((p for key, p in law.items() if key[1] == 'timeout'), F())
    if not 0 < zero < 1 or not 0 < idle < 1 or not 0 < timeout < 1:
        raise ArithmeticError('the controlled full-history counter witnesses vanished')
    return {'profile': profile, 'full_configuration_derivatives': rows,
            'complete_controlled_history': {'epsilon': F(1, 8), 'initial_load': initial,
                'target_executed_repairs': 1, 'slot_cap': 2,
                'nonzero_record_labels': len(law), 'law_sha256': law_digest(law),
                'some_zero_current_execution_probability': zero,
                'some_idle_slot_probability': idle, 'timeout_probability': timeout,
                'timeout_fields_are_terminal_fields': True},
            'exact_stopped_likelihood_control': likelihood_control(profile, F(1, 8)),
            'continuum_normalization_and_resources': [continuum_response_budget(level) for level in (0, 1, 2)],
            'continuum_mean_derivative': '2*(I-exp(-tau*A))*f',
            'finite_clock_mean_derivative': '2*chi*K*(I-K^r)*h',
            'finite_slot_mean_derivative': '2*chi*(I-((I+K)/2)^q)*h',
            'weak_force_if_f_in_form_domain': '2*A*f, in the dual of the closed form domain',
            'natural_physical_force_or_atomic_clock_identified': False}


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def law_digest(law):
    rows = [[_json(key), str(value)] for key, value in sorted(law.items())]
    return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()


def preparation_run():
    states = source().states
    left = next(s for s in states if s[0] and not s[1])
    right = (left[1], left[0])+left[2:]
    # Far from equilibrium, but locally swap-balanced and with nonzero pulse response.
    balanced = tuple(F(int(s in (left, right)), 2) for s in states)
    invariant = preparation_defect(balanced)
    if not invariant['exact_marked_identity'] or invariant['source_probability_TV_to_uniform'] == 0:
        raise ArithmeticError('local balance did not remove global stationarity')
    local = local_balancing_preparation(tuple(F(int(s == left)) for s in states))
    if local['prior'] != balanced:
        raise ArithmeticError('the native one-repair preparation differs from endpoint averaging')
    epsilon = F(1, 132)
    sharp = tuple(F(1, 66)+epsilon*(int(s == left)-int(s == right)) for s in states)
    tight = preparation_defect(sharp)
    actual = stopped_law(2, prior=sharp)
    forecast = stopped_law(2, prior=sharp, reverse_forecast=True)
    l1 = sum((abs(actual.get(k, F())-forecast.get(k, F()))
              for k in actual.keys() | forecast.keys()), F())
    timeout_defect = sum((actual.get(k, F())-forecast.get(k, F())
                         for k in actual.keys() | forecast.keys() if k[3] == 'timeout'), F())
    if (l1 != epsilon or tight['source_probability_TV_to_uniform'] != epsilon
            or sum(forecast.values()) != 1-epsilon or timeout_defect != epsilon*F(5, 6)**2):
        raise ArithmeticError('the tight complete-law preparation control failed')
    def upper(value):
        scaled = value*10**12
        return F(-(-scaled.numerator//scaled.denominator), 10**12)
    rows = []
    for n in (0, 8, 64, 221):
        prior = warmup_prior(left, n)
        defect, budget = preparation_defect(prior), preparation_budget(n)
        if defect['source_probability_TV_to_uniform'] > budget['preparation_half_trace_distance_upper']:
            raise ArithmeticError('native warm-up exceeds its inherited error budget')
        rows.append({'warmup_attempts': n,
            'measured_source_TV_upper': upper(defect['source_probability_TV_to_uniform']),
            'measured_weighted_forecast_L1_upper': upper(defect['marked_forecast_L1_error']),
            'uniform_source_and_forecast_error_upper': upper(budget['preparation_half_trace_distance_upper']),
            'prior_sha256': hashlib.sha256(str(prior).encode()).hexdigest()})
    return {'identity_requires_only_pulse_endpoint_swap_invariance': True,
            'exact_local_preparation': {k: v for k, v in local.items() if k != 'prior'},
            'nonstationary_swap_invariant_control': invariant,
            'exact_defect': 'L1(A_mu-W_mu)=TV(mu,swap(mu))/2 <= TV(mu,uniform)',
            'sharp_complete_stopped_control': {**tight, 'attempt_cap': 2,
                'actual_law_sha256': law_digest(actual), 'weighted_forecast_sha256': law_digest(forecast),
                'complete_stopped_L1_defect': l1, 'timeout_event_defect': timeout_defect,
                'forecast_renormalized': False},
            'finite_warmup': rows,
            'contract': {'binary_source_density_required': True,
                'any_admissible_current_record_law_covered': True,
                'warmup_seam_schedule_exogenous_with_positive_coverage': True,
                'warmup_tape_and_controller_unavailable_to_readout': True,
                'all_pulse_marks_zero_currents_and_timeouts_retained': True,
                'later_classical_record_continuation_common_and_complete': True,
                'physical_duration_or_record_law_selection_claimed': False}}


def run():
    actual, predicted = pulse_distribution(), pulse_distribution(reverse_forecast=True)
    if actual != predicted:
        raise ArithmeticError('native pulse and passive-record forecast disagree')
    rows = []
    for cap in (1, 2):
        for clock_b in (1, 2):
            law = stopped_law(cap, clock_b=clock_b)
            if law != stopped_law(cap, clock_b=clock_b, reverse_forecast=True):
                raise ArithmeticError('the complete clock/field response forecast failed')
            mark_mass = tuple(sum((v for key, v in law.items() if key[0] == mark), F()) for mark in MARKS)
            timeout = sum((v for key, v in law.items() if key[3] == 'timeout'), F())
            if timeout != F(5, 6)**cap:
                raise ArithmeticError('the watchdog tail differs from native incidence')
            signal = sum((v for key, v in law.items() if key[3] == 'tick'
                          and key[-1][1] == 1 and key[-2][1] >= 1), F())
            baseline = stopped_law(cap, clock_b=clock_b, intervention=False)
            def future_signal(key):
                return (key[3] == 'tick' and key[-1][1] == 1
                        and key[-2][1]-key[0][0]*int(clock_b in PULSE) >= 1)
            response = sum((v for key, v in law.items() if future_signal(key)), F())
            passive = sum((v for key, v in baseline.items() if future_signal(key)), F())
            covariance = sum((v*F(key[1][1]-key[1][0], 2) for key, v in baseline.items()
                              if future_signal(key)), F())
            if response-passive != covariance or covariance < 0:
                raise ArithmeticError('the held-out joint event violates its passive response forecast')
            rows.append({'attempt_cap': cap, 'second_clock_port': clock_b,
                'nonzero_complete_labels': len(law), 'law_sha256': law_digest(law),
                'pulse_mark_probabilities': mark_mass, 'timeout_probability': timeout,
                'joint_signal': 'completed port-zero tick, final port-one field one, second clock at least one',
                'joint_signal_probability': signal,
                'matched_continuation_signal': 'completed tick, final port-one field one, at least one subsequent B tick',
                'passive_signal_probability': passive, 'feedback_signal_probability': response,
                'signal_increase_predicted_from_passive_covariance': covariance,
                'all_current_zero_events_retained': True,
                'pulse_is_charged_to_each_incident_clock': True})
    chosen = next(i for i, state in enumerate(source().states) if state[0] == 1 and state[1] == 0)
    prior = tuple(F(int(i == chosen)) for i in range(66))
    nonstationary = pulse_distribution(prior=prior)
    invalid_forecast = pulse_distribution(prior=prior, reverse_forecast=True)
    if nonstationary == invalid_forecast:
        raise ArithmeticError('the nonstationary premise control failed')
    return _json({'schema': 'oph.native_repair_response.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'pulse': {'seam': PULSE, 'policy': 'repair iff disclosed initial port-zero occupation equals one',
            'complete_marks_trigger_coin_current': MARKS,
            'no_pulse_coin_sentinel': -1,
            'sum_density': '1+(n_1-n_0)/2',
            'passive_future_event_response': 'Delta Pr(F)=E[(n_1-n_0) 1_F]/2'},
        'stopped_joint_forecasts': rows,
        'mean_field_responses': {str(n): density_response(n) for n in (0, 1, 2, 4)},
        'discounted_spatial_response': {str(z): discounted_response(z) for z in (1, 2)},
        'controlled_continuum_response': controlled_response_run(),
        'finite_preparation_response': preparation_run(),
        'controls': {'nonstationary_prior_rejects_stationary_weight_formula': True,
            'nonstationary_actual_total_mass': sum(nonstationary.values()),
            'nonstationary_invalid_forecast_total_mass': sum(invalid_forecast.values()),
            'nonstationary_invalid_forecast_is_not_renormalized': True,
            'future_law_includes_trigger_coin_current': True},
        'scope': {'stationary_uniform_fixed_number_reference_law': True,
            'exact_stationary_preparation_required_for_finite_error_forecast': False,
            'source_readouts_and_future_continuation_fixed_before_comparison': True,
            'response_parameter_fitted': False, 'pulse_has_measured_physical_duration': False,
            'unmarked_response_requires_continuation_independent_of_pulse_mark': True,
            'marked_response_may_use_a_separately_declared_continuation_for_each_mark': True,
            'arbitrary_continuation_theorem_requires_identical_continuation_in_both_arms': True,
            'matched_public_intervention_dataset_available': False,
            'new_discriminator_against_SM_GR_established': False,
            'Lorentzian_or_CFT_limit_established': False}, 'all_pass': True})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(result, strict_json.load(args.check))
        print('Native passive-record response and complete clock/field forecasts replayed.')
    else:
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
