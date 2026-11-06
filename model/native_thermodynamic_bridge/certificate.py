"""Exact native heat baths, common reference and entropy identities.

All transition probabilities come from source endpoint execution. Entropy
identities are checked as rational linear combinations of logarithms of
primes, without rounded logarithms or inferred positivity from floating data.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from tools import strict_json
from ..native_born_bridge import certificate as born
from ..native_equilibration import cut_completion as cut

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASIS = born.native.states(binary=True)
INDEX = {state: i for i, state in enumerate(BASIS)}
N = len(BASIS)
require = born.require


def law(value, *, positive=False, size=N):
    require(type(value) is tuple and len(value) == size, 'complete exact probability tuple required')
    result = tuple(born.rational(x) for x in value)
    require(sum(result) == 1 and all(x > 0 if positive else x >= 0 for x in result),
            'normalized positive law required' if positive else 'normalized nonnegative law required')
    return result


def kernel(edge):
    """Actual source, input rows and output columns, with both placement coins."""
    born.seam(edge)
    return cut.seam_population(edge)


def sparse_rows(K):
    require(type(K) is tuple and K and all(type(row) is tuple and len(row) == len(K) for row in K),
            'square exact kernel required')
    exact = tuple(tuple(born.rational(x) for x in row) for row in K)
    rows = tuple(tuple((j, x) for j, x in enumerate(row) if x) for row in exact)
    require(all(all(w >= 0 for _, w in row) and sum(w for _, w in row) == 1 for row in rows),
            'nonnegative row-stochastic kernel required')
    return rows


def push(p, K):
    p = law(p, size=len(K))
    result = [Q() for _ in K]
    for i, row in enumerate(sparse_rows(K)):
        for j, w in row:
            result[j] += p[i]*w
    return tuple(result)


def multiply(P, QK):
    require(len(P) == len(QK), 'kernel dimensions differ')
    a, b = sparse_rows(P), sparse_rows(QK)
    result = [[Q() for _ in P] for _ in P]
    for i, row in enumerate(a):
        for k, first in row:
            for j, second in b[k]:
                result[i][j] += first*second
    return tuple(map(tuple, result))


def random_scan(weights=None):
    edges = born.carrier_seams()
    weights = tuple(Q(1, len(edges)) for _ in edges) if weights is None else law(weights, size=len(edges))
    result = [[Q() for _ in BASIS] for _ in BASIS]
    for w, edge in zip(weights, edges):
        for i, row in enumerate(sparse_rows(kernel(edge))):
            for j, value in row:
                result[i][j] += w*value
    return tuple(map(tuple, result))


def log_rational(value):
    """Formal exact log(value), expanded at positive rational prime factors."""
    value = born.rational(value)
    require(value > 0, 'a positive rational logarithm argument is required')
    return _log_rational(value)


@lru_cache(maxsize=2048)
def _log_rational(value):
    coefficients = defaultdict(Q)
    for sign, integer in ((1, value.numerator), (-1, value.denominator)):
        for prime, exponent in s.factorint(integer).items():
            coefficients[int(prime)] += sign*int(exponent)
    return tuple(sorted((p, c) for p, c in coefficients.items() if c))


def combine(*terms):
    """Each term is (rational coefficient, formal logarithmic expression)."""
    result = defaultdict(Q)
    for coefficient, expression in terms:
        coefficient = born.rational(coefficient)
        for prime, value in expression:
            result[prime] += coefficient*value
    return tuple(sorted((p, c) for p, c in result.items() if c))


def entropy(p):
    p = law(p, size=len(p))
    return combine(*[(-x, log_rational(x)) for x in p if x])


def relative_entropy(p, q):
    p, q = law(p, size=len(p)), law(q, size=len(p))
    require(all(not x or y for x, y in zip(p, q)), 'relative entropy is infinite outside reference support')
    return combine(*[(x, log_rational(x/y)) for x, y in zip(p, q) if x])


def projection_entropy(p, edge):
    p = law(p)
    q = push(p, kernel(edge))
    change = combine((1, entropy(q)), (-1, entropy(p)))
    divergence = relative_entropy(p, q)
    require(change == divergence, 'native projection entropy identity failed')
    return change


def random_scan_entropy(p, weights=None):
    p = law(p)
    edges = born.carrier_seams()
    weights = tuple(Q(1, 30) for _ in edges) if weights is None else law(weights, size=30)
    images = [push(p, kernel(edge)) for edge in edges]
    q = push(p, random_scan(weights))
    first = combine(*[(w, relative_entropy(p, image)) for w, image in zip(weights, images) if w])
    mixing = combine(*[(w, relative_entropy(image, q)) for w, image in zip(weights, images) if w])
    change = combine((1, entropy(q)), (-1, entropy(p)))
    require(change == combine((1, first), (1, mixing)), 'random-scan entropy split lost its mixing term')
    return {'entropy_change': change, 'mean_local_dissipation': first, 'mixing_information': mixing}


def fluctuation_average(p, K, *, require_positive=True):
    """Exact exp(-sigma) average; sigma=log(p_x/q_y) for uniform reference.

    No logarithm of zero is evaluated. For support-deficient p, the returned
    supported forward average can be less than one and is explicitly separate.
    """
    require(type(require_positive) is bool, 'explicit Boolean support contract required')
    p = law(p, positive=require_positive, size=len(K))
    q = push(p, K)
    require(all(sum(row[j] for row in K) == 1 for j in range(len(K))),
            'uniform-reference fluctuation identity requires column sums one')
    value = sum((p[i]*w*q[j]/p[i] for i, row in enumerate(sparse_rows(K)) if p[i]
                 for j, w in row), Q())
    if require_positive:
        require(all(y > 0 for y in q) and value == 1, 'positive-support fluctuation equality failed')
    require(0 <= value <= 1, 'supported fluctuation average outside its bound')
    return value


def heat_bath_certificate():
    uniform = tuple(Q(1, N) for _ in BASIS)
    adjacency = [set() for _ in BASIS]
    rows = []
    for edge in born.carrier_seams():
        a, b = edge
        fibres = tuple(tuple(x[i] for i in range(12) if i not in edge)+(x[a]+x[b],) for x in BASIS)
        actual = kernel(edge)
        expected = cut.population_heat_bath(uniform, fibres)
        require(actual == expected, 'actual endpoint kernel is not its uniform conditional resampling')
        swap = []
        for i, state in enumerate(BASIS):
            moved = list(state)
            moved[a], moved[b] = moved[b], moved[a]
            j = INDEX[tuple(moved)]
            swap.append(j)
            adjacency[i].add(j)
        require(all(swap[swap[i]] == i for i in range(N)), 'seam exchange is not an involution')
        require(all(actual[i][j] == Q(int(i == j)+int(swap[i] == j), 2)
                    for i in range(N) for j in range(N)), 'hold/swap bridge failed')
        require(multiply(actual, actual) == actual, 'individual seam is not a conditional projection')
        require(push(uniform, actual) == uniform and all(actual[i][j] == actual[j][i]
                    for i in range(N) for j in range(N)), 'common reference or detailed balance failed')
        sizes = {label: fibres.count(label) for label in set(fibres)}
        rows.append({'seam': edge, 'singleton_fibres': sum(n == 1 for n in sizes.values()),
                     'two_point_fibres': sum(n == 2 for n in sizes.values()),
                     'exact_entries_checked': N*N})
    distances = [-1]*N
    distances[0] = 0
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j in adjacency[i]:
            if distances[j] < 0:
                distances[j] = distances[i]+1
                queue.append(j)
    require(min(distances) >= 0, 'combined native token graph is disconnected')
    scan = random_scan()
    square = multiply(scan, scan)
    witness = next((i, j) for i in range(N) for j in range(N) if square[i][j] != scan[i][j])
    return {'states': N, 'seams': len(rows), 'source_coin_branches': N*30*2,
            'local_heat_baths': rows, 'connected_token_vertices': sum(x >= 0 for x in distances),
            'token_root_eccentricity': max(distances), 'common_stationary_weight': '1/66',
            'random_scan_idempotent': False,
            'nonidempotence_witness': {'indices': witness, 'P': str(scan[witness[0]][witness[1]]),
                                     'P_squared': str(square[witness[0]][witness[1]])}}


def startup_control():
    """The full 78-state startup is not a bistochastic binary heat bath."""
    p = tuple(Q(1, 78) for _ in born.BASIS)
    q = [Q() for _ in born.BASIS]
    for i, state in enumerate(born.BASIS):
        for coin in (0, 1):
            final, _ = born.endpoint(state, born.native.E, coin)
            q[born.INDEX[final]] += p[i]/2
    q = tuple(q)
    delta = combine((1, entropy(q)), (-1, entropy(p)))
    require(delta == ((3, Q(-1, 26)),), 'full-carrier entropy-drop control changed')
    # A covering scan removes each collision when a seam first touches its port.
    covered = p
    for edge in born.carrier_seams():
        next_law = [Q() for _ in born.BASIS]
        for i, state in enumerate(born.BASIS):
            for coin in (0, 1):
                final, _ = born.endpoint(state, edge, coin)
                next_law[born.INDEX[final]] += covered[i]/2
        covered = tuple(next_law)
    require(all(w == 0 for state, w in zip(born.BASIS, covered) if max(state) == 2),
            'covering scan did not remove all collision support')
    return {'uniform_78_is_stationary': False, 'one_seam_entropy_change': '-log(3)/26',
            'one_covering_scan_collision_probability': '0',
            'binary_sector_reference_cannot_be_used_before_attachment': True}


def encoded(expression):
    return {str(prime): str(coefficient) for prime, coefficient in expression}


def report():
    require(strict_json.load(ROOT/'model/native_born_bridge/results.json')
            == json.loads(json.dumps(born.report())), 'upstream native record bridge receipt is stale')
    heat = heat_bath_certificate()
    p = tuple(Q(i+1, 2211) for i in range(N))
    u = tuple(Q(1, N) for _ in BASIS)
    require(entropy(p) == combine((1, log_rational(N)), (-1, relative_entropy(p, u))),
            'same-reference state entropy identity failed')
    chosen = born.carrier_seams()[:3]
    current, terms = p, []
    for edge in chosen:
        terms.append(projection_entropy(current, edge))
        current = push(current, kernel(edge))
    require(combine(*[(1, term) for term in terms]) == combine((1, entropy(current)), (-1, entropy(p))),
            'fixed schedule entropy telescoping failed')
    random = random_scan_entropy(p)
    scan = random_scan()
    require(fluctuation_average(p, scan) == 1, 'source fluctuation equality failed')
    edge = born.carrier_seams()[0]
    P = kernel(edge)
    moving = next(i for i in range(N) if P[i][i] == Q(1, 2))
    point = tuple(Q(i == moving) for i in range(N))
    support_control = fluctuation_average(point, P, require_positive=False)
    require(support_control == Q(1, 2), 'missing-support fluctuation control changed')
    tree = born.source_tree()
    full = tuple(p[INDEX[state]] if state in INDEX else Q() for state in born.BASIS)
    joint = tuple(born.execute(tree, full).values())
    tape_increment = combine((1, entropy(joint)), (-1, entropy(p)))
    require(tape_increment == ((2, Q(2)),), 'retained native tape entropy omitted a branch')
    inherited_path = 'model/native_equilibration/record_results.json'
    inherited = strict_json.load(ROOT/inherited_path)
    for path, expected in inherited['source_sha256'].items():
        require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected,
                f'stale inherited quantum thermodynamics source: {path}')
    paths = ['model/native_thermodynamic_bridge/certificate.py',
             'model/native_thermodynamic_bridge/test_certificate.py',
             'model/native_thermodynamic_bridge/PROOF.md',
             'model/native_born_bridge/certificate.py', 'model/native_born_bridge/results.json',
             'model/observer/native.py', 'model/observer/quantum_mixing.py',
             'model/native_equilibration/cut_completion.py',
             'model/native_equilibration/record_thermodynamics.py', inherited_path,
             'Lean/Thermodynamics/FluctuationTheorems.lean', 'Lean/Thermodynamics/CapFirstLaw.lean']
    return {'schema': 'oph.native_thermodynamic_bridge.v1', 'all_pass': True,
            'native_heat_bath': heat,
            'same_reference_optimization': {'state_optimizer': 'uniform on the supplied connected binary total-two sector',
                                          'transition_optimizer': 'uniform on each actual native seam fibre',
                                          'state_entropy_gap': 'D(p || uniform66)',
                                          'row_entropy_gap': 'D(row || uniform_fibre)',
                                          'full_78_state_claim': False},
            'entropy_checks': {'exact_logarithmic_arithmetic': True,
                               'fixed_schedule_telescoping_checked': True,
                               'random_scan': {k: encoded(v) for k, v in random.items()},
                               'strictly_positive_fluctuation_average': '1',
                               'zero_support_fluctuation_control': str(support_control),
                               'complete_two_attempt_tape_entropy_increment': '2 log(2)'},
            'startup_control': startup_control(),
            'inherited_quantum_record_thermodynamics': {'receipt': inherited_path,
                 'source_pins_verified': len(inherited['source_sha256']),
                 'fresh_full_quantum_execution_or_new_quantum_entropy_proof': False},
            'scope': {'native_classical_heat_bath_attachment_derived': True,
                      'binary_sector_and_fair_coin_measure_supplied': True,
                      'exogenous_fixed_schedule_or_independent_random_scan': True,
                      'arbitrary_record_feedback_has_same_unread_law': False,
                      'uniform_reference_modular_energy_is_constant': True,
                      'physical_energy_temperature_clock_or_cooling_derived': False,
                      'zero_temperature_or_refinement_limit_claim': False},
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = report()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:
        require(strict_json.load(path) == json.loads(json.dumps(result)), 'native thermodynamic bridge receipt stale')
    print(json.dumps({'all_pass': True, 'native_seams': 30, 'binary_states': N,
                      'physical_temperature_derived': False}))


if __name__ == '__main__':
    main()
