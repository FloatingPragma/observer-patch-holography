"""Exact finite preparation of the bounded local family's faithful reference.

Let E be the source channel after discarding the NEW preparation records, and
R_pi(X)=Tr(X) pi. For any preexisting untouched R, the half trace distance
between (E**n tensor id)(rho) and pi tensor rho_R is at most

    min(1, (65/2) r**n + (1-beta)**floor(n/m)).

The first term bounds surviving off-diagonal blocks, with the complete tape
retained before contractively discarding it. The second bounds the classical
channel by a certified Doeblin minorization of P**m. This is a half-diamond
bound on the source replacement channel. It is NOT a replacement bound on the
new preparation archive: that archive can remain correlated with the source.
Fixed source/record-independent selectors are essential. No source pulses,
physical service errors, or adaptive-feedback stationary law are supplied.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from math import gcd, lcm
from pathlib import Path

import numpy as np

from . import local_bias as local
from ..observer import native, record_recovery as records

ROOT = Path(__file__).resolve().parents[2]
DIMENSION = 66
MAX_POWER = 64
MAX_ATTEMPTS = 100_000


def _fraction(value, *, floor=0):
    """Validate original slots before arithmetic and detach mutable Fractions."""
    if type(value) is not F:
        raise ValueError('an original canonical bounded Fraction is required')
    n, d = value.numerator, value.denominator
    if (type(n) is not int or type(d) is not int or d <= 0
            or gcd(n, d) != 1 or not floor*d <= n <= d or d > 256):
        raise ValueError('an original canonical bounded Fraction is required')
    return F(n, d)


def _parameters(parameters):
    if type(parameters) is not local.Parameters:
        raise ValueError('the frozen local Parameters type is required')
    try:
        ports, seams = parameters.port_weights, parameters.seam_weights
        c = _fraction(parameters.c_weight)
        repair = _fraction(parameters.repair_probability, floor=F(1, 16))
    except AttributeError as error:
        raise ValueError('all original parameter fields are required') from error
    if (type(ports) is not tuple or len(ports) != 12
            or any(type(x) is not int or not 1 <= x <= 32 for x in ports)
            or gcd(*ports) != 1 or max(ports) > 4*min(ports)):
        raise ValueError('twelve primitive bounded integer port weights required')
    if (type(seams) is not tuple or len(seams) != 30
            or any(type(x) is not int or not 1 <= x <= 16 for x in seams)):
        raise ValueError('thirty bounded positive integer seam weights required')
    return ports, seams, c.numerator, c.denominator, repair.numerator, repair.denominator


def canonical_parameters():
    """Fresh C endpoint; the canonical certificate also covers E and mixtures."""
    return local.Parameters(port_weights=(25, 25)+(8,)*10,
                            c_weight=F(1), repair_probability=F(1))


def _count(value, maximum, label):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(f'{label} must be an integer in [0,{maximum}]')
    return value


def _native_payload():
    """Read the actual endpoint law before entering numerical-result caches.

    Keys include this immutable payload, so changed native maps cannot reuse
    previously cached integer calculations. Every semantic zero slot remains.
    Reserved physical raw addresses are the continuing instrument's concern.
    """
    states, seams = native.states(binary=True), records.seams()
    if len(states) != DIMENSION or len(set(states)) != DIMENSION or len(seams) != 30:
        raise ArithmeticError('the frozen binary carrier or seam menu changed')
    index = {x: i for i, x in enumerate(states)}
    rows = []
    for law in ('C', 'E'):
        for edge_index, edge in enumerate(seams):
            for coin in (0, 1):
                slots = {label: [] for label in (range(-1, 2) if law == 'C' else range(-2, 3))}
                for i, state in enumerate(states):
                    output, record = records.recorded_step(law, state, edge, coin)
                    if output not in index or record[:3] != (*edge, coin) or record[-1] not in slots:
                        raise ArithmeticError('native output or complete record left the frozen carrier')
                    slots[record[-1]].append((i, index[output]))
                for label, pairs in slots.items():
                    if len({j for _, j in pairs}) != len(pairs):
                        raise ArithmeticError('a complete record is not a partial injection')
                    rows.append((law, edge_index, coin, label, tuple(pairs)))
    return tuple(seams), tuple(rows)


@lru_cache(maxsize=32)
def _dynamics_payload(key, source):
    ports, seam_weights, cn, cd, rn, rd = key
    c, repair = F(cn, cd), F(rn, rd)
    seams, rows = source
    branches = []
    for law, edge_index, coin, _, pairs in rows:
        a, b = seams[edge_index]
        weight = (repair * F(seam_weights[edge_index], sum(seam_weights))
                  * (c if law == 'C' else 1-c)
                  * F(ports[a if coin else b], ports[a]+ports[b]))
        branches.append((weight, pairs))
    branches.append((1-repair, tuple((i, i) for i in range(DIMENSION))))
    denominator = lcm(*(w.denominator for w, _ in branches))
    population = [[0]*DIMENSION for _ in range(DIMENSION)]
    agreement = [[0]*DIMENSION for _ in range(DIMENSION)]
    for weight, pairs in branches:
        scaled = weight.numerator*(denominator//weight.denominator)
        for i, j in pairs:
            population[i][j] += scaled
            for other, _ in pairs:
                agreement[i][other] += scaled
    if (any(sum(row) != denominator for row in population)
            or any(agreement[i][i] != denominator for i in range(DIMENSION))):
        raise ArithmeticError('the complete native instrument is not trace preserving')
    return denominator, tuple(map(tuple, population)), tuple(map(tuple, agreement))


def probability_kernel(parameters):
    """Exact source population kernel; all returned scalars are detached."""
    denominator, population, _ = _dynamics_payload(_parameters(parameters), _native_payload())
    return tuple(tuple(F(x, denominator) for x in row) for row in population)


def coherence_rate(parameters):
    """Maximum common-record probability for two distinct input addresses."""
    denominator, _, agreement = _dynamics_payload(_parameters(parameters), _native_payload())
    return F(max(agreement[i][j] for i in range(DIMENSION)
                 for j in range(DIMENSION) if i != j), denominator)


@lru_cache(maxsize=16)
def _power_payload(population, denominator, block):
    # Exact integer matrix multiplication; no floating eigensolver or power.
    matrix = np.linalg.matrix_power(np.array(population, dtype=object), block)
    scale = denominator**block
    rows = tuple(tuple(int(x) for x in row) for row in matrix)
    if any(sum(row) != scale for row in rows):
        raise ArithmeticError('exact population power lost normalization')
    return rows, scale


def population_power(parameters, attempts):
    """Reconstructible integer numerator matrix and common denominator."""
    _count(attempts, MAX_POWER, 'exact matrix-power horizon')
    denominator, population, _ = _dynamics_payload(_parameters(parameters), _native_payload())
    return _power_payload(population, denominator, attempts)


def population_minorization(parameters, block=20):
    """P**block >= beta times the row-independent column-minimum law."""
    _count(block, MAX_POWER, 'minorization block')
    if block == 0:
        raise ValueError('a positive minorization block is required')
    rows, denominator = population_power(parameters, block)
    minima = tuple(min(rows[i][j] for i in range(DIMENSION)) for j in range(DIMENSION))
    return {'block': block, 'denominator': denominator, 'column_minima': minima,
            'beta': F(sum(minima), denominator)}


def mixing_bound(parameters, attempts, block=20):
    """Accessible source/preexisting-reference half-diamond upper bound.

    Exported preparation records are traced out in this comparison. They are
    not asserted independent, erased for free, or included in a reset target.
    """
    _count(attempts, MAX_ATTEMPTS, 'attempt count')
    minor = population_minorization(parameters, block)
    rate = coherence_rate(parameters)
    return min(F(1), F(DIMENSION-1, 2)*rate**attempts
               +(1-minor['beta'])**(attempts//block))


def _upper(value, denominator=10**12):
    """A small exact outward rational, not a rounded floating diagnostic."""
    return F(-(-value.numerator*denominator//value.denominator), denominator)


def mixing_bound_summary(parameters, attempts, block=20):
    bound = mixing_bound(parameters, attempts, block)
    minor = population_minorization(parameters, block)
    return {'attempts': attempts, 'population_block': block,
            'population_minorization': str(minor['beta']),
            'off_diagonal_rate': str(coherence_rate(parameters)),
            'half_diamond_upper': str(_upper(bound)),
            'new_preparation_archive_in_replacement_target': False}


def general_certificate():
    """Uniform, deliberately conservative constants for the whole class."""
    frozen = local.contract()
    boundary = frozen['parameter_boundary']
    if boundary != {'law_tags': ['C', 'E'], 'port_weights': {
            'count': 12, 'integer_max': 32, 'integer_min': 1,
            'max_to_min_ratio': 4, 'primitive_gcd': 1},
            'probability_max_denominator': 256, 'repair_probability_min': '1/16',
            'seam_weights': {'count': 30, 'integer_max': 16, 'integer_min': 1},
            'source_and_record_independent_selection': True}:
        raise ArithmeticError('uniform proof constants do not cover the changed class')
    seams, rows = _native_payload()
    adjacency = [set() for _ in range(DIMENSION)]
    for law, _, _, _, pairs in rows:
        if law == 'C':
            for i, j in pairs:
                if i != j:
                    adjacency[i].add(j)
    diameter = 0
    for start in range(DIMENSION):
        distances, todo = {start: 0}, deque((start,))
        while todo:
            i = todo.popleft()
            for j in adjacency[i]:
                if j not in distances:
                    distances[j] = distances[i]+1
                    todo.append(j)
        if len(distances) != DIMENSION:
            raise ArithmeticError('the binary configuration graph disconnected')
        diameter = max(diameter, max(distances.values()))
    if diameter != 4:
        raise ArithmeticError('the certified four-step population paths changed')
    uniform_c = local.Parameters(c_weight=F(1), repair_probability=F(1))
    uniform_e = local.Parameters(c_weight=F(0), repair_probability=F(1))
    if coherence_rate(uniform_c) != F(5, 6) or coherence_rate(uniform_e) != F(7, 10):
        raise ArithmeticError('the fair uniform comparison rates changed')
    # Minimum moving probability: repair * seam * endpoint coin. Self loops
    # have >=1/5, so paths of length <=4 can be padded to exactly four.
    step = F(1, 16)*F(1, 480)*F(1, 5)
    beta = DIMENSION*step**diameter
    # Every complete event contains this fraction of a fair/uniform tagged
    # instrument. The remainder has common-record probability at most one.
    comparison_fraction = F(1, 16)*F(1, 16)*F(2, 5)
    return {'dimension': DIMENSION, 'configuration_graph_diameter': diameter,
            'moving_transition_lower': str(step), 'self_loop_lower': '1/5',
            'population_block': diameter, 'population_minorization_lower': str(beta),
            'fair_uniform_instrument_fraction': str(comparison_fraction),
            'off_diagonal_rate_upper': str(1-comparison_fraction*F(1, 6)),
            'E_only_off_diagonal_rate_upper': str(1-comparison_fraction*F(3, 10)),
            'bound': 'min(1,65*r^n/2+(1-beta)^floor(n/4))',
            'unique_full_quantum_source_stationary_state': True,
            'fixed_source_and_record_independent_selector_required': True,
            'all_preexisting_untouched_references': True,
            'new_preparation_archive_in_replacement_target': False}


def canonical_certificate():
    parameters = canonical_parameters()
    minor = population_minorization(parameters)
    rate = coherence_rate(parameters)
    endpoint_e = local.Parameters(port_weights=parameters.port_weights,
                                  c_weight=F(0), repair_probability=F(1))
    e_rate = coherence_rate(endpoint_e)
    pi, kernel = local.reference(parameters), probability_kernel(parameters)
    if probability_kernel(endpoint_e) != kernel:
        raise ArithmeticError('C/E endpoint population kernels disagree')
    if (sum(pi) != 1 or min(pi) <= 0
            or tuple(sum(pi[i]*kernel[i][j] for i in range(DIMENSION))
                     for j in range(DIMENSION)) != pi):
        raise ArithmeticError('the canonical product reference is not faithful and stationary')
    if minor['beta'] < F(33, 200) or rate != F(893, 990) or e_rate != F(7, 10):
        raise ArithmeticError('the exact canonical minorization or coherence control changed')
    bound = F(65, 2)*rate**1600+F(167, 200)**80
    if not bound < F(1, 10**6):
        raise ArithmeticError('the canonical finite preparation bound failed')
    return {'parameters': parameters.description(), 'population_block': 20,
            'population_integer_denominator': minor['denominator'],
            'population_column_minima': list(minor['column_minima']),
            'population_minorization': str(minor['beta']),
            'population_minorization_lower': '33/200',
            'C_off_diagonal_rate': str(rate), 'E_off_diagonal_rate': str(e_rate),
            'any_fixed_tagged_C_E_mixture_covered': True,
            'attempts': 1600, 'half_diamond_upper': str(_upper(bound)),
            'target_error': '1/1000000',
            'minimum_reference_weight': str(min(pi)),
            'new_preparation_archive_in_replacement_target': False}


def source_bindings():
    paths = set(local.contract()['source_sha256']) | {
        'model/local_repair_family/mixing.py', 'model/local_repair_family/test_mixing.py',
        'model/local_repair_family/local_bias.py', 'model/local_repair_family/campaign.json'}
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def report():
    return {'schema': 'oph.local_repair_family.mixing.v1',
            'source_sha256': source_bindings(),
            'general_family': general_certificate(), 'canonical': canonical_certificate(),
            'proof': [
                'Every committed native branch is a partial injection. A distinct pair surviving a common history remains distinct; conditioning on the next event bounds its common-record weight by r.',
                'For positive source/reference blocks, sum of off-diagonal trace norms is at most65. Thus the complete-tape dephasing comparison has half trace error at most65*r^n/2; discarding the new tape contracts this distance.',
                'Column minima certify P^m >= beta*nu for a probability row nu. Dobrushin contraction gives max_i TV(P^n(i),pi) <= (1-beta)^floor(n/m), since pi is stationary.',
                'On a dephased source/reference state the diagonal reference blocks are positive, so that population bound is uniform over arbitrary preexisting untouched references. Triangle inequality proves the stated half-diamond replacement bound.',
                'The bound tends to zero for each fixed admitted member. Its product reference is therefore the unique quantum source stationary state; no stationary or independent full exported archive is inferred.'],
            'scope': {'native_source_pulses_executed': False,
                      'physical_service_errors_bounded': False,
                      'adaptive_feedback_stationarity_claimed': False,
                      'new_preparation_archive_in_replacement_target': False}}
