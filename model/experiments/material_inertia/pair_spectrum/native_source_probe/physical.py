"""Independent five-qutrit standards and exact physical routing of two targets."""

from fractions import Fraction
from itertools import product

import numpy as np
from scipy.linalg import expm

from model.experiments.interval_echo_consensus.signed_effects import fixed_pointer as pointer

WORDS = np.array(list(product(range(3), repeat=5)))


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def endpoint_word(path):
    forward = list(zip(path, path[1:]))
    return forward + forward[-2::-1]


def validate_route(board, path, source_a, source_b):
    check(not set(path) & {source_a, source_b}, 'route touches the source')
    check(all(b in board.neighbours[a] for a, b in zip(path, path[1:])), 'non-native route edge')
    labels = np.arange(board.n)
    for a, b in endpoint_word(path):
        labels[a], labels[b] = labels[b], labels[a]
    expected = np.arange(board.n)
    a, b = path[0], path[-1]
    expected[a], expected[b] = expected[b], expected[a]
    check(np.array_equal(labels, expected), 'unknown-data endpoint transposition fails')
    return len(endpoint_word(path))


def primitive(word, state):
    # Logical axes: target1,target2,balancing_reference,source_park,source_cross.
    out = np.array(state, complex, copy=True)
    for kind, a, b, ratio in word:
        time = float(ratio) * np.pi
        if kind == 'reference':
            out *= np.exp(-8j * time * WORDS[:, a]**2)[:, None]
            continue
        h = np.zeros((9, 9))
        for n, (left, right) in enumerate(product(range(3), repeat=2)):
            if left != right:
                h[n, n] = 8 * (left - right)**2 + 0.5
                h[3 * right + left, n] = -0.5
        u = expm(-1j * time * h)
        rest = [i for i in range(5) if i not in (a, b)]
        axes = [a, b, *rest, 5]
        shaped = out.reshape(*([3] * 5), out.shape[1]).transpose(axes)
        shaped = (u @ shaped.reshape(9, -1)).reshape(3, 3, *([3] * 3), out.shape[1])
        out = shaped.transpose(np.argsort(axes)).reshape(243, out.shape[1])
    return out


def audit(board, profiles, meters):
    route_count = 0
    for profile in profiles:
        app = profile['apparatus']
        for path in app['endpoint_transposition_paths']:
            route_count += validate_route(board, path, app['source_park'], app['source_cross'])
    profile = next(p for p in profiles if p['center'] == 0 and p['radius'] == 1)
    app = profile['apparatus']
    p1, p2 = app['target_ports'][:2]
    outside = min(p for p in range(board.n) if p // 12 not in profile['owners']
                  and p not in (app['source_park'], app['source_cross'], app['interaction_port']))
    physical_ports = [p1, p2, outside, app['source_park'], app['source_cross']]
    inputs = [[0, 1, 2, 2, 0], [1, 1, 1, 2, 0], [2, 1, 0, 2, 0],
              [2, 0, 1, 2, 0], [0, 2, 1, 2, 0]]
    region_charges = [-1, 0, 1, 0, 0]
    source = np.zeros((243, len(inputs)), complex)
    for i, row in enumerate(inputs):
        source[int(np.asarray(row) @ [81, 27, 9, 3, 1]), i] = 1
    H = pointer.neutral_basis(3, 4)
    prepared = primitive(H, source)
    maximum_error, records = 0.0, []
    for multiplier in (1, 2):
        field = []
        for target in (0, 1):
            field += pointer.diagonal_word(target, 4, Fraction(8 * multiplier, 31))
        driven = primitive(field, prepared)
        for axis in ('X', 'Y'):
            analyzer = ([] if axis == 'X' else [pointer.pulse('reference', 4, 4, Fraction(1, 64))])
            analyzer += pointer.inverse(H)
            final = primitive(analyzer, driven)
            probabilities = []
            for col, q in enumerate(region_charges):
                p = np.bincount(3 * WORDS[:, 3] + WORDS[:, 4], weights=abs(final[:, col])**2, minlength=9)
                error = float(max(abs(p - meters[multiplier][q][axis])))
                maximum_error = max(maximum_error, error)
                probabilities.append(p.tolist())
            records.append(dict(multiplier=multiplier, axis=axis, complete_meter_probabilities=probabilities))
    check(maximum_error < 2e-11, 'physical compensated standards disagree with native meter calibration')
    return dict(profile=dict(center=0, radius=1), physical_ports=physical_ports,
                physical_inputs=inputs, region_charges=region_charges,
                logical_order='target1,target2,untouched outside balancing rail,source park,source cross',
                routed_endpoint_swap_count=route_count,
                all_routes_checked=len([path for p in profiles for path in p['apparatus']['endpoint_transposition_paths']]),
                maximum_physical_standard_error=maximum_error, standards=records,
                scope='Five-qutrit quantum pulse replay plus exact all-qutrit route permutations; other target hosts add one common phase')
