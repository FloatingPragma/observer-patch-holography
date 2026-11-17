"""Finite, budgeted phase readout of the existing row-1 spatial field.

This is an additive decoder amendment, not a new source law or an official
row decision. It replaces the old |F_e| <= 3 sampling restriction with a
per-edge resolution proved sufficient for the declared interpolating field.
"""
import cmath
import math


def read_phase(geometry, current, region, *, deformed=True, sample_budget=100000):
    if len(current) != len(geometry['seams']) or any(type(j) is not int for j in current):
        raise ValueError('one actual integer current per seam is required')
    if type(sample_budget) is not int or sample_budget < 0:
        raise ValueError('sample_budget must be a nonnegative integer')
    region = tuple(region)
    if len(set(region)) != len(region) or any(type(v) is not int or not 0 <= v < 12 for v in region):
        raise ValueError('region must be a subset of the twelve ports')
    traversals = [pair for v in region for pair in geometry['dual_port_boundaries'][v]]
    edge_ids = sorted(set(e for e, _ in traversals))
    # n >= 8|j|+8 and |g_b-g_a| <= .62, |a| <= .23 imply
    # sup |phi'|/n <= (2*pi*|j|+.62+2*pi*.23)/(8*|j|+8) < pi/4.
    # Every sampled increment is therefore its unique principal increment.
    resolutions = {e: 8 * abs(current[e]) + 8 for e in edge_ids}
    required = sum(n + 1 for n in resolutions.values())
    if required > sample_budget:
        return {'status': 'resource_stop', 'phase_degree': None,
                'required_samples': required, 'sample_budget': sample_budget,
                'reason': 'complete readout exceeds the declared sample budget'}
    increments = {}
    max_increment = 0.0
    for e in edge_ids:
        a, b = geometry['dual_edges'][e]
        ga = .31 * math.sin(a + 1) if deformed else 0.0
        gb = .31 * math.sin(b + 1) if deformed else 0.0
        amp = .23 * math.cos(a + b + 1) if deformed else 0.0
        n = resolutions[e]
        values = [cmath.exp(1j * (ga + (gb - ga) * k/n +
                  2 * math.pi * current[e] * k/n + amp * math.sin(2 * math.pi * k/n)))
                  for k in range(n + 1)]
        angles = [cmath.phase(b * a.conjugate()) for a, b in zip(values, values[1:])]
        max_increment = max(max_increment, max(map(abs, angles), default=0.0))
        increments[e] = math.fsum(angles)
    if max_increment >= math.pi / 2:
        return {'status': 'numerical_failure', 'phase_degree': None,
                'required_samples': required, 'sample_budget': sample_budget,
                'reason': 'finite precision violated the certified sampling margin'}
    degree = math.fsum(sign * increments[e] for e, sign in traversals) / (2 * math.pi)
    return {'status': 'read', 'phase_degree': degree,
            'required_samples': required, 'sample_budget': sample_budget,
            'max_principal_increment': max_increment,
            'analytic_increment_bound': math.pi / 4,
            'resolutions': resolutions}
