"""Two full-graph modes and complete ideal arc records, without producer imports.

Run from the repository root: python3 model/openwave_energy_clock_join/check.py
Only graph.json is consumed. No receipt, expected spectrum or producer is read.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

try:
    from . import bounds as iv
except ImportError:
    import bounds as iv

ROOT = Path(__file__).resolve().parents[2]
GRAPH = ROOT/'code/openwave_static_signed_source/evidence/openwave/shared/localized_rest_join/evidence/graph.json'
V, K, P, J = 240, 1260, 256, 2048
D = P*J
SUCCESS = F(1, V*K)
Q = sp.sqrt(5)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def current(d, c):
    return (d-(2*c-1)*(abs(d) % 2))//2


def laplacian(edges, size):
    out = sp.zeros(size)
    for a, b in edges:
        out[a, a] += 1
        out[b, b] += 1
        out[a, b] -= 1
        out[b, a] -= 1
    return out


def zero(vector):
    return all(sp.simplify(x) == 0 for x in vector)


def admitted(state):
    state = np.asarray(state)
    return (state.shape == (D,) and np.isfinite(state).all()
            and not np.any(state.reshape(P, J)[V:])
            and not np.any(state.reshape(P, J)[:V, K:]))


@lru_cache(maxsize=1)
def source(path=GRAPH):
    raw = path.read_bytes()
    data = json.loads(raw)
    edges = tuple(tuple(e) for e in data['seams'])
    require(data['ports'] == V and len(edges) == 630 and len(set(edges)) == 630,
            '240 ports and 630 distinct seams required')
    require(all(0 <= a < b < V for a, b in edges), 'seam domain')
    cell = tuple(e for e in edges if e[1] < 12)
    require(len(cell) == 30, 'thirty internal cell seams')
    for c in range(20):
        local = tuple((a % 12, b % 12) for a, b in edges if a//12 == b//12 == c)
        require(local == cell, 'identical ordered internal cells required')
    L = laplacian(edges, V)
    # All 302400 native event columns are covered: off-seam d=0 holds;
    # the four endpoint/coin cases below are the only nonidentity arcs.
    counts = K*np.eye(V, dtype=np.int64)
    reverse = np.arange(D, dtype=np.int64).reshape(P, J)
    for ei, (a, b) in enumerate(edges):
        for c in (0, 1):
            for p in (a, b):
                j = current(int(p == a)-int(p == b), c)
                if j:
                    q = b if j == 1 else a
                    reverse[p, 2*ei+c] = q*J+2*ei+1-c
                    counts[p, p] -= 1
                    counts[q, p] += 1
    reverse = reverse.ravel()
    require(current(0, 0) == current(0, 1) == 0, 'all off-seam holds')
    require(np.array_equal(reverse[reverse], np.arange(D)), 'S²=I on every padded column')
    require(sp.Matrix(counts) == K*sp.eye(V)-L, 'all-column compression A*SA=I-L/1260')
    require(max(L[i, i] for i in range(V)) == 6, 'positive spectrum bound needs dmax=6')
    require(all(sum(counts[:, i]) == K for i in range(V)), 'holds retained in every column')
    return edges, L, reverse, hashlib.sha256(raw).hexdigest()


@lru_cache(maxsize=1)
def modes(path=GRAPH):
    edges, L, _, _ = source(path)
    cell_L = laplacian(tuple(e for e in edges if e[1] < 12), 12)
    f = sp.zeros(V, 1)
    for p, value in zip((4, 5, 6, 7), (1, -1, -1, 1)):
        f[p] = value
    require(zero(L*f-6*f), 'actual compact anchor eigen-equation')
    B = sp.Matrix.hstack(*(cell_L-(5+Q)*sp.eye(12)).nullspace())
    require(B.shape == (12, 3), 'actual upper cell band dimension')
    external = tuple((a, b) for a, b in edges if a//12 != b//12)
    g, profiles = sp.zeros(V, 1), {}
    for c in range(20):
        boundary = tuple(sorted({p % 12 for e in external for p in e if p//12 == c}))
        require(len(boundary) == 3, 'three intercell ports')
        if boundary not in profiles:
            C = B[list(boundary), :]
            require(sp.simplify(C.det()) != 0, 'boundary interpolation is unique')
            profiles[boundary] = (B*C.inv()*sp.ones(3, 1)).applyfunc(sp.simplify)
        g[12*c:12*c+12, 0] = profiles[boundary]
    require(all(g[a] == g[b] == 1 for a, b in external), 'external seam differences vanish')
    require(zero(L*g-(5+Q)*g), 'second mode belongs to the full 240-port source')
    require(sp.simplify(sum(g)) == 0 and sum(f) == 0, 'zero-reference orthogonality')
    return f, g, {str(k): list(map(str, v)) for k, v in profiles.items()}


def embed(profile):
    values = np.array([float(x) for x in profile])
    values /= np.linalg.norm(values)*math.sqrt(K)
    out = np.zeros((P, J))
    out[:V, :K] = values[:, None]
    return out.ravel()


def reference():
    edges, _, reverse, _ = source()
    out = np.zeros(D)
    for a, b in ((0, 1), (1, 5), (5, 0)):
        ei = edges.index(tuple(sorted((a, b))))
        c = int(a > b)
        arc = a*J+2*ei+c
        require(reverse[arc]//J == b, 'triangle uses an actual moving arc')
        out[arc], out[reverse[arc]] = 1, -1
    require(np.array_equal(out[reverse], -out), 'Sv=-v')
    require(not np.any(out.reshape(P, J).sum(axis=1)), 'A*v=0')
    require(out@out == 6, 'six distinct reference arcs')
    return out/math.sqrt(6)


def walk(state):
    reverse = source()[2]
    out = -state.copy().reshape(P, J)
    a = state.reshape(P, J)
    out[:V, :K] += 2*a[:V, :K].mean(axis=1)[:, None]
    # Published padding extension: identity on invalid position addresses.
    out[V:] = a[V:]
    return out.ravel()[reverse]


def law(state):
    require(admitted(state), 'preparation outside physical arc sector')
    require(abs(np.vdot(state, state).real-1) < 1e-12, 'normalized successful preparation')
    return np.r_[1-float(SUCCESS), float(SUCCESS)*abs(state)**2]


def same_law(a, b):
    return (a.shape == b.shape == (D+1,) and np.isfinite(a).all() and np.isfinite(b).all()
            and min(a.min(), b.min()) >= 0 and abs(a.sum()-1) < 1e-12
            and abs(b.sum()-1) < 1e-12 and np.max(abs(a-b)) < 1e-14)


def same_conversion(first, second):
    # A scale is one mathematical parameter, shared before either comparison.
    return first == second


def algebraic_interval(value):
    value = sp.expand(value)
    b = value.coeff(Q)
    a = sp.simplify(value-b*Q)
    a, b = F(int(sp.numer(a)), int(sp.denom(a))), F(int(sp.numer(b)), int(sp.denom(b)))
    return iv.add((a, a), iv.mul((b, b), iv.sqrt_point(F(5))))


def one_word_observable_gap(profile, transfer, predicted_phase):
    """Certified full-protocol mass at coordinate(5,3), actual minus prediction."""
    require(source()[2][5*J+3]//J == 0, 'fixed witness reversal ends at port zero')
    norm = algebraic_interval((profile.T*profile)[0])
    denom = iv.sqrt(iv.mul((F(K), F(K)), norm))
    fp, fq = algebraic_interval(profile[5]), algebraic_interval(profile[0])
    u = iv.div(fp, denom)
    su = iv.div(fq, denom)
    sine = iv.sqrt(iv.sub((F(1), F(1)), iv.mul(transfer, transfer)))
    w = iv.div(iv.sub(su, iv.mul(transfer, u)), sine)
    v = iv.div((F(1), F(1)), iv.sqrt_point(F(6)))
    actual = iv.add(su, v)
    predicted = iv.add(iv.add(iv.mul(iv.cos(predicted_phase), u),
                              iv.mul(iv.sin(predicted_phase), w)), v)
    require(0 < actual[0] and 0 < predicted[0], 'positive witness amplitudes')
    gap = iv.mul(iv.sub(iv.mul(actual, actual), iv.mul(predicted, predicted)),
                 (SUCCESS/2, SUCCESS/2))
    require(gap[0] > 0, 'disjoint complete-protocol witness probabilities')
    return gap


def observed(profile, transfer, predicted_phase):
    u, v = embed(profile), reference()
    su = u[source()[2]]
    s = math.sqrt(1-transfer*transfer)
    w = (su-transfer*u)/s
    require(max(abs(u@v), abs(w@v), abs(u@w), abs(w@w-1)) < 2e-12, 'orthonormal reducing sector')
    require(np.max(abs(walk(v)-v)) < 1e-15, 'stationary reference')
    theta = math.acos(transfer)
    actual = (u+v)/math.sqrt(2)
    largest = 0.
    for n in range(17):
        predicted = (math.cos(n*theta)*u+math.sin(n*theta)*w+v)/math.sqrt(2)
        a, b = law(actual), law(predicted)
        require(same_law(a, b), 'complete original-weight coordinate law')
        largest = max(largest, float(np.max(abs(a-b))))
        actual = walk(actual)
    n = 3
    moving = math.cos(n*theta)*u+math.sin(n*theta)*w
    baseline = law((moving+v)/math.sqrt(2))
    double = law((math.cos(2*n*theta)*u+math.sin(2*n*theta)*w+v)/math.sqrt(2))
    reverse_phase = law((math.cos(n*theta)*u-math.sin(n*theta)*w+v)/math.sqrt(2))
    shifted_reference = law((moving+1j*v)/math.sqrt(2))
    dephased = (law((moving+v)/math.sqrt(2))+law((moving-v)/math.sqrt(2)))/2
    controls = {'double_frequency': double, 'reversed_phase': reverse_phase,
                'reference_quarter_turn': shifted_reference,
                'reference_dephased': dephased, 'reference_removed': law(moving)}
    controls['anchored_log_phase_prediction'] = law((math.cos(n*predicted_phase)*u
                                                   +math.sin(n*predicted_phase)*w+v)/math.sqrt(2))
    # On the calibration mode this prediction is the positive control.
    if abs(predicted_phase-theta) < 1e-12:
        require(same_law(baseline, controls.pop('anchored_log_phase_prediction')),
                'anchor predicts its own complete law')
    controls['failures_removed'] = np.r_[0., abs((moving+v)/math.sqrt(2))**2]
    deltas = {}
    for name, changed in controls.items():
        require(not same_law(baseline, changed), name+' must change the same complete-law predicate')
        deltas[name] = float(np.max(abs(baseline-changed)))
    aliased = law((math.cos(n*(theta+2*math.pi))*u+math.sin(n*(theta+2*math.pi))*w+v)/math.sqrt(2))
    require(same_law(baseline, aliased), 'integer-word 2pi alias is non-discriminant')
    witness = 5*J+3
    c1 = float(v[witness]*u[witness])
    require(abs(c1) > 1e-6, 'theta harmonic must be visible at the same arc')
    return {'theta_display': theta, 'witness_arc': witness, 'theta_cosine_coefficient_display': c1,
            'complete_labels': D+1, 'sampled_word_counts': [0, 16],
            'complete_probability_residual_display': largest,
            'negative_control_max_coordinate_differences_display': deltas,
            'integer_word_2pi_alias_non_discriminant': True}


def result():
    edges, _, _, digest = source()
    f, g, profiles = modes()
    sq5 = iv.sqrt_point(F(5))
    t1 = (F(209, 210), F(209, 210))
    t2 = iv.sub((F(1), F(1)), iv.div(iv.add((F(5), F(5)), sq5), (F(K), F(K))))
    e1, phase1 = iv.spectral(t1)
    e2, phase2 = iv.spectral(t2)
    k = iv.div(e1, phase1)
    prediction = iv.mul(k, phase2)
    residual = iv.sub(e2, prediction)
    require(iv.disjoint(e2, prediction) and residual[0] > 0, 'certified scalar-join failure')
    hidden_k = iv.div(e2, phase2)
    require(same_conversion(k, k) and not same_conversion(k, hidden_k), 'hidden recalibration control')
    # Positive source control: two distinct compact eigenrays share exactly the
    # same eigenvalue, so one calibration applies identically to both.
    other = sp.zeros(V, 1)
    for p, value in zip((16, 17, 18, 19), (1, -1, -1, 1)):
        other[p] = value
    require(zero(source()[1]*other-6*other) and zero(f.T*other), 'degenerate positive control')
    require(iv.spectral(t1) == (e1, phase1), 'identical spectral data accepted')
    predicted_phase2 = iv.div(e2, k)
    probability_gap = one_word_observable_gap(g, t2, predicted_phase2)
    records = [observed(f, float(t1[0]), float(sum(phase1)/2)),
               observed(g, float(sum(t2)/2), float(sum(predicted_phase2)/2))]
    bad = embed(f)
    bad[240*J] = 1
    require(not admitted(bad), 'outside-sector control')
    for branch in (-1, 1):
        try:
            iv.minus_log_point(t1[0], branch)
        except ValueError:
            pass
        else:
            raise ValueError('complex log branch admitted as real transfer energy')
    return {'diagnostic': 'static-two-mode-arc-diagnostic-v1', 'graph_sha256': digest,
            'source_edges': len(edges), 'L_eigenvalues': ['6', '5+sqrt(5)'],
            'second_mode_cell_profiles': profiles,
            'transfer_intervals': [iv.strings(t1), iv.strings(t2)],
            'energy_intervals': [iv.strings(e1), iv.strings(e2)],
            'phase_intervals': [iv.strings(phase1), iv.strings(phase2)],
            'anchor_scale_interval': iv.strings(k), 'second_energy_prediction_interval': iv.strings(prediction),
            'second_energy_residual_interval': iv.strings(residual),
            'second_energy_residual_display': list(map(float, residual)),
            'second_phase_prediction_interval': iv.strings(predicted_phase2),
            'one_word_coordinate_probability_gap_interval': iv.strings(probability_gap),
            'one_word_coordinate_probability_gap_display': list(map(float, probability_gap)),
            'single_conversion_refuted_on_two_modes_and_zero_reference': True,
            'positive_degenerate_mode_control': True, 'preparation_success_exact': str(SUCCESS),
            'preparation_fail_exact': str(1-SUCCESS), 'records': records,
            'native_attempt_duration_of_reflection_word': None,
            'second_reading': 'provisional; internal different-method crosscheck only',
            'physical_energy_clock_identification': False}


if __name__ == '__main__':
    print(json.dumps(result(), indent=2, allow_nan=False))
