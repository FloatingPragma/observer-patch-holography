"""Exact all-lag correction and controlled reduced native-memory dynamics.

The semigroup comparison uses stationary product preparation and discards
memory records. Complete native instruments define the underlying process.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import scipy.linalg as la
import sympy as s

from ..common_source import native_family_record_memory as memory

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
POLICIES = ('always', 'C', 'E', 'CE')
FLAGS = tuple(itertools.product((False, True), repeat=2))
SUPPORT = (2, 3, 9, 10)
RAYS = (2, 3, 10)


def clean(matrix):
    return s.Matrix(matrix).applyfunc(lambda x: s.factor(s.cancel(s.expand(x))))


def zero(matrix):
    return all(s.cancel(s.expand(x)) == 0 for x in matrix)


def entries(matrix):
    return [[str(x) for x in row] for row in s.Matrix(matrix).tolist()]


def number(matrix):
    return np.asarray(s.Matrix(matrix).evalf(30), complex)


def enabled(policy, c, e):
    if policy not in POLICIES:
        raise ValueError('all-record pulse duty must be one of the four declared policies')
    return policy == 'always' or policy == 'C' and c or policy == 'E' and e or policy == 'CE' and c and e


def source_transition(p, q, policy):
    data = memory.architecture()
    T, E = s.zeros(12), s.zeros(12)
    weights = ((1-p)*(1-q), (1-p)*q, p*(1-q), p*q)
    for old in range(12):
        for (c, e), w in zip(FLAGS, weights):
            new = memory._destination(old, c, e, 2, data['gC'], data['gE'])
            T[new, old] += w
            if enabled(policy, c, e):
                E[new, old] += w
    return clean(T), clean(E)


@lru_cache(maxsize=1)
def symbolic_response():
    p, q = s.symbols('p q', real=True)
    a, b = p*(1-q), p*q
    delta = 1-a*a+a*b
    T = s.Matrix([[1-a, 1-b, 1-a], [a, 0, 0], [0, b, a]])
    pi = s.Matrix([1-a, a*(1-a), a*b])/delta
    Pi = pi*s.ones(1, 3)
    R = (s.eye(3)+T-2*Pi)/delta
    assert zero(T*pi-pi) and s.cancel(sum(pi)-1) == 0
    assert zero((s.eye(3)-T)*R-(s.eye(3)-Pi))
    assert zero(R*(s.eye(3)-T)-(s.eye(3)-Pi))
    assert zero(R*pi) and zero(s.ones(1, 3)*R)
    transitions = {
        'always': T,
        'C': s.Matrix([[b, a, b], [a, 0, 0], [0, b, a]]),
        'E': s.Matrix([[q, q-b, q], [0, 0, 0], [0, b, 0]]),
        'CE': s.Matrix([[b, 0, b], [0, 0, 0], [0, b, 0]])}
    native = memory.architecture()
    P = [native['projectors'][k] for k in RAYS]
    assert zero(native['projectors'][9]-native['projectors'][10])
    C = s.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1]])
    policies = {}
    for policy, E in transitions.items():
        T12, E12 = source_transition(p, q, policy)
        T4, E4 = T12.extract(SUPPORT, SUPPORT), E12.extract(SUPPORT, SUPPORT)
        assert zero(C*T4-T*C) and zero(C*E4-E*C)
        w = clean(E*pi)
        K = clean(s.diag(*w)/2-w*w.T/2+E*R*s.diag(*w))
        D = clean(K+K.T)
        A = clean(-s.I/2*sum((K[i, j]*(P[i]*P[j]-P[j]*P[i])
                              for i in range(3) for j in range(3)), s.zeros(3)))
        Y = clean(sum((w[i]*P[i] for i in range(3)), s.zeros(3)))
        assert zero(A-A.H) and zero(A+A.T) and zero(A+A.conjugate())
        assert zero(Y-Y.T) and zero(Y-Y.conjugate())
        policies[policy] = dict(E=E, w=w, K=K, D=D, A=A, Y=Y)
    compact = s.I*a*b*(1-a)/(10*delta**2)*s.Matrix(
        [[0, 3+s.sqrt(5), 0], [-3-s.sqrt(5), 0, -2], [0, 2, 0]])
    assert zero(policies['always']['A']-compact)
    return dict(p=p, q=q, T=T, pi=pi, Pi=Pi, R=R, policies=policies,
                compact_potential=compact, denominator=delta)


def exact_parameters(p, q):
    p, q = s.sympify(p), s.sympify(q)
    if (any(x.has(s.Float) or x.free_symbols or x.is_real is not True for x in (p, q))
            or p.is_nonnegative is not True or (1-p).is_positive is not True
            or q.is_nonnegative is not True or (1-q).is_nonnegative is not True):
        raise ValueError('exact 0<=p<1 and 0<=q<=1 required')
    return p, q


@lru_cache(maxsize=64)
def response(p=s.Rational(1, 12), q=s.Rational(1, 6), policy='always'):
    p, q = exact_parameters(p, q)
    if policy not in POLICIES:
        raise ValueError('unknown pulse duty')
    data = symbolic_response()
    sub = {data['p']: p, data['q']: q}
    values = {key: clean(data[key].subs(sub)) for key in ('T', 'pi', 'R')}
    values.update({key: clean(v.subs(sub)) for key, v in data['policies'][policy].items()})
    return values


def exact_rational(value):
    if isinstance(value, (float, np.floating, bool)):
        raise ValueError('exact rational bound parameters required')
    try:
        return F(value)
    except (TypeError, ValueError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('exact rational bound parameters required') from exc


def finite_time_bound(steps, time, p, duty, gain=1):
    """Exact bound conditional on supplied upper bounds for reset p, duty and gain.

    The caller must establish that these dominate its actual source parameters.
    The native producer supplies the exact admitted rates and duties. These
    are not inferred from a measured spectrum or from a numerical norm test.
    """
    if type(steps) is not int or steps <= 0:
        raise ValueError('positive integer refinement count required')
    time, p, duty, g = map(exact_rational, (time, p, duty, gain))
    if not time >= 0 or not 0 <= p < 1 or not 0 <= duty <= 1 or not g >= 0:
        raise ValueError('bound parameters outside their declared domain')
    h, u = time/steps, 1-p
    alpha, beta = 2*g*duty, 2*g*g*duty
    c = 2*alpha/u
    b = beta+alpha*c
    gamma = b+alpha*alpha/2
    d2 = (2*beta+3*alpha*c)/u
    mu = 4*g**3*duty/3
    e3 = alpha*gamma+h*gamma*gamma/2+(alpha+h*gamma)**3/6
    r3 = 2*alpha*d2+(beta+b)*c+h*(beta+b)*d2+(mu+e3)*(1+h*c+h*h*d2)
    initial = alpha*h*h*(c+h*d2)/u
    propagated = steps*h**3*r3
    return {'unhalved_diamond_bound_exact': str(min(F(2), initial+propagated)),
            'unclipped_bound_exact': str(initial+propagated),
            'initial_slip_bound_exact': str(initial),
            'intertwining_defect_bound_exact': str(propagated),
            'h_exact': str(h), 'p_upper_exact': str(p), 'duty_upper_exact': str(duty),
            'gain_upper_exact': str(g),
            'scope': 'Unread field with arbitrary reference, stationary product preparation; not a full-record or conditioned channel bound.'}


def coefficient_tail_bound(steps, p, duty, gain=1):
    if type(steps) is not int or steps <= 0:
        raise ValueError('positive integer horizon required')
    p, duty, g = map(exact_rational, (p, duty, gain))
    if not 0 <= p < 1 or not 0 <= duty <= 1 or g < 0:
        raise ValueError('coefficient bound outside declared domain')
    return 8*g*g*duty*duty*(p**(steps-1)/(1-p)+1/(steps*(1-p)**2))


def lindblad_generator(values):
    P = [number(memory.architecture()['projectors'][k]) for k in RAYS]
    X = np.array([[0, 1], [1, 0]], complex)
    H = [np.kron(X, v) for v in P]
    L = np.asarray([-1j*(np.kron(v, np.eye(6))-np.kron(np.eye(6), v.T)) for v in H])
    w, K, D = number(values['w']).real[:, 0], number(values['K']).real, number(values['D']).real
    G0 = sum((w[i]*L[i] for i in range(3)), np.zeros((36, 36), complex))
    G1 = sum((K[i, j]*L[i]@L[j] for i in range(3) for j in range(3)), np.zeros((36, 36), complex))
    potential = np.kron(np.eye(2), number(values['A']))
    coherent = -1j*(np.kron(potential, np.eye(6))-np.kron(np.eye(6), potential.T))
    dissipative = sum((D[i, j]*L[i]@L[j]/2 for i in range(3) for j in range(3)), np.zeros((36, 36), complex))
    assert la.norm(G1-coherent-dissipative, 2) < 2e-12
    return G0, G1, L, H


def transfer(T, E, H, h):
    n, d = T.shape[0], H[0].shape[0]**2
    pulses = [la.expm(-1j*h*v) for v in H]
    return np.block([[(T[i, j]-E[i, j])*np.eye(d)+E[i, j]*np.kron(pulses[i], pulses[i].conj())
                      for j in range(n)] for i in range(n)])


def normalized_choi(superoperator):
    return superoperator.reshape(6, 6, 6, 6).transpose(2, 0, 3, 1).reshape(36, 36)/6


def choi_error(a, b):
    value = normalized_choi(a-b)
    assert la.norm(value-value.conj().T, 2) < 5e-12
    return float(np.sum(abs(la.eigvalsh((value+value.conj().T)/2))))


def source_join():
    native = memory.architecture()
    p, q = native['p'], native['q']
    rows = {}
    C = s.zeros(3, 12)
    C[0, 2], C[1, 3], C[2, 9], C[2, 10] = 1, 1, 1, 1
    for policy in POLICIES:
        T, E = source_transition(p, q, policy)
        assert T == native['T']
        enabled_raw = s.zeros(12)
        for record in native['records']:
            _, _, ec, er, c, weight = record
            e = ec == 0 and er in (4, 5)
            for old in range(12):
                new = memory._destination(old, c, e, 2, native['gC'], native['gE'])
                if enabled(policy, c, e):
                    enabled_raw[new, old] += weight
        assert E == enabled_raw
        Pi = native['stationary']*s.ones(1, 12)
        R = (s.eye(12)-T+Pi).inv()-Pi
        w = E*native['stationary']
        K = s.diag(*w)/2-w*w.T/2+E*R*s.diag(*w)
        compact = response(p, q, policy)
        assert zero(C*K*C.T-compact['K'])
        assert zero(C*native['stationary']-compact['pi'])
        rows[policy] = {'T_source_exact': True, 'E_raw_records_exact': True,
                        'full12_lag_sum_equals_three_ray_exact': True}
    # Direct actual source channel, independent of symbolic correlation formulas.
    h = float(memory.READOUT_STEPS[0])
    P = [number(v) for v in native['projectors']]
    H12 = [np.kron(np.array([[0, 1], [1, 0]]), v) for v in P]
    full = transfer(number(native['T']).real, number(native['T']).real, H12, h)
    residual = float(la.norm(full-memory.joint_channel(h), 2))
    assert residual < 3e-12
    return {'native_raw_labels': 128, 'all_old_memory_inputs': 12,
            'full_CQ_dimension': 432, 'policies': rows,
            'actual_source_transfer_residual_diagnostic': residual,
            'lumping_scope': 'Stationary reduced experiment with identical P9=P10; full old-memory tape is not identified by this three-ray reduction.'}


def convergence_cases():
    rows = []
    for policy in POLICIES:
        values = response(policy=policy)
        T, E = number(values['T']).real, number(values['E']).real
        pi = number(values['pi']).real[:, 0]
        J, L = np.kron(pi[:, None], np.eye(36)), np.kron(np.ones((1, 3)), np.eye(36))
        G0, G1, _, H = lindblad_generator(values)
        d = {'always': F(1), 'C': F(1, 12), 'E': F(1, 6), 'CE': F(1, 72)}[policy]
        covariance_min = float(np.linalg.eigvalsh(number(values['D']).real).min())
        assert covariance_min > -2e-13
        for n in (8, 16, 32, 64, 128):
            h = 1/n
            exact = L@np.linalg.matrix_power(transfer(T, E, H, h), n)@J
            corrected, mean = la.expm(G0+h*G1), la.expm(G0)
            error, original = choi_error(exact, corrected), choi_error(exact, mean)
            bound = finite_time_bound(n, F(1), F(1, 12), d)
            assert error <= float(F(bound['unhalved_diamond_bound_exact']))+3e-12
            assert np.linalg.eigvalsh((normalized_choi(corrected)+normalized_choi(corrected).conj().T)/2).min() > -3e-12
            rows.append({'policy': policy, 'steps': n,
                         'mean_unitary_normalized_Choi_trace_error': original,
                         'corrected_semigroup_normalized_Choi_trace_error': error,
                         'covariance_min_eigenvalue_diagnostic': covariance_min,
                         **bound})
    return rows


def build():
    data = symbolic_response()
    native = memory.architecture()
    sub = {data['p']: native['p'], data['q']: native['q']}
    prefactor = s.simplify(data['compact_potential'][0, 1].subs(sub)/(s.I*(3+s.sqrt(5))))
    assert prefactor == s.Rational(603, 6666724)
    symbolic = {policy: {key: entries(value) for key, value in row.items()}
                for policy, row in data['policies'].items()}
    exact_native = {policy: {key: entries(value.subs(sub)) for key, value in row.items()}
                    for policy, row in data['policies'].items()}
    paths = ['model/common_source/native_family_record_memory.py',
             'model/common_source/native_repair_mass.py',
             'model/common_source/native_gauge_family_records.py',
             'model/common_source/native_gauge_object.py',
             'model/particle_memory_response/certificate.py']
    return {'schema': 'oph.particle_memory_response.v1',
            'domain': '0<=p<1, 0<=q<=1, all four declared duties; h>=0 and stationary product initial preparation for the reduced semigroup bound.',
            'T_symbolic': entries(data['T']), 'pi_symbolic': entries(data['pi']),
            'group_inverse_symbolic': entries(data['R']),
            'all_lag_sum_formula': 'K=diag(w)/2-ww^T/2+E R diag(w), R=sum_(k>=0)(T^k-Pi)',
            'symbolic_policies': symbolic, 'native_policies_exact': exact_native,
            'native_always_ordering_prefactor_exact': str(prefactor),
            'native_always_ordering_norm_exact': str(prefactor*s.sqrt(18+6*s.sqrt(5))),
            'native_always_ordering_norm_diagnostic': float(prefactor*s.sqrt(18+6*s.sqrt(5))),
            'source_join': source_join(), 'convergence_diagnostics': convergence_cases(),
            'scope': {'all_memory_lags_summed': True, 'positive_diffusion_retained': True,
                      'finite_time_complete_reference_error_bound': True,
                      'full_record_semigroup_equivalence': False,
                      'physical_mass_or_CP_identification': False,
                      'finite_h_complete_channel_logarithm_asserted': False},
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = build(); out = HERE/'results.json'
    if args.check:
        def same(a, b):
            if type(a) is float:
                assert np.isclose(a, b, atol=2e-12, rtol=1e-9)
            elif isinstance(a, dict):
                assert set(a) == set(b)
                for k in a: same(a[k], b[k])
            elif isinstance(a, list):
                assert len(a) == len(b)
                for x, y in zip(a, b): same(x, y)
            else: assert a == b
        same(result, json.loads(out.read_text()))
        print('particle memory response: receipt verified')
    else:
        out.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n');print(out)


if __name__ == '__main__':
    main()
