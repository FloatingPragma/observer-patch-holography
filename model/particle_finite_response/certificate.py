"""Exact finite identities, native source joins, and numerical response controls.

The symbolic statements concern declared finite instruments. No measured
particle parameter, calibrated mass producer, or continuum identification is
used. See PROOF.md for domains and proofs of the universal statements.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TOLERANCE = 4e-12


def clean(matrix):
    return s.Matrix(matrix).applyfunc(s.simplify)


def entries(matrix):
    return [[str(v) for v in row] for row in s.Matrix(matrix).tolist()]


def _rational_matrix(matrix):
    matrix = s.Matrix(matrix)
    if any(v.is_Rational is not True for v in matrix):
        raise ValueError('exact rational matrix entries required')
    return matrix


def _instrument(maps):
    if not maps:
        raise ValueError('at least one complete outcome required')
    matrices = {label: _rational_matrix(value) for label, value in maps.items()}
    n = next(iter(matrices.values())).rows
    if any(m.shape != (n, n) or any(v < 0 for v in m) for m in matrices.values()):
        raise ValueError('nonnegative square outcome maps of one dimension required')
    if s.ones(1, n)*sum(matrices.values(), s.zeros(n)) != s.ones(1, n):
        raise ValueError('all outcomes must sum to a column-stochastic channel')
    return matrices


def observation_intertwining(coarse, large, small):
    """Certify complete-history comparison, retaining each outcome separately.

    Returns the exact combined induced-l1 defect. If it vanishes, compatible
    initial distributions give identical joint records at every finite length.
    Otherwise their total-variation distance is at most (initial_l1+n*defect)/2.
    """
    large, small = _instrument(large), _instrument(small)
    if set(large) != set(small):
        raise ValueError('both instruments must retain the same outcome alphabet')
    coarse = _rational_matrix(coarse)
    nl = next(iter(large.values())).rows
    ns = next(iter(small.values())).rows
    if coarse.shape != (ns, nl) or any(v < 0 for v in coarse):
        raise ValueError('a nonnegative coarse-graining map is required')
    if s.ones(1, ns)*coarse != s.ones(1, nl):
        raise ValueError('coarse graining must preserve total probability')
    defects = [coarse*large[y]-small[y]*coarse for y in large]
    defect = max(sum(abs(m[i, j]) for m in defects for i in range(ns)) for j in range(nl))
    return {'exact': defect == 0, 'defect_l1': defect,
            'outcome_count': len(large), 'large_dimension': nl, 'small_dimension': ns}


def word_probability(instrument, initial, word):
    instrument = _instrument(instrument)
    vector = _rational_matrix(initial)
    n = next(iter(instrument.values())).rows
    if vector.shape != (n, 1) or any(v < 0 for v in vector) or sum(vector) != 1:
        raise ValueError('a normalized nonnegative initial column is required')
    for y in word:
        vector = instrument[y]*vector
    return sum(vector)


def total_variation_bound(initial_l1, defect_l1, steps):
    if type(steps) is not int or steps < 0:
        raise ValueError('a nonnegative integer history length is required')
    initial_l1, defect_l1 = s.sympify(initial_l1), s.sympify(defect_l1)
    if any(x.is_Rational is not True or x < 0 for x in (initial_l1, defect_l1)):
        raise ValueError('nonnegative rational defects required')
    return min(s.Integer(1), (initial_l1+steps*defect_l1)/2)


def observation_certificate():
    small = {'0': s.Matrix([[0, 0], [1, 0]]), '1': s.Matrix([[0, 1], [0, 0]])}
    # One hidden bit changes at each step. Its capacity changes no record.
    flip = s.Matrix([[0, 1], [1, 0]])
    large = {y: s.kronecker_product(m, flip) for y, m in small.items()}
    coarse = s.kronecker_product(s.eye(2), s.ones(1, 2))
    result = observation_intertwining(coarse, large, small)
    assert result['exact']
    initial = s.Matrix([s.Rational(1, 3), 0, 0, s.Rational(2, 3)])
    checked = 0
    for n in range(7):
        for word in itertools.product(small, repeat=n):
            assert word_probability(large, initial, word) == word_probability(small, coarse*initial, word)
            checked += 1
    # The unread channel is identical while the retained outcome is opposite.
    first = {'0': s.eye(1), '1': s.zeros(1)}
    second = {'0': s.zeros(1), '1': s.eye(1)}
    defect = observation_intertwining(s.eye(1), first, second)
    assert sum(first.values(), s.zeros(1)) == sum(second.values(), s.zeros(1))
    assert defect['defect_l1'] == 2 and not defect['exact']
    assert total_variation_bound(0, 2, 1) == 1
    return {'all_history_theorem': 'outcome-wise intertwining plus compatible initial state',
            'exact_refinement': {k: str(v) if isinstance(v, s.Basic) else v for k, v in result.items()},
            'finite_words_checked': checked,
            'unread_channel_equality_counterexample_tv_at_one_step': 1,
            'determinism_alone_implies_equal_records': False,
            'error_bound': 'TV(records_n) <= min(1,(initial_l1+n*defect_l1)/2)'}


def phase_bound():
    a = (1-s.Rational(95, 96)*s.sqrt(s.Rational(653, 848)))/2
    bound = s.simplify(a*a/(1-2*a))
    return {'domain': '0<=p<=1/12, 0<h<=1/2; unshifted reference extension',
            'omega_upper_exact': str(a), 'omega_upper': float(a),
            'normal_adjacent_squared_gap_ratio_upper_exact': str(bound),
            'normal_adjacent_squared_gap_ratio_upper': float(bound),
            'inverted_squared_gap_ratio_lower': float(1-(a/(1-a))**2),
            'attachment': '(m1,m2,m3)=M*(0,omega,1-omega), M>0, normal ordering',
            'physical_particle_mass_theorem': False}


def neutral_roots(p, h):
    p, h = float(p), float(h)
    if not math.isfinite(p) or not math.isfinite(h) or not 0 <= p <= 1/12 or not 0 < h <= 1/2:
        raise ValueError('the certified domain is 0<=p<=1/12 and 0<h<=1/2')
    angle = h/2
    y = (1-p)*math.cos(angle)+1j*math.sqrt((1-p)**2+4*p/5)*math.sin(angle)
    slow = (y+np.sqrt(y*y+4*p))/2
    fast = -p/slow
    return np.exp(-1j*angle)*np.array([slow, np.conj(slow), fast, np.conj(fast)])


def neutral_certificate():
    from model.common_source import native_family_record_memory as memory
    native = memory.architecture()
    assert (native['p'], native['q'], native['coincidence']) == (s.Rational(1, 12), s.Rational(1, 6), s.Rational(1, 72))
    p2, p3 = native['projectors'][2:4]
    assert clean(p2*p2-p2) == s.zeros(3) and clean(p3*p3-p3) == s.zeros(3)
    assert s.simplify(s.trace(p2*p3)) == s.Rational(1, 5)
    assert clean(p2+p3).rank() == 2
    neutral = memory.architecture_from_factors(s.Rational(1, 12), 0)
    transition = neutral['T']
    recurrent, transient = [2, 3], [i for i in range(12) if i not in (2, 3)]
    assert transition.extract(recurrent, recurrent) == s.Matrix([[s.Rational(11, 12), 1], [s.Rational(1, 12), 0]])
    assert transition.extract(transient, recurrent) == s.zeros(10, 2)
    assert all(sum(transition[i, j] for i in transient) <= s.Rational(1, 12) for j in transient)
    p, z, lam = s.symbols('p z lambda')
    P, Q = s.diag(1, 0), s.Matrix([[1, 2], [2, 4]])/5
    U, V = s.eye(2)+(z-1)*P, s.eye(2)+(z-1)*Q
    F = ((1-p)*U).row_join(U).col_join((p*V).row_join(s.zeros(2)))
    polynomial = s.expand(F.charpoly(lam).as_expr())
    claimed = (lam**4-(1-p)*(1+z)*lam**3
               +((1-p)**2*z-p*(2*z+(z-1)**2/5))*lam**2
               +p*(1-p)*z*(1+z)*lam+p*p*z*z)
    assert s.expand(polynomial-claimed) == 0
    assert s.factor(polynomial.subs(z, 1)) == (lam-1)**2*(lam+p)**2
    # The determinant and trace determine the static bright spectrum.
    mean = (P+p*Q)/(1+p)
    assert s.simplify(mean.det()-4*p/(5*(1+p)**2)) == 0
    assert s.simplify(s.trace(mean)-1) == 0
    u, v = s.symbols('u v')
    C, B = (1-p)**2+4*p/5, (1+p)**2+4*p*v/5
    reduced = 4*p*u*u-B*u+C*v
    assert s.factor(reduced.subs(u, v)) == 16*p*v*(v-1)/5
    # Rational endpoint inequalities in the proof, separate from sampling.
    assert C.subs(p, s.Rational(1, 12)) == s.Rational(653, 720)
    assert B.subs({p: s.Rational(1, 12), v: s.Rational(1, 16)}) == s.Rational(53, 45)
    controls, maximum = [], 0.
    for pv in (0., 1/120, 1/24, 1/12):
        for h in (1/16, 1/8, 1/4, 1/2):
            roots = neutral_roots(pv, h)
            concrete = np.array(F.subs({p: pv, z: np.exp(-1j*h)}), complex)
            numerical = np.linalg.eigvals(concrete)
            error = max(min(abs(root-numerical)) for root in roots)
            assert error < TOLERANCE
            slow = roots[:2]
            omega = min(-np.angle(slow)/h)
            assert -TOLERANCE <= omega <= phase_bound()['omega_upper']+TOLERANCE
            assert max(abs(slow)) <= 1+TOLERANCE and min(abs(slow)) > .5
            assert max(abs(roots[2:])) < .5
            maximum = max(maximum, error)
            controls.append({'p': pv, 'h': h, 'omega': float(omega), 'slow_modulus': float(abs(slow[0]))})
    return {'source_event_probabilities': ['1/12', '1/6'], 'native_projector_overlap_exact': '1/5',
            'bright_characteristic_polynomial': str(claimed),
            'static_bright_trace': '1', 'static_bright_determinant': '4*p/(5*(1+p)^2)',
            'bound': phase_bound(), 'controls': controls, 'max_eigenvalue_residual': maximum}


def detector_certificate():
    # Canonical bright plane plus its exact dark line.
    P, Q = s.diag(1, 0, 0), s.Matrix([[1, 2, 0], [2, 4, 0], [0, 0, 0]])/5
    J = s.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
    B, D = s.diag(1, 1, 0), s.diag(0, 0, 1)
    X = s.Matrix([[0, 1], [1, 0]])
    plus, minus = s.Matrix([1, 1])/s.sqrt(2), s.Matrix([1, -1])/s.sqrt(2)
    A = s.kronecker_product(plus*minus.T, J)
    O = A+A.H
    PB, PD = s.kronecker_product(s.eye(2), B), s.kronecker_product(s.eye(2), D)
    assert O*O == PB
    for R in (P, Q):
        assert R*J+J*R == J
    effects = [(PB+O)/2, (PB-O)/2, PD]
    assert sum(effects, s.zeros(6)) == s.eye(6)
    assert all(E*E == E for E in effects)
    assert all(a*b == s.zeros(6) for i, a in enumerate(effects) for b in effects[:i])
    n, h = 7, .125
    p = 1/12
    blocks = [np.array(PB/4, complex)/(1+p), p*np.array(PB/4, complex)/(1+p)]
    matrices = [np.array(R, complex) for R in (P, Q)]
    unitaries = [np.eye(6)+(math.cos(h)-1)*np.kron(np.eye(2), R)-1j*math.sin(h)*np.kron(np.array(X, complex), R) for R in matrices]
    effects_n = list(map(lambda m: np.array(m, complex), effects))
    correlation = 0.
    mass = 0.
    for value, E in zip((1, -1, 0), effects_n):
        inserted = np.array([E@rho@E for rho in blocks])
        for _ in range(n):
            a, b = inserted
            inserted = np.array([unitaries[0]@((1-p)*a+b)@unitaries[0].conj().T,
                                 unitaries[1]@(p*a)@unitaries[1].conj().T])
        for next_value, next_E in zip((1, -1, 0), effects_n):
            probability = float(np.trace(next_E@sum(inserted)).real)
            assert probability > -TOLERANCE
            mass += probability
            correlation += value*next_value*probability
    assert abs(mass-1) < TOLERANCE and abs(correlation-math.cos(n*h)) < TOLERANCE
    # A common pulse phase cancels from density histories but changes an
    # added reference-to-response coherence; it is not a mass construction.
    phase = np.exp(-1j*h/3)
    rho = np.array(PB/4, complex)
    U = unitaries[0]
    assert np.max(abs(U@rho@U.conj().T-(phase*U)@rho@(phase*U).conj().T)) < TOLERANCE
    assert np.linalg.norm(U-phase*U) > 0
    return {'bright_correlation': 'cos(n*h)', 'dark_correlation': '0',
            'tracial_correlation': '(2/3)*cos(n*h)', 'retained_joint_outcomes': 9,
            'sample_correlation_residual': abs(correlation-math.cos(n*h)),
            'reference_extension_density_invariant': True,
            'absolute_amplitude_frequency_selected_by_density_instrument': False}


def passive_certificate():
    from model.common_source.full_repair_response import stationary_selection as native
    weights = native.maximal_stationary_weights()
    data, geometry = native.schur_family(weights), native.family_geometry()
    B, G, U = geometry['basis'], geometry['gram'], data['U']
    S, Z, x = data['S0'], data['Z'], geometry['conserved_family_coordinate']
    assert clean(Z-s.eye(3)-(U-B).T*(U-B)/G) == s.zeros(3)
    assert clean(S*x) == s.zeros(3, 1)
    assert data['rank'] == 2
    assert clean(native.port_generator(weights)*geometry['conserved_port_vector']) == s.zeros(12, 1)
    for edge in native.source()['admitted']:
        assert clean(native.port_generator({edge: 1})*geometry['conserved_port_vector']) == s.zeros(12, 1)
    conservation = native.conservation_certificate()
    return {'source_seams': len(weights), 'retained_timeout_probability': str(1-sum(weights.values())),
            'rank': data['rank'], 'positive_metric_identity': 'Z=I+(U-B)^T*(U-B)/G',
            'conserved_coordinate': entries(x), 'S0': entries(S), 'Z': entries(Z),
            'generalized_coefficients': [str(v) for v in data['canonical'].eigenvals()],
            'conservation': conservation, 'positive_normalization_lifts_zero': False}


def algebra_certificate():
    X, Y, Z = s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)
    alpha = [s.kronecker_product(Z, a, s.eye(2)) for a in (X, Y, Z)]
    M, A = s.Matrix([[2, 1], [1, 3]]), s.Matrix([[0, -s.I], [s.I, 0]])
    betaM = s.kronecker_product(X, s.eye(2), M)
    potential = s.kronecker_product(s.eye(4), A)
    assert all(clean(a*betaM+betaM*a) == s.zeros(8) for a in alpha)
    assert clean((potential-alpha[0]*potential*alpha[0])/2) == s.zeros(8)
    h, k = s.Rational(1, 8), s.Rational(1, 10)
    H0, Hk = betaM+h*potential, betaM+h*potential+k*alpha[2]
    witness = s.trace((Hk*Hk-k*k*s.eye(8))**2)-s.trace(H0**4)
    assert s.simplify(witness-16*h*h*k*k*s.trace(A*A)) == 0
    # Cyclic trace proof tested on noncommuting finite matrices.
    H, P, K = s.Matrix([[2, s.I], [-s.I, 1]]), X, s.Matrix([[1, 2], [3, 4]])
    ET, EB = P*H-H*P, K*P-P*K-s.I*H
    assert s.expand(s.trace(K*ET)-s.trace(EB*H)-s.I*s.trace(H*H)) == 0
    # Exact quantum-link occupation counterexample with complete Fock space.
    annihilation = [s.kronecker_product(s.Matrix([[0, 1], [0, 0]]), s.eye(2)),
                    s.kronecker_product(Z, s.Matrix([[0, 1], [0, 0]]))]
    number = s.zeros(8)
    for i in range(2):
        for j in range(2):
            link_entry = s.zeros(2); link_entry[i, j] = s.Rational(1, 2)
            number += s.kronecker_product(annihilation[i].H*annihilation[j], link_entry)
    eigenvalues = number.eigenvals()
    assert eigenvalues == {0: 5, 1: 1, s.Rational(1, 2): 2}
    return {'scalar_projection': '(B-alpha_1*B*alpha_1)/2 for spin-invariant B',
            'potential_projection_zero_exact': True,
            'dispersion_trace_identity': '16*h^2*|k|^2*Tr(A^2)', 'dispersion_control_exact': str(witness),
            'finite_translation_identity': 'i*Tr(H^2)=Tr(K*[P,H])-Tr(([K,P]-i*H)*H)',
            'finite_translation_result': '[P,H]=0, [K,P]=i*H, H=H* imply H=0',
            'quantum_link_naive_number_eigenvalues': {str(k): v for k, v in eigenvalues.items()},
            'quantum_link_projector_implies_integer_occupation': False}


def run():
    result = {'schema': 'oph.particle_finite_response.v1',
              'observation_comparison': observation_certificate(),
              'neutral_response': neutral_certificate(), 'detector_response': detector_certificate(),
              'passive_repair': passive_certificate(), 'finite_algebra': algebra_certificate(),
              'numerical_admission': TOLERANCE,
              'coverage': {
                  'exact_recomputed': ['neutral bright characteristic polynomial and native overlap',
                      'native passive harmonic metric and conserved zero',
                      'Clifford scalar projection and dispersion trace control',
                      'finite translation trace identity', 'quantum-link noninteger occupation counterexample',
                      'complete-outcome capacity comparison certificate'],
                  'proved_with_numerical_controls': ['neutral finite-step hierarchy inequality',
                      'stationary detector correlations', 'reference-phase density equivalence'],
                  'inherited_not_reexecuted_here': ['full 343-site gauge-covariant band and polar transport',
                      'full 60-state coordinate-link completion', 'all-lag memory correction producer'],
                  'not_inferred': ['laboratory particle identities', 'selected vacuum and SI clock',
                      'source-to-observer capacity map', 'equality with measured spectrum']},
              'scope': 'Exact finite identities with declared attachments; no laboratory particle spectrum or continuum theorem.'}
    paths = {HERE/'certificate.py', HERE/'PROOF.md', HERE/'test_certificate.py'}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.is_relative_to(ROOT/'model') and path.suffix == '.py':
                paths.add(path)
    result['source_bindings'] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in sorted(paths)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--output', type=Path, default=HERE/'results.json')
    args = parser.parse_args()
    result = run()
    if args.check:
        if json.loads(args.output.read_text()) != result:
            raise SystemExit('finite-response receipt differs from independent replay')
        print('finite-response receipt verified')
    else:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
        print(f'wrote {args.output.relative_to(ROOT) if args.output.is_relative_to(ROOT) else args.output}')


if __name__ == '__main__':
    main()
