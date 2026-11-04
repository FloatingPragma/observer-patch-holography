"""Exact history symmetries and retained finite-memory nonunitarity."""
from functools import lru_cache

import sympy as s

from . import algebra, intervals as iv

X = s.Matrix([[0, 1], [1, 0]])
Z = s.diag(1, -1)
X6 = s.kronecker_product(X, s.eye(3))
Z6 = s.kronecker_product(Z, s.eye(3))


def _exact(value):
    value = s.sympify(value)
    if value.has(s.Float) or value.free_symbols or value.is_real is not True or value.is_finite is not True:
        raise ValueError('a finite exact real number is required')
    return value


def _parameters(p, q, cosine):
    p, q, cosine = map(_exact, (p, q, cosine))
    if (p.is_nonnegative is not True or (1-p).is_positive is not True
            or q.is_nonnegative is not True or (1-q).is_nonnegative is not True
            or (1-cosine).is_nonnegative is not True or (1+cosine).is_nonnegative is not True):
        raise ValueError('require 0<=p<1, 0<=q<=1 and -1<=cosine<=1')
    return p, q, cosine


@lru_cache(maxsize=1)
def formulas():
    data = algebra.symbolic_memory()
    p, q = data['p'], data['q']
    c = s.Symbol('c', real=True)
    a, b, D = p*(1-q), p*q, data['denominator']
    pi2 = (1-a)/D
    weights = {'always': (0, pi2, a*pi2, a*b/D),
               'C': (1-p, pi2-(1-p), a*pi2, a*b/D),
               'E': (1-q, q-a*b*pi2, 0, a*b*pi2),
               'CE': (1-b, b-a*b*pi2, 0, a*b*pi2)}
    defects = {}
    t0, t1 = (2+c)/3, (7+8*c)/15
    for policy, values in weights.items():
        w0, w2, w3, w10 = map(s.sympify, values)
        actual = data['weights'][policy]
        expected = (1-sum(actual), actual[2], actual[3], actual[9]+actual[10])
        if any(s.cancel(a-b) != 0 for a, b in zip(values, expected, strict=True)):
            raise ArithmeticError('stationary grouping lost a source record')
        active = 1-w0
        defects[policy] = (2*w0*active*(1-t0*t0)
                           +(active*active-w2*w2-w3*w3-w10*w10)*(1-t1*t1))
    return {'p': p, 'q': q, 'c': c, 'weights': weights, 'defects': defects,
            'identity_overlap': t0, 'distinct_projector_overlap': t1}


def stationary_defect(p, q, cosine, policy='always'):
    if policy not in algebra.POLICIES:
        raise ValueError('unknown source duty')
    p, q, cosine = _parameters(p, q, cosine)
    f = formulas()
    return s.cancel(f['defects'][policy].subs({f['p']: p, f['q']: q, f['c']: cosine}))


def finite_history_bound(p, q, cosine, steps, policy='always'):
    """Choi purity lower bound; arbitrary independent initial memory, N>=2."""
    if policy not in algebra.POLICIES or type(steps) is not int or steps < 2:
        raise ValueError('a source duty and integer N>=2 are required')
    p, q, cosine = _parameters(p, q, cosine)
    a, b = p*(1-q), p*q
    if policy == 'always':
        return s.cancel(2*a*(1-p)**(2*steps-1)*(1-((7+8*cosine)/15)**2))
    inactive = (1-p)*(1-q)
    final = {'C': a, 'E': (1-p)*q, 'CE': b}[policy]
    return s.cancel(2*final*inactive**(2*steps-1)*(1-((2+cosine)/3)**2))


def pulse(projector, cosine, sine):
    cosine, sine = _exact(cosine), _exact(sine)
    if s.simplify(cosine*cosine+sine*sine-1) != 0:
        raise ValueError('an exact point of the trigonometric unit circle is required')
    if (not isinstance(projector, s.MatrixBase) or projector.shape != (3, 3)
            or not algebra.zero(projector*projector-projector)
            or not algebra.zero(projector.H-projector)):
        raise ValueError('a Hermitian three-dimensional projector is required')
    return (s.eye(6)+(cosine-1)*s.kronecker_product(s.eye(2), projector)
            -s.I*sine*s.kronecker_product(X, projector))


def hermitian_grades(matrix):
    """Recognize the exact principal-generator symmetry class, not a logarithm."""
    if (not isinstance(matrix, s.MatrixBase) or matrix.shape != (6, 6)
            or not algebra.zero(matrix.H-matrix)
            or not algebra.zero(matrix*X6-X6*matrix)
            or not algebra.zero(Z6*matrix.conjugate()*Z6+matrix)):
        raise ValueError('the Hermitian X/antiunitary grading is required')
    A = algebra.clean((matrix[:3, :3]+matrix[3:, 3:])/(2*s.I))
    S = algebra.clean((matrix[:3, 3:]+matrix[3:, :3])/2)
    if (not algebra.zero(A+A.T) or not algebra.zero(S-S.T)
            or any(s.im(x) != 0 for x in A) or any(s.im(x) != 0 for x in S)):
        raise ArithmeticError('the real generator grades changed')
    return A, S


def coherent_projection(weighted_unitaries):
    """Hilbert--Schmidt projection onto traceless Hamiltonian commutators.

    This dimensionless tangent diagnostic is not an exact channel generator.
    """
    result = s.zeros(6)
    total = s.Integer(0)
    for weight, U in weighted_unitaries:
        weight = _exact(weight)
        if weight.is_nonnegative is not True or U.shape != (6, 6) or not algebra.zero(U.H*U-s.eye(6)):
            raise ValueError('nonnegative weights and exact six-dimensional unitaries required')
        total += weight
        trace = s.trace(U)
        result += s.I*weight*(s.conjugate(trace)*U-trace*U.H)/12
    if s.simplify(total-1) != 0:
        raise ValueError('the complete history weights must sum to one')
    return algebra.clean(result)


@lru_cache(maxsize=1)
def exact_controls():
    projectors = algebra.symbolic_memory()['projectors']
    cosine, sine = s.Rational(99, 101), s.Rational(20, 101)
    pulses = [pulse(P, cosine, sine) for P in projectors]
    for U in pulses:
        if (not algebra.zero(U.H*U-s.eye(6)) or not algebra.zero(U*X6-X6*U)
                or not algebra.zero(Z6*U.conjugate()*Z6-U)):
            raise ArithmeticError('an actual pulse lost its complete history symmetry')
    for i in (2, 3, 10):
        for j in (2, 3, 10):
            if i != j and s.simplify(s.trace(projectors[i]*projectors[j])-s.Rational(1, 5)) != 0:
                raise ArithmeticError('source projector overlap changed')
    a, b = s.Rational(5, 72), s.Rational(1, 72)
    histories = [((2, 2), (1-a)**2), ((2, 3), a*(1-a)),
                 ((3, 2), a*(1-b)), ((3, 10), a*b)]
    H = coherent_projection([(w, pulses[path[1]]*pulses[path[0]]) for path, w in histories])
    A, S = hermitian_grades(H)
    td = (10+16*cosine+4*cosine*cosine)/5
    commutator = lambda x, y: x*y-y*x
    expected = td*cosine*(cosine-1)/6*(a*(b-a)*commutator(projectors[3], projectors[2])
                                       +a*b*commutator(projectors[10], projectors[3]))
    if not algebra.zero(A-expected) or A[0, 2].is_positive is not True:
        raise ArithmeticError('the retained two-step diagonal correction changed')
    first = s.diag(1, 2, 3)
    second = s.Matrix([[1, 1, s.I], [1, 2, 1], [-s.I, 1, 3]])
    for K in (first, second):
        imaginary = s.im(K)
        real = s.re(K)
        hermitian_grades(s.kronecker_product(s.eye(2), s.I*imaginary)+s.kronecker_product(X, real))
    invariant = algebra.cp_invariant(first, second)
    if invariant != -34560*s.I:
        raise ArithmeticError('the symmetry-class CP countercontrol changed')
    skew = s.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
    alias = s.kronecker_product(X, s.I*skew)
    lam = s.Symbol('lambda')
    if s.factor(alias.charpoly(lam).as_expr()-lam**2*(lam*lam-1)**2) != 0:
        raise ArithmeticError('the nonprincipal identity alias changed')
    return {'source_projectors_checked': 12,
            'two_step_diagonal_A02_exact': str(A[0, 2]),
            'two_step_point': {'p': '1/12', 'q': '1/6', 'cosine': '99/101', 'sine': '20/101'},
            'two_step_all_four_histories_retained': True,
            'two_step_offdiagonal_S_real_symmetric': True,
            'general_Hermitian_plus_block_CP_counterexample': str(invariant),
            'counterexample_is_not_a_physical_Yukawa_assignment': True,
            'identity_alias_spectrum_in_units_2pi': [-1, -1, 0, 0, 1, 1],
            'nonprincipal_alias_has_imaginary_mass_block': True}


def memory_log_certificate():
    data = algebra.symbolic_memory()
    p, q = data['p'], data['q']
    a, b = p*(1-q), p*q
    T = data['transition']
    if any(T[i, j] != 0 for i in (0, 1) for j in range(12)):
        raise ArithmeticError('a forbidden core memory output became nonzero')
    indices = [2, 3, 9, 10]
    recurrent = T.extract(indices, indices)
    expected = s.Matrix([[1-a, 1-b, 1-a, 1-a], [a, 0, 0, 0],
                         [0, 0, 0, a], [0, b, a, 0]])
    lam = s.Symbol('lambda')
    polynomial = (lam-1)*(lam+a)*(lam*lam-a*a+a*b)
    if not algebra.zero(recurrent-expected) or s.expand(recurrent.charpoly(lam).as_expr()-polynomial) != 0:
        raise ArithmeticError('the recurrent memory spectrum changed')
    return {'complete_CQ_operator_dimension': 432, 'rank_upper_bound': 360,
            'forbidden_core_memory_output_blocks': [0, 1],
            'complete_CQ_channel_has_any_matrix_logarithm': False,
            'recurrent_ports': indices, 'recurrent_characteristic_polynomial': str(polynomial),
            'recurrent_identity_sector_eigenvalue': '-p*(1-q)',
            'recurrent_principal_log_off_negative_real_axis_exists_at_positive_gates': False,
            'other_logarithm_branches_or_embeddability_classified': False}


def certificate():
    f = formulas()
    rows = []
    cosine = 1-2*iv.sine(iv.I(s.Rational(1, 8)))**2
    for name, gate in algebra.source_gates().items():
        substitutions = {f['p']: iv.expression(gate['p']), f['q']: iv.expression(gate['q']), f['c']: cosine}
        for policy in algebra.POLICIES:
            defect = iv.expression(f['defects'][policy], substitutions)
            if defect.lo <= 0:
                raise ArithmeticError('a positive-gate stationary Choi witness lost positivity')
            rows.append({'channel': name, 'duty': policy,
                         'normalized_Choi_purity_defect_interval': iv.rounded_out(defect, 12).strings()})
    return {'history_symmetry': '[U,X]=0; (Z tensor I) conjugate(U) (Z tensor I)=U',
            'principal_generator_domain': 'unitary history with no eigenvalue -1; sum(abs(pulse steps))<2 is sufficient',
            'principal_generator_form': 'I tensor iA + X tensor S; A real antisymmetric, S real symmetric',
            'continuous_unitary_limit_has_same_form_modulo_scalar': True,
            'symmetry_is_physical_CP_identification': False,
            'whole_fine_current_output_antiunitary_claimed': False,
            'stationary_product_memory_weights': {k: list(map(str, v)) for k, v in f['weights'].items()},
            'normalized_trace_overlaps': {'idle_active': str(f['identity_overlap']),
                                           'distinct_active': str(f['distinct_projector_overlap'])},
            'stationary_Choi_defect': {k: str(v) for k, v in f['defects'].items()},
            'stationary_nonunitarity_domain': '0<p,q<1; epsilon not in 2pi Z; all four duties',
            'diamond_norm_convention': 'unhalved',
            'distance_from_every_unitary_lower_bound': '2*(1-sqrt(1-delta)) >= delta',
            'finite_history_domain': 'N>=2, any initial memory independent of feedback/reference, 0<p,q<1, epsilon not in 2pi Z',
            'finite_history_bounds': {'always': '2*p*(1-q)*(1-p)^(2N-1)*(1-((7+8*cos(epsilon))/15)^2)',
                                     'C_E_CE': '2*d_D*((1-p)*(1-q))^(2N-1)*(1-((2+cos(epsilon))/3)^2)',
                                     'd_C_d_E_d_CE': ['p*(1-q)', '(1-p)*q', 'p*q']},
            'finite_bounds_give_positive_asymptotic_noise_floor': False,
            'pulse_angle_for_stationary_intervals': '1/4',
            'stationary_interval_witnesses': rows,
            'channel_logarithms': memory_log_certificate(),
            'exact_controls': exact_controls(),
            'coherent_projection_is_an_exact_channel_generator': False,
            'full_correlated_channel_identified_with_static_mean_action': False,
            'physical_particle_or_interacting_continuum_claim': False}
