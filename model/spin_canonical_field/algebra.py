"""Exact spin identities; finite controls do not replace the limiting proof."""
import sympy as s


def spin():
    return (s.Matrix([[0, 1], [1, 0]])/2,
            s.Matrix([[0, -s.I], [s.I, 0]])/2,
            s.diag(1, -1)/2)


def comm(a, b):
    return a*b-b*a


def exact_checks():
    v = spin()
    identity = s.eye(2)
    h = sum((s.kronecker_product(a, a) for a in v), s.zeros(4))-s.eye(4)/4
    x, y = s.symbols('x y', real=True)
    for axis in range(3):
        left = s.kronecker_product(v[axis], identity)
        right = s.kronecker_product(identity, v[axis])
        a = x*left+y*right
        transverse = sum((s.kronecker_product(v[b], v[b])
                          for b in range(3) if b != axis), s.zeros(4))
        residual = comm(a, comm(h, a))+(x-y)**2*transverse
        if residual.applyfunc(s.simplify) != s.zeros(4):
            raise ArithmeticError('double commutator identity')
        if max(abs(e) for e in transverse.eigenvals()) != s.Rational(1, 2):
            raise ArithmeticError('transverse norm')
    # Q=p^2 eps S3 and P=-p S2 have [Q,P]=i p^3 eps S1.
    if comm(v[2], -v[1]) != s.I*v[0]:
        raise ArithmeticError('canonical orientation')
    # All three staggered macroscopic coordinates commute in the limit;
    # at finite scale the transverse momenta are NOT a second free pair.
    if comm(v[1], v[2]) != s.I*v[0]:
        raise ArithmeticError('uncontrolled second-pair commutator')
    return {'all_three_axis_double_commutators': True,
            'transverse_bond_norm': '1/2', 'canonical_sign': '[Q,P]=+iM',
            'second_transverse_momentum_commutator_is_nonzero': True,
            'finite_controls_are_ground_scaling_evidence': False}
