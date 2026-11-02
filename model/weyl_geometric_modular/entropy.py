"""Finite, dimensionless bounds for the vacuum-relative first law."""
from fractions import Fraction as F
import sympy as s


def relative_entropy_bound(epsilon, radius):
    if type(epsilon) not in (int, F) or type(radius) not in (int, F):
        raise ValueError('exact rational perturbation and radius required')
    epsilon, radius = F(epsilon), F(radius)
    if not 0 <= abs(epsilon) <= radius < 1:
        raise ValueError('require |epsilon|<=radius<1')
    return epsilon**2/(2*(1-radius**2))


def certificate():
    e = s.symbols('epsilon', real=True)
    relative = ((1+e)*s.log(1+e)+(1-e)*s.log(1-e))/2
    if s.simplify(s.diff(relative, e, 2)-1/(1-e*e)) != 0:
        raise ArithmeticError('binary relative entropy Hessian changed')
    if relative.subs(e, 0) != 0 or s.diff(relative, e).subs(e, 0) != 0:
        raise ArithmeticError('relative entropy has an erroneous first-order term')
    return {'perturbations': 'binary normal-state mixtures with faithful vacuum as midpoint',
            'quantum_bound_input': 'Araki relative entropy decreases under normal state preparation channels',
            'binary_relative_entropy_second_derivative': '1/(1-epsilon^2)',
            'relative_entropy_bound': 'epsilon^2/[2*(1-radius^2)]',
            'domain': '|epsilon|<=radius<1',
            'example_epsilon_1_10_radius_1_2': str(relative_entropy_bound(F(1, 10), F(1, 2))),
            'relative_entropy_first_derivative_at_vacuum': '0',
            'vacuum_subtracted_entropy_first_law': 'd(DeltaS)/d(epsilon)=d(DeltaK_B)/d(epsilon) at zero',
            'absolute_type_III_von_Neumann_entropy_assumed_finite': False,
            'relative_entropy_definition_identified_with_a_specific_lattice_subtraction': False,
            'generalized_entropy_stationarity_or_area_law_inferred': False}
