"""Complete local-even CAR classification, including quartic parity terms."""
from functools import lru_cache

import sympy as s

from ..native_mass_channels.algebra import SIGMA
from ..radial_interacting_core.algebra import source_matrix


def annihilators(modes):
    if type(modes) is not int or not 1 <= modes <= 6:
        raise ValueError('one through six complete CAR modes required')
    result = []
    for mode in range(modes):
        c = s.zeros(2**modes)
        for mask in range(2**modes):
            if mask & (1 << mode):
                c[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
        result.append(s.ImmutableSparseMatrix(c))
    return tuple(result)


def onsite_basis(c, site):
    """I, n-1, three spin currents, pair quadratures, and quartic parity."""
    a, b = c[2*site:2*site+2]
    eye = s.eye(a.rows)
    n = a.H*a+b.H*b
    spin = tuple(sum((p[i, j]*c[2*site+i].H*c[2*site+j]
                      for i in range(2) for j in range(2)), s.zeros(a.rows)) for p in SIGMA)
    pair = a.H*b.H
    parity = (eye-2*a.H*a)*(eye-2*b.H*b)
    return (eye, n-eye, *spin, pair+pair.H, -s.I*(pair-pair.H), parity)


@lru_cache(maxsize=3)
def edge_classification(axis):
    """Solve all 256 exact matrix equations for all sixteen local coefficients."""
    if type(axis) is not int or axis not in range(3):
        raise ValueError('axis zero, one or two required')
    c = annihilators(4)
    basis = onsite_basis(c, 0)+onsite_basis(c, 1)
    forward = sum((s.I*SIGMA[axis][i, j]*c[i].H*c[2+j]/2
                   for i in range(2) for j in range(2)), s.zeros(16))
    h = forward+forward.H
    columns = [s.Matrix(q*h-h*q).reshape(256, 1) for q in basis]
    equations = s.Matrix.hstack(*columns)
    expected = []
    for position in (0, 8):
        v = s.zeros(16, 1); v[position] = 1; expected.append(v)
    for position in range(1, 7):
        v = s.zeros(16, 1); v[position] = 1
        sign = 1 if position == 1 or position == axis+2 else -1
        v[8+position] = sign
        expected.append(v)
    kernel = s.Matrix.hstack(*expected)
    if equations.rank() != 8 or kernel.rank() != 8 or equations*kernel != s.zeros(256, 8):
        raise ArithmeticError('complete even edge-current classification failed')
    return {'axis': axis, 'local_even_dimension_each_site': 8,
            'unknown_coefficients': 16, 'complete_Fock_dimension': 16,
            'equation_count': 256, 'equation_rank': 8, 'kernel_dimension': 8,
            'kernel_columns': [[int(x) for x in kernel.col(j)] for j in range(8)],
            'both_quartic_parity_coefficients_forced_zero': True,
            'number_coefficient_transport': 'a_y=a_x',
            'spin_coefficient_transport': 'A_y=sigma_j A_x sigma_j',
            'pair_coefficient_transport': 'f_y=-f_x',
            'site_constants_unconstrained': True}


def local_completeness():
    c = annihilators(2); basis = onsite_basis(c, 0)
    parity = basis[-1]
    flattened = s.Matrix.hstack(*(s.Matrix(q).reshape(16, 1) for q in basis))
    if flattened.rank() != 8 or any(q != q.H or q*parity != parity*q for q in basis):
        raise ArithmeticError('the full local even Hermitian algebra was truncated')
    plus = (basis[5]+s.I*basis[6])/2
    z = basis[1]/2
    if (plus*plus.H-plus.H*plus != 2*z or z*plus-plus*z != plus
            or plus*parity != parity*plus):
        raise ArithmeticError('onsite eta algebra or centered contact commutation failed')
    for sign in (-1, 1):
        raising = sign*plus
        for up, down in ((c[0].H, sign*c[1]), (c[1].H, -sign*c[0])):
            if (z*up-up*z != up/2 or z*down-down*z != -down/2
                    or raising*down-down*raising != up
                    or raising.H*up-up*raising.H != down
                    or raising*up-up*raising != s.zeros(4)
                    or raising.H*down-down*raising.H != s.zeros(4)):
                raise ArithmeticError('the complete CAR eta-Nambu current action failed')
    return {'all_even_Hermitian_dimension': 8,
            'basis': ['I', 'n-1', 'Sx', 'Sy', 'Sz', 'pair+pair*', '-i(pair-pair*)', '(-1)^n'],
            'quartic_local_term_retained': True,
            'eta_commutators_exact': True, 'centered_contact_commutes_with_all_even_local_currents': True,
            'full_CAR_Nambu_doublets': ['(c_up_dagger,D_x c_down)', '(c_down_dagger,-D_x c_up)'],
            'all_eta_Nambu_ladder_equations_exact': True}


def texture_constraints():
    """Two transported texture values determine all spin coefficients."""
    x, y, z = s.symbols('x y z', real=True)
    spin = x*SIGMA[0]+y*SIGMA[1]+z*SIGMA[2]
    columns = []
    for axis in (0, 2):
        commutator = spin*SIGMA[axis]-SIGMA[axis]*spin
        columns.extend(list(commutator))
    matrix, rhs = s.linear_eq_to_matrix(columns, (x, y, z))
    if matrix.rank() != 3 or rhs != s.zeros(8, 1):
        raise ArithmeticError('two nonparallel texture directions did not remove the spin current')
    return {'transported_nonparallel_directions': ['e1', 'e3'],
            'spin_unknowns': 3, 'constraint_rank': 3,
            'nonconstant_local_conserved_even_current_dimension': 3,
            'Lie_algebra': 'su(2)_eta',
            'number_preserving_subalgebra': 'u(1) generated by N-V',
            'massless_or_collinear_texture_subject_to_separate_larger_commutant': True}


@lru_cache(maxsize=1)
def literal_source_pairing():
    """Check the pairing equation on every row of the original 54-row source."""
    h, points = source_matrix()
    J = s.I*SIGMA[1]
    pair = s.diag(*((-1)**(sum(x) % 2) for x in points))
    pairing = s.kronecker_product(pair, J)
    if h*pairing+pairing*h.T != s.zeros(54):
        raise ArithmeticError('the literal full radial source breaks eta pairing')
    if pairing.H*pairing != s.eye(54):
        raise ArithmeticError('pairing must retain every source mode')
    return {'complete_one_body_dimension': 54,
            'pairing_matrix': 'diag((-1)^sum(x)) tensor i sigma_2',
            'h Delta + Delta h^T_exact_zero': True,
            'all_source_modes_and_edges_retained': True}


def certificate():
    return {'complete_local_algebra': local_completeness(),
            'all_three_edge_classifications': [edge_classification(j) for j in range(3)],
            'texture_commutant': texture_constraints(),
            'literal_radial_pairing': literal_source_pairing()}
