"""Derive the LH invariant from all twelve admitted native currents."""
from functools import lru_cache

import sympy as s

from ..common_source import native_matter_attachment as native
from ..common_source import calibrated_matter as calibrated
from ..common_source.chiral_candidate import action


@lru_cache(maxsize=1)
def dictionary():
    r = native.representation()
    selected = list(r['selected'])
    q = r['operators']['q'].extract(selected, selected)
    indices = [i for i in range(15) if q[i, i] == -3]
    if len(indices) != 2:
        raise ArithmeticError('the native lepton doublet changed')
    variables = s.symbols('t0:4')
    candidate = s.zeros(15, 2)
    for a, i in enumerate(indices):
        for j in range(2):
            candidate[i, j] = variables[2*a+j]
    equations = []
    for g, lifted in zip(r['matter_currents'], r['lifted'], strict=True):
        equations.extend(g.T*candidate+candidate*lifted[3:, 3:])
    matrix, rhs = s.linear_eq_to_matrix(equations, variables)
    kernel = matrix.nullspace()
    if rhs != s.zeros(matrix.rows, 1) or len(kernel) != 1:
        raise ArithmeticError('the native LH invariant is not one dimensional')
    tensor = candidate.subs(dict(zip(variables, kernel[0], strict=True)))
    tensor /= s.sqrt((tensor[:, 1].H*tensor[:, 1])[0])
    tensor = s.ImmutableSparseMatrix(tensor)
    old = action._representation()
    native_subsets = [tuple(j for j in range(5) if mask >> j & 1) for mask in selected]
    permutation = s.zeros(15)
    for j, subset in enumerate(old[0]):
        permutation[j, native_subsets.index(subset)] = 1
    for j, g in enumerate(r['matter_currents']):
        old_g = old[2][r['native']['native_to_fixture_port'][j]]
        if not native.zero(old_g*permutation-permutation*g):
            raise ArithmeticError('native and retained current dictionaries disagree')
    transported = permutation*tensor
    retained = calibrated._LH_tensor()
    row, column = next(iter(retained.todok()))
    scale = s.simplify(transported[row, column]/retained[row, column])
    if scale not in (-1, 1) or not native.zero(transported-scale*retained):
        raise ArithmeticError('the native LH line differs from the retained Weinberg line')
    return {'native': r, 'tensor': tensor, 'permutation': s.ImmutableSparseMatrix(permutation),
            'retained_tensor': retained, 'line_scale': scale,
            'equation_rank': matrix.rank(), 'equation_count': matrix.rows}


def certificate():
    data = dictionary()
    h0, h1 = s.symbols('H0 H1')
    coefficients = s.symbols('k00 k01 k02 k11 k12 k22')
    k00, k01, k02, k11, k12, k22 = coefficients
    kappa = s.Matrix([[k00, k01, k02], [k01, k11, k12], [k02, k12, k22]])
    higgs = s.Matrix([h0, h1])
    ell = data['tensor']*higgs
    square = ell*ell.T
    count = 0
    for g, lifted in zip(data['native']['matter_currents'], data['native']['lifted'], strict=True):
        derivative = data['tensor']*lifted[3:, 3:]*higgs
        if not native.zero(g.T*ell+derivative):
            raise ArithmeticError('the arbitrary-Higgs LH current identity failed')
        if not native.zero(g.T*square+square*g+derivative*ell.T+ell*derivative.T):
            raise ArithmeticError('the quadratic Higgs/current identity failed')
        count += 1
    p = data['permutation']
    full = s.kronecker_product(p, s.eye(3))
    source = s.kronecker_product(square, kappa)
    target_ell = data['retained_tensor']*higgs
    target = s.kronecker_product(target_ell*target_ell.T, kappa)
    if not native.zero(full*source*full.T-target) or source.T != source:
        raise ArithmeticError('the six-complex-parameter K45 tensor transport failed')
    neutral = target.subs({h0: 0, h1: 1})
    charge = action.electric_charge()
    if not native.zero(charge.T*neutral+neutral*charge):
        raise ArithmeticError('the actual Majorana slot is electrically charged')
    return {'admitted_source': 'native_sm_current with its declared five-mode spin lift and nonvacuum even exterior sector',
            'native_current_count': count, 'LH_coefficient_unknowns': 4,
            'LH_equation_count': data['equation_count'], 'LH_equation_rank': data['equation_rank'],
            'LH_invariant_line_dimension': 1,
            'native_LH_tensor': action._serialize(data['tensor']),
            'retained_LH_tensor': action._serialize(data['retained_tensor']),
            'native_to_retained_permutation': action._serialize(p),
            'LH_line_scale': str(data['line_scale']),
            'all_Higgs_components_and_six_complex_family_coefficients_retained': True,
            'symmetric_family_coefficient_complex_dimension': 6,
            'quadratic_current_covariance_exact': True,
            'complete_K45_transport_exact': True, 'neutral_electric_charge_exact': True,
            'internal_dimension': 45, 'independent_right_neutrinos_added': 0,
            'dimension_five_coefficient_selected_by_native_source': False,
            'native_control_word_for_this_interaction_executed': False}
