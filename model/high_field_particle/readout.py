"""Complete local preparation and exact speed bounds on the declared spin source."""
from fractions import Fraction as F

import sympy as s

from model.native_ordered_ground import source
from model.spin_twist_closure import field_gap


def creation_certificate():
    """Retain both pointer results; distinguish spin charge from CAR parity."""
    identity = s.eye(2)
    x = s.Matrix([[0, 1], [1, 0]])
    z = s.diag(1, -1)
    h = s.Matrix([[1, 1], [1, -1]]) / s.sqrt(2)
    p_plus, p_minus = (identity + x)/2, (identity - x)/2
    success, failure = z*p_plus, p_minus
    if success.T*success + failure.T*failure != identity:
        raise ArithmeticError('complete local preparation lost a pointer outcome')
    lowering, occupation = s.Matrix([[0, 0], [1, 0]]), s.diag(0, 1)
    if s.simplify(h*success*h) != lowering or s.simplify(h*failure*h) != occupation:
        raise ArithmeticError('source pointer and feedback have a different frame image')
    charges = []
    for sign in (1, -1):
        # Relative to the product vacuum, physical transformed S1 is -sign*n.
        q = -sign*occupation
        if q*lowering - lowering*q != -sign*lowering:
            raise ArithmeticError('successful local preparation has the wrong charge')
        if q*occupation != occupation*q:
            raise ArithmeticError('failure changed charge')
        charges.append(-sign)
    # Full onsite Fock order: empty, up, down, double. Failure also retains
    # both states outside the single-occupancy code.
    full_success, full_failure = s.zeros(4), s.diag(1, 0, 0, 1)
    full_success[1:3, 1:3] = success
    full_failure[1:3, 1:3] = failure
    grading = s.diag(1, -1, -1, 1)
    if full_success.T*full_success + full_failure.T*full_failure != s.eye(4):
        raise ArithmeticError('onsite Fock failure completion is incomplete')
    if grading*full_success != full_success*grading:
        raise ArithmeticError('spin creator changed fermion grading')
    left = s.kronecker_product(full_success, s.eye(4))
    right = s.kronecker_product(s.eye(4), full_success)
    if left*right != right*left:
        raise ArithmeticError('disjoint even creators failed to commute')
    # The complete instrument is an isometry, including arbitrary untouched
    # references; this checks all source columns rather than one input state.
    isometry = full_success.col_join(full_failure)
    if isometry.T*isometry != s.eye(4):
        raise ArithmeticError('complete preparation is not trace preserving')
    return {
        'source_frame': 'Marshall coordinates of native_ordered_ground.source.matrix',
        'source_success_Kraus': 'Z (I+X)/2',
        'source_failure_Kraus_on_spin_code': '(I-X)/2',
        'implementation_on_protected_spin_code': 'complete X pointer, followed by logical Z on its + outcome',
        'abstract_full_Fock_completion_equals_literal_native_extension': False,
        'source_physical_frame_on_odd_sites': 'relabel the X pointer using the Marshall onsite Z frame',
        'transformed_success_Kraus': '|1><0|',
        'transformed_failure_Kraus': '|1><1|',
        'relative_charge_convention': '-sum_x epsilon_x n_x',
        'successful_charge_by_sublattice': charges,
        'full_onsite_Fock_dimension': 4,
        'complete_isometry_rows': 8,
        'full_Fock_completeness_and_grading_exact': True,
        'disjoint_creators_commute_exactly': True,
        'untouched_reference_preserved_by_complete_instrument': True,
        'failure_and_out_of_code_records_retained': True,
        'success_renormalized': False,
        'finite_coupling_success_is_a_pure_particle_band_projection': False,
        'uniform_numeric_particle_residue_floor_supplied': False,
        'single_physical_fermion_under_grading_preserving_dictionary': False,
        'native_compiler_execution_performed_in_this_certificate': False,
    }


def speed_bound(spacing, ratio_ceiling, rest_energy_ceiling):
    """Exact consequence of gap>=B/8 and source Lieb--Robinson speed<=72J.

    Inputs are positive rationals. The formula uses the SAME ground-relative
    Hamiltonian and a common positive clock factor; it is not a speed fit.
    """
    values = (spacing, ratio_ceiling, rest_energy_ceiling)
    if any(type(v) is not F or v <= 0 for v in values):
        raise ValueError('positive exact spacing, coupling-ratio and rest-energy ceilings required')
    return F(576)*spacing*ratio_ceiling*rest_energy_ceiling


def finite_preparation_example():
    """Reuse the existing complete finite-ground planner, without executing its tape."""
    edges = source.cube()[0]
    return {'coupling': '1', 'selection_field': '100',
            'requested_complete_error': '1/100',
            'plan': field_gap.resource_plan(8, edges, F(1), F(100), F(1, 100)),
            'finite_ground_bound_applies': True,
            'value_1_over_100_certified_below_Yarotsky_threshold': False,
            'native_tape_executed': False,
            'infinite_volume_preparation_error_supplied': False}
