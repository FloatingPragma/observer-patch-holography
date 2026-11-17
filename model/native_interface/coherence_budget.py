"""Sharp preparation floors for the committed native-only observer grammar.

The C/E Kraus maps are checked partial injections, not assumed incoherent.
Classical adaptive control, fresh independent diagonal ancillas and heralding
preserve diagonality of the joint accessible work marginal. Diagonal marginals
of correlated accessible registers alone are insufficient. Coherent access to
an uncommitted record dilation, or a coherent work encoding, is outside this
grammar.

All reported bounds are rational. Numerical source guards enter only as
explicit conditional allowances; this producer does not reexecute their
source words or turn those guards into physical fidelity certificates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import sympy as sp

from tools import strict_json
from . import construction as native_interface
from ..common_source.slow_geometry import source
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..predictions import inverse_clock as preparation
from ..predictions import prepared_clock, prepared_clock_source
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*native_interface.SOURCE_PATHS,
    *preparation.SOURCE_PATHS, 'model/native_interface/coherence_budget.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/predictions/prepared_clock.py',
    'model/predictions/prepared_clock_source.py'))))


def _fraction(value, *, maximum=None):
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or math.gcd(value.numerator, value.denominator) != 1
            or value < 0 or maximum is not None and value > maximum):
        raise ValueError('an original canonical nonnegative Fraction in range is required')
    return value


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original integer in the declared finite range is required')
    return value


def herald_floor(probability, conditional_floor, conditional_target_error=F(0)):
    """A complete-output floor, including any change of success probability.

    For a target success block p*rho with distance c from normalized diagonal
    states, D(p*rho,q*delta)+abs(p-q)/2 >= p*c. The rest of the normalized
    output supplies the second term. An epsilon approximation to rho reduces
    c by at most epsilon. This is not an error bound conditioned on a rare event.
    """
    p = _fraction(probability, maximum=F(1))
    c = _fraction(conditional_floor, maximum=F(1))
    epsilon = _fraction(conditional_target_error, maximum=F(1))
    return p*max(F(0), c-epsilon)


def robust_floor(ideal_complete_floor, *, target_error=F(0), input_error=F(0),
                 per_step_error=F(0), attempt_cap=0):
    """Triangle bound for a finite approximate native-only protocol.

    Each per-step number is a uniform half-diamond error of a complete CPTP
    instrument, including its classical output. The cap is a maximum over
    adaptive paths, not an average stopping time. Initial error and every
    service outside the ideal grammar must be charged. No physical value for
    these errors is inferred here.
    """
    floor = _fraction(ideal_complete_floor, maximum=F(1))
    target = _fraction(target_error, maximum=F(1))
    initial = _fraction(input_error, maximum=F(1))
    step = _fraction(per_step_error, maximum=F(1))
    cap = _integer(attempt_cap, 0, 10**9)
    native_error = min(F(1), initial+cap*step)
    return max(F(0), floor-target-native_error)


def sqrt5_interval(digits=9):
    """Strict rational enclosure certified by integer squares."""
    _integer(digits, 1, 30)
    grid = 10**digits
    numerator = math.isqrt(5*grid*grid)
    lower, upper = F(numerator, grid), F(numerator+1, grid)
    if not lower*lower < 5 < upper*upper:
        raise ArithmeticError('the rational square-root enclosure failed')
    return lower, upper


@lru_cache(maxsize=1)
def _native_geometry():
    seams = native_interface.quantum.mixing.seams()
    laplacian = sp.zeros(12)
    for a, b in seams:
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    variable = sp.Symbol('x')
    expected = variable*(variable-6)**5*(variable**2-10*variable+20)**3
    if laplacian.charpoly(variable).as_expr().expand() != expected.expand():
        raise ArithmeticError('the native Laplacian spectral multiplicities changed')
    rotations = tuple(tuple(permutation) for permutation in source().rotations())
    if len(rotations) != 60 or len(set(rotations)) != 60:
        raise ArithmeticError('the supplied native rotations are not sixty distinct permutations')
    for permutation in rotations:
        if sorted(permutation) != list(range(12)) or any(
                laplacian[permutation[a], permutation[b]] != laplacian[a, b]
                for a in range(12) for b in range(12)):
            raise ArithmeticError('a native rotation fails exact Laplacian invariance')
    # Uniform averaging must send every diagonal density to I/12. Checking
    # all orbit counts proves precisely the twirling property used below.
    if any(sum(permutation[a] == b for permutation in rotations) != 5
           for a in range(12) for b in range(12)):
        raise ArithmeticError('the native rotations do not give uniform transitive twirling')
    return sp.ImmutableMatrix(laplacian), rotations


def prepared_reference_floor(digits=9):
    """Sharp diagonal-state floor for the actual capped sampler polynomial.

    Vertex transitivity fixes rho and twirls every diagonal density to I/12.
    Contractivity therefore proves that I/12 is a closest diagonal state.
    Exact spectral multiplicities and rational mass enclosures identify the
    positive eigenvalues of rho-I/12 without floating eigenvalue decisions.
    """
    lower, upper = sqrt5_interval(digits)
    laplacian, rotations = _native_geometry()
    weights = preparation.subordinated_filter_weights()
    if type(weights) is not tuple or len(weights) != 5:
        raise ArithmeticError('the declared capped positive sampler changed')
    for weight in weights:
        _fraction(weight, maximum=F(1))
    if not 0 < weights[0] <= sum(weights) <= 1 or not any(weights[1:]):
        raise ArithmeticError('the declared capped positive sampler changed')
    lazy = sp.eye(12)-laplacian/120
    polynomial = sum((sp.Rational(weight)*lazy**(2*n)
                      for n, weight in enumerate(weights)), sp.zeros(12))
    success = F(sp.trace(polynomial))/12
    if not 0 < success < 1:
        raise ArithmeticError('the successful preparation mass is outside its domain')

    def mass(eigenvalue):
        return sum((weight*(1-eigenvalue/120)**(2*n)
                    for n, weight in enumerate(weights)), F(0))

    m0, m6 = mass(F(0)), mass(F(6))
    low = (mass(5-lower), mass(5-upper))
    high = (mass(5+upper), mass(5+lower))
    # The positive-weight polynomial decreases on 0 <= lambda < 120.
    if not (0 < high[0] <= high[1] < success and m6 < success
            < low[0] <= low[1] < m0):
        raise ArithmeticError('the rational enclosures do not resolve all spectral signs')
    complete = tuple((m0+3*m-4*success)/12 for m in low)
    conditional = tuple(value/success for value in complete)
    if not 0 < conditional[0] <= conditional[1] < 1:
        raise ArithmeticError('the prepared-reference distance enclosure failed')
    return {'weights': weights, 'sqrt5_interval': (lower, upper),
        'spectral_mass_intervals': {'0': (m0, m0), '5-sqrt5': low,
                                   '6': (m6, m6), '5+sqrt5': high},
        'success_probability': success,
        'conditional_distance_to_diagonal_interval': conditional,
        'complete_success_failure_distance_floor_interval': complete,
        'positive_bands': ('0', '5-sqrt5'), 'positive_multiplicity': 4,
        'negative_multiplicity': 8,
        'closest_diagonal_state': 'I12/12',
        'complete_floor_formula': '(m0+3*m_(5-sqrt5)-4*Psuccess)/12',
        'exact_rotation_count': len(rotations),
        'exact_transitive_orbit_count_per_address': len(rotations)//12,
        'sharpness_scope': 'closest diagonal successful-work state and coarse success/failure '
                           'output with failed quantum systems discarded; no sharpness claim '
                           'for the actual failed-work or retained reference outputs'}


def native_phase_floor():
    """Existing native sign-herald targets, with an exact attainable floor.

    Each target has ten disjoint balanced coherent pairs. Its nearest
    diagonal state is the equal mixture of both orientations. One further
    unread fixed-seam repair prepares that mixture natively after either
    sign herald. The comparison retains the original sign and failure; the
    additional repair's record is explicitly forgotten by this access map.
    """
    audits = {law: native_interface.exact_native_audit(law) for law in ('C', 'E')}
    pairs = native_interface.phase_pairs()
    dimension = len(native_interface.quantum.STATES)
    if (not pairs or len({index for pair in pairs for index in pair}) != 2*len(pairs)
            or any(not row['native_outcomes_commute_with_work_dephasing']
                   or not row['native_classical_feedback_cannot_create_coherence_from_diagonal_inputs']
                   for row in audits.values())):
        raise ArithmeticError('the actual native partial-injection premise failed')
    per_sign = F(len(pairs), 2*dimension)
    for law, audit in audits.items():
        if F(audit['exact_herald_probability_per_sign']) != per_sign:
            raise ArithmeticError('the source sign-herald probability changed')
        # Check the saturating native preparation by exact population counts.
        branches = native_interface.fixed_maps(law)
        for sign in (-1, 1):
            first = [F(0)]*dimension
            for branch in branches:
                if branch.current == sign:
                    for _, output in branch.mapping:
                        first[output] += F(1, 2*dimension)
            second = [F(0)]*dimension
            for branch in branches:
                for initial, output in branch.mapping:
                    second[output] += first[initial]/2
            support = {index for pair in pairs for index in pair}
            if any(value != (per_sign/(2*len(pairs)) if i in support else 0)
                   for i, value in enumerate(second)):
                raise ArithmeticError('one unread native repair does not saturate the diagonal floor')
    conditional = len(pairs)*F(1, 2*len(pairs))
    complete = 2*per_sign*conditional
    return {'carrier_dimension': dimension, 'orthogonal_orientation_pairs': len(pairs),
        'per_sign_success_probability': per_sign,
        'total_success_probability': 2*per_sign,
        'failure_probability': 1-2*per_sign,
        'conditional_distance_to_diagonal': conditional,
        'complete_herald_work_distance_floor': complete,
        'sharp_saturation_additional_unread_native_repairs_on_success': 1,
        'saturation_extra_record_forgotten_by_declared_access': True,
        'fine_record_or_reference_sharpness_claimed': False,
        'native_partial_injection_audits': {law: {
            'recorded_maps': row['recorded_partial_injection_maps_checked'],
            'source_addresses': row['native_branch_addresses_checked'],
            'dephasing_intertwining': row['native_outcomes_commute_with_work_dephasing']}
            for law, row in audits.items()}}


def prepared_source_allowance(weights, success):
    """Reproduce the existing guard arithmetic without executing its words.

    Its full allowance is conservatively larger than the stopped-preparation
    allowance, so it also bounds the raw preparation checkpoint under the
    same guards. The reference floor is not a claim about the subsequently
    measured detector output, whose work state has been dephased.
    """
    if type(weights) is not tuple or weights != preparation.subordinated_filter_weights():
        raise ValueError('the original capped preparation weights are required')
    for weight in weights:
        _fraction(weight, maximum=F(1))
    p = _fraction(success, maximum=F(1))
    if p != prepared_reference_floor()['success_probability']:
        raise ValueError('the exact successful polynomial mass is required')
    basis = _fraction(prepared_clock_source.BASIS_GUARD, maximum=F(1))
    phase = _fraction(prepared_clock_source.PHASE_GUARD, maximum=F(1))
    pair = _fraction(prepared_clock_source.PAIR_GUARD, maximum=F(1))
    filtering = _fraction(prepared_clock.FILTER_GUARD, maximum=F(1))
    arithmetic = _fraction(prepared_clock.ARITHMETIC_ALLOWANCE, maximum=F(1))
    clock_gain = (1+basis)**2*(1+phase)
    kick_gain = (1+basis)**2*(1+pair)
    detector_gain = (1+basis)*(1+pair)
    continuation_gain = (clock_gain**2*kick_gain*detector_gain)**2
    prep_error = sum((weight*((1+filtering)**(2*n)-1)/2
                      for n, weight in enumerate(weights)), F(0))
    exact = continuation_gain*prep_error+p*(continuation_gain-1)/2+arithmetic
    grid = 10**12
    return F(math.ceil(exact*grid), grid)


def _json(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: _json(item) for key, item in value.items()}
    if type(value) in (tuple, list):
        return [_json(item) for item in value]
    return value


def run():
    native = native_phase_floor()
    prepared = prepared_reference_floor()
    floor = prepared['complete_success_failure_distance_floor_interval'][0]
    allowance = prepared_source_allowance(prepared['weights'], prepared['success_probability'])
    after_allowance = robust_floor(floor, target_error=allowance)
    if after_allowance <= 0:
        raise ArithmeticError('the existing target allowance fails to resolve the native-only floor')
    return _json({'schema': 'oph.native_interface.coherence_budget.v1',
        'source_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in SOURCE_PATHS}),
        'native_sign_herald_targets': native,
        'prepared_reference_target': prepared,
        'existing_prepared_source_comparison': {
            'comparison_stage': 'preparation checkpoint, before clock and detector',
            'complete_target_error_allowance_under_existing_numerical_guards': allowance,
            'distance_from_native_only_outputs_lower_under_those_guards': after_allowance,
            'source_word_numerical_guards_reexecuted_here': False,
            'physical_fidelity_or_directed_rounding_certificate': False},
        'robustness': {
            'formula': 'max(0, ideal_complete_floor-target_error-input_error-N*step_error)',
            'step_error': 'uniform half-diamond error of each complete CPTP instrument',
            'N': 'maximum finite adaptive attempt count, not mean stopping time',
            'rare_herald_formula': 'Psuccess*max(0, conditional_floor-conditional_target_error)',
            'changed_success_probability_and_all_failure_mass_included': True},
        'scope': {
            'native_C_E_partial_injections_checked_from_actual_maps': True,
            'classical_adaptation_and_heralding_included': True,
            'joint_input_marginal_of_all_accessible_quantum_registers_diagonal_required': True,
            'fresh_ancillas_independent_and_diagonal_required': True,
            'diagonal_marginals_of_correlated_accessible_registers_suffice': False,
            'work_basis_and_any_extraction_must_preserve_its_dephasing': True,
            'uncommitted_coherent_record_dilations_are_incoherent': False,
            'coherent_encoders_or_source_pulses_inside_native_only_grammar': False,
            'preparation_floor_applied_to_terminal_measured_detector_output': False,
            'all_OPH_observer_interfaces_ruled_out': False,
            'new_physical_data_or_physical_service_errors_inferred': False}})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args(argv)
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native coherence-budget receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
