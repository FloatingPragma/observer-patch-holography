"""Infer a configured coherent controller from its finite source spectrum.

This is parameter observability on a specified source and clock. It does
not select a preferred controller or erase target-data calibration ancestry.
Only small exact matrices and the inherited spectral inverse are used.
"""
import json

import sympy as s

from model.particle_candidate_atlas import attenuation
from model.particle_calibration_selection import phase_attenuation as moments
from . import bare_effective
from . import charged_higgs_records as charged


SOURCE_INPUTS = (
    'model/particle_candidate_atlas/attenuation.py',
    'model/particle_calibration_selection/phase_attenuation.py',
    'model/standard_model_source_map/koide_source.py',
    'model/particle_observer_outputs/bare_effective.py',
    'model/particle_observer_outputs/BARE_EFFECTIVE_BOUNDARY.md',
    'model/particle_observer_outputs/charged_higgs_records.py',
    'Lean/ObserverPatchHolography/KoideCirculant.lean',
)


def _clean(value):
    return s.simplify(s.sqrtdenest(s.factor(s.radsimp(s.expand(value)))))


def _real(value):
    return bare_effective._real(value)


def _sign(value):
    return bare_effective._sign(value)


def inverse_moments(trace_y, trace_y_squared):
    """Exact inverse, in the literal source coefficient normalization.

    This checks the moment-domain constraints. Full-spectrum positive-root
    admission belongs to controller_from_coefficients below.
    """
    first, second = map(_real, (trace_y, trace_y_squared))
    if _sign(first) <= 0 or _sign(second) <= 0:
        raise ValueError('strictly positive source spectral moments are required')
    r2 = _clean(36*first/5-1)
    if _sign(r2) < 0 or _sign(1-r2) < 0:
        raise ValueError('the reconstructed attenuation must satisfy 0<=r^2<=1')
    r = _clean(s.sqrt(r2))
    if r2 == 0:
        if _clean(second-s.Rational(25, 3888)) != 0:
            raise ValueError('zero attenuation requires its phase-independent second moment')
        return {'r_squared': r2, 'r': r,
                'cos_three_theta': None, 'principal_abs_theta': None,
                'phase_observable': False,
                'phase_equivalence': 'every unit phase produces the same source coefficient at r=0'}
    residual = second-s.Rational(25, 3888)-s.Rational(25, 648)*r2-s.Rational(25, 2592)*r2*r2
    cosine = _clean(1944*residual/(25*s.sqrt(2)*r**3))
    if _sign(1-cosine) < 0 or _sign(1+cosine) < 0:
        raise ValueError('the reconstructed cos(3theta) must lie in [-1,1]')
    return {'r_squared': r2, 'r': r, 'cos_three_theta': cosine,
            'principal_abs_theta': _clean(s.acos(cosine)/3),
            'phase_observable': True,
            'phase_equivalence': 'theta=+/-phi+2*pi*k/3, k=0,1,2; phi=acos(cos(3theta))/3 in[0,pi/3]'}


def controller_from_coefficients(coefficients, coefficient_scale=1):
    """Consume all three positive sigma(Y), in one known source-clock scale.

    The scale converts recovered energy coefficients to the source's fixed
    C and Y convention. It is not a fitted particle-mass normalization.
    """
    values = tuple(map(_real, coefficients))
    scale = _real(coefficient_scale)
    if len(values) != 3 or _sign(scale) <= 0 or any(_sign(value) <= 0 for value in values):
        raise ValueError('all three positive coefficients and a positive declared scale are required')
    values = tuple(_clean(value/scale) for value in values)
    first, second = _clean(sum(values)), _clean(sum(value*value for value in values))
    result = inverse_moments(first, second)
    a, b = moments.source_coefficients()
    roots = tuple(_clean(s.sqrt(3*value)) for value in values)
    if _clean(sum(roots)-3*a) != 0:
        raise ValueError('the complete spectrum is outside the admitted positive-C-root source branch')
    if result['r_squared'] == 0 and any(_clean(value-a*a/3) != 0 for value in values):
        raise ValueError('zero attenuation must have the complete triply degenerate spectrum')
    return dict(result, sigma_y=values, trace_Y=first, trace_Y_squared=second,
                source_a=a, source_b=b, positive_root_trace=_clean(sum(roots)),
                positive_root_branch_verified=True,
                coefficient_scale=scale,
                common_source_clock_and_normalization_required=True,
                unordered_family_permutations_unresolved=True,
                controller_uniquely_selected_by_source_dynamics=False,
                target_calibration_ancestry_erased=False,
                physical_particle_mass_claimed=False)


def controller_from_reconstructed(rows, channel='left_doublet', coefficient_scale=1):
    """Join the owning charged spectral inverse without losing degeneracies."""
    copies = {'left_doublet': 4, 'charged_singlet': 2}.get(channel)
    if copies is None:
        raise ValueError('an actual source charged channel is required')
    values = []
    for row in rows:
        if row['internal_channel'] != channel:
            continue
        rank = row['multiplicity']
        if type(rank) is not int or rank <= 0 or rank % copies:
            raise ValueError('the full source channel multiplicity must be retained')
        if row['physical_rest_mass_claimed']:
            raise ValueError('this inverse consumes finite source coefficients, not physical mass assignments')
        probabilities = row['original_addition_probabilities']
        if len(probabilities) != 2 or any(_sign(p) <= 0 for p in probabilities):
            raise ValueError('both original addition-line probabilities must be strictly positive')
        values.extend([row['reconstructed_sigma']]*(rank//copies))
    result = controller_from_coefficients(values, coefficient_scale)
    result.update({'internal_channel': channel, 'full_family_multiplicity': len(values),
                   'upstream_complete_probe_and_residue_contract_required': True,
                   'conditioning_or_renormalization_of_original_records': False})
    return result


def controller_from_ratios(proportional_coefficients):
    """Remove a common scale using the admitted positive-root trace identity.

    This assumes the coefficient spectrum belongs to this source branch.
    It is not a procedure for turning arbitrary physical masses into source
    data or for removing the ancestry of a calibrated controller.
    """
    values = tuple(map(_real, proportional_coefficients))
    if len(values) != 3 or any(_sign(value) <= 0 for value in values):
        raise ValueError('three positive proportional source coefficients are required')
    a, _ = moments.source_coefficients()
    root_sum = _clean(sum(s.sqrt(value) for value in values))
    scale = _clean(root_sum**2/(3*a*a))
    result = controller_from_coefficients(values, scale)
    koide_ratio = _clean(sum(values)/root_sum**2)
    if _clean(result['r_squared']-(3*koide_ratio-1)) != 0:
        raise ArithmeticError('the positive-root ratio inverse lost scale invariance')
    result.update({'ratio_only_source_input': True,
                   'positive_root_trace_assumed_for_scale_reconstruction': True,
                   'common_coefficient_scale_reconstructed': scale,
                   'source_spectral_Koide_ratio': koide_ratio,
                   'r_squared_from_ratio': _clean(3*koide_ratio-1),
                   'physical_energy_calibration_required': False})
    return result


def source_example(r, phase, label):
    """Verify source controls with their literal 3x3 matrices and all flags."""
    data = attenuation.family(r, phase)
    z = data['phase']
    roots, eigenvalues = moments.roots_and_masses(s.re(z), s.im(z), r)
    if any(_sign(value) <= 0 for value in roots):
        raise ValueError('this source example requires the positive-C-root chamber')
    x = s.symbols('x')
    expected_poly = s.expand(s.prod(x-value for value in eigenvalues))
    actual_poly = data['Y'].charpoly(x).as_expr()
    if _clean(expected_poly-actual_poly) != 0:
        raise ArithmeticError('the complete literal family spectrum differs from the source eigenvalues')
    if (_clean(sum(eigenvalues)-s.trace(data['Y'])) != 0
            or _clean(sum(v*v for v in eigenvalues)-s.trace(data['Y']**2)) != 0):
        raise ArithmeticError('source spectral moments differ from the literal matrix')
    original_records = data['original_tracial_records']
    if len(original_records) != 16 or _clean(sum(original_records)-1) != 0:
        raise ArithmeticError('the coherent source analyzer lost a complementary flag')
    result = controller_from_coefficients(eigenvalues)
    if _clean(result['r']-data['r']) != 0:
        raise ArithmeticError('the source attenuation was not reconstructed')
    if result['phase_observable'] and _clean(result['cos_three_theta']-s.re(z**3)) != 0:
        raise ArithmeticError('the reconstructed spectral phase class differs')
    return {'label': label, 'declared_r': data['r'], 'declared_unit_phase': z,
            'literal_Y': data['Y'], 'full_spectrum': eigenvalues,
            'all16_original_tracial_flag_probabilities': original_records,
            'controller_readout': result}


def spectral_inverse_join():
    """Evaluate the inherited exact line formula; no large carrier simulation.

    The charged source's full-CAR Gram theorem supplies this line formula.
    Only its three-dimensional Y input and scalar spectral data are needed
    for this executable composition into the owning inverse.
    """
    y = attenuation.family('1', 1)['Y']
    g2 = _clean(4*s.trace(y.H*y))
    gap = _clean(s.sqrt(1+4*g2))
    small_input = {'Y': y, 'pair_cost': s.S.One, 'fermion_energy': s.S.Zero,
                   'neutral': {'g_squared': g2, 'gap': gap, 'pair_cost': s.S.One,
                               'ground_energy': (1-gap)/2}}
    original_rows = charged.response_lines(small_input)
    reconstructed = [bare_effective.channel_reconstruction(gap, row) for row in original_rows]
    results = {channel: controller_from_reconstructed(reconstructed, channel)
               for channel in ('left_doublet', 'charged_singlet')}
    for result in results.values():
        if result['r'] != 1 or result['cos_three_theta'] != 1 or result['principal_abs_theta'] != 0:
            raise ArithmeticError('complete charged-line inversion changed the source controller class')
    return {'source_control': 'r=1, z=1 (theta=0)',
            'complete_original_spectral_rows': original_rows,
            'reconstructed_channels': results,
            'same_full_CAR_Gram_theorem_consumed': True,
            'large_carrier_simulation_reexecuted': False,
            'repeated_family_eigenvalue_retained': True}


def _serial(value):
    if isinstance(value, s.MatrixBase):
        return [[str(item) for item in row] for row in value.tolist()]
    if isinstance(value, s.Basic):
        return str(value)
    if isinstance(value, dict):
        return {key: _serial(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serial(item) for item in value]
    return value


def certificate():
    moments.moment_certificate()
    z = (63+16*s.I)/65
    controls = [source_example('1', 1, 'unattenuated zero phase'),
                source_example('1/2', z, 'half attenuation and z=(63+16i)/65'),
                source_example('0', z, 'zero attenuation, phase invisible')]
    ratio_control = controller_from_ratios([7*value for value in controls[1]['full_spectrum']])
    return _serial({'schema': 'oph.observer_outputs.controller_readout.v1',
        'source_definition': 'C=aI+rb(zR+conjugate(z)R^2), Y=C^2/3',
        'source_constants': {'a_squared': '5/36', 'b_squared': '5/72'},
        'inverse_r_squared': '36 Tr(Y)/5-1',
        'inverse_cos_three_theta': '1944[Tr(Y^2)-25/3888-25r^2/648-25r^4/2592]/(25sqrt(2)r^3)',
        'principal_phase_magnitude': 'acos(cos(3theta))/3 in[0,pi/3]',
        'unresolved_equivalence': 'theta -> +/-theta+2pi*k/3 with family permutation',
        'positive_C_root_domain_required_by_full_spectrum_admission': True,
        'complete_spectral_inverse_join': spectral_inverse_join(),
        'ratio_only_control': ratio_control,
        'source_controls': controls,
        'eliminated_unknown': 'the configured attenuation and principal phase class are reconstructible from the admitted same-clock finite source spectrum',
        'source_dynamics_selects_a_target_controller': False,
        'empirical_calibration_ancestry_is_retained': True,
        'measured_particle_inputs': [], 'physical_mass_predictions': [],
        'source_inputs': list(SOURCE_INPUTS)})


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
