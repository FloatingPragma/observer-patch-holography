"""Three-input inverse controller; no held-out mass or source selection."""
from fractions import Fraction as F
import json
import math
from pathlib import Path

from model.native_flavor_orbit.intervals import I, sqrt_rational

SOURCE_INPUTS = (
    'model/particle_candidate_atlas/reference_values.json',
    'model/particle_candidate_atlas/attenuation.py',
    'model/particle_candidate_atlas/BARE_EFFECTIVE_CONTROLLER.md',
    'model/particle_candidate_atlas/CHARGED_FIELD_ADMISSION.md',
    'model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md',
)


def root(x):
    x = x if isinstance(x, I) else I(x)
    return I(sqrt_rational(x.lo).lo, sqrt_rational(x.hi).hi)


def midpoint_rational(x, digits=40):
    denominator = 10**digits
    return F(round((x.lo+x.hi)*denominator/2), denominator)


def fit(masses):
    """Calibrate an ordered positive triplet in the admitted attenuator range."""
    if len(masses) != 3 or any(type(x) not in (str, int, F) for x in masses):
        raise ValueError('three exact mass inputs required')
    masses = tuple(map(F, masses))
    if not 0 < masses[0] < masses[1] < masses[2]:
        raise ValueError('strictly ordered positive charged masses required')
    light, middle, heavy = [root(x) for x in masses]
    mean = (light+middle+heavy)/3
    real = (heavy-(light+middle)/2)/3
    imag = root(3)*(middle-light)/6
    modulus = root(real**2+imag**2)
    amplitude = root(2)*modulus/mean
    if not 0 < amplitude.lo <= amplitude.hi < 1:
        raise ValueError('target is outside the certified attenuation-only range')
    tangent = imag/(modulus+real)
    scale = F(108, 5)*mean**2
    r, t, gain = map(midpoint_rational, (amplitude, tangent, scale))
    x, y = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    if x*x+y*y != 1 or not 0 < r < 1:
        raise ArithmeticError('rational complete-controller conditions failed')
    a, b, c = root(5)/6, root(10)/12, root(30)/12
    roots = (a-r*b*x-r*c*y, a-r*b*x+r*c*y, a+2*r*b*x)
    if roots[0].lo <= 0 or any(u.hi >= v.lo for u, v in zip(roots, roots[1:])):
        raise ArithmeticError('controller left the ordered positive-root chamber')
    outputs = [gain*q**2/3 for q in roots]
    errors = [max(abs(q.lo-m), abs(q.hi-m)) for q, m in zip(outputs, masses)]
    if max(errors) >= F(1, 10**36):
        raise ArithmeticError('rounded controller did not meet the certified target budget')
    q_target = I(sum(masses))/(light+middle+heavy)**2
    return {
        'input_mass_energies_GeV': list(map(str, masses)),
        'all_three_masses_used_as_inputs': True,
        'held_out_mass_predictions': 0,
        'source_selected_controller': False,
        'physical_mass_matching_proved': False,
        'exact_algebraic_fit_proved_by_inverse_DFT': True,
        'implemented_controller_uses_rational_approximations': True,
        'native_tape_for_this_fitted_point_executed': False,
        'pass_amplitude': str(r),
        'pass_amplitude_algebraic_interval': amplitude.strings(),
        'phase_half_angle_tangent': str(t),
        'phase_real': str(x), 'phase_imaginary': str(y),
        'phase_radians_approx': 2*math.atan(float(t)),
        'common_mass_scale_GeV': str(gain),
        'scale_is_selected_physical_Higgs_vev': False,
        'complete_retained_flags': 16,
        'loss_probability': str(F(5, 36)*(1-r*r)),
        'target_Koide_quotient_interval': q_target.strings(),
        'controller_Koide_quotient': str((1+r*r)/3),
        'controller_output_mass_energies_GeV': [q.strings() for q in outputs],
        'absolute_mass_error_bounds_GeV': list(map(str, errors)),
        'all_output_errors_less_than_1e_minus_36_GeV': True,
        'roots_are_strictly_positive_and_ordered': True,
        'phase_in_balanced_positive_chamber': (imag/real).hi < (2-root(3)).lo,
        'correction_computed_from_repair': False,
    }


def certificate():
    reference = json.loads(Path(__file__).with_name('reference_values.json').read_text())
    rows = {x['target']: x for x in reference['masses']}
    if any(rows[k]['unit'] != 'GeV' for k in ('electron', 'muon', 'tau')):
        raise ValueError('reference mass units changed')
    result = fit([rows[k]['value'] for k in ('electron', 'muon', 'tau')])
    result['schema'] = 'oph.calibrated_three_mass_controller.v1'
    result['input_reference_sources'] = {
        k: {'source': rows[k]['source'], 'url': rows[k]['source_url']}
        for k in ('electron', 'muon', 'tau')}
    result['reference_intervals_or_covariance_fitted'] = False
    result['reference_central_values_are_exact_true_masses'] = False
    return result
