"""Read family singular coefficients from complete finite interacting spectra.

The neutral gap and both originally weighted charged lines determine a
coefficient without the cross-sector additive offset. These are coefficients
of the specified finite vertex, not measured relativistic particle masses.
"""
from fractions import Fraction as F
from functools import cmp_to_key
import json

import sympy as s

from model.native_flavor_orbit import intervals as iv
from . import charged_higgs_records as charged
from . import higgs_records as neutral


SOURCE_INPUTS = (
    'model/particle_observer_outputs/charged_higgs_records.py',
    'model/particle_observer_outputs/higgs_records.py',
    'model/particle_observer_outputs/thermal_clock.py',
    'model/particle_calibration_selection/phase_attenuation.py',
    'model/particle_candidate_atlas/attenuation.py',
    'model/native_flavor_orbit/intervals.py',
)


def _real(value):
    value = neutral._clean(neutral._exact(value))
    if value.is_real is not True:
        raise ValueError('an exact real spectral coordinate is required')
    return value


def _sign(value):
    value = _real(value)
    if value == 0:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    enclosure = iv.expression(value)
    if enclosure.lo > 0:
        return 1
    if enclosure.hi < 0:
        return -1
    raise ValueError('the exact sign is unresolved by the admitted arithmetic')


def channel_reconstruction(neutral_gap, row):
    """A measured neutral gap and two resolved nonzero addition residues.

    row is one channel exported by the literal source response_lines.
    The energies may contain any common additive charged-sector offset.
    No branch or original event weight is renormalized.
    """
    channel = row['internal_channel']
    if channel not in ('left_doublet', 'charged_singlet'):
        raise ValueError('the actual charged source channel is required')
    kappa = 1 if channel == 'left_doublet' else 2
    lines = row['addition_lines']
    if len(lines) != 2 or {r['branch'] for r in lines} != {'lower', 'upper'}:
        raise ValueError('both complete addition branches are required')
    by_branch = {r['branch']: r for r in lines}
    probabilities = tuple(_real(by_branch[name]['unconditional_probe_probability'])
                          for name in ('lower', 'upper'))
    if any(_sign(p) <= 0 or _sign(1-p) < 0 for p in probabilities):
        raise ValueError('both original addition probabilities must be strictly positive')
    gap = _real(neutral_gap)
    lower, upper = (_real(by_branch[name]['ground_referenced_energy'])
                    for name in ('lower', 'upper'))
    split = neutral._clean(upper-lower)
    if _sign(gap) <= 0 or _sign(split) <= 0:
        raise ValueError('positive neutral gap and resolved charged separation required')
    squared = neutral._clean((gap*gap-split*split)/(4*kappa))
    if _sign(squared) < 0:
        raise ValueError('the spectral coordinates reconstruct a negative squared coefficient')
    return {'internal_channel': channel, 'kappa': kappa,
            'multiplicity': row['multiplicity'],
            'neutral_gap': gap, 'charged_line_separation': split,
            'original_addition_probabilities': probabilities,
            'reconstructed_sigma_squared': squared,
            'reconstructed_sigma': neutral._clean(s.sqrtdenest(s.sqrt(squared))),
            'charged_additive_offset_used': False,
            'physical_rest_mass_claimed': False}


def reconstruct_source(data):
    """Replay actual source spectra and compare against its family operator."""
    rows = charged.response_lines(data)
    recovered = [channel_reconstruction(data['neutral']['gap'], row) for row in rows]
    for actual, source_row in zip(recovered, rows, strict=True):
        if neutral._clean(actual['reconstructed_sigma_squared']-source_row['sigma_squared']) != 0:
            raise ArithmeticError('complete spectral inversion differs from the source family coefficient')
    return recovered


def coefficient_ratios(reconstructed, channel='left_doublet'):
    copies = 4 if channel == 'left_doublet' else 2 if channel == 'charged_singlet' else None
    if copies is None:
        raise ValueError('an actual source channel is required')
    values = []
    for row in reconstructed:
        if row['internal_channel'] == channel:
            if type(row['multiplicity']) is not int or row['multiplicity'] <= 0 or row['multiplicity'] % copies:
                raise ValueError('the complete source channel multiplicity differs')
            values.extend([row['reconstructed_sigma']]*(row['multiplicity']//copies))
    if len(values) != 3:
        raise ValueError('all three declared family coordinates must be retained')
    ordered = sorted(values, key=cmp_to_key(lambda a, b: _sign(a-b)))
    if _sign(ordered[0]) <= 0:
        raise ValueError('a positive reference coefficient is required for ratios')
    return tuple(neutral._clean(x/ordered[0]) for x in ordered)


def _interval(bounds):
    if isinstance(bounds, iv.I):
        return bounds
    if (not isinstance(bounds, (tuple, list)) or len(bounds) != 2
            or any(type(x) not in (int, F, str) for x in bounds)):
        raise ValueError('two exact rational interval endpoints are required')
    return iv.I(*bounds)


def interval_reconstruction(neutral_gap, lower_line, upper_line, *, kappa,
                            original_weight_lower_bounds):
    """Conservative rectangular data errors, with no hidden correlations.

    Experimental/source/detector errors must be included in the supplied
    intervals. An unresolved line or nonpositive original-weight lower
    bound is rejected; the method does not infer visibility from a fit.
    """
    if type(kappa) is not int or kappa not in (1, 2):
        raise ValueError('the source density factor must be one or two')
    if (not isinstance(original_weight_lower_bounds, (tuple, list))
            or len(original_weight_lower_bounds) != 2
            or any(type(x) not in (int, F, str) for x in original_weight_lower_bounds)):
        raise ValueError('both exact original-weight lower bounds are required')
    weights = tuple(F(x) for x in original_weight_lower_bounds)
    if any(not 0 < x <= 1 for x in weights):
        raise ValueError('both original lines must be certified visible')
    gap, lower, upper = map(_interval, (neutral_gap, lower_line, upper_line))
    split = upper-lower
    if gap.lo <= 0 or split.lo <= 0:
        raise ValueError('the complete intervals do not resolve positive gaps')
    squared = iv.I((gap.lo**2-split.hi**2)/(4*kappa),
                   (gap.hi**2-split.lo**2)/(4*kappa))
    if squared.lo < 0:
        raise ValueError('the uncertainty does not certify a nonnegative squared coefficient')
    coefficient = iv.I(iv.sqrt_rational(squared.lo).lo, iv.sqrt_rational(squared.hi).hi)
    return {'squared_coefficient': squared, 'coefficient': coefficient,
            'charged_line_separation': split,
            'original_weight_lower_bounds': weights,
            'physical_matching_error_included': False}


def ratio_interval(numerator, denominator):
    numerator, denominator = map(_interval, (numerator, denominator))
    if numerator.lo < 0 or denominator.lo <= 0:
        raise ValueError('nonnegative numerator and strictly positive denominator required')
    return numerator/denominator


def ambient_ground_control():
    """Actual Q*=2 variational witness, not a physical vacuum identification."""
    y = s.diag(1, 2, 3)
    vacuum_image = neutral.source.vertex_action(y, 0, 0)
    u = {(0, mask): value for (scalar, mask), value in vacuum_image.items() if scalar == 1}
    def apply_vertex(vector):
        result = {}
        for state, value in vector.items():
            for out, amplitude in neutral.source.vertex_action(y, *state).items():
                result[out] = result.get(out, s.S.Zero)+value*amplitude
        return {key: neutral._clean(value) for key, value in result.items() if value != 0}
    wu = apply_vertex(u)
    norm_u = sum(s.conjugate(x)*x for x in u.values())
    norm_w = sum(s.conjugate(x)*x for x in wu.values())
    if (norm_u, norm_w) != (28, 1764):
        raise ArithmeticError('the full-CAR trial-vector norms differ')
    if any(mask.bit_count()-2*(scalar != 0) != 2 for scalar, mask in (*u, *wu)):
        raise ArithmeticError('the trial leaves its stated Q*=2 domain')
    coupling_squared = neutral._clean(norm_w/norm_u)
    bound = (1-s.sqrt(1+4*coupling_squared))/2
    if _sign(bound+7) >= 0:
        raise ArithmeticError('the exact ambient trial does not lower the neutral-sector energy')
    twice = apply_vertex(wu)
    invariant = all(neutral._clean(twice.get(key, 0)-coupling_squared*u.get(key, 0)) == 0
                    for key in twice.keys() | u.keys())
    return {'fixture': 'Y=diag(1,2,3), epsilon_F=0, epsilon_H=1',
            'trial_Q_star': 2, 'unnormalized_u_norm_squared': norm_u,
            'unnormalized_Wu_norm_squared': norm_w, 'trial_coupling_squared': coupling_squared,
            'ambient_ground_energy_upper_bound': bound,
            'neutral_sector_ground_energy': s.Integer(-7),
            'trial_span_is_invariant': invariant,
            'neutral_sector_ground_is_full_ambient_ground': False,
            'trial_gauge_singlet_claimed': False,
            'physical_Gauss_vacuum_minimum_claimed': False,
            'compensating_reference_and_energy_would_require_a_separate_domain': True}


def _serial(value):
    if isinstance(value, dict):
        return {k: _serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(v) for v in value]
    if isinstance(value, s.Basic):
        return str(value)
    if isinstance(value, F):
        return str(value)
    return value


def certificate():
    fixture = reconstruct_source(charged.charged_data(s.diag(1, 2, 3), 1, 0))
    coherent = reconstruct_source(charged.coherent_charged_data())
    fixture_ratios = coefficient_ratios(fixture)
    coherent_ratios = coefficient_ratios(coherent)
    if fixture_ratios != (1, 2, 3) or any(neutral._clean(a-b) != 0 for a, b in
            zip(coherent_ratios, (1, 1, 34+24*s.sqrt(2)), strict=True)):
        raise ArithmeticError('spectral inversion changed a declared source family pattern')
    return _serial({'schema': 'oph.finite_spectral_coefficient_inverse.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'inverse': 'sigma_j(Y)^2=(Delta^2-(omega_j,+-omega_j,-)^2)/(4*kappa)',
        'line_visibility': 'both original nonzero addition probabilities required',
        'charged_additive_offset_used': False,
        'common_positive_clock_scale_cancels_from_coefficient_ratios': True,
        'fixture': {'reconstructed_channels': fixture, 'coefficient_ratios': fixture_ratios},
        'coherent_source': {'z': '1', 'phase_angle': '0', 'pass_amplitude': '1',
                            'reconstructed_channels': coherent, 'coefficient_ratios': coherent_ratios,
                            'actual_repeated_pair_retained': True},
        'ambient_ground_control': ambient_ground_control(),
        'chosen_phase_or_attenuation_inferred_as_unique': False,
        'target_fitted_phase_or_attenuator_derived': False,
        'measured_particle_inputs': [], 'physical_mass_predictions': []})


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
