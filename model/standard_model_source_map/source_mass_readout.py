"""Exact ratios of the existing source-selected static action coefficients."""
from fractions import Fraction as F

from model.native_flavor_orbit import packet as flavor, intervals as iv
from tools import strict_json


def _box(pair):
    if type(pair) is not list or len(pair) != 2 or any(type(x) is not str for x in pair):
        raise ValueError('two exact source interval endpoints required')
    lo, hi = map(F, pair)
    if lo < 0 or hi < lo:
        raise ValueError('nonnegative ordered source interval required')
    return lo, hi


def koide_interval(boxes):
    """Physical-positive-root Koide ratio for the stated action-mass intervals."""
    values = [iv.I(lo, hi) for lo, hi in boxes]
    if len(values) != 3 or sum((x.lo for x in values), F()) <= 0:
        raise ValueError('three nonnegative masses with positive total lower bound required')
    roots = [iv.I(iv.sqrt_rational(x.lo).lo, iv.sqrt_rational(x.hi).hi) for x in values]
    result = sum(values, iv.I(0))/sum(roots, iv.I(0))**2
    return {'positive_root_Koide_ratio_interval': result.strings(),
            'balanced_value_two_thirds_excluded': not result.lo <= F(2, 3) <= result.hi,
            'this_checks_reported_action_masses_not_an_identified_physical_spectrum': True}


def summarize(source):
    """Transform a source receipt; caller owns authentication of its input."""
    if source['schema'] != 'oph.native_flavor_orbit.v1':
        raise ValueError('the owning native flavor source is required')
    scope = source['scope']
    if (scope['physical_pole_masses_or_absolute_units_derived'] is not False
            or scope['static_mean_action_equals_finite_correlated_memory_channel'] is not False):
        raise ValueError('source scope changed; physical interpretation requires review')
    policies = []
    for policy in source['current_gated_joint_spectrum']['policies']:
        channels = []
        for channel in policy['channels']:
            boxes = [_box(pair) for pair in channel['dimensionless_singular_values']]
            if len(boxes) != 3 or any(a[1] >= b[0] for a, b in zip(boxes, boxes[1:])):
                raise ValueError('three separated ascending source coefficients required')
            if boxes[0][0] > 0:
                ratios = [['1', '1']]+[[str(lo/boxes[0][1]), str(hi/boxes[0][0])]
                                      for lo, hi in boxes[1:]]
                ratio_status = 'finite, with propagated source intervals'
            elif boxes[0] == (0, 0):
                ratios = None
                ratio_status = 'undefined: the smallest action coefficient is exactly zero'
            else:
                raise ValueError('the first source coefficient is not resolved from zero')
            channels.append({'channel': channel['channel'],
                'ascending_dimensionless_action_mass_intervals':
                    [[str(lo), str(hi)] for lo, hi in boxes],
                'ratios_to_smallest_interval': ratios,
                'ratio_status': ratio_status,
                'Koide_diagnostic': koide_interval(boxes),
                'K45_multiplicity_per_family': channel['K45_multiplicity_per_family']})
        if len(channels) != 3 or {c['channel'] for c in channels} != {'up', 'down', 'charged_lepton'}:
            raise ValueError('incomplete charged source channels')
        policies.append({'policy': policy['policy'], 'channels': channels,
            'minimal_neutrino_action_coefficients': ['0', '0', '0'],
            'charged_current_CP_invariant': policy['charged_current_CP_commutator_trace_cube_exact'],
            'units': policy['units']})
    if [p['policy'] for p in policies] != ['always', 'C', 'E', 'CE']:
        raise ValueError('all four configured feedback policies must be retained')
    return {'source': 'native_flavor_orbit static mean-feedback action',
        'policies': policies,
        'same_Higgs_and_canonical_kinetic_quadratic_mass_rule': 'm_f=w_H*sigma_f(Y)',
        'reported_neutral_Higgs_norm_w_H': '1',
        'physical_particle_labels_for_family_indices': [],
        'physical_pole_mass_table': [],
        'empirical_mass_inputs_consumed': [],
        'source_parameters_optimized_by_this_producer': False,
        'runtime_empirical_separation_proves_target_independent_design_history': False,
        'source_static_action_replaces_full_correlated_channel': False,
        'quadratic_action_mass_is_interacting_pole_mass': False,
        'GeV_conversion_supplied': False,
        'massless_neutrinos_identified_in_nature': False}


def certificate(*, cold_replay=False):
    """Authenticate the owner; an explicit cold mode replays its full experiment."""
    if type(cold_replay) is not bool:
        raise ValueError('the replay policy must be a boolean')
    pins = flavor.source_bindings()
    source = flavor.report() if cold_replay else strict_json.load(flavor.OUTPUT)
    if source.get('source_bindings') != pins or pins != flavor.source_bindings():
        raise ValueError('native flavor receipt is stale or incompletely bound')
    result = summarize(source)
    result['authenticated_source_bindings'] = source['source_bindings']
    result['owning_source_receipt_authenticated'] = True
    result['owning_source_cold_scientific_replay_performed_here'] = cold_replay
    return result
