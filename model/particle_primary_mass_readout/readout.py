"""Cooling transports original transition weights to final excitation energy.

Momentum, a relativistic shell and physical time units require their own
same-source attachments. The primary report supplies none of those premises.
"""
from fractions import Fraction as F

import sympy as s

from model.particle_common_contract import contract
from model.source_mass_shell import units
from model.relativistic_sector_audit.sector import sqrt_interval

CLAIM = 'same-source transfer weights to ground-relative excitation intervals'
SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_filtered_probe/PROOF.md',
    'model/particle_positive_window/PROOF.md',
    'model/particle_primary_survival/survival.py',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_joint_reconstruction/POSITIVE_ENERGY.md',
    'model/source_mass_shell/units.py',
    'model/source_mass_shell/CALIBRATION_PROOF.md',
    'model/relativistic_sector_audit/sector.py',
    'model/spin_shell_readout/PROOF.md',
    'model/particle_joint_translation_readout/PROOF.md',
    'model/particle_gauge_characteristics/PROOF.md',
)


def nonnegative(value, name):
    value = units.exact(value)
    if value < 0:
        raise ValueError(name+' must be nonnegative')
    return value


def transport_bound(mean_excess, norm_squared=1, filter_half_width=None):
    """An unnormalized coupling; no division by probe success probability."""
    epsilon = nonnegative(mean_excess, 'initial ground excess')
    cap = nonnegative(norm_squared, 'squared probe norm')
    cost = cap*epsilon
    result = {'coupling_transport_cost_upper': cost,
              'Lipschitz_test_error_per_unit_Lipschitz_constant': cost,
              'original_total_mass_upper': cap,
              'equal_original_marginal_masses': True,
              'zero_response_allowed': True,
              'success_probability_division': False}
    if filter_half_width is not None:
        width = units.exact(filter_half_width)
        if width <= 0:
            raise ValueError('positive source-time filter half width required')
        result['triangular_filter_half_width'] = width
        result['transfer_first_absolute_moment_upper'] = 2*cap/width
        result['final_excitation_first_moment_upper'] = 2*cap/width+cost
    return result


def excitation_window(transfer_interval, original_weight_lower,
                      mean_excess, initial_energy_cut, norm_squared=1):
    """Infer final E-E0 weight from an ORIGINAL positive-transfer window.

The input weight must include the spectral readout's complete error budget.
The inferred interval uses the certified cooling premise, not an executed
absolute-energy measurement or a numerical determination of E0.
"""
    lower, upper = units.interval(transfer_interval, nonnegative=True)
    cap = nonnegative(norm_squared, 'squared probe norm')
    weight = nonnegative(original_weight_lower, 'original weight lower bound')
    cut = units.exact(initial_energy_cut)
    if weight > cap or cut <= 0:
        raise ValueError('weight at most the norm cap and positive initial cut required')
    cost = transport_bound(mean_excess, cap)['coupling_transport_cost_upper']
    loss = min(cap, cost/cut)
    return {'transfer_energy_interval': (lower, upper),
            'initial_energy_cut': cut,
            'initial_excitation_tail_weight_upper': loss,
            'inferred_final_excitation_interval': (lower, upper+cut),
            'original_final_excitation_weight_lower': max(F(0), weight-loss),
            'ground_origin': 'true ground of the same physical source carrier',
            'numerical_ground_energy_required_for_this_inference': False,
            'absolute_final_energy_measured': False,
            'rest_mass_identified': False}


def absolute_to_excitation(final_energy_interval, ground_interval=None):
    """A different readout route: subtract a supplied certified ground interval."""
    lower, upper = units.interval(final_energy_interval)
    if ground_interval is None:
        return {'status': 'GROUND_ENERGY_INTERVAL_REQUIRED',
                'final_excitation_interval': None,
                'numerical_ground_energy_computed': False}
    floor, ceiling = units.interval(ground_interval)
    if upper < floor:
        raise ValueError('final energy and ground floor have no physical realization')
    return {'status': 'CONDITIONAL_ON_CERTIFIED_GROUND_INTERVAL',
            'final_excitation_interval': (max(F(0), lower-ceiling), upper-floor),
            'supplied_ground_interval': (floor, ceiling),
            'numerical_ground_energy_computed': False}


def momentum_from_transfer(transfer_norm_interval, initial_radius,
                           initial_tail_probability, norm_squared=1):
    """Conditional on actual joint translation spectra, not Fourier labels."""
    lower, upper = units.interval(transfer_norm_interval, nonnegative=True)
    radius = nonnegative(initial_radius, 'initial momentum radius about vacuum')
    tail = nonnegative(initial_tail_probability, 'initial momentum tail')
    cap = nonnegative(norm_squared, 'squared probe norm')
    if tail > 1:
        raise ValueError('initial momentum tail is a probability')
    return {'final_momentum_norm_interval': (max(F(0), lower-radius), upper+radius),
            'additional_original_weight_loss_upper': cap*tail,
            'joint_translation_representation_assumed': True,
            'initial_momentum_control_derived_from_cooling_alone': False}


def conditional_rest_interval(excitation_interval, momentum_norm_interval,
                              speed_interval=(1, 1), squared_shell_error=0):
    """Enclose one admitted branch's M from |E^2-v^2|p|^2-M^2|<=eta.

This arithmetic does not identify that branch, prove a pole, or construct
commuting translation generators. Input intervals refer to the same event.
"""
    elo, ehi = units.interval(excitation_interval, nonnegative=True)
    plo, phi = units.interval(momentum_norm_interval, nonnegative=True)
    vlo, vhi = units.interval(speed_interval, positive=True)
    defect = nonnegative(squared_shell_error, 'squared shell error')
    lower = max(F(0), elo*elo-vhi*vhi*phi*phi-defect)
    upper = ehi*ehi-vlo*vlo*plo*plo+defect
    if upper < lower:
        raise ValueError('energy, momentum and nonnegative shell premises are inconsistent')
    rest = sqrt_interval(lower)[0], sqrt_interval(upper)[1]
    return {'status': 'CONDITIONAL_SAME_BRANCH_REST_ENERGY_ENCLOSURE',
            'rest_energy_squared_interval': (lower, upper),
            'rest_energy_interval': rest,
            'physical_momentum_or_shell_constructed': False,
            'positive_original_weight_required': True,
            'isolated_pole_or_particle_identity_proved': False}


def architectural_energy_units(model_energy_interval, gev_per_source_energy,
                               declaration):
    """An explicit unit assignment has no implied physical clock calibration."""
    lower, upper = units.interval(model_energy_interval, nonnegative=True)
    slo, shi = units.interval(gev_per_source_energy, positive=True)
    if not isinstance(declaration, str) or not declaration.strip():
        raise ValueError('an explicit architecture declaration is required')
    return {'status': 'DECLARED_ARCHITECTURAL_ENERGY_UNIT',
            'energy_GeV_interval_under_declaration': (lower*slo, upper*shi),
            'declaration': declaration,
            'physical_SI_clock_identification_proved': False,
            'source_derived_physical_mass_GeV': None}


def declared_SI_comparison(model_energy_interval, *, metres_per_source_coordinate=None,
                           seconds_per_source_time=None, model_speed=1):
    """Choose one absolute scale and match the declared signal speed to c.

This permitted comparison convention supplies units immediately. Agreement
of the model signal and excitation with observed physics is a separate claim.
"""
    lower, upper = units.interval(model_energy_interval, nonnegative=True)
    velocity = units.exact(model_speed)
    if velocity <= 0:
        raise ValueError('a positive declared source signal speed is required')
    if (metres_per_source_coordinate is None) == (seconds_per_source_time is None):
        raise ValueError('declare exactly one absolute source length or time scale')
    if metres_per_source_coordinate is not None:
        length = units.interval(metres_per_source_coordinate, positive=True)
        time = tuple(x*velocity/units.LIGHT_SPEED_M_S for x in length)
        origin = {'metres_per_source_coordinate': length}
    else:
        time = units.interval(seconds_per_source_time, positive=True)
        origin = {'seconds_per_source_time': time}
    conversion = units.c_normalized_units(time, velocity)
    elo, ehi = map(F, conversion['energy_unit_GeV_interval'])
    return contract.serial({
        'status': 'DECLARED_SI_COMPARISON_SCALE', 'absolute_scale_declaration': origin,
        'model_signal_speed_declared_as_c': velocity,
        'seconds_per_source_time_interval': time,
        'metres_per_source_coordinate_interval': conversion['length_unit_m_interval'],
        'GeV_per_source_energy_interval': (elo, ehi),
        'energy_GeV_interval': (lower*elo, upper*ehi),
        'mass_unit_for_a_valid_rest_energy': 'GeV/c^2',
        'c_alone_fixes_absolute_scale': False,
        'empirical_agreement_established_by_unit_declaration': False,
        'mass_or_shell_identified_by_unit_declaration': False})


def primary_excitation_bound(member, transfer_interval, original_weight_lower,
                             initial_energy_cut, *, evidence_binding,
                             provenance, norm_squared=1):
    """Authenticate the primary premise, retaining supplied evidence as supplied."""
    contract.verify_member(member)
    required = {'program': contract.PROGRAM, 'member_id': member['member_id'], **member['ids']}
    if evidence_binding != required:
        raise ValueError('spectral evidence must use this source, state, reference and clock')
    if (not isinstance(provenance, (tuple, list)) or not provenance
            or any(not isinstance(item, str) or not item.strip() for item in provenance)):
        raise ValueError('explicit spectral-weight provenance required')
    result = excitation_window(transfer_interval, original_weight_lower,
                               member['state']['target_ground_excess'],
                               initial_energy_cut, norm_squared)
    return contract.serial({**required, 'result': result,
                            'spectral_weight_provenance': provenance,
                            'supplied_spectral_weight_truth_verified_by_arithmetic': False,
                            'physical_rest_mass_GeV': None})


def common_contract_payload(member):
    contract.verify_member(member)
    n = member['source']['level']
    epsilon = F(member['state']['target_ground_excess'])
    cut = F(1, 2**(n//2))
    if (not member['state']['source_stationary']
            or member['reference']['idle_H'] != 0
            or member['clock']['GeV_or_proper_time_conversion'] is not None):
        raise ValueError('the selected source stationarity, idle reference or unit domain changed')
    return contract.serial({
        'program': contract.PROGRAM, 'member_id': member['member_id'], **member['ids'],
        'claim': CLAIM, 'source_level': n,
        'transport_for_contraction': transport_bound(epsilon, 1, F(1, 16)),
        'shrinking_initial_excess_cut': cut,
        'window_weight_allowance_for_contraction': epsilon/cut,
        'selected_reference_counter_evolved': False,
        'source_energy_zero_changed': False,
        'true_ground_used_only_in_cooling_and_inferred_excitation_definition': True,
        'numerically_computed_primary_ground_energy': None,
        'measured_primary_transition_window_consumed_here': None,
        'numerically_evaluated_primary_original_window_weight': None,
        'actual_primary_translation_momentum_supplied': False,
        'quadratic_characteristic_match_is_Lorentz_shell': False,
        'primary_rest_energy_interval': None,
        'seconds_per_source_time': None,
        'primary_physical_mass_GeV': None,
        'unit_relation': 'E_GeV=hbar_GeV_s E_source/T; T=seconds per source time',
        'physical_unit_proof': 'model/source_mass_shell/CALIBRATION_PROOF.md',
        'raw_transition_is_rest_energy': False,
        'broad_positive_energy_window_is_particle_pole': False,
        'static_Yukawa_is_measured_mass': False,
        'all_original_response_weights_retained': True,
    })


def common_contract_bound(member):
    payload = common_contract_payload(member)
    return {'schema': 'oph.primary_mass_readout_bound.v1', 'program': contract.PROGRAM,
            'member_id': member['member_id'], **member['ids'], 'claim': CLAIM,
            'result': payload, 'result_sha256': contract.digest(payload)}


def verify_bound(member, receipt):
    if receipt != common_contract_bound(member):
        raise ValueError('the excitation-readout bound differs from this primary member')


def primary_survival_excitation(member):
    """Consume the actual factored finite band without dividing its small weight."""
    from model.particle_primary_survival import survival
    contract.verify_member(member)
    witness = survival.selected_witness(member)
    width = F(member['state']['spectral_upper'])-F(member['state']['spectral_lower'])
    if F(witness['original_positive_band']['upper']) != width:
        raise ValueError('the original band must use this full source width')
    return {'program': contract.PROGRAM, 'member_id': member['member_id'], **member['ids'],
            'inherited_witness_sha256': contract.digest(witness),
            'inferred_final_excitation_interval': witness['original_positive_band'],
            'original_final_excitation_weight_lower': witness['original_positive_band_weight_lower'],
            'finite_filter_final_excitation_weight_lower': witness[
                'finite_apparatus_original_positive_band_weight_lower'],
            'proof': 'E_final-E0 >= E_final-E_initial; 0 <= E_final-E0 <= U-L',
            'finite_filter': 'actual algebraic Cayley/quadrature blank operator; original H spectral measure',
            'additional_weight_loss': '0',
            'same_original_weight_retained': True,
            'factored_probability_or_lower_energy_materialized': False,
            'positive_finite_excitation_weight_proved': True,
            'finite_filter_positive_excitation_weight_proved': True,
            'fixed_band_noncollapse_or_rest_mass_proved': False,
            'GeV_rest_mass': None}


def finite_control():
    """Separate exact mixed-state control; no selected primary matrix is substituted."""
    h = s.diag(0, 3, 4)
    rho = s.diag(s.Rational(15, 16), s.Rational(1, 16), 0)
    b = s.Matrix([[0, 0, 0], [0, 0, 0], [0, 1, 0]])
    failure = s.diag(1, 0, 1)
    if b.H*b+failure.H*failure != s.eye(3) or h*rho != rho*h:
        raise ArithmeticError('the independent complete stationary control failed')
    rows = []
    for initial in range(3):
        for final in range(3):
            weight = rho[initial, initial]*abs(b[final, initial])**2
            rows.append({'initial_energy': int(h[initial, initial]),
                         'final_energy': int(h[final, final]),
                         'transfer': int(h[final, final]-h[initial, initial]),
                         'weight': F(int(s.numer(weight)), int(s.denom(weight)))})
    cost = sum((row['weight']*row['initial_energy'] for row in rows), F(0))
    if cost != F(3, 16):
        raise ArithmeticError('the sharp initial-excess transport cost changed')
    return {'H': h, 'rho': rho, 'probe': b, 'failure': failure,
            'joint_energy_rows': rows,
            'receipt': {'is_primary_source': False, 'source_dimension': 3,
                'reference_dimension': 3, 'original_response_mass': F(1, 16),
                'original_failure_mass': F(15, 16), 'mean_initial_excess': F(3, 16),
                'transfer_energy': 1, 'final_ground_relative_energy': 4,
                'exact_transport_cost': cost, 'transport_norm_bound_saturated': True,
                'initial_and_final_momentum_control': (1, 1),
                'momentum_transfer_control': 0,
                'naive_transfer_invariant_squared': 1,
                'final_invariant_squared_at_speed_one': 15,
                'these_control_momenta_define_a_primary_translation': False}}


def certificate():
    return contract.serial({
        'schema': 'oph.primary_mass_readout.v1',
        'selected_members': [common_contract_bound(contract.member(n)) for n in (0, 2, 4)],
        'actual_primary_finite_excitation_band': primary_survival_excitation(contract.member(1)),
        'independent_mixed_control': finite_control()['receipt'],
        'direct_absolute_readout_without_ground': absolute_to_excitation((4, 4)),
        'conditional_floor_example': absolute_to_excitation((4, 4), ('-1/10', '1/10')),
        'conditional_rest_example': conditional_rest_interval((5, 5), (3, 3)),
        'declared_one_metre_comparison_example': declared_SI_comparison(
            (4, 4), metres_per_source_coordinate=(1, 1)),
        'unresolved_unit_example': units.energy_conversion((4, 4),
            generator_id='separate conditional shell example; not a primary observation'),
        'physical_primary_mass_table': [],
        'accepted_physical_primary_clock_inputs': [],
    })
