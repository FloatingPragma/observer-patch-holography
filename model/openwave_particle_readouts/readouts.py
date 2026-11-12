"""Rational localization, persistence and rotation bounds on the native source."""
from functools import lru_cache
from pathlib import Path
import hashlib
import json

import sympy as s

from model.coherent_defect_wave.chiral_joint.interacting_charged_response.transfer import data as charged_data
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import operator_certificate
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import outward, sqrt_interval
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.particle_common_contract.provenance import import_closure

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANCHOR_GEV = s.Rational('0.00051099895069')
HORIZONS = (0, 1, 1000, 100000, 1000000)


def rational(x):
    value = s.sympify(x)
    if value.is_Rational is not True:
        raise ValueError('an exact rational value is required')
    return value


def up(x, bits=80):
    return s.Rational(outward((0, rational(x)), bits)[1])


def down(x, bits=80):
    return s.Rational(outward((rational(x), rational(x)), bits)[0])


def interval(values):
    return outward(tuple(map(rational, values)), 80)


def decimal(x):
    return str(s.N(x, 16))


def prepared_line_bounds(weight, error, perturbation, time):
    """Bounds for chi close to an atomic eigenvector supported at the origin.

    This controls a specific complete preparation; an arbitrary vector in
    the same energy bin need not satisfy its localization bound.
    """
    lo, hi = map(rational, weight)
    eps, v, t = map(rational, (error, perturbation, time))
    if not (0 < lo <= hi <= 1 and eps >= 0 and v >= 0 and t >= 0):
        raise ValueError('ordered positive weight and nonnegative errors required')
    root_lo, root_hi = sqrt_interval((lo, hi), 100)
    if root_lo <= eps:
        raise ValueError('the complete line has no proved positive probability')
    probability_low = (root_lo-eps)**2
    probability_high = min(s.Integer(1), (root_hi+eps)**2)
    outside = min(s.Integer(1), (eps+t*v*root_hi)**2/probability_low)
    # Compare the evolved normalized prepared state to its own initial state.
    # The two epsilon terms retain both preparation endpoints.
    phase_distance = min(s.Integer(2), (2*eps+t*v*root_hi)/(root_lo-eps))
    survival = max(s.Integer(0), 1-phase_distance**2/2)**2
    return {'probability': (probability_low, probability_high),
            'outside_origin_including_screened_cloud_upper': outside,
            'origin_probability_lower': 1-outside,
            'norm_distance_to_initial_ray_upper': phase_distance,
            'initial_state_return_probability_lower': survival}


def disturbance_bound(outside_upper, state_norm_error, extra_generator_norm, time):
    """Complete localization after a state disturbance and a bounded drive."""
    probability, delta, extra, t = map(rational, (outside_upper, state_norm_error, extra_generator_norm, time))
    if not (0 <= probability <= 1 and 0 <= delta <= 2 and extra >= 0 and t >= 0):
        raise ValueError('probability, state distance and nonnegative drive required')
    root_hi = sqrt_interval((probability, probability), 100)[1]
    return 1-min(s.Integer(1), (root_hi+delta+t*extra)**2)


@lru_cache(maxsize=2)
def inputs(family):
    if family not in ('ec', 'L'):
        raise ValueError('only source-derived charged lepton sectors admitted')
    source = charged_data(family)
    prep = source['preparation']['finite_energy_agreement_filter']
    trace = s.Rational(prep['conditional_one_source_trace_distance_upper'])
    fidelity = s.Rational(prep['conditional_ground_pair_fidelity_lower'])
    if not trace**2 >= 1-fidelity:
        raise ArithmeticError('purification error requires the source fidelity bound')
    # The actual prepared purification is within 2 trace of g_lambda tensor
    # an environment vector. The same contraction acts on both vectors.
    error = up(source['window_error']+2*trace)
    local = source['local']
    weight = tuple(s.Rational(v) for v in interval(local['rows'][0]['weight_interval']))
    energy = local['rows'][0]['energy_interval']
    ground = local['e_interval']
    atomic_gap = (energy[0]-ground[1], energy[1]-ground[0])
    return {'source': source, 'weight': weight, 'error': error,
            'preparation_probability_lower': s.Rational(prep['success_probability_lower']),
            'prepared_purification_error_upper': 2*trace,
            'atomic_gap': tuple(map(s.Rational, interval(atomic_gap)))}


def native_rotation_controls():
    """Native spin phases on both spin rows, with their quantum reference."""
    gauges = representation()[2]
    sz = s.diag(*([s.Rational(1, 2)]*15+[-s.Rational(1, 2)]*15))
    assert all(sz*s.kronecker_product(s.eye(2), g)
               == s.kronecker_product(s.eye(2), g)*sz for g in gauges)
    quarter = s.diag(*([s.exp(-s.I*s.pi/4)]*15+[s.exp(s.I*s.pi/4)]*15))
    assert quarter**4 == -s.eye(30) and quarter**8 == s.eye(30)
    # Every source Fock basis vector is covered by the occupation identity.
    # Counts, rather than 2^240 bitstrings, exhaust the phase possibilities.
    count_controls = 0
    for n_up in range(121):
        for n_down in range(121):
            rotation_sign = (-1)**((n_up-n_down) % 2)
            parity_sign = (-1)**(n_up+n_down)
            if rotation_sign != parity_sign:
                raise ArithmeticError('native spin rotation is not fermion parity')
            count_controls += 1
    apparatus = operator_certificate()
    if not all(r['full_spin_commutators'] == 3 for r in apparatus['all_five_native_elementary_probe_controls']):
        raise ArithmeticError('the actual charged pump must intertwine source/reference spin')
    p = s.symbols('p', nonnegative=True)
    # A pointer starts in |+>. On the odd source, controlled central
    # rotation is Z_pointer. Joint source/reference rotation is identity.
    plus = s.Matrix([1, 1])/s.sqrt(2)
    hadamard = s.Matrix([[1, 1], [1, -1]])/s.sqrt(2)
    full_turn = hadamard*s.diag(1, -1)*plus
    double_turn = hadamard*plus
    assert full_turn == s.Matrix([0, 1]) and double_turn == s.Matrix([1, 0])
    # A broken detector that also rotates the compensating odd reference
    # returns the opposite record. Dropping the coherent pointer fails too.
    joint_turn = hadamard*s.diag(1, (-1)*(-1))*plus
    assert joint_turn == double_turn and joint_turn != full_turn
    pointer_density=plus*plus.H
    dephased=s.diag(pointer_density[0,0],pointer_density[1,1])
    dephased_output=hadamard*s.diag(1,-1)*dephased*s.diag(1,-1)*hadamard.H
    assert dephased_output==s.eye(2)/2
    return {'native_modes_per_cell': 30, 'source_cells': 8,
            'source_spin_z_generator': '(N_up-N_down)/2',
            'native_one_particle_quarter_turn_power_four': '-I_30',
            'native_one_particle_quarter_turn_power_eight': 'I_30',
            'all_occupation_count_controls': count_controls,
            'full_source_central_rotation': 'exp[-2pi i S_z]=(-1)^N',
            'initial_neutral_ground_parity': 'even, continued from eight odd local ground doublets',
            'successful_native_removal_and_energy_filter_parity': 'odd',
            'source_rotation_commutes_with_all_12_Gauss_generators': True,
            'parity_preserved_by_every_H0_and_hopping_term': True,
            'same_source_clock_times': 'all real t',
            'two_pi_pointer_probabilities_on_line': ['0', '1'],
            'four_pi_pointer_probabilities_on_line': ['1', '0'],
            'joint_source_reference_rotation_pointer_probabilities': ['1', '0'],
            'source_ray_period': '2pi', 'source_amplitude_period': '4pi',
            'frame_coordinate_convention': 'phi=theta/2: ray period pi, amplitude period 2pi',
            'coordinate_relabelling_is_not_an_apolar_local_tensor_derivation': True,
            'native_rotation_is_declared_spin_action_on_the_actual_CAR_source': True,
            'intermediate_continuous_spin_rotation_is_a_symmetry_of_H_claimed': False,
            'reference_odd_parity_compensation_retained': True,
            'apparatus': 'controlled finite native spin rotation and two pointer Hadamards; finite coherent apparatus service',
            'passive_CE_compilation_or_laboratory_spin_identification_claimed': False,
            'broken_joint_rotation_rejected': True,
            'broken_pointer_dephasing_gives_half_half_instead_of_deterministic_record': True}


def certificate():
    from .motion import certificate as motion_certificate
    from .magnetic import certificate as magnetic_certificate

    ec = inputs('ec')
    # This rational centre is a unit convention, not a second observed mass.
    anchor_native = sum(ec['atomic_gap'])/2
    gev_per_native = ANCHOR_GEV/anchor_native
    families = []
    for family in ('ec', 'L'):
        d = inputs(family); source = d['source']
        initial = prepared_line_bounds(d['weight'], d['error'], source['v'], 0)
        rows = []
        for t in HORIZONS:
            b = prepared_line_bounds(d['weight'], d['error'], source['v'], t)
            # pi < 22/7 gives a rational lower bound on atomic phase cycles.
            cycles = s.Rational(t)*d['atomic_gap'][0]*7/44
            row = {'source_time': str(t), 'atomic_source_phase_cycles_lower': decimal(down(cycles)),
                   'origin_probability_lower': str(down(b['origin_probability_lower'])),
                   'origin_probability_lower_decimal': decimal(down(b['origin_probability_lower'])),
                   'all_nonorigin_including_cloud_probability_upper': str(up(b['outside_origin_including_screened_cloud_upper'])),
                   'initial_state_return_probability_lower': str(down(b['initial_state_return_probability_lower'])),
                   'initial_state_return_probability_lower_decimal': decimal(down(b['initial_state_return_probability_lower'])),
                   'norm_distance_to_initial_atomic_phase_upper': str(up(b['norm_distance_to_initial_ray_upper']))}
            disturbed=disturbance_bound(b['outside_origin_including_screened_cloud_upper'],
                                       s.Rational(1,1000), source['v']/100, t)
            row['origin_probability_lower_with_state_and_generator_disturbance']=str(down(disturbed))
            row['origin_probability_lower_with_disturbance_decimal']=decimal(down(disturbed))
            rows.append(row)
        physical_bin = tuple(x*gev_per_native for x in source['source_excitation'])
        q = d['preparation_probability_lower']
        p_low = initial['probability'][0]
        families.append({'family': family,
            'source_hypercharge_on_removal': '-1' if family == 'ec' else '+1/2',
            'native_electroweak_representation': 'weak and colour singlet' if family=='ec' else 'weak doublet, colour singlet',
            'atom_weight_interval': interval(d['weight']),
            'actual_preparation_to_atomic_line_vector_error_upper': str(d['error']),
            'prepared_purification_vector_error_upper': str(up(d['prepared_purification_error_upper'])),
            'complete_preparation_success_lower': str(down(q)),
            'line_probability_conditional_on_preparation': interval(initial['probability']),
            'line_probability_from_original_input_lower': str(down(q*p_low)),
            'atomic_source_gap_interval': interval(d['atomic_gap']),
            'atomic_source_gap_decimal': [decimal(x) for x in d['atomic_gap']],
            'source_energy_bin_after_explicit_reference_subtraction': interval(source['source_excitation']),
            'conditional_GeV_energy_bin': interval(physical_bin),
            'conditional_GeV_energy_bin_decimal': [decimal(x) for x in physical_bin],
            'conditional_atomic_GeV_gap_decimal': decimal(sum(d['atomic_gap'])/2*gev_per_native),
            'persistence': rows,
            'complete_records': {
                'preparation_failure': '1-q', 'pump_failure': 'q p_failure',
                'energy_bin_complement': 'q p_outside',
                'line_and_origin': 'q p_line P(origin|line,t)',
                'line_and_nonorigin_including_cloud': 'q p_line [1-P(origin|line,t)]',
                'normalization': '(1-q)+q(p_failure+p_outside+p_line)=1',
                'retained_source_environment_and_static_reference': True,
                'cloud_discarded_or_assumed_absent': False},
            'spin_pointer_records': {'preparation_failure':'1-q', 'pump_failure':'q p_failure',
                'energy_bin_complement':'q p_outside',
                'line_pointer_zero_two_pi':'0', 'line_pointer_one_two_pi':'q p_line',
                'line_pointer_zero_four_pi':'q p_line', 'line_pointer_one_four_pi':'0',
                'all_source_and_reference_outputs_retained':True}})
    v = ec['source']['v']
    if any(s.Rational(row['persistence'][-1]['origin_probability_lower']) <= s.Rational(999, 1000)
           for row in families):
        raise ArithmeticError('million-time-unit localization target failed')
    seed_paths = [str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')]
    seed_paths += ['model/coherent_defect_wave/chiral_joint/interacting_charged_response/PROOF.md',
                  'model/coherent_defect_wave/chiral_joint/interacting_charge_localization/PROOF.md',
                  'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json',
                  'model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/results.json']
    pins = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in sorted(import_closure(seed_paths))}
    return {'schema':'oph.prepared_particle_readouts.v1', 'input_sha256':pins,
            'source':{'Hamiltonian':'same full N2 H0+2^-34 T+4 N_ref',
                'full_hopping_perturbation_norm_upper':str(v),
                'source_clock':'exp[-it(H_lambda+4 N_ref)]',
                'fixed_reference_cell':0, 'reference_cost_subtracted_from_gap': '4',
                'energy_radius':'1/1024',
                'physical_ground_preparation_and_all_failures_retained':True},
            'calibration':{'electron_mass_anchor_GeV':str(ANCHOR_GEV),
                'native_ec_atomic_gap_centre':str(anchor_native),
                'GeV_per_native_energy':str(gev_per_native),
                'GeV_per_native_energy_decimal':decimal(gev_per_native),
                'scale_fitted_once_to_ec_channel':True,
                'electron_identification_or_GeV_mass_prediction_claimed':False,
                'L_transition_identified_as_neutrino_or_muon_mass':False},
            'prepared_lines':families,
            'prepared_line_motion':motion_certificate(),
            'same_source_hypercharge_background':magnetic_certificate(ec['source']['neutral']['coupling']),
            'native_spin_rotation':native_rotation_controls(),
            'clock':{'measured_object':'Ramsey characteristic function of the actual prepared line under its same total Hamiltonian',
                'native_atomic_phase':'a_selected+7 e_ground+4',
                'subtract_ground_and_reference_phase_explicitly':True,
                'phase_distance_formula':'(2 epsilon + |t| v sqrt(w_upper))/sqrt(p_line_lower)',
                'bounded_clock_is_self_starting_frequency_rigid_time_crystal':False,
                'retained_reference_interferometry_is_declared_apparatus':True},
            'disturbance_control':{'normalized_prepared_state_norm_distance_upper':'1/1000',
                'additional_Hermitian_generator_norm_upper':str(v/100),
                'state_and_generator_must_preserve_physical_source_constraints':True,
                'localization_formula':'1-min(1,[sqrt(original_outside_bound)+delta+|t| eta_H]^2)',
                'changed_generator_is_a_control_not_silently_the_original_source':True},
            'scope':{'prepared_line_screened_cloud_localization_control':True,
                'uniform_all_states_in_bin_origin_localization_claimed':False,
                'finite_horizon_localization_and_state_persistence':True,
                'infinite_time_localization_or_asymptotic_free_particle':False,
                'fixed_compensator_binding_removed':False,
                'source_spin_double_cover_on_actual_prepared_odd_state':True,
                'relativistic_mass_shell_or_physical_electron_dictionary':False},
            'verdict':'PASS_PREPARED_LINE_MOTION_SPIN_AND_STATIC_HYPERCHARGE_RESPONSE'}


def validate_pins(value):
    for name, expected in value['input_sha256'].items():
        p = Path(name)
        if p.is_absolute() or '..' in p.parts or hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=expected:
            raise ValueError('source pin mismatch: '+name)
