"""A dimension-independent frequency cut and complete radial source clocks."""
from fractions import Fraction as F
from math import comb

from model.particle_common_contract import contract, provenance
from model.particle_candidate_scalar_gauge_response import model as radial
from model.particle_primary_joint_spectroscopy import spectroscopy

SOURCE_INPUTS = tuple('model/particle_candidate_scalar_spectrum/' + name for name in
    ('__init__.py', '__main__.py', 'model.py', 'test_model.py', 'PROOF.md', 'README.md')) + (
    'model/particle_candidate_scalar_gauge_response/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/spin_shell_readout/PROOF.md',
    *radial.SOURCE_INPUTS, *spectroscopy.SOURCE_INPUTS,
)


def ceil_log2(value):
    value = F(value)
    if value <= 0:
        raise ValueError('positive dyadic comparison required')
    n, d = value.numerator, value.denominator
    k = n.bit_length() - d.bit_length()
    below = (n <= d << k) if k >= 0 else (n << (-k) <= d)
    return k if below else k + 1


def dyadic(coefficient, exponent):
    coefficient = F(coefficient)
    if coefficient <= 0 or not isinstance(exponent, int) or exponent < 0:
        raise ValueError('positive rational and nonnegative exponent required')
    return {'coefficient': coefficient, 'base': 2, 'negative_exponent': exponent}


def plateau(frequency, cutoff):
    """Exact multiplier of 2 f_(2 lambda)-f_lambda."""
    x, cutoff = abs(F(frequency)), F(cutoff)
    if cutoff <= 0:
        raise ValueError('positive frequency cutoff required')
    return 2 * max(F(), 1 - x/(2*cutoff)) - max(F(), 1 - x/cutoff)


def derivative_certificate():
    # pi >= 3. On |s|<=1, |g'|<=1/(6 pi); on |s|>=1,
    # |g'| <=1/(pi s^2)+4/(pi |s|^3).
    central = F(1, 9)
    tail = F(2)
    total = central + tail
    if total >= 3:
        raise ArithmeticError('elementary kernel derivative bound failed')
    return {'kernel': 'f_lambda(t)=lambda sinc^2(lambda t/2)/(2pi)',
            'integral': 1, 'Fourier_transform': '(1-|omega|/lambda)_+',
            'unit_kernel_central_derivative_L1_upper': central,
            'unit_kernel_tail_derivative_L1_upper': tail,
            'unit_kernel_derivative_L1_upper': total,
            'simple_unit_kernel_derivative_L1_upper': 3,
            'plateau_kernel': '2 f_(2lambda)-f_lambda',
            'plateau_L1_upper': 3,
            'plateau_derivative_L1_coefficient_upper': 15,
            'plateau_multiplier_on_abs_omega_le_lambda': 1,
            'plateau_multiplier_on_abs_omega_ge_2lambda': 0,
            'high_multiplier_range': [0, 1],
            'auxiliary_filter_is_executed': False}


def clock_plan(norm, cutoff, weight_coefficient, floor_exponent):
    """Compile huge clocks using exponents, without constructing 2^N integers."""
    R, lam, q = F(norm), F(cutoff), F(weight_coefficient)
    N = floor_exponent
    if min(R, lam, q) <= 0 or q >= 1 or not isinstance(N, int) or N < 1:
        raise ValueError('positive norm/cutoff, weight coefficient below one and positive exponent required')
    eta = lam/8
    tau = 1/(4*R)
    if tau*eta > 1:
        raise ValueError('the actual signed energy chart is not admitted')
    margin = tau*eta/7
    k = max(0, ceil_log2(16/(margin*q)))
    bits = N + k
    h = max(0, ceil_log2(1/(48*q)))
    substep_bits = (bits + N + h + 1)//2
    # N cancels from each exact rational check. Never form a source-sized bitset.
    if F(1 << k)*margin*q/16 < 1:
        raise ArithmeticError('complete Fejer clock tail budget failed')
    if F(1, 768*(1 << (2*substep_bits-bits-N))) > q/16:
        raise ArithmeticError('complete Cayley clock error budget failed')
    return {'energy_resolution': eta, 'clock_step': tau,
            'turn_margin': margin, 'clock_qubits_each': bits,
            'clock_bins': {'base': 2, 'exponent': bits},
            'Cayley_substeps_per_base_tick': {'base': 2, 'exponent': substep_bits},
            'maximum_source_duration_each': {'coefficient': tau, 'power_two_minus_one': bits},
            'one_clock_tail_allowance': dyadic(q/16, N),
            'one_complete_clock_implementation_allowance': dyadic(q/16, N),
            'total_two_clock_event_loss_upper': dyadic(q/4, N),
            'source_energy_zero_changed': False, 'source_clock_changed': False,
            'signed_energy_centers': '2pi*(j if 2j<J else j-J)/(J tau)',
            'source_phase_radius_upper': F(1, 4), 'alias_probability_upper': 0,
            'huge_register_or_histogram_allocated': False,
            'two_bit_compact_apparatus_bound_applies': False}


def degree_census(side):
    if not isinstance(side, int) or side < 2:
        raise ValueError('actual cubic source side at least two required')
    rows = [{'degree': 3+j, 'cells': comb(3,j)*2**(3-j)*(side-2)**j}
            for j in range(4)]
    rows = [r for r in rows if r['cells']]
    if sum(r['cells'] for r in rows) != side**3:
        raise ArithmeticError('source degree census lost cells')
    return rows


def bind_member(member, cutoff='1/40'):
    contract.verify_member(member)
    lam = F(cutoff)
    a = F(member['source']['spacing'])
    source = member['preparation_plan']['source']
    R = max(abs(F(source['spectral_lower'])), abs(F(source['spectral_upper'])))
    D = F(source['spectral_upper'])-F(source['spectral_lower'])
    dynamics = radial.scalar_dynamics(member)
    local = radial.radial_control()['summary']
    records = radial.scalar_records(member)
    floor = dynamics['source_density_floor']
    N = floor['negative_exponent']
    slope_positive = F(dynamics['selected_modular_slope']['modular_slope_lower']) > 0
    if not slope_positive:
        raise ArithmeticError('the original selected density does not certify passivity')
    witness_norm = F(3, 4)/a
    if lam <= 0 or 15*lam > witness_norm/2:
        raise ValueError('frequency cutoff exceeds the proved radial commutator margin')
    gamma = F(9, 16)/a**2
    if F(dynamics['physical_witness']['A_first_norm_squared_coefficient']) != F(9,16):
        raise ArithmeticError('the actual radial commutator witness differs')
    degree = degree_census(member['source']['side'])
    if sum(r['degree']*r['cells'] for r in degree) != 2*member['source']['interior_hopping_bonds']:
        raise ArithmeticError('source edge/degree census differs')
    q = gamma/(8*D**2)
    improved = (witness_norm-15*lam)**2/(2*D**2)
    if improved < q or not 0 < q < 1:
        raise ArithmeticError('positive original frequency-cut bound failed')
    plan = clock_plan(R, lam, q, N)
    V = member['source']['side']**3
    elastic = F(records['Jensen_effect_lower'])
    if elastic < F(3481,160000):
        raise ArithmeticError('the inherited radial mean lost its elastic floor')
    if member['source']['level'] == 0 and elastic <= F(7,10):
        raise ArithmeticError('the member-zero elastic margin differs')
    result = contract.serial({'member_id': member['member_id'], **member['ids'],
        'same_source_state_reference_clock': True,
        'source_level': member['source']['level'], 'source_spacing': a,
        'source_width': D, 'source_H_norm_upper': R,
        'local_radial_spectrum': local['spectrum'],
        'source_density_floor': floor, 'actual_degree_census': degree,
        'actual_corner_witness_norm_squared': gamma,
        'all_site_witness_norm_squared_lower': gamma,
        'each_incident_bond_norm_squared_coefficient': F(3,16),
        'selected_slope_is_positive': slope_positive,
        'frequency_cut': {'cutoff': lam, 'kernel_certificate': derivative_certificate(),
            'low_filtered_commutator_norm_upper': 15*lam,
            'high_filtered_commutator_norm_lower': witness_norm-15*lam,
            'high_filtered_operator_norm_lower': (witness_norm-15*lam)/D,
            'positive_spectral_weight_coefficient': q,
            'sharper_positive_spectral_weight_coefficient': improved,
            'dimension_factor_used': False,
            'witness_replaces_selected_state': False},
        'original_ideal_positive_transfer_interval': [lam,D],
        'original_ideal_positive_transfer_weight_lower': dyadic(q,N),
        'complete_clock_plan': plan,
        'reported_positive_transfer_interval': [3*lam/4,D+lam/4],
        'original_reported_positive_transfer_weight_lower': dyadic(3*q/4,N),
        'true_transfer_lower_on_good_clock_records': lam/2,
        'complete_probe': {'success_rows': 'A_x/sqrt(V), A_x=sum_a Phi_(x,a)^dagger Phi_(x,a)/2',
            'failure': 'sqrt(I-sum_x A_x^2/V)', 'probe_outcomes': V+1,
            'fine_outcomes': {'coefficient': V+1, 'base': 2,
                              'exponent': 2*plan['clock_qubits_each']},
            'two_energy_registers_retained': True, 'all_addresses_and_failures_retained': True,
            'quantum_reference_retained': True, 'initial_source_marginal_preserved': True,
            'initial_reference_joint_state_claimed_unchanged': False,
            'original_large_histogram_executed': False},
        'original_elastic_weight_lower': elastic,
        'uniform_elastic_weight_lower': F(3481,160000),
        'member_zero_simple_elastic_strict_lower': F(7,10) if member['source']['level']==0 else None,
        'uniform_total_detection_is_uniform_inelastic_weight': False,
        'positive_cutoff_is_capacity_independent_in_source_clock': lam==F(1,40),
        'upper_endpoint_is_capacity_independent': False,
        'uniform_compact_positive_energy_weight_proved': False,
        'isolated_scalar_pole_or_mass_gap_proved': False,
        'physical_Higgs_field_identified': False,
        'physical_particle_map': [], 'physical_mass_table_GeV': [],
        'next_needed': 'populated compact spectral branch and momentum dependence, a physical scalar/Higgs-phase dictionary, and controlled mass and capacity errors'})
    return {'schema': 'oph.particle_candidate_scalar_spectrum_bound.v1',
            'result': result, 'result_sha256': contract.digest(result)}


def verify_bound(member, value, cutoff='1/40'):
    if value != bind_member(member,cutoff):
        raise ValueError('scalar spectral window, complete clock or physical domain differs')


def certificate():
    paths = provenance.definition_inputs(SOURCE_INPUTS)
    pins = provenance.definition_hashes(paths)
    value = {'schema': 'oph.particle_candidate_scalar_spectrum.v1',
             'members': [bind_member(contract.member(n)) for n in (0,1,2,8)],
             'source_sha256': pins}
    if provenance.definition_hashes(paths) != pins:
        raise ValueError('scalar spectral source inputs changed during replay')
    return value


build = certificate
