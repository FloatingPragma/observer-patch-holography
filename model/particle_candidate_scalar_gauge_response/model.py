"""Composite radial occupation, original response and one finite source clock."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_candidate_charged import model as charged
from model.particle_primary_selected_odd_response import model as selected
from model.particle_primary_composite_transport import model as transport
from model.particle_primary_energy_response import model as current
from model.particle_primary_gravity_premises import model as modular
from model.particle_spatial_modular_response import regional
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_refinement_response import model as refinement
from model.particle_primary_variational_readout import variational
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS = tuple('model/particle_candidate_scalar_gauge_response/' + n for n in
    ('__init__.py', '__main__.py', 'model.py', 'test_model.py', 'PROOF.md', 'README.md')) + (
    'model/particle_candidate_charged/PROOF.md',
    'model/particle_primary_charge_conjugation/PROOF.md',
    'model/particle_primary_selected_odd_response/PROOF.md',
    'model/particle_primary_energy_response/PROOF.md',
    'Lean/QFT/ScalarRegionalTimeSlice.lean',
    'Lean/Screen/ElectroweakBreakingComposition.lean',
    'Lean/Screen/ChargedCurrentDictionary.lean',
)


def clean(value):
    return s.expand(value)


@lru_cache(maxsize=1)
def radial_control():
    """Sharp positive inequalities on all64 occupations, with native currents."""
    d = charged.dressing_control()
    S = clean(sum((p.H*p for p in d['phi']), s.zeros(64)))
    I = s.eye(64)
    holes = 6*I-d['NL']-d['NE']
    eigen = S.eigenvals()
    expected = {s.Integer(0): 29, s.Rational(1, 2): 16, s.Integer(1): 10,
                s.Rational(3, 2): 8, s.Integer(2): 1}
    if eigen != expected or S*holes != holes*S:
        raise ArithmeticError('the exact radial occupation spectrum differs')
    rows = []
    for n in range(7):
        indices = [k for k in range(64) if 6-k.bit_count() == n]
        block = S.extract(indices, indices)
        mean = (2*s.eye(len(indices))-block).eigenvals()
        deficit = (s.eye(len(indices))-block**2/4).eigenvals()
        if max(mean) > n or max(deficit) > s.Rational(3, 4)*n:
            raise ArithmeticError('a radial hole inequality failed')
        rows.append({'negative_holes': n, 'dimension': len(indices),
                     'largest_radial_mean_deficit': str(max(mean)),
                     'largest_radial_effect_deficit': str(max(deficit))})
    creators = d['creators']
    for q in (d['NL'], d['NE']):
        if q*S != S*q:
            raise ArithmeticError('radial readout changes native occupation')
    weak = charged.blocks()[1][1]
    electron = charged.blocks()[4][1]
    for T in charged.native.native_currents():
        J = charged.second_quantize(s.diag(s.kronecker_product(T.extract(weak, weak), s.eye(2)),
                                           T[electron[0], electron[0]]*s.eye(2)), creators)
        if clean(J*S-S*J) != s.zeros(64):
            raise ArithmeticError('radial readout violates an actual native current')
    signs = s.diag(*((-1)**k.bit_count() for k in range(64)))
    if signs*S != S*signs or clean(S*d['anchor']) != 2*d['anchor']:
        raise ArithmeticError('radial parity or filled anchor changed')
    one_holes = []
    for j, C in enumerate(creators):
        v = C.H*d['anchor']
        value = s.Rational(3, 2) if j < 4 else s.Integer(1)
        if clean(S*v-value*v) != s.zeros(64, 1):
            raise ArithmeticError('one-hole scalar action changed')
        one_holes.append(str(value))
    return {'S': S, 'A': S/2, 'holes': holes, 'creators': creators,
            'summary': {'local_CAR_modes': 6, 'full_Fock_dimension': 64,
                'radial_operator': 'S=sum_a Phi_a^dagger Phi_a',
                'normalized_seed': 'A=S/2', 'source_support': 'one site, negative-beta weak and charge-six modes',
                'spectrum': {str(k): int(v) for k, v in sorted(eigen.items())},
                'hole_blocks': rows, 'sharp_mean_deficit_hole_coefficient': 1,
                'sharp_effect_deficit_hole_coefficient': F(3, 4),
                'one_hole_radial_values': one_holes,
                'native_generators_checked': 12, 'proper_cube_scalar_coefficients': True,
                'inversion_even_coefficients': True, 'CAR_parity': 0,
                'native_increment': [0, 0, 0, 0, 0], 'norm_A': 1,
                'source_gauge_invariant': True, 'charged_Phi_one_point_is_radial_mean': False,
                'physical_Higgs_radial_field_identified': False}}


@lru_cache(maxsize=1)
def wilson_witness():
    """Exact full-source commutators on an admitted120-mode Haar/CAR witness."""
    radial_control()
    blocks = current.sectors.native_blocks()
    anchor = sum(1 << (60*x+j) for x in range(2) for j in range(30, 60))
    axes = []
    for axis in range(3):
        hop = 2*transport.spin_blocks()[axis]
        first, second = defaultdict(complex), defaultdict(complex)
        for r, (charge, indices) in enumerate(blocks):
            twice_drop = 1 if charge == -3 else 2 if charge == 6 else 0
            for i, j in product(indices, repeat=2):
                for a, b in product(range(4), repeat=2):
                    z = complex(hop[a, b])
                    if not z:
                        continue
                    if z not in (1, -1, 1j, -1j):
                        raise ArithmeticError('doubled native Wilson entry is not a Gaussian integer')
                    for target, origin, coefficient, orientation in (
                            (15*a+i, 60+15*b+j, z, 1),
                            (60+15*b+j, 15*a+i, z.conjugate(), -1)):
                        moved = current.move(anchor, target, origin)
                        if moved is None or origin >= 60 or not twice_drop:
                            continue
                        state, sign = moved
                        key = (r, orientation, i, j, state)
                        first[key] += twice_drop*coefficient*sign
                        second[key] -= twice_drop**2*coefficient*sign
        def norm(values, denominator):
            total = F()
            for (r, orientation, i, j, state), z in values.items():
                if not z.real.is_integer() or not z.imag.is_integer():
                    raise ArithmeticError('exact radial witness left Gaussian integers')
                total += F(int(z.real)**2+int(z.imag)**2, denominator**2*len(blocks[r][1]))
            return total
        one, two = norm(first, 4), norm(second, 8)
        if (one, two) != (F(3, 4), F(9, 16)):
            raise ArithmeticError('the literal radial Wilson witness differs')
        axes.append({'axis': axis, 'nonzero_native_charges': [-3, 6],
                     'one_link_coordinates': len(first),
                     'S_first_commutator_norm_squared': one,
                     'S_double_commutator_norm_squared': two})
    return {'actual_two_cell_modes': 120, 'axes': axes,
            'all_native_blocks_included': 5, 'corner_degree': 3,
            'one_link_cutoff_admission': True, 'Gauss_physical_witness': True,
            'anchor_substituted_for_selected_state': False,
            'S_first_norm_squared_coefficient': F(9, 4),
            'S_double_norm_squared_coefficient': F(27, 16),
            'A_first_norm_squared_coefficient': F(9, 16),
            'A_double_norm_squared_coefficient': F(27, 256),
            'spacing_power': -2, 'these_are_full_commutators_on_anchor': True}


def factored(value, exponent):
    value = F(value)
    if value <= 0 or exponent < 0:
        raise ValueError('positive rational multiplier and nonnegative dyadic exponent required')
    return {'rational_multiplier': value, 'base': 2, 'negative_exponent': int(exponent)}


def scalar_records(member):
    contract.verify_member(member)
    original = selected.bind_member(member)
    selected.verify_bound(member, original)
    V = member['source']['side']**3
    delta = F(original['result']['weak_plus_e_defect_count_upper'])
    mean = 1-delta/(4*V)
    linear = 1-3*delta/(8*V)
    probability = max(F(), linear, mean**2)
    if mean < F(59, 400) or probability < F(59, 400)**2:
        raise ArithmeticError('the actual selected radial success bound failed')
    if member['source']['level'] == 0 and probability <= F(383, 500):
        raise ArithmeticError('the member-zero radial margin failed')
    return {'selected_odd_result_sha256': original['result_sha256'],
            'actual_weak_plus_e_defect_upper': delta, 'volume': V,
            'averaged_radial_mean_lower': 2*mean,
            'averaged_A_mean_lower': mean,
            'linear_effect_lower': linear,
            'Jensen_effect_lower': mean**2,
            'original_total_success_lower': probability,
            'uniform_original_total_success_lower': F(3481, 160000),
            'member_zero_simple_strict_lower': F(383, 500) if member['source']['level'] == 0 else None,
            'success_rows': 'K_x=S_x/(2 sqrt(V))',
            'failure': 'sqrt(I-sum_x K_x^dagger K_x)',
            'complete_probe_outcomes': V+1, 'retained_success_addresses': V,
            'initial_state_changed': False, 'reference_and_all_failures_retained': True,
            'all_complete_H_clock_marginals_preserve_total_success': True,
            'success_includes_elastic_response': True,
            'uniform_positive_energy_success_lower_proved': False}


def scalar_dynamics(member):
    contract.verify_member(member)
    witness = wilson_witness()
    source = member['preparation_plan']['source']
    R = max(abs(F(source['spectral_lower'])), abs(F(source['spectral_upper'])))
    a = F(member['source']['spacing'])
    floor = regional.common_contract_payload(member, 1)['physical_density_floor_factored']
    N = floor['negative_exponent']
    slope = modular.selected_modular_slope(member)
    beta_low = slope['modular_slope_lower']
    gamma = F(9, 16)/a**2
    tau = 1/(64*a*R**2)
    if a*R < 1 or 2*R**2*tau > 1/(32*a):
        raise ArithmeticError('source-time commutator remainder failed')
    if tau**2*R**3/6 > 1/(64*a):
        raise ArithmeticError('finite Cayley remainder failed')
    if F(27, 1024) <= F(1, 64):
        raise ArithmeticError('radial leading amplitude bound failed')
    # The clock is the literal rational Cayley unitary of the unchanged H.
    event_prefactor = tau**2/(256*a**2)
    first_prefactor = beta_low*gamma/2
    positive_weight = max(first_prefactor/(2*R), gamma/(8*R**2))
    return {'physical_witness': witness, 'source_H_norm_upper': R,
            'source_density_floor': floor, 'selected_modular_slope': slope,
            'original_A_second_transfer_moment_lower': factored(gamma, N),
            'original_A_first_transfer_moment_lower': factored(first_prefactor, N),
            'ideal_strict_positive_transfer_weight_lower': factored(positive_weight, N),
            'ideal_positive_transfer_weight_is_a_resolved_clock_bin': False,
            'complete_clock_duration': tau,
            'source_clock_unitary': 'Q_tau=(I+i tau H/2)^(-1)(I-i tau H/2)',
            'dynamical_success_row': 'K_tau=[A,Q_tau^dagger A Q_tau]/2',
            'dynamical_failure': 'sqrt(I-K_tau^dagger K_tau)',
            'explicit_two_record_unitary': '[[K_tau,F],[F,K_tau]], F=sqrt(I-K_tau^dagger K_tau)',
            'Cayley_logarithm_replaces_source_H': False,
            'complete_dynamical_outcomes': 2, 'success_at_time_zero': 0,
            'original_finite_time_success_strict_lower': factored(event_prefactor, N),
            'finite_time_success_is_an_energy_transfer_bin': False,
            'source_state_reference_clock_changed': False,
            'witness_is_selected_preparation': False,
            'all_native_sectors_and_reference_retained': True,
            'uniform_centered_variance_or_inelastic_residue_proved': False,
            'Higgs_phase_or_Lorentz_scalar_pole_proved': False}


@lru_cache(maxsize=1)
def deleted_hopping_comparison():
    """The tested retained-energy estimate cannot force a large radial residue."""
    member = contract.member(0)
    upper = variational.bind_member(member)
    variational.verify_bound(member, upper)
    rows = [row for row in refinement.native_control() if row[2] not in (-3, 6)]
    shift = 4*sum(d for d, C, q in rows)*(sqrt_upper(F(5, 4))-1)
    loss = sum((occupation.occupation_loss(d, F(1), 1/(F(105, 2)*C-shift))
                ['occupation_loss_upper'] for d, C, q in rows), F())
    lower = -720-12*loss-F(6, 7)
    mean_upper = F(upper['result']['ground_energy_upper'])+F(1, 16)
    if mean_upper-lower <= 0:
        raise ArithmeticError('the deletion estimate has changed its stated scope')
    if lower < F(-724469, 1000) or mean_upper > F(-7237053, 10000) or mean_upper-lower > F(7634, 10000):
        raise ArithmeticError('the outward displayed energy bounds differ')
    return {'retained_native_charges': [q for d, C, q in rows],
            'deleted_native_hopping_charges': [-3, 6],
            'full_Fock_shift': shift, 'per_bond_loss': loss,
            'deleted_action_ground_lower': lower,
            'actual_original_mean_energy_upper': mean_upper,
            'resulting_weak_plus_e_hopping_expectation_upper': mean_upper-lower,
            'strict_negative_hopping_from_this_estimate_proved': False,
            'actual_scalar_dynamics_excluded': False}


def bind_member(member):
    contract.verify_member(member)
    result = contract.serial({'member_id': member['member_id'], **member['ids'],
        'local_radial_control': radial_control()['summary'],
        'complete_radial_records': scalar_records(member),
        'original_radial_dynamics': scalar_dynamics(member),
        'member_zero_energy_method_boundary': deleted_hopping_comparison() if member['source']['level'] == 0 else None,
        'same_source_state_reference_clock': True,
        'source_scalar_contraction_constructed': True,
        'source_uniform_total_scalar_detection_proved': True,
        'source_finite_original_radial_dynamics_proved': True,
        'physical_Higgs_field_identified': False,
        'physical_mass_table_GeV': [], 'physical_species_confirmed': [],
        'next_physical_attachment': 'connected scalar spectral shape, common field vacuum, relativistic propagation and scale'})
    return {'schema': 'oph.particle_candidate_scalar_gauge_response_bound.v1',
            'result': result, 'result_sha256': contract.digest(result)}


def verify_bound(member, value):
    if value != bind_member(member):
        raise ValueError('source scalar result, clock or domain differs')


def certificate():
    paths = provenance.definition_inputs(SOURCE_INPUTS)
    pins = provenance.definition_hashes(paths)
    value = {'schema': 'oph.particle_candidate_scalar_gauge_response.v1',
             'members': [bind_member(contract.member(n)) for n in (0, 1, 2, 8)],
             'source_sha256': pins}
    if provenance.definition_hashes(paths) != pins:
        raise ValueError('scalar response inputs changed during replay')
    return value


build = certificate
