"""Actual eta multiplets, source-current comparison and complete-Fock leakage."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from ..common_source import native_matter_attachment as target
from ..radial_interacting_core.core import core_bounds, encode


def low_representation():
    """The basis is C_minus,G1,G2,C_plus; C_plus=eta_plus C_minus."""
    z = s.diag(-s.Rational(1, 2), 0, 0, s.Rational(1, 2))
    plus = s.zeros(4); plus[3, 0] = 1
    x, y = (plus+plus.H)/2, (plus-plus.H)/(2*s.I)
    casimir = x*x+y*y+z*z
    if (plus*plus.H-plus.H*plus != 2*z or z*plus-plus*z != plus
            or casimir != s.diag(s.Rational(3, 4), 0, 0, s.Rational(3, 4))):
        raise ArithmeticError('the actual isolated eta representation is inconsistent')
    row = core_bounds()
    if row['actual_charged_gap_upper'] >= row['other_charged_excitations_threshold_lower']:
        raise ArithmeticError('the complete four-state isolation is required')
    return {'basis': ['C_minus', 'G1', 'G2', 'C_plus=eta_plus C_minus'],
            'particle_numbers': [row['sites']-1, row['sites'], row['sites'], row['sites']+1],
            'eta_z': [str(z[i, i]) for i in range(4)],
            'eta_plus': [[str(v) for v in plus.row(i)] for i in range(4)],
            'eta_Casimir': [str(casimir[i, i]) for i in range(4)],
            'charged_ladder_norm_squared': '1',
            'ground_eta_singlet_multiplicity': 2,
            'charged_eta_doublet_multiplicity': 1,
            'Hamiltonian_and_eta_intertwining_exact': True,
            'eta_current_leakage_from_complete_low_quartet': '0',
            'actual_charge_addition_removal_amplitude_identity': 'A_plus(x)=D_x i sigma_2 A_minus(x)',
            'amplitude_convention': 'A_plus(x,s,a)=<C_plus|c_xs_dagger|G_a>; A_minus(x,s,a)=<C_minus|c_xs|G_a>',
            'amplitude_identity_uses_actual_interacting_states': True,
            'no_same_energy_states_in_other_number_sectors_required': True,
            'source_core_bounds': encode(row),
            'charge_doublet_is_a_selected_Standard_Model_multiplet': False}


def leakage_bounds():
    row = core_bounds()
    b_squared = F(3, 2)*row['velocity']**2/row['spacing']**2
    alpha_upper = 1/(1+b_squared/(row['mass']+row['hopping_bound'])**2)
    angle = row['projection_angle_bound']
    charged_leak_squared = (1-alpha_upper**2)/2-angle
    alpha_lower = row['origin_weight_lower']
    one_body_leak_squared = alpha_lower*(1-alpha_upper)
    if alpha_upper != F(2904, 2905) or charged_leak_squared <= F(1, 4096):
        raise ArithmeticError('the full-source leakage witness lost its strict lower bound')
    return encode({'origin_zero_weight_lower': alpha_lower,
                   'origin_zero_weight_upper': alpha_upper,
                   'interacting_projector_angle': angle,
                   'each_actual_charged_center_Pauli_leakage_squared_lower': charged_leak_squared,
                   'each_actual_charged_center_Pauli_leakage_norm_strict_lower': F(1, 64),
                   'one_body_zero_frame_center_Pauli_leakage_squared_lower': one_body_leak_squared,
                   'one_body_compressed_Pauli': 'alpha sigma_j',
                   'one_body_compressed_bracket_defect': '2i alpha(1-alpha) epsilon_jkl sigma_l',
                   'actual_charged_full_low_projection_of_center_Pauli': '0',
                   'actual_charged_leakage_is_outside_entire_low_quartet': True,
                   'selected_zero_orbital_projection_used_for_interacting_bound': False})


@lru_cache(maxsize=1)
def native_current_comparison():
    """Retain the actual twelve-port current and its declared exterior lift."""
    native = target.canonical_currents()
    r = target.representation()
    selected = list(r['selected'])
    lifted = [target.exterior_derivative(target.spin_lift(g)).extract(selected, selected)
              for g in native]
    color_rank = s.Matrix.hstack(*(s.Matrix(g).reshape(225, 1) for g in lifted[:8])).rank()
    weak_rank = s.Matrix.hstack(*(s.Matrix(g).reshape(225, 1) for g in lifted[8:11])).rank()
    if color_rank != 8 or weak_rank != 3 or len(native) != 12:
        raise ArithmeticError('the complete declared native current target changed')
    if any(a*b-b*a != s.zeros(15) for a in lifted[:8] for b in lifted[8:]):
        raise ArithmeticError('the target color action fails its commuting-factor identity')
    return {'native_ports': 12, 'declared_exterior_matter_dimension': 15,
            'native_color_current_rank': color_rank, 'native_weak_current_rank': weak_rank,
            'target_multiplets_from_source_Casimirs': target.decoded_multiplets(),
            'every_homomorphism_su3_to_source_local_eta_su2_is_zero': True,
            'weak_eta_identification_leaves_no_nonconstant_commuting_local_hypercharge': True,
            'faithful_full_native_current_action_in_source_local_conserved_algebra': False,
            'auxiliary_five_mode_or_fifteen_label_register_added_to_radial_source': False,
            'broken_gauge_currents_with_dynamical_links_excluded_by_this_theorem': False,
            'nonlocal_or_emergent_current_algebras_excluded_by_this_theorem': False}


def certificate():
    return {'actual_interacting_low_representation': low_representation(),
            'strict_complete_Fock_leakage': leakage_bounds(),
            'same_native_current_comparison': native_current_comparison()}
