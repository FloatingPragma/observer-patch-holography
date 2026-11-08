"""Exact local composite-force and original low-energy probability bounds."""
from fractions import Fraction as F

import sympy as s

from model.source_family_compatibility.local_matter import annihilator
from model.field_symmetry import symmetry
from model.full_fock_field.field import FieldGeometry


def _fraction(value, *, positive=False):
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError('exact rational parameter required')
    result = F(value)
    if result < 0 or (positive and not result):
        raise ValueError('parameter outside its nonnegative domain')
    return result


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('integer outside its admitted domain')
    return value


def upper_admissible(alpha_squared, force_squared, energy, probability):
    """Exact squared test for alpha*sqrt(w)-beta*sqrt(1-w)>=E."""
    a2 = _fraction(alpha_squared, positive=True)
    b2, e, w = map(_fraction, (force_squared, energy, probability))
    if w > 1:
        raise ValueError('probability must be at most one')
    if e*e >= a2:
        return w == 1
    difference = (a2+b2)*w-b2-e*e
    return difference >= 0 and difference*difference >= 4*e*e*b2*(1-w)


def low_window_upper(alpha_squared, force_squared, energy, bits=48):
    """Outward rational enclosure; no conditional-success normalization."""
    a2 = _fraction(alpha_squared, positive=True)
    b2, e = map(_fraction, (force_squared, energy))
    _integer(bits, 1)
    if e*e >= a2:
        return F(1)
    if not e:
        return b2/(a2+b2)
    lower, upper = F(0), F(1)
    for _ in range(bits):
        middle = (lower+upper)/2
        if upper_admissible(a2, b2, e, middle):
            upper = middle
        else:
            lower = middle
    if not upper_admissible(a2, b2, e, upper):
        raise ArithmeticError('outward spectral probability enclosure failed')
    return upper


def local_budget(coupling, contact, degree, energy, bits=48):
    c, u = (_fraction(value, positive=True) for value in (coupling, contact))
    _integer(degree)
    e = _fraction(energy)
    r2 = c*c*degree/4
    upper = low_window_upper(u*u/4, 9*r2, e, bits)
    return {'alpha_squared': u*u/4, 'neighbor_orbital_norm_squared': r2,
            'force_operator_norm_squared': 5*r2,
            'force_retarded_norm_squared_upper': 9*r2,
            'original_signed_low_window_probability_upper': upper,
            'original_signed_outside_probability_lower': 1-upper,
            'each_local_PH_charge_low_window_probability_upper': upper/2,
            'each_local_PH_charge_outside_probability_lower': (1-upper)/2,
            'full_many_body_charge_gap_proved': False,
            'arbitrary_packet_bound_from_local_estimate': False}


def cubic_packet_budget(coupling, contact, energy, bits=48):
    """Cubic-invariant infinite-source force symbol, with its spatial premise."""
    c, u = (_fraction(value, positive=True) for value in (coupling, contact))
    g0_upper, g2_absolute_upper = 27*c*c/2, 3*c*c/4
    force_squared = g0_upper+6*g2_absolute_upper
    upper = low_window_upper(u*u/4, force_squared, energy, bits)
    return {'cubic_ground_symmetry_and_range_two_force_symbol_required': True,
            'local_force_coefficient_upper': g0_upper,
            'second_neighbor_coefficient_absolute_upper': g2_absolute_upper,
            'all_packet_force_squared_upper': force_squared,
            'original_signed_low_window_probability_upper': upper,
            'original_signed_outside_probability_lower': 1-upper,
            'individual_smeared_charge_has_half_weight': False,
            'full_many_body_charge_gap_proved': False}


def time_average_budget(coupling, contact, halfwidth, bits=48):
    """Complete retarded score averaged over [-delta,delta], with unit mass."""
    delta = _fraction(halfwidth, positive=True)
    row = cubic_packet_budget(coupling, contact, 2/delta, bits)
    deficit = row['original_signed_outside_probability_lower']/2
    return {'time_halfwidth': delta, 'physical_energy_cutoff': 2/delta,
            'one_minus_real_time_average_lower': deficit,
            'real_time_average_upper': 1-deficit,
            'complete_signed_response_mass': F(1),
            'filtered_response_excluded': False}


def _comm(a, b):
    return a*b-b*a


def _anti(a, b):
    return a*b+b*a


def star_algebra():
    """Four actual CAR modes after the normalized neighbor-orbital transform."""
    a = tuple(annihilator(2, j) for j in range(4))
    identity = s.eye(16)
    t = sum((a[j].H*a[j+2]+a[j+2].H*a[j] for j in (0, 1)), s.zeros(16))
    b = (2*a[1].H*a[1]-identity)*a[0]
    y = _comm(b, t)
    gram, squared = _anti(y.H, y), y.H*y
    if gram.eigenvals() != {s.Integer(1): 10, s.Integer(5): 4, s.Integer(9): 2}:
        raise ArithmeticError('exact composite-force anticommutator spectrum changed')
    if squared.eigenvals() != {s.Integer(0): 8, s.Integer(1): 4, s.Integer(5): 4}:
        raise ArithmeticError('exact composite-force operator norm changed')
    # The interaction part contributes alpha*a exactly, without a projection.
    contact = sum((a[2*x].H*a[2*x]*a[2*x+1].H*a[2*x+1]
                   -(a[2*x].H*a[2*x]+a[2*x+1].H*a[2*x+1])/2
                   for x in (0, 1)), s.zeros(16))
    if _comm(b, contact) != a[0]/2:
        raise ArithmeticError('centered composite contact equation failed')
    return {'full_Fock_dimension': 16, 'central_and_neighbor_CAR_modes': 4,
            'force_anticommutator_eigenvalues': {'1': 10, '5': 4, '9': 2},
            'force_squared_eigenvalues': {'0': 8, '1': 4, '5': 4},
            'retarded_force_factor': 9, 'operator_force_factor': 5,
            'full_contact_equation_exact': True}


def literal_weyl_controls():
    """Literal three-site Weyl sources, all Pauli-link axes and both spins."""
    rows = []
    identity = s.eye(64)
    for axis in range(3):
        sides = tuple(3 if j == axis else 1 for j in range(3))
        source = symmetry.operators(FieldGeometry(sides, 1., 8.))
        clock = (s.Matrix(source.K.real.toarray().tolist())
                 +s.I*s.Matrix(source.K.imag.toarray().tolist()))/source.scale
        hopping = clock-s.diag(*clock.diagonal())
        for spin in (0, 1):
            a = annihilator(3, 2+spin)
            opposite = annihilator(3, 2+(spin ^ 1))
            b = (2*opposite.H*opposite-identity)*a
            y = _comm(b, hopping)
            if _comm(b, clock)-y != 4*a:
                raise ArithmeticError('literal Weyl composite equation changed')
            g = _anti(y.H, y)
            r2 = s.Rational(1, 2)
            if (g-r2*identity)*(g-5*r2*identity)*(g-9*r2*identity) != s.zeros(64):
                raise ArithmeticError('literal source fails the exact star-force polynomial')
            if g.eigenvals() != {r2: 40, 5*r2: 16, 9*r2: 8}:
                raise ArithmeticError('literal source force bound or sharpness changed')
        rows.append({'axis': axis, 'full_Fock_dimension': 64, 'probe_site_degree': 2,
                     'c': F(1), 'U': F(8), 'neighbor_orbital_norm_squared': F(1, 2),
                     'both_original_flavours_checked': True,
                     'original_Weyl_action_used': True})
    return rows


def response_countercontrol():
    """Sharp projection control; no assertion that this matrix is the source."""
    alpha, beta = s.Integer(4), s.Integer(3)
    clock = s.Matrix([[0, alpha, 0], [alpha, 0, beta], [0, beta, 0]])
    v, u, force = s.eye(3)[:, 0], s.eye(3)[:, 1], beta*s.eye(3)[:, 2]
    zero = s.eye(3)-clock*clock/25
    weight = (v.H*zero*v)[0]
    if (clock*u != alpha*v+force or (v.H*force)[0] != 0
            or zero*zero != zero or weight != s.Rational(9, 25)):
        raise ArithmeticError('orthogonal response saturation control failed')
    if low_window_upper(16, 9, 0) != F(9, 25):
        raise ArithmeticError('exact zero-window bound is not sharp')
    return {'dimension': 3, 'signed_energies': [-5, 0, 5],
            'original_probabilities': [F(8, 25), F(9, 25), F(8, 25)],
            'force_orthogonal_to_original_vector': True,
            'zero_window_bound_attained': True,
            'actual_Hubbard_source_realization_claimed': False,
            'gap_follows_from_response_bounds': False}


def certificate():
    return {'exact_local_force': star_algebra(),
            'literal_Weyl_source_controls': literal_weyl_controls(),
            'default_ratio': F(8),
            'degree_six_zero_scaled_window': local_budget(1, 8, 6, 0),
            'degree_six_fixed_window_at_c_1024': local_budget(1024, 8192, 6, 1),
            'cubic_ground_all_packet_zero_scaled_window': cubic_packet_budget(1, 8, 0),
            'cubic_ground_time_average_at_c_1024': time_average_budget(1024, 8192, 1),
            'default_ballistic_time_average_deficit_liminf': F(4, 17),
            'unfiltered_unit_mass_regular_limit_excluded_on_stated_cubic_family': True,
            'sharp_orthogonal_projection_countercontrol': response_countercontrol(),
            'physical_particle_species': [], 'physical_particle_masses': []}
