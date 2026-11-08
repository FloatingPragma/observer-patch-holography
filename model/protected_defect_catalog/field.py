"""Source charge action and odd preparation in the inherited field."""
from fractions import Fraction as F

from model.spin_canonical_field.algebra import spin, comm
import sympy as s


def configuration_charge():
    """Audit the finite identities used by the configuration-subnet proof."""
    x, y, z = spin()
    identity = s.eye(2)
    raising, lowering = z-s.I*y, z+s.I*y
    if comm(x, raising) != raising or comm(x, lowering) != -lowering:
        raise ArithmeticError('transverse charge orientation changed')
    rotation = (2*identity+2*s.I*x)/s.sqrt(5)
    if (rotation*z*rotation.H != s.Rational(3,5)*z+s.Rational(4,5)*y
            or rotation*y*rotation.H != s.Rational(3,5)*y-s.Rational(4,5)*z):
        raise ArithmeticError('native rotation fails the configuration-field action')
    left = [s.kronecker_product(a,identity) for a in (x,y,z)]
    right = [s.kronecker_product(identity,a) for a in (x,y,z)]
    h = sum((a*b for a,b in zip(left,right)), s.zeros(4))-s.eye(4)/4
    b = s.Rational(1,4)
    h -= b*(left[0]-right[0])
    swap = s.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    conjugation = swap*s.kronecker_product(2*y,2*y)
    q2, q3 = left[1]-right[1], left[2]-right[2]
    total = left[0]+right[0]
    if (conjugation.H*conjugation != s.eye(4)
            or conjugation*h*conjugation != h
            or conjugation*total*conjugation != -total
            or conjugation*q2*conjugation != -q2
            or conjugation*q3*conjugation != q3):
        raise ArithmeticError('odd translation and source rotation lose charge conjugation')
    return {'source_transverse_charge_orientation_checked':True,
            'source_combined_translation_rotation_checked':True,
            'field_symmetry_attached':'U1 on the transverse configuration-time-history subnet',
            'configuration_fields':['Q3(F)','Q2(F)'],
            'charged_field_vectors':{'plus':'(Q3-i Q2)Omega','minus':'(Q3+i Q2)Omega'},
            'equal_time_configuration_groups_strongly_commute':True,
            'joint_energy_momentum_measure_each':'2 times the actual Q3(F)Omega measure',
            'charge_conjugation_preserves_source_state_and_scaled_clock':True,
            'positive_energy_and_gapless_annulus_weights_inherited':True,
            'U1_action_on_entire_original_Q_P_history_net_proved':False,
            'second_canonical_momentum_pair_attached':False,
            'physical_antiparticle_or_Lorentz_spin_identified':False,
            'finite_cube_local_branches_identified_with_continuum_poles':False}


def configuration_commutator_bound(scale, absolute_test_product_sum):
    """Exact ||[Q2_p(F),Q3_p(G)]|| bound, before a Riemann limit."""
    if any(type(x) not in (int,F) for x in (scale,absolute_test_product_sum)):
        raise TypeError('exact rational parameters required')
    scale, total = map(F,(scale,absolute_test_product_sum))
    if scale <= 0 or total < 0:
        raise ValueError('positive scale and nonnegative test sum required')
    return scale**4*total/2


def preparation_bounds(m0=F(1,1000), regularizer=F(1), resolution=F(1)):
    if any(type(x) not in (int,F) for x in (m0,regularizer,resolution)):
        raise TypeError('exact rational parameters required')
    m0, regularizer, resolution = map(F,(m0,regularizer,resolution))
    if not 0 < m0 <= F(1,2) or regularizer <= 0 or resolution <= 0:
        raise ValueError('positive field order, regularizer and resolution required')
    n0 = 25*m0*m0/2304
    h0 = 12/(m0*m0)
    limiting = n0**3/(h0+regularizer*n0)**2
    uniform = (n0/2)**3/(2*h0+regularizer*n0/2)**2
    return {'order_lower_bound':str(m0), 'odd_sine_amplitude':str(2/m0),
            'conjugate_test_amplitude':str(m0/4),
            'liminf_original_odd_preparation_probability':str(n0),
            'limsup_unnormalized_scaled_energy':str(h0),
            'regularizer_in_the_scaled_clock':str(regularizer),
            'liminf_odd_filtered_original_probability':str(regularizer**2*limiting),
            'sufficiently_small_scale_filtered_probability_lower':str(regularizer**2*uniform),
            'effective_scale_threshold_supplied':False,
            'all_zero_energy_atoms_excluded_for_fixed_sine_amplitude':False,
            'positive_energy_odd_Stone_field_vectors_inherited':True,
            'positive_energy_sine_channel_for_some_rational_amplitude_proved':True,
            'field_symmetry_attached':'Z2 native spin parity',
            'second_canonical_pair_attached':False,
            'physical_mass_shell_or_species_identified':False}
