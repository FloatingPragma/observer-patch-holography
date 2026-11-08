"""Quantitative odd-sector preparation in the inherited canonical field."""
from fractions import Fraction as F


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
            'full_U1_or_second_canonical_pair_attached':False,
            'physical_mass_shell_or_species_identified':False}
