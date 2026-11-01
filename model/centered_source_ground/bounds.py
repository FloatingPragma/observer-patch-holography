"""Rational centered-ground estimates and source-selected clock coefficients."""
from fractions import Fraction as F


def pull_through(ratio):
    if type(ratio) not in (int,F) or ratio<0:
        raise ValueError('nonnegative exact ratio U/c required')
    return min(F(2),F(11,8)*ratio*ratio)


def original_pull_through(ratio):
    if type(ratio) not in (int,F) or ratio<0:
        raise ValueError('nonnegative exact ratio U/c required')
    return min(F(2),F(11,4)*ratio*ratio)


def certificate():
    return {'clock_coefficient':'U/2',
            'centered_density_per_site':'1',
            'centered_density_per_spin':'1/2',
            'kinetic_occupation_defect_upper':'min(2,11(U/c)^2/8)',
            'original_uncentered_occupation_defect_upper':'min(2,11(U/c)^2/4)',
            'inverse_sine_square_integral_upper':'pi sqrt(3)/2 < 11/4',
            'centered_first_spectral_moment':'h',
            'centered_second_spectral_moment':'h^2+U^2/4',
            'continuum_first_spectral_moment':'h_W',
            'continuum_second_moment_bound':'||h_W f||^2+(U0^2/4)||f||^2',
            'occupation_examples':[{ 'U_over_c':str(F(8,2**power)),
                                    'D_upper':str(pull_through(F(8,2**power)))}
                                   for power in (8,16,24)],
            'centered_thermodynamic_ground_constructed':True,
            'old_ground_reused_without_check':False,
            'interacting_Lorentz_continuum_proved':False,
            'physical_mass_or_species_identified':False}
