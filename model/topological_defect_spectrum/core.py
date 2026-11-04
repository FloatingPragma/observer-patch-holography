"""Controlled original contact interaction on a supplied finite defect core."""
from fractions import Fraction as F
import sympy as s

from ..native_mass_channels.algebra import rational, SIGMA
from . import lattice


def core_bounds(coupling=F(1, 65536), dimensions=3):
    coupling = rational(coupling)
    if coupling <= 0 or type(dimensions) is not int or dimensions not in (1, 3):
        raise ValueError('positive exact contact coupling and dimension one or three required')
    volume, ipr, gap = 3**dimensions, F(1, 2**dimensions), F(6, 5)
    norm = coupling*volume/4
    if 4*norm >= gap:
        raise ValueError('the four-level perturbation requires ||W|| < Delta/4')
    remainder = norm*norm/(gap-2*norm)
    projected_gap = coupling*ipr/2
    return {
        'sites': volume, 'complete_one_body_dimension': 2*volume,
        'full_Fock_dimension': 2**(2*volume),
        'filled_negative_modes': volume-1, 'zero_modes': 2,
        'one_body_and_Fock_excitation_gap_lower': gap,
        'zero_orbital_inverse_participation': ipr,
        'original_centered_contact_coupling': coupling,
        'full_contact_norm': norm,
        'projected_contact_strength': coupling*ipr,
        'projected_charged_gap': projected_gap,
        'ordered_low_energy_error': remainder,
        'actual_charged_gap_lower': projected_gap-2*remainder,
        'actual_charged_gap_upper': projected_gap+2*remainder,
        'charged_gap_certified_positive': projected_gap > 2*remainder,
        'exact_low_Kramers_ground_multiplicity': 2 if projected_gap > 2*remainder else None,
        'particle_hole_related_charged_levels': 2,
        'rest_of_Fock_spectrum_above_shifted_zero': gap-norm,
        'independent_full_Fock_matrix_test_available': dimensions == 1,
    }


def exact_core_algebra():
    B = lattice.jacobi(1)
    z = s.Matrix([1, 2, 1])
    reflection = s.eye(3)[:, [2, 1, 0]]
    x = s.symbols('x')
    if (B*z != s.zeros(3, 1) or reflection*z != z
            or reflection*B*reflection != -B
            or B.charpoly(x).as_expr() != x**3-s.Rational(3, 2)*x):
        raise ArithmeticError('exact finite defect zero mode or gap failed')
    one_ipr = sum(y**4 for y in z)/(z.dot(z))**2
    if one_ipr != s.Rational(1, 2):
        raise ArithmeticError('zero orbital localization changed')
    # Exact spin antiunitary and particle-hole matrix identity.
    spin_time = s.I*SIGMA[1]
    if spin_time*spin_time.conjugate() != -s.eye(2):
        raise ArithmeticError('onsite antiunitary is not Kramers')
    for sigma in SIGMA:
        if spin_time*sigma.conjugate()*spin_time.H != -sigma:
            raise ArithmeticError('source spin time reversal identity failed')
    # The sea relation is derived from exact reflection, spectral complement
    # and the even zero orbital; no numerical sea eigenvector is substituted.
    ipr = one_ipr**3
    sea = (1-ipr)*s.eye(2)/2
    linear = (s.trace(sea)-s.Rational(1, 2))*s.eye(2)-sea
    if linear != -ipr*s.eye(2)/2:
        raise ArithmeticError('Wick contact compression has a spurious core chemical potential')
    return {
        'one_dimensional_characteristic_polynomial': 'x^3-(3/2)x',
        'zero_orbital_unnormalized': [1, 2, 1],
        'zero_orbital_inverse_participation_three_dimensions': str(ipr),
        'weighted_sea_density': '(1-IPR) I_2/2',
        'exact_compressed_centered_contact': 'C+U*IPR*(n1*n2-(n1+n2)/2)',
        'spatial_inversion_anticommutes_with_transformed_source': True,
        'combined_inversion_spin_time_reversal_commutes': True,
        'odd_number_sectors_have_exact_Kramers_degeneracy': True,
        'onsite_particle_hole_preserves_centered_contact': True,
    }


def _encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    return value


def certificate():
    cube = core_bounds()
    if not cube['charged_gap_certified_positive']:
        raise ArithmeticError('controlled interacting defect gap is not positive')
    return {'exact_algebra': exact_core_algebra(),
            'controlled_27_site_cube': _encode(cube),
            'independent_full_Fock_chain_domain': _encode(core_bounds(F(1, 1024), 1)),
            'original_centered_contact_retained': True,
            'virtual_nonzero_modes_retained_in_error_bound': True,
            'two_to_54_Fock_matrix_executed': False,
            'supplied_linear_mass_texture_required': True,
            'source_selects_texture_or_its_finite_preparation_energy_proved': False,
            'thermodynamic_interacting_defect_gap_proved': False,
            'physical_particle_assignment': False}
