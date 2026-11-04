"""Exact odd-box zero modes and complete-Fock contact bounds."""
from fractions import Fraction as F
from itertools import product

import sympy as s

from ..native_mass_channels.algebra import rational, SIGMA


def parameters(length, dimensions):
    if (type(length) is not int or not 3 <= length <= 127 or length % 2 != 1
            or type(dimensions) is not int or dimensions not in (1, 3)):
        raise ValueError('odd length 3..127 and dimension one or three required')
    if length**dimensions > 4096:
        raise ValueError('explicit source geometry permits at most 4096 sites')


def core_bounds(coupling=F(1, 65536), length=3, speed=F(1), dimensions=3):
    parameters(length, dimensions)
    coupling, speed = rational(coupling), rational(speed)
    if coupling <= 0 or speed <= 0:
        raise ValueError('positive exact source speed and contact required')
    volume = length**dimensions
    ipr = F(2, length+1)**dimensions
    # sin(pi/4)=1/sqrt(2)>2/3. Other lengths use sin(x)>=2x/pi.
    gap = speed*(F(2, 3) if length == 3 else F(2, length+1))
    norm = coupling*volume/4
    if 4*norm >= gap:
        raise ValueError('complete-Fock perturbation requires ||W|| < Delta/4')
    epsilon = norm*norm/(gap-2*norm)
    projected_gap = coupling*ipr/2
    return {
        'length': length, 'dimensions': dimensions, 'sites': volume,
        'source_speed': speed,
        'complete_one_body_dimension': 2*volume,
        'full_Fock_dimension': 2**(2*volume),
        'filled_negative_modes': volume-1, 'zero_modes': 2,
        'one_body_and_Fock_excitation_gap_lower': gap,
        'zero_orbital_inverse_participation': ipr,
        'original_centered_contact_coupling': coupling,
        'full_contact_norm': norm,
        'projected_contact_strength': coupling*ipr,
        'projected_charged_gap': projected_gap,
        'ordered_low_energy_error': epsilon,
        'projection_angle_bound': norm/(gap-2*norm),
        'actual_charged_gap_lower': projected_gap-2*epsilon,
        'actual_charged_gap_upper': projected_gap+2*epsilon,
        'charged_gap_certified_positive': projected_gap > 2*epsilon,
        'exact_low_Kramers_ground_multiplicity': 2 if projected_gap > 2*epsilon else None,
        'ground_particle_number': volume if projected_gap > 2*epsilon else None,
        'charged_particle_numbers': [volume-1, volume+1],
        'particle_hole_related_charged_levels': 2,
        'rest_of_Fock_spectrum_above_shifted_zero': gap-norm,
        'supplied_mass_texture_required': False,
        'independent_full_Fock_matrix_test_available': dimensions == 1,
    }


def chain(length=3, speed=F(1)):
    parameters(length, 1)
    speed = rational(speed)
    if speed <= 0:
        raise ValueError('positive exact source speed required')
    result = s.zeros(length)
    for j in range(length-1):
        result[j, j+1] = s.I*s.Rational(speed.numerator, 2*speed.denominator)
        result[j+1, j] = -result[j, j+1]
    return s.ImmutableMatrix(result)


def zero_vector(length=3):
    parameters(length, 1)
    return s.ImmutableMatrix([1 if j % 2 == 0 else 0 for j in range(length)])


def exact_source_anchor():
    """Check every original c=1 open-cube hopping entry exactly."""
    import numpy as np
    from scipy import sparse
    from ..full_fock_field.field import FieldGeometry

    geometry = FieldGeometry((3, 3, 3), speed=1., interaction=1/65536)
    raw = geometry.single_h()
    if (type(raw) is not sparse.csr_matrix or raw.shape != (54, 54)
            or type(raw.data) is not np.ndarray or raw.data.dtype != np.dtype('complex128')
            or not np.isfinite(raw.data).all()
            or raw.indptr[0] != 0 or raw.indptr[-1] != len(raw.data)):
        raise ValueError('the declared canonical source CSR matrix is required')
    actual = s.zeros(54)
    for row in range(54):
        indices = raw.indices[int(raw.indptr[row]):int(raw.indptr[row+1])]
        if any(int(b) <= int(a) for a, b in zip(indices, indices[1:])):
            raise ValueError('strictly ordered canonical source CSR rows required')
        for index in range(int(raw.indptr[row]), int(raw.indptr[row+1])):
            value = raw.data[index]
            actual[row, int(raw.indices[index])] = (
                s.Rational(float(value.real))+s.I*s.Rational(float(value.imag)))
    sites = tuple(product(range(3), repeat=3))
    expected = s.zeros(54)
    for x, position in enumerate(sites):
        for axis in range(3):
            neighbor = tuple(p+(axis == j) for j, p in enumerate(position))
            if neighbor not in sites:
                continue
            y = sites.index(neighbor)
            for a in range(2):
                for b in range(2):
                    expected[2*x+a, 2*y+b] = s.I*SIGMA[axis][a, b]/2
                    expected[2*y+b, 2*x+a] = s.conjugate(expected[2*x+a, 2*y+b])
    if actual != expected or actual != actual.H or geometry.sites != sites:
        raise ArithmeticError('unmodified original Weyl links changed')
    one = chain()
    separated = sum((s.kronecker_product(
        *(one if j == axis else s.eye(3) for j in range(3)), SIGMA[axis])
        for axis in range(3)), s.zeros(54))
    if actual != separated:
        raise ArithmeticError('the source does not separate into the stated axes')
    zero = s.kronecker_product(zero_vector(), zero_vector(), zero_vector(), s.eye(2))
    if actual*zero != s.zeros(54, 2) or zero.H*zero != 8*s.eye(2):
        raise ArithmeticError('the declared native zero pair failed')
    return {'original_hopping_entries_checked': 54*54,
            'source_geometry': [3, 3, 3], 'source_speed': '1',
            'onsite_mass_entries_added': 0,
            'zero_orbital_unnormalized_one_dimension': [1, 0, 1],
            'zero_orbital_corner_probability': '1/8',
            'zero_orbital_central_probability': '0',
            'complete_one_body_dimension': 54}


def exact_algebra():
    one = chain()
    reversal = s.eye(3)[:, [2, 1, 0]]
    stagger = s.diag(1, -1, 1)
    theta = s.I*SIGMA[1]
    x = s.symbols('x')
    if (one*zero_vector() != s.zeros(3, 1)
            or reversal*one*reversal != -one
            or stagger*one*stagger != -one
            or reversal*zero_vector() != zero_vector()
            or one.charpoly(x).as_expr() != x**3-x/2
            or theta*theta.conjugate() != -s.eye(2)):
        raise ArithmeticError('odd-box zero, gap or reflection identity failed')
    for sigma in SIGMA:
        link = s.kronecker_product(one, sigma)
        time = s.kronecker_product(s.eye(3), theta)
        ph = s.kronecker_product(stagger, theta)
        if (time*link.conjugate()*time.H != link
                or ph*link.T*ph.H != -link):
            raise ArithmeticError('source Kramers or particle-hole identity failed')
    ipr = s.Rational(1, 8)
    sea = (1-ipr)*s.eye(2)/2
    if (s.trace(sea)-s.Rational(1, 2))*s.eye(2)-sea != -ipr*s.eye(2)/2:
        raise ArithmeticError('filled-sea centered contact linear term failed')
    return {'one_dimensional_characteristic_polynomial': 'x^3-x/2',
            'exact_one_body_gap': 'c*sin(pi/(L+1))',
            'exact_zero_orbital_IPR_dimension_d': '(2/(L+1))^d',
            'weighted_sea_density': '(1-IPR) I_2/2',
            'exact_compressed_centered_contact': 'C+U*IPR*(n1*n2-(n1+n2)/2)',
            'source_spin_time_reversal_commutes': True,
            'staggered_spin_particle_hole_preserves_source': True,
            'exact_odd_half_filled_Kramers_degeneracy': True}


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def certificate():
    cube = core_bounds()
    family = []
    for length in (3, 5, 7, 15):
        row = core_bounds(F(1, 2**16*(length+1)**10), length)
        if (not row['charged_gap_certified_positive']
                or 2*row['ordered_low_energy_error'] > row['projected_charged_gap']/2**21):
            raise ArithmeticError('the controlled odd-box diagonal family failed')
        family.append(row)
    if not cube['charged_gap_certified_positive']:
        raise ArithmeticError('finite unmodified-source charged gap failed')
    return {'source_anchor': exact_source_anchor(), 'algebra': exact_algebra(),
            'controlled_27_site_cube': encode(cube),
            'controlled_odd_length_family': encode(family),
            'family_coupling': 'U_L=c/(2^16*(L+1)^10)',
            'family_projected_charged_gap': '4*c/(2^16*(L+1)^13)',
            'full_Fock_virtual_states_retained': True,
            'numerical_eigenvalues_used_for_proof_admission': False,
            'two_to_54_Fock_matrix_executed': False,
            'boundary_standing_wave_zero_mode': True,
            'topological_defect_or_autonomous_mass_field_supplied': False,
            'volume_uniform_nonzero_charged_gap_proved': False}
