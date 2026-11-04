"""Original Weyl source with a supplied linear spin texture and its Jacobi map."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial

import mpmath as mp
import sympy as s

from ..native_mass_channels import algebra as source


def parameters(v=1, kappa=1, spacing=1):
    result = tuple(source.rational(x) for x in (v, kappa, spacing))
    if any(x <= 0 for x in result):
        raise ValueError('positive exact speed, texture slope and spacing required')
    return result


def _radius(radius):
    if type(radius) is not int or radius < 1:
        raise ValueError('positive integer symmetric cutoff radius required')
    return radius


def jacobi(radius, v=1, kappa=1, spacing=1):
    """B=epsilon*(v D+ kappa*a*n), ordered n=-R,...,R, open ends."""
    radius = _radius(radius)
    speed, slope, a = (s.Rational(x) for x in parameters(v, kappa, spacing))
    entries = {}
    for index, n in enumerate(range(-radius, radius+1)):
        parity = (-1)**(n % 2)
        entries[index, index] = slope*a*n*parity
        if n < radius:
            entries[index, index+1] = entries[index+1, index] = parity*speed/(2*a)
    return s.ImmutableSparseMatrix(2*radius+1, 2*radius+1, entries)


def first_order(radius, v=1, kappa=1, spacing=1):
    radius = _radius(radius)
    speed, slope, a = (s.Rational(x) for x in parameters(v, kappa, spacing))
    entries = {}
    for index, n in enumerate(range(-radius, radius+1)):
        entries[index, index] = slope*a*n
        if n < radius:
            entries[index, index+1] = speed/(2*a)
            entries[index+1, index] = -speed/(2*a)
    return s.ImmutableSparseMatrix(2*radius+1, 2*radius+1, entries)


@lru_cache(maxsize=8)
def _parity_phase(parities):
    result = s.prod(1+s.I*(-1)**j for j in parities)/s.sqrt(8)
    return s.simplify(s.expand_complex(result))


def phase(site):
    if type(site) is not tuple or len(site) != 3 or any(type(x) is not int for x in site):
        raise ValueError('three integer source site coordinates required')
    return _parity_phase(tuple(x % 2 for x in site))


def source_box(radius, v=1, kappa=1, spacing=1):
    """Full two-spin original nearest-neighbor Weyl links and onsite texture."""
    radius = _radius(radius)
    speed, slope, a = (s.Rational(x) for x in parameters(v, kappa, spacing))
    sites = tuple(product(range(-radius, radius+1), repeat=3))
    addresses = {site: j for j, site in enumerate(sites)}
    entries = {}
    for site, address in addresses.items():
        onsite = sum((slope*a*site[j]*(-1)**(site[j] % 2)*source.SIGMA[j]
                      for j in range(3)), s.zeros(2))
        for alpha in range(2):
            for beta in range(2):
                if onsite[alpha, beta]:
                    entries[2*address+alpha, 2*address+beta] = onsite[alpha, beta]
        for axis in range(3):
            other = tuple(n+int(j == axis) for j, n in enumerate(site))
            if other not in addresses:
                continue
            target = addresses[other]
            link = s.I*speed/(2*a)*source.SIGMA[axis]
            for alpha in range(2):
                for beta in range(2):
                    if link[alpha, beta]:
                        entries[2*address+alpha, 2*target+beta] = link[alpha, beta]
                        entries[2*target+beta, 2*address+alpha] = s.conjugate(link[alpha, beta])
    return s.ImmutableSparseMatrix(2*len(sites), 2*len(sites), entries)


def transformed_box(radius, v=1, kappa=1, spacing=1):
    radius = _radius(radius)
    one = jacobi(radius, v, kappa, spacing)
    identity = s.eye(one.rows)
    result = s.zeros(2*one.rows**3)
    for axis in range(3):
        factors = [one if j == axis else identity for j in range(3)]
        result += s.kronecker_product(*factors, source.SIGMA[axis])
    return s.ImmutableSparseMatrix(result)


def bessel_interval(order, z=1, terms=40):
    """Positive-series rational enclosure of I_order(z), with geometric tail."""
    if type(order) is not int or type(terms) is not int or terms < 0:
        raise ValueError('integer Bessel order and nonnegative series cutoff required')
    z = source.rational(z)
    if z <= 0:
        raise ValueError('positive exact Bessel argument required')
    n = abs(order)
    term = (z/2)**n/factorial(n)
    total = term
    for j in range(1, terms+1):
        term *= z*z/(4*j*(n+j))
        total += term
    following = term*z*z/(4*(terms+1)*(n+terms+1))
    ratio = z*z/(4*(terms+2)*(n+terms+2))
    if ratio >= 1:
        raise ValueError('increase the series cutoff to certify its tail')
    return total, total+following/(1-ratio)


def bessel_zero(site, v=1, kappa=1, spacing=1):
    """High-precision diagnostic amplitude, normalized on the infinite lattice."""
    phase(site)
    speed, slope, a = parameters(v, kappa, spacing)
    z = speed/(slope*a*a)
    value = mp.mpf(z.numerator)/z.denominator
    p = mp.fprod(mp.mpc(1, (-1)**(n % 2)) for n in site)/mp.sqrt(8)
    return p*mp.fprod(mp.besseli(abs(n), value) for n in site)/mp.besseli(0, 2*value)**mp.mpf('1.5')


def certificate():
    for site in product(range(-2, 3), repeat=3):
        value = phase(site)
        if s.simplify(s.conjugate(value)*value) != 1:
            raise ArithmeticError('source phase map is not unitary')
        for j in range(3):
            for direction in (-1, 1):
                other = tuple(n+direction*int(k == j) for k, n in enumerate(site))
                if s.simplify(phase(other)+s.I*(-1)**(site[j] % 2)*value) != 0:
                    raise ArithmeticError('source nearest-neighbor phase identity failed')
    finite = jacobi(3)
    A = first_order(3)
    reflection = s.eye(7)[:, list(reversed(range(7)))]
    if finite != finite.H or finite*finite != A.H*A or reflection*finite*reflection != -finite:
        raise ArithmeticError('source Jacobi square or reflection changed')
    return {'source_onsite_texture': 'kappa*a*sum_j x_j*(-1)^x_j*sigma_j',
            'source_kinetic_coefficient': 'c=v/a',
            'supplied_texture_unbounded_at_spatial_infinity': True,
            'phase': 'product_j(1+i*(-1)^x_j)/sqrt(8)',
            'transformed_operator': 'sum_j sigma_j*B_j, B=epsilon*(v*D+kappa*a*n)',
            'squared_operator': 'sum_j A_j^* A_j',
            'Jacobi_reflection_anticommutes': True,
            'one_dimensional_zero_mode': 'I_n(z)/sqrt(I_0(2z)), z=v/(kappa*a^2)',
            'full_exact_lattice_zero_mode_dimension': 2,
            'complete_lattice_spectrum': '0 and +/-sqrt(b_nx^2+b_ny^2+b_nz^2), b_0=0',
            'nonzero_signed_multiplicity_per_ordered_positive_index_triple': '2^(number of nonzero indices); coincident energies add',
            'one_dimensional_eigenvalues_simple': True,
            'uniform_low_spectrum_continuum_rate_asserted': False,
            'contact_selected_texture_or_physical_species_proved': False,
            'I0_at_one_enclosure': [str(x) for x in bessel_interval(0, 1)],
            'I1_at_one_enclosure': [str(x) for x in bessel_interval(1, 1)]}
