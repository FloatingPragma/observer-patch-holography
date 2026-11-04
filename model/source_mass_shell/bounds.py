"""Exact rational energy enclosures and radial exterior-packet bounds."""
from fractions import Fraction as F

from model.contact_scattering.onshell import pi_interval
from model.native_mass_channels.algebra import rational, vector
from model.relativistic_sector_audit.sector import sine_interval, sqrt_interval


def _square(interval):
    lo, hi = interval
    return (F(0) if lo <= 0 <= hi else min(lo*lo, hi*hi), max(lo*lo, hi*hi))


def _root(interval):
    return sqrt_interval(interval[0])[0], sqrt_interval(interval[1])[1]


def _strings(interval):
    return [str(x) for x in interval]


def parameters(spacing, velocity, coefficients):
    a, v = rational(spacing), rational(velocity)
    mass = vector(coefficients, 4)
    if a <= 0 or v <= 0:
        raise ValueError('positive exact spacing and velocity required')
    return a, v, mass


def energy_rows(momentum, spacing=F(1, 16), velocity=1, coefficients=(0, 0, 0, 1)):
    """All signed energies; intervals enclose irrational mass directions too.

    The rows label the two commuting eight-dimensional mass sectors when the
    spin mass is nonzero. Coincident signed energies add their multiplicities.
    With zero spin mass there is a single sixteen-dimensional sector.
    """
    a, v, mass = parameters(spacing, velocity, coefficients)
    p = vector(momentum, 3)
    if max(abs(a*x) for x in p) > 1:
        raise ValueError('every absolute physical-model momentum times spacing must be at most one')
    m0, *spin = mass
    r2 = sum((x*x for x in spin), F())
    rlo, rhi = sqrt_interval(r2)
    sine = [sine_interval(abs(a*x)) for x in p]
    kinetic = tuple(v*v/a**2*sum((q[j]**2 for q in sine), F()) for j in (0, 1))
    linear = v*v*sum((x*x for x in p), F())
    shell_error = v*v*a*a*sum((x**4 for x in p), F())/3
    energy_error = v*a*a*sum((abs(x)**3 for x in p), F())/6
    branches = (-1, 1) if r2 else (0,)
    rows = []
    for branch in branches:
        signed_mass = ((m0-rhi, m0-rlo) if branch == -1 else
                       (m0+rlo, m0+rhi) if branch == 1 else (m0, m0))
        mass_squared = _square(signed_mass)
        if m0*m0 == r2 and (branch == 0 or m0*branch <= 0):
            mass_squared = (F(0), F(0))
        rest = _root(mass_squared)
        continuum_squared = (linear+mass_squared[0], linear+mass_squared[1])
        lattice_squared = (kinetic[0]+mass_squared[0], kinetic[1]+mass_squared[1])
        continuum, lattice = _root(continuum_squared), _root(lattice_squared)
        zero = not any(p) and m0*m0 == r2 and (branch == 0 or m0*branch <= 0)
        rank = 8 if branch == 0 else 4
        signs = (0,) if zero else (-1, 1)
        for sign in signs:
            rows.append({
                'mass_sector': branch, 'energy_sign': sign,
                'rank': 2*rank if zero else rank,
                'rest_energy_interval': _strings(rest),
                'rest_mass_in_model_units_interval': _strings(tuple(x/(v*v) for x in rest)),
                'lattice_energy_interval': _strings((-lattice[1], -lattice[0]) if sign == -1 else lattice),
                'continuum_energy_interval': _strings((-continuum[1], -continuum[0]) if sign == -1 else continuum),
                'mass_shell_defect_interval': [str(-shell_error), '0'],
                'absolute_energy_error_upper': str(energy_error),
            })
    if sum(row['rank'] for row in rows) != 16:
        raise ArithmeticError('energy rows do not retain the complete source fiber')
    return {'momentum': [str(x) for x in p], 'spacing': str(a), 'velocity': str(v),
            'mass_coefficients': [str(x) for x in mass], 'full_fiber_dimension': 16,
            'all_valleys': 8, 'rows': rows,
            'energy_formula': 'E=+/-sqrt((v/a)^2 sum_j sin(a*p_j)^2+(m0+s*sqrt(sum_j m_j^2))^2)',
            'continuum_formula': 'E=+/-sqrt(v^2 |p|^2+(m0+s*sqrt(sum_j m_j^2))^2)',
            'interacting_pole_energy': False}


def exterior_packet_bound(spacing, velocity, mass, cells, displacement, momentum=(0, 0, 0)):
    """Normalized compact packets at anchor (0,0,R), with L cells per axis.

    A cell has eight physical sites and two CAR rows per site. The envelope
    is the product of normalized discrete Dirichlet sine vectors. R is even
    so the original staggered mass has the declared cell frame.
    """
    a, v, masses = parameters(spacing, velocity, (0, 0, 0, mass))
    m = masses[3]
    p = vector(momentum, 3)
    if m <= 0 or max(abs(a*x) for x in p) > 1:
        raise ValueError('positive radial mass and all |a*p_j|<=1 required')
    if type(cells) is not int or cells < 1:
        raise ValueError('positive integer cell width required')
    if type(displacement) is not int or displacement < 2 or displacement % 2:
        raise ValueError('positive even lattice displacement required')
    pi_upper = pi_interval()[1]
    sqrt3_upper = sqrt_interval(F(3))[1]
    kinetic = 3*pi_upper*v/(a*(cells+1))
    texture = 4*sqrt3_upper*m*cells/displacement
    lattice_error = kinetic+texture
    shell_linearization = min(v*a*a*sum((abs(x)**3 for x in p), F())/6,
                              v*v*a*a*sum((x**4 for x in p), F())/(6*m))
    return {'spacing': str(a), 'velocity': str(v), 'radial_mass': str(m),
            'cells_per_axis': cells, 'anchor_lattice_coordinates': [0, 0, displacement],
            'support_sites': (2*cells)**3, 'full_packet_rows': 2*(2*cells)**3,
            'kinetic_envelope_error_upper': str(kinetic),
            'radial_texture_error_upper': str(texture),
            'lattice_eigenvector_residual_upper': str(lattice_error),
            'continuum_energy_linearization_upper': str(shell_linearization),
            'continuum_eigenvector_residual_upper': str(lattice_error+shell_linearization),
            'negative_and_positive_branches': True,
            'eigensolver_or_large_packet_execution_used_for_bound': False,
            'interacting_charged_residue_inferred': False}


def certificate():
    return {
        'uniform_cube_bounds': {
            'domain': '|p_j|<=Lambda, a*Lambda<=1',
            'mass_shell_defect': '-v^2*a^2*Lambda^4 <= E_a^2-v^2*|p|^2-M_s^2 <=0',
            'absolute_energy_error': '<=v*a^2*Lambda^3/2',
            'group_velocity_norm': '<=v',
            'mass_sign_and_zero_branches_retained': True},
        'bulk_samples': [energy_rows(p, a, 1, coefficients)
                         for coefficients in ((0, 0, 0, 1), (1, 0, 0, 1), (2, 1, 1, 0))
                         for a in (F(1, 4), F(1, 16), F(1, 64))
                         for p in ((0, 0, 0), (1, 0, 0), (1, 1, 1))],
        'radial_bulk_energy_table': [energy_rows(p, a, 1, (0, 0, 0, 64))
                                    for a in (F(1), F(1, 4), F(1, 16), F(1, 64))
                                    for p in ((0, 0, 0), (1, 0, 0), (1, 1, 1))],
        'same_radial_source_exterior_witness': exterior_packet_bound(
            F(1, 64), 1, 64, 10**6, 10**15, (1, 1, 1)),
        'radial_continuum_sequence': 'a=1/n,L=n^3,R=2*n^8; n positive integer',
        'radial_continuum_residual_order': 'O(v/n^2)+O(m/n^5)',
    }
