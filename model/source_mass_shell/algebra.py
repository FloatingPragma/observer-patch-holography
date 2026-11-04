"""Exact full-source Fourier join and mass-sector spectral projectors."""
from functools import lru_cache
from itertools import product

import sympy as s

from model.native_mass_channels import algebra as source
from model.periodic_momentum.momentum import PeriodicGeometry, one_particle


def spectral_projectors(kinetic, coefficients):
    """Exact rational kinetic coordinates; algebraic mass and energy roots."""
    k, masses = source.vector(kinetic, 3), source.vector(coefficients, 4)
    h = source.folded_hamiltonian(k, masses)
    _, beta = source.folded_generators()
    m0, *spin = tuple(s.Rational(x) for x in masses)
    r2 = sum(x*x for x in spin)
    r = s.sqrt(r2)
    if r:
        q = beta[0]*sum((x*b for x, b in zip(spin, beta[1:], strict=True)), s.zeros(16))/r
        sectors = [(sign, (s.eye(16)+sign*q)/2, m0+sign*r) for sign in (-1, 1)]
    else:
        sectors = [(0, s.eye(16), m0)]
    rows = []
    for label, sector, mass in sectors:
        energy = s.sqrt(sum(s.Rational(x)**2 for x in k)+mass*mass)
        for sign in ((-1, 1) if energy != 0 else (0,)):
            projector = sector if sign == 0 else sector*(s.eye(16)+sign*h/energy)/2
            rows.append((label, sign*energy, s.ImmutableMatrix(projector.applyfunc(s.simplify))))
    return tuple(rows)


@lru_cache(maxsize=1)
def source_join():
    """All128 physical rows on a periodic4^3 cube, using original source links."""
    geometry = PeriodicGeometry((4, 4, 4))
    kinetic = one_particle(geometry)
    coefficients = (s.Rational(2, 7), s.Rational(1, 3), s.Rational(2, 5), s.Rational(-3, 11))
    rational_coefficients = tuple(source.rational(int(z.p))/int(z.q) for z in coefficients)
    h = s.SparseMatrix(kinetic)+s.SparseMatrix(s.diag(*(
        source.onsite(x, rational_coefficients) for x in geometry.sites)))
    sine = (0, 1, 0, -1)
    labels = set()
    dft = s.Matrix(4, 4, lambda x, k: s.I**((x*k) % 4)/2)
    if dft.H*dft != s.eye(4):
        raise ArithmeticError('exact one-axis Fourier basis is not unitary')
    for base in product((0, 1), repeat=3):
        embed = s.SparseMatrix(128, 16, {
            (2*x+spin, 2*j+spin): s.I**(sum(n*(q+2*b) for n, q, b in zip(site, base, bits, strict=True)) % 4)/8
            for x, site in enumerate(geometry.sites)
            for j, bits in enumerate(source.SITES) for spin in range(2)})
        block = s.SparseMatrix(source.folded_hamiltonian(
            tuple(geometry.speed*sine[q] for q in base),
            rational_coefficients))
        if any(s.expand(z) != 0 for z in h*embed-embed*block):
            raise ArithmeticError('full-source folded Fourier identity failed')
        if embed.H*embed != s.eye(16):
            raise ArithmeticError('folded Fourier columns are not orthonormal')
        labels.update(tuple(q+2*b for q, b in zip(base, bits, strict=True)) for bits in source.SITES)
    # The tensor product of the exact one-axis DFT is unitary. Its full64
    # momentum labels each carry the two independent physical spin columns.
    if labels != set(product(range(4), repeat=3)):
        raise ArithmeticError('folded source basis omitted or repeated a valley')
    return {'periodic_sides': [4, 4, 4], 'source_rows': 128,
            'reduced_momenta': 8, 'rows_per_reduced_momentum': 16,
            'all_four_mass_channels_present': True,
            'all_source_Fourier_columns_exact': True,
            'supplied_profile_is_period_two': True}


def clifford_certificate():
    gamma, beta = source.folded_generators()
    one, zero = s.eye(16), s.zeros(16)
    lorentz = (beta[0], *(beta[0]*g for g in gamma))
    metric = (1, -1, -1, -1)
    for i, j in product(range(4), repeat=2):
        target = 2*metric[i]*one if i == j else zero
        if lorentz[i]*lorentz[j]+lorentz[j]*lorentz[i] != target:
            raise ArithmeticError('massive full-row Lorentz Clifford identity failed')
    q = beta[0]*(beta[1]+2*beta[2]+2*beta[3])/3
    if q*q != one or s.trace(q) != 0:
        raise ArithmeticError('exact rational-direction mass-sector grading failed')
    if any(q*g-g*q != zero for g in lorentz):
        raise ArithmeticError('mass sector is not preserved by the Dirac algebra')
    for label, energy, projector in spectral_projectors((1, 2, 3), (2, 1, 2, 2)):
        h = source.folded_hamiltonian((1, 2, 3), (2, 1, 2, 2))
        if ((projector*projector-projector).applyfunc(s.simplify) != zero or
                (h*projector-energy*projector).applyfunc(s.simplify) != zero or s.trace(projector) != 4):
            raise ArithmeticError('exact signed mass-sector projector failed')
    return {'Clifford_checks': 16, 'metric': [1, -1, -1, -1],
            'full_dimension': 16, 'spin_nonzero_mass_sector_dimensions': [8, 8],
            'Dirac_copies_per_mass_sector': 2,
            'generic_each_signed_energy_rank': 4,
            'source_join': source_join(),
            'physical_Lorentz_spin_identification': False}
