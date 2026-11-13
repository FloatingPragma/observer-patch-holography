"""Actual spin and quantum-link operators; a declared intraband source."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from model.openwave_transverse_completion.source import flux
import model.openwave_common_mass_source.source as parent

ROOT = Path(parent.__file__).resolve().parents[2]
ELECTRIC = F(1, 32)
MAGNETIC = F(1, 8)
# The quadratic Maxwell symbol has c²=2 E B. A centered momentum has
# symbol 2 h sin(k), so h²=E B/2. Its square is exactly rational.
MOMENTUM_SQUARED = ELECTRIC*MAGNETIC/2
SPEED_SQUARED = 4*MOMENTUM_SQUARED
COORDINATES = ((0, 0), (1, 0), (1, 1), (0, 1))
SIGMA = (np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex),
         np.diag([1., -1.]).astype(complex))


def charged_masses():
    receipt = json.loads((ROOT/'model/openwave_common_mass_source/results.json').read_text())
    return [(r['channel'], F(r['energy_electron_units_exact']))
            for r in receipt['same_source_free_Dirac_branches'][:3]]


def orbital_operators(cap, charge, theta=0.):
    if type(cap) is not int or cap < 1 or type(charge) is not int:
        raise ValueError('integer cap and native charge required')
    basis = [(x, n) for x in range(4) for n in range(-cap, cap+1)]
    index = {b: i for i, b in enumerate(basis)}
    shape = (len(basis),)*2
    axes = [sparse.lil_matrix(shape, dtype=complex) for _ in range(2)]
    first = [sparse.lil_matrix(shape, dtype=complex) for _ in range(2)]
    second = [sparse.lil_matrix(shape, dtype=complex) for _ in range(2)]
    lap = 2*sparse.eye(len(basis), dtype=complex, format='lil')
    lap_first = sparse.lil_matrix(shape, dtype=complex)
    lap_second = sparse.lil_matrix(shape, dtype=complex)
    for (x, n), i in index.items():
        y = (x+1) % 4
        nn = n-charge if x == 3 else n
        if (y, nn) not in index:
            continue
        j = index[(y, nn)]
        axis = x % 2
        spatial_sign = 1 if x in (0, 1) else -1
        q = charge if x == 0 else 0
        phase = np.exp(1j*q*theta)
        value = -1j*spatial_sign*phase
        for target, derivative in ((axes, value), (first, 1j*q*value),
                                   (second, -q*q*value)):
            target[axis][j, i] = derivative
            target[axis][i, j] = derivative.conjugate()
        for target, value in ((lap, -phase), (lap_first, -1j*q*phase),
                              (lap_second, q*q*phase)):
            target[j, i] = value
            target[i, j] = value.conjugate()
    matrices = [[x.tocsr() for x in row] for row in (axes, first, second)]
    return basis, index, matrices, [a.tocsr() for a in (lap, lap_first, lap_second)]


def pauli_operator(cap, charge, theta=0.):
    basis, index, ops, lap = orbital_operators(cap, charge, theta)
    matrices = [sum((sparse.kron(SIGMA[i], row[i], format='csr')
                     for i in range(2)), sparse.csr_matrix((2*len(basis),)*2, dtype=complex))
                for row in ops]
    d, dp, dpp = matrices
    scale = float(MOMENTUM_SQUARED)
    a = scale*(d@d)
    ap = scale*(dp@d+d@dp)
    app = scale*(dpp@d+2*(dp@dp)+d@dpp)
    px, py = ops[0]
    magnetic = 1j*scale*sparse.kron(SIGMA[2], px@py-py@px, format='csr')
    scalar = scale*sparse.kron(sparse.eye(2), px@px+py@py, format='csr')
    residual = a-scalar-magnetic
    assert not residual.nnz or np.max(abs(residual.data)) < 2e-17
    return {'orbital_basis': basis, 'orbital_index': index, 'D': np.sqrt(scale)*d,
            'D2': a, 'D2_first': ap, 'D2_second': app,
            'magnetic_commutator': magnetic, 'scalar_square': scalar,
            'covariant_laplacian_derivatives': [sparse.kron(sparse.eye(2), z, format='csr') for z in lap]}


def binomial_half(k):
    value = F(1)
    for j in range(k):
        value *= (F(1, 2)-j)/(j+1)
    return value


def root_tail(mass, degree):
    mass = F(mass)
    if mass < 1 or degree < 1:
        raise ValueError('this charged-family expansion requires m>=1 and degree>=1')
    rho = SPEED_SQUARED/(mass*mass)
    return mass*rho**(degree+1)/(1-rho)


def polynomial_member(cap, charge, mass, degree=6, theta=0., regulator=F(1, 16)):
    """P_cap times the infinite Pauli polynomial, not a squared clipped Dirac matrix."""
    mass = F(mass)
    regulator = F(regulator)
    if regulator < 0:
        raise ValueError('the declared Wilson energy coefficient is nonnegative')
    tail = root_tail(mass, degree)
    # At most two elementary charged hops per D² factor; all return paths remain.
    pad = 2*degree*abs(charge)+2
    op = pauli_operator(cap+pad, charge, theta)
    a, ap, app = op['D2'], op['D2_first'], op['D2_second']
    size = a.shape[0]
    power = sparse.eye(size, dtype=complex, format='csr')
    first = sparse.csr_matrix((size, size), dtype=complex)
    second = first.copy()
    result, current, contact = first.copy(), first.copy(), first.copy()
    for k in range(1, degree+1):
        second = app@power+2*ap@first+a@second
        first = ap@power+a@first
        power = a@power
        coefficient = float(binomial_half(k)*mass**(1-2*k))
        result += coefficient*power
        current += coefficient*first
        contact += coefficient*second
    # Positive covariant graph energy; its free continuum effect is O(a).
    # This is a declared intraband regulator, not a canonical-QED FW transfer.
    lap, lap_first, lap_second = op['covariant_laplacian_derivatives']
    result += float(regulator/2)*lap
    current += float(regulator/2)*lap_first
    contact += float(regulator/2)*lap_second
    orbital = op['orbital_basis']
    kept = [s*len(orbital)+i for s in range(2)
            for i, (_, n) in enumerate(orbital) if abs(n) <= cap]
    crop = lambda matrix: matrix[kept][:, kept].toarray()
    basis = [(s, x, n) for s in range(2) for x in range(4) for n in range(-cap, cap+1)]
    idx = {b: i for i, b in enumerate(basis)}
    field = np.zeros((len(basis),)*2, complex)
    for (spin, x, n), i in idx.items():
        field[i, i] = float(ELECTRIC*sum(e*e for e in flux(x, n, charge))+MAGNETIC)
        for step in (-1, 1):
            if (spin, x, n+step) in idx:
                field[idx[(spin, x, n+step)], i] = -float(MAGNETIC/2)
    kinetic = crop(result)
    h = field+kinetic
    assert np.max(abs(h-h.conj().T)) < 1e-14
    # The compensating reference sits at (0,0). This relative dipole is a
    # gauge-invariant charge observable on the same reduced Gauss carrier.
    x_coordinate = [x for x, _ in COORDINATES]
    dipole = np.array([charge*x_coordinate[x] for _, x, _ in basis], dtype=float)
    difference = dipole[None, :]-dipole[:, None]
    dipole_current = 1j*difference*h
    dipole_contact = -(difference*difference)*h
    return {'basis': basis, 'index': idx, 'H_minus_mass': h,
            'kinetic': kinetic, 'first_gauge_derivative': crop(current),
            'second_gauge_derivative': crop(contact),
            'magnetic_commutator': crop(op['magnetic_commutator']),
            'relative_x_dipole_diagonal': dipole,
            'relative_dipole_current': dipole_current,
            'relative_dipole_second_derivative': dipole_contact,
            'relative_dipole_charge_unit': 'native third-charge q=3Q',
            'electron_charge_relative_x_dipole_diagonal': dipole/3,
            'electron_charge_relative_dipole_current': dipole_current/3,
            'electron_charge_relative_dipole_second_derivative': dipole_contact/9,
            'plaquette_coordinates': COORDINATES,
            'reference_coordinate': COORDINATES[0],
            'Maxwell_electric_coefficient_exact': str(ELECTRIC),
            'Maxwell_magnetic_coefficient_exact': str(MAGNETIC),
            'common_free_speed_squared_exact': str(SPEED_SQUARED),
            'positive_Wilson_energy_coefficient_exact': str(regulator),
            'mass_exact': str(mass), 'polynomial_degree': degree,
            'uniform_root_operator_error_exact': str(tail),
            'external_probe': 'classical link0 phase in matter covariant momentum; field Hamiltonian fixed',
            'full_compact_flux_source': True, 'source_is_number_preserving': True,
            'is_full_signed_Dirac_plus_field_unitary_transfer': False}
