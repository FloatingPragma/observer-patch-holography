"""Original-H doubled Weyl limit, with every outgoing contact state retained.

The continuum embedding is an initial-state comparison, not a modification
of the source Hamiltonian. Infinite-lattice Haar fractions and finite-torus
checks are kept separate. See README.md for the operator-domain proof.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial, isqrt

import numpy as np
import sympy as s
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..contact_scattering.onshell import pi_interval
from ..dilute_momentum_ward.ward import sqrt_upper
from ..interacting_field.lattice import TwoFermions, node_table
from ..periodic_momentum.momentum import PeriodicGeometry, one_particle

PAULI = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))


def exact(value, name):
    if type(value) not in (int, F):
        raise ValueError(f'exact rational {name} required')
    return F(value)


def sqrt_interval(value, bits=80):
    value = exact(value, 'radicand')
    if value < 0 or type(bits) is not int or not 1 <= bits <= 256:
        raise ValueError('nonnegative radicand and 1..256 bits required')
    scale = 1 << bits
    numerator = value.numerator*scale*scale
    root = isqrt(numerator//value.denominator)
    lower = F(root, scale)
    upper = lower if root*root*value.denominator == numerator else F(root+1, scale)
    return lower, upper


def sine_interval(value):
    """Exact alternating enclosure on [-1,1], with no floating sine input."""
    value = exact(value, 'angle')
    if abs(value) > 1:
        raise ValueError('absolute angle at most one required')
    if value < 0:
        lo, hi = sine_interval(-value)
        return -hi, -lo
    lower = sum(((-1)**j*value**(2*j+1)/factorial(2*j+1) for j in range(8)), F())
    return lower, lower+value**17/factorial(17)


def continuum_budget(spacing, cutoff, maximum_number, duration, *, speed=F(1), interaction=F(1)):
    """Outward rational full-vector bound for the infinite-lattice embedding.

    The calculator uses a*Lambda<=1, a sufficient strict nonoverlap condition.
    The proof needs only a*Lambda<pi/2. The physical momentum cube, particle
    number, v and time are fixed when taking a to zero. U>=0 is the source domain.
    """
    a, cutoff, t, v, u = (exact(x, name) for x, name in
                         ((spacing, 'spacing'), (cutoff, 'cutoff'), (duration, 'time'),
                          (speed, 'speed'), (interaction, 'interaction')))
    if (a <= 0 or cutoff <= 0 or a*cutoff > 1 or v <= 0 or u < 0 or
            type(maximum_number) is not int or not 0 <= maximum_number <= 4096):
        raise ValueError('positive a,Lambda,v, a*Lambda<=1, U>=0 and Nmax in 0..4096 required')
    pi_lo, pi_hi = pi_interval()
    fraction_lower = 8*(a*cutoff/pi_hi)**3
    fraction_upper = 8*(a*cutoff/pi_lo)**3
    d = maximum_number//2
    collision = abs(u)*d*sqrt_upper(fraction_upper, bits=80)
    linearization = maximum_number*v*a*a*cutoff**3/2
    total = min(F(2), abs(t)*(collision+linearization))
    return {'spacing': a, 'cutoff': cutoff, 'maximum_number': maximum_number,
            'duration': t, 'speed': v, 'interaction': u,
            'band_fraction_lower': fraction_lower, 'band_fraction_upper': fraction_upper,
            'contact_on_band_norm_upper': collision,
            'free_generator_error_upper': linearization,
            'complete_vector_error_upper': total,
            'half_trace_distance_upper': min(F(1), total),
            'one_particle_projector_error_upper': a*a*cutoff*cutoff/6,
            'one_particle_mass_shell_defect_magnitude_upper': v*v*a*a*cutoff**4}


def energy_row(momentum, spacing=F(1, 64), speed=F(1)):
    """Rational enclosures for both original-H branches, shared by all valleys."""
    if type(momentum) not in (tuple, list) or len(momentum) != 3:
        raise ValueError('three exact momentum components required')
    p = tuple(exact(x, 'momentum') for x in momentum)
    a, v = exact(spacing, 'spacing'), exact(speed, 'speed')
    if a <= 0 or v <= 0 or max(abs(a*x) for x in p) > 1:
        raise ValueError('positive spacing and speed, all |a*p_j|<=1 required')
    intervals = [sine_interval(abs(a*x)) for x in p]
    lattice_squared = tuple(v*v/a**2*sum((interval[j]**2 for interval in intervals), F())
                            for j in (0, 1))
    lo = sqrt_interval(lattice_squared[0])[0]
    hi = sqrt_interval(lattice_squared[1])[1]
    continuum_squared = v*v*sum((x*x for x in p), F())
    clo, chi = sqrt_interval(continuum_squared)
    dispersion_bound = v*a*a*sum((abs(x)**3 for x in p), F())/6
    shell_bound = v*v*a*a*sum((x**4 for x in p), F())/3
    return {'momentum': [str(x) for x in p], 'spacing': str(a), 'speed': str(v),
            'lattice_positive_energy_interval': [str(lo), str(hi)],
            'lattice_negative_energy_interval': [str(-hi), str(-lo)],
            'continuum_positive_energy_interval': [str(clo), str(chi)],
            'continuum_negative_energy_interval': [str(-chi), str(-clo)],
            'continuum_energy_squared': str(continuum_squared),
            'absolute_dispersion_error_upper': str(dispersion_bound),
            'lattice_mass_shell_defect_interval': [str(-shell_bound), '0'],
            'continuum_mass': '0', 'all_eight_valleys': True,
            'lattice_positive_energy_display': float((lo+hi)/2),
            'continuum_positive_energy_display': float((clo+chi)/2),
            'zero_energy_degenerate_spin_space': continuum_squared == 0}


def energy_table():
    return [energy_row(p) for p in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1))]


def relativistic_algebra_certificate():
    """Exact Clifford algebra of an opposite-chirality pair; no added mass."""
    identity = s.eye(2)
    zero = s.zeros(2)
    beta = zero.row_join(identity).col_join(identity.row_join(zero))
    alpha = tuple(s.diag(-sigma, sigma) for sigma in PAULI)
    gamma = (beta, *(beta*a for a in alpha))
    metric = (1, -1, -1, -1)
    for i, j in product(range(4), repeat=2):
        target = 2*metric[i]*s.eye(4) if i == j else s.zeros(4)
        if gamma[i]*gamma[j]+gamma[j]*gamma[i] != target:
            raise ArithmeticError('massless opposite-chirality Clifford pair failed')
    coefficients = s.symbols('m0:4', real=True)
    mass = coefficients[0]*identity+sum((coefficients[j+1]*PAULI[j] for j in range(3)), zero)
    equations = [z for sigma in PAULI for z in mass*sigma+sigma*mass]
    matrix, rhs = s.linear_eq_to_matrix(equations, coefficients)
    if matrix.rank() != 4 or rhs != s.zeros(len(equations), 1):
        raise ArithmeticError('single two-component constant mass obstruction failed')
    return {'opposite_chirality_Dirac_pair_Clifford_checks': 16,
            'metric_signature': [1, -1, -1, -1],
            'single_two_component_constant_anticommuting_mass_dimension': 0,
            'massive_term_added': False,
            'valley_pairing_selects_physical_Dirac_species': False,
            'exclusion_scope': 'a constant 2x2 Hermitian mass anticommuting with all three kinetic Pauli matrices; not a many-body mass theorem'}


@lru_cache(maxsize=1)
def valley_certificate():
    """Exact Pauli algebra and all 128 source Fourier columns on a 4^3 torus."""
    geometry = PeriodicGeometry((4, 4, 4))
    h = one_particle(geometry)
    rows = []
    for entry in node_table():
        bits = tuple(entry['pi_coordinates'])
        signs = tuple((-1)**b for b in bits)
        chirality = -np.prod(signs).item()
        if chirality != entry['jacobian_chirality']:
            raise ArithmeticError('source chirality convention changed')
        frame = tuple(-chirality*z for z in signs)
        rotation = s.eye(2) if frame == (1, 1, 1) else s.I*PAULI[frame.index(1)]
        if rotation.conjugate().T*rotation != s.eye(2) or rotation.det() != 1:
            raise ArithmeticError('valley rotation is not SU2')
        if any(rotation*PAULI[j]*rotation.conjugate().T != frame[j]*PAULI[j] for j in range(3)):
            raise ArithmeticError('valley frame lacks its exact SU2 lift')
        spin = tuple(frame[j]*PAULI[j]/2 for j in range(3))
        for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            if spin[i]*spin[j]-spin[j]*spin[i] != s.I*spin[k]:
                raise ArithmeticError('valley spin frame is not a proper SU2 frame')
        if any(2*chirality*spin[j] != -signs[j]*PAULI[j] for j in range(3)):
            raise ArithmeticError('valley helicity frame does not match source')
        rows.append({'pi_coordinates': list(bits), 'linear_symbol_axis_signs': [-z for z in signs],
                     'jacobian_chirality': chirality, 'proper_spin_frame_signs': list(frame),
                     'positive_energy_helicity': str(F(chirality, 2)),
                     'negative_energy_helicity': str(F(-chirality, 2)),
                     'continuum_mass': '0'})
    # Exact 4-point DFT identity implies normalization of its 3D tensor product.
    dft = s.Matrix(4, 4, lambda x, k: s.I**((x*k) % 4)/2)
    if dft.conjugate().T*dft != s.eye(4):
        raise ArithmeticError('exact discrete Fourier normalization failed')
    sine = (0, 1, 0, -1)
    for labels in product(range(4), repeat=3):
        orbital = [s.I**(sum(x*k for x, k in zip(site, labels, strict=True)) % 4)/8
                   for site in geometry.sites]
        embed = s.SparseMatrix(128, 2, {(2*x+spin, spin): z
                                      for x, z in enumerate(orbital) for spin in range(2)})
        symbol = -s.Rational(geometry.speed)*sum((sine[labels[j]]*PAULI[j] for j in range(3)), s.zeros(2))
        residual = h*embed-embed*symbol
        if any(s.expand(z) != 0 for z in residual):
            raise ArithmeticError('source hopping does not match original-H Fourier symbol')
    return {'valleys': rows, 'positive_chirality_count': 4, 'negative_chirality_count': 4,
            'finite_torus_sides': [4, 4, 4], 'CAR_modes': 128,
            'exact_Fourier_columns_checked': 128,
            'all_discrete_momenta_retained': True,
            'original_H_no_scalar_shift_or_time_reversal': True,
            'projector_at_exact_node': 'undefined; rank-two energy degeneracy retained'}


def finite_contact_check():
    """All eight zero nodes on the actual 4^3 source; complete N2 propagation.

    Integer identities prove the contact norm and escape. The time evolution
    is an independently labeled floating diagnostic on all 8128 N2 states.
    No native source-word execution is asserted by this check.
    """
    geometry = PeriodicGeometry((4, 4, 4))
    exact_h = one_particle(geometry)
    h = np.asarray(exact_h, dtype=complex)
    real = np.asarray(s.re(16*exact_h), dtype=np.int64)
    imag = np.asarray(s.im(16*exact_h), dtype=np.int64)
    same_parity = np.array([[all((a-b) % 2 == 0 for a, b in zip(x, y, strict=True))
                            for y in geometry.sites] for x in geometry.sites], dtype=np.int64)
    projector_times_eight = np.kron(same_parity, np.eye(2, dtype=np.int64))
    if not np.array_equal(projector_times_eight@projector_times_eight, 8*projector_times_eight):
        raise ArithmeticError('eight-valley projection lost idempotence')
    up = np.array([int(spin == 0 and all(x % 2 == 0 for x in site))
                   for site in geometry.sites for spin in range(2)], dtype=np.int64)
    down = np.array([int(spin == 1 and all(x % 2 == 0 for x in site))
                     for site in geometry.sites for spin in range(2)], dtype=np.int64)
    initial_times_eight = np.outer(up, down)-np.outer(down, up)
    if any(np.any(part@initial_times_eight+initial_times_eight@part.T) for part in (real, imag)):
        raise ArithmeticError('eight-valley initial pair is not a free zero mode')
    contact_times_eight = np.zeros_like(initial_times_eight)
    for site in range(64):
        contact_times_eight[2*site, 2*site+1] = initial_times_eight[2*site, 2*site+1]
        contact_times_eight[2*site+1, 2*site] = initial_times_eight[2*site+1, 2*site]
    projected = projector_times_eight@contact_times_eight@projector_times_eight.T
    if not np.array_equal(projected, 8*initial_times_eight):
        raise ArithmeticError('exact contact principal angle is not 1/8')
    if int(np.sum(initial_times_eight**2)) != 128 or int(np.sum(contact_times_eight**2)) != 16:
        raise ArithmeticError('exact initial or contact norm changed')
    pair = TwoFermions(128)
    initial = initial_times_eight[pair.first, pair.second].astype(complex)/8
    kinetic = pair.lift_one_body(sparse.csr_matrix(h))
    diagonal = np.zeros(pair.dimension)
    for site in range(64):
        diagonal[pair.index(2*site, 2*site+1)] = 1
    H = kinetic+sparse.diags(diagonal)
    duration = F(1, 8)
    evolved = expm_multiply(-1j*float(duration)*H, initial)
    difference = float(np.linalg.norm(evolved-initial))
    bound = float(duration*sqrt_upper(F(1, 8)))
    if not 1e-4 < difference <= bound+1e-14:
        raise ArithmeticError('full outgoing evolution violates contact comparison')
    # Negative original-H energy exists even when U is arbitrarily large:
    # at k=(pi/2,0,0), the +sigma_x spinor has energy -c on N1.
    return {'sides': [4, 4, 4], 'all_N2_states_retained': pair.dimension,
            'initial_spatial_band_modes': 8, 'exact_band_fraction': '1/8',
            'initial_state_norm_squared': '1', 'contact_norm_squared': '1/8',
            'outside_band_contact_norm_squared': '7/64',
            'contact_bound_is_sharp_on_this_state': True,
            'initial_free_energy': '0', 'duration': str(duration),
            'full_N2_evolution_error_diagnostic': difference,
            'rigorous_vector_error_upper': str(duration*sqrt_upper(F(1, 8))),
            'normalization_error_diagnostic': float(abs(np.vdot(evolved, evolved)-1)),
            'floating_diagnostic_tolerance': 1e-12,
            'vacuum_original_energy': '0', 'one_particle_original_energy_below_vacuum': '-1/8',
            'contact_on_N0_and_N1': 'exactly zero for every U',
            'vacuum_lowering_changes_number_and_parity': True,
            'compiled_native_word_executed': False}
