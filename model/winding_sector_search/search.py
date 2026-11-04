"""Complete finite one-body spectra for a declared rational-map family.

These diagnostics do not determine interacting particles or continuum limits.
Every finite site, spin row and eigenvalue is retained.
"""
import numpy as np
from scipy import linalg
from scipy.sparse import coo_matrix

from ..bounded_defect_index.finite import sites
from ..native_mass_channels.algebra import SIGMA, rational
from .algebra import winding


def direction(point, degree):
    """Normalize the homogeneous rational map of degree abs(degree)."""
    degree = winding(degree)
    if type(point) is not tuple or len(point) != 3 or any(type(x) is not int for x in point):
        raise ValueError('three exact integer coordinates required')
    if point == (0, 0, 0):
        return np.zeros(3)
    x, y, z = point
    if x == 0 and y == 0:
        return np.array([0., 0., float(np.sign(z)*np.sign(degree))])
    radius = np.linalg.norm(point)
    a = complex(x, y)
    # The second expression avoids cancellation near the north pole.
    b = (x*x+y*y)/(radius+z) if z > 0 else radius-z
    scale = max(abs(a), abs(b))
    a, b = (a/scale)**abs(degree), (b/scale)**abs(degree)
    denominator = abs(a)**2+b*b
    result = np.array([2*(a*b).real, 2*(a*b).imag, abs(a)**2-b*b])/denominator
    return np.sign(degree)*result


def source_matrix(radius, degree, mass=1, spacing=1, velocity=1):
    degree = winding(degree)
    if type(radius) is not int or radius < 1:
        raise ValueError('positive integer cube radius required')
    m, a, v = (rational(x) for x in (mass, spacing, velocity))
    if min(m, a, v) <= 0:
        raise ValueError('positive exact source parameters required')
    points = sites(radius)
    address = {point: j for j, point in enumerate(points)}
    sigma = [np.array(p, complex) for p in SIGMA]
    rows, cols, values = [], [], []

    def add(i, j, block):
        for spin in range(2):
            for other in range(2):
                if block[spin, other] != 0:
                    rows.append(2*i+spin)
                    cols.append(2*j+other)
                    values.append(block[spin, other])

    for point, i in address.items():
        n = direction(point, degree)
        add(i, i, float(m)*sum(n[j]*(-1)**(point[j] % 2)*sigma[j] for j in range(3)))
        for axis in range(3):
            target = tuple(x+int(j == axis) for j, x in enumerate(point))
            if target in address:
                hop = 1j*float(v/a)*sigma[axis]/2
                add(i, address[target], hop)
                add(address[target], i, hop.conj().T)
    return coo_matrix((values, (rows, cols)), shape=(2*len(points),)*2).tocsr()


def finite_spectrum(radius, degree, mass=1):
    h = source_matrix(radius, degree, mass).toarray()
    energies, columns = linalg.eigh(h, driver='evr')
    residual = np.linalg.norm(h@columns-columns*energies, 'fro')
    orthogonal = np.linalg.norm(columns.conj().T@columns-np.eye(len(energies)), 'fro')
    if max(residual, orthogonal) > 1e-8:
        raise ArithmeticError('finite spectral calculation failed its residual admission')
    points = sites(radius)
    origin = points.index((0, 0, 0))
    # Fixed, declared energy window. Projector traces are invariant under
    # arbitrary basis rotations of a degenerate eigenspace.
    threshold = float(rational(mass))/4
    if np.min(np.abs(np.abs(energies)-threshold)) < 1e-7:
        raise ArithmeticError('search window is numerically unresolved')
    mask = np.abs(energies) < threshold
    density = np.sum(np.abs(columns[:, mask])**2, axis=1).reshape(-1, 2).sum(axis=1)
    boundary = np.array([max(map(abs, point)) == radius for point in points])
    distance2 = np.array([sum(x*x for x in point) for point in points])
    count = int(np.count_nonzero(mask))
    zero_like = int(np.count_nonzero(np.abs(energies) < 1e-8))
    if degree % 2 and int(mass) == 4:
        if zero_like != 2 or np.min(np.abs(energies[np.abs(energies) >= 1e-8])) < 2-1e-8:
            raise ArithmeticError('the finite control contradicts the analytic strong-mass theorem')
    return {
        'radius': radius, 'degree': degree, 'mass': str(mass), 'spacing': '1', 'velocity': '1',
        'dimension': len(energies), 'all_eigenvalues': [float(x) for x in energies],
        'complete_eigenbasis_residual_admission_upper': 1e-8,
        'numerical_zero_threshold': 1e-8, 'numerical_zero_count': zero_like,
        'window_abs_energy_below': threshold, 'window_rank': count,
        'window_origin_projector_trace': float(density[origin]),
        'window_boundary_projector_trace': float(density[boundary].sum()),
        'window_total_radius_squared_trace': float(density@distance2),
        'window_states_are_certified_infinite_volume_bound_particles': False,
        'raw_roundoff_used_as_analytic_error_bound': False,
    }


def search_report():
    cases = [finite_spectrum(radius, degree, mass)
             for mass in (1, 4) for degree in (1, 2, 3, 4) for radius in (1, 2, 3)]
    followup = [finite_spectrum(radius, degree, 1)
                for degree in (1, 3) for radius in (4, 5)]
    return {
        'declared_search_grid': {'mass': [1, 4], 'degree': [1, 2, 3, 4], 'cube_radius': [1, 2, 3]},
        'initial_case_count': len(cases), 'all_initial_cases': cases,
        'volume_followup': {
            'reason': 'initial degree-three window counts varied with volume; compare degree one and three at larger radii',
            'chosen_after_initial_search': True,
            'grid': {'mass': [1], 'degree': [1, 3], 'cube_radius': [4, 5]},
            'all_cases': followup,
        },
        'completed_case_count': len(cases)+len(followup),
        'negative_degree_coverage': 'H[-d]=-Ds H[d] Ds; particle-hole symmetry gives the same spectral multiset',
        'continuum_index_controls_lattice_zero_count_at_fixed_spacing': False,
        'finite_matrix_spectra_are_numerical_diagnostics': True,
        'source_selected_degree_or_mass': False,
        'observed_species_or_mass_ratios_assigned': False,
    }
