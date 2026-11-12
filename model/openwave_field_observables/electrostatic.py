"""Coulomb shape and force from the cochain source's Gauss minimizer."""
import numpy as np
from . import lattice as l


def pair_energy(green, displacement, coupling=1.0):
    if not np.isfinite(coupling) or coupling <= 0:
        raise ValueError('positive declared coupling required')
    return float(coupling**2*(green[0, 0, 0]-green[tuple(displacement)]))


def radial_force(green, position):
    """Central difference of the source Green function, not a 1/r ansatz."""
    n = len(green)
    position = np.asarray(position)
    radius = np.linalg.norm(position)
    force = []
    for j in range(3):
        left, right = position.copy(), position.copy()
        left[j] -= 1; right[j] += 1
        force.append((green[tuple(left % n)]-green[tuple(right % n)])/2)
    return float(np.dot(force, position/radius))


def run():
    small = l.green(12)
    rho = np.zeros_like(small); rho[0, 0, 0] = 1; rho[3, 0, 0] = -1
    # Translation superposes the actual +1 and -1 source Green functions.
    phi = small-np.roll(small, 3, axis=0)
    electric = l.gradient(phi)
    gauss_error = float(np.max(np.abs(l.divergence_adjoint(electric)-rho)))
    field_energy = float(np.vdot(electric, electric).real/2)
    green_energy = float(np.vdot(rho, phi).real/2)
    assert gauss_error < 1e-12
    assert abs(field_energy-green_energy) < 1e-12
    assert abs(green_energy-pair_energy(small, (3, 0, 0))) < 1e-12
    rng = np.random.default_rng(897)
    perturb = l.curl_adjoint(rng.normal(size=(12, 12, 12, 3)))
    orthogonality = float(abs(np.vdot(electric, perturb)))
    assert orthogonality < 1e-12
    assert np.linalg.norm(electric+perturb)**2 > np.linalg.norm(electric)**2
    rows = []
    for n in (64, 96, 128):
        green = l.green(n)
        for r in (4, 6, 8):
            force = radial_force(green, (r, 0, 0))
            potential_difference = float(green[r, 0, 0]-green[2*r, 0, 0])
            factor = 4*np.pi*r*r*force
            rows.append({'grid': n, 'radius_sites': r,
                         'Green_difference_r_2r': round(potential_difference, 13),
                         'potential_amplitude_ratio': round(float(8*np.pi*r*potential_difference), 10),
                         'radial_force': round(force, 13),
                         'force_amplitude_ratio': round(float(factor), 10)})
        assert all(0.9 < row['force_amplitude_ratio'] < 1.16 for row in rows if row['grid'] == n)
    final = l.green(128)
    force_values = [radial_force(final, (r, 0, 0)) for r in (6, 8, 10, 12)]
    slopes = [-np.log(force_values[i+1]/force_values[i])/np.log((8, 10, 12)[i]/(6, 8, 10)[i])
              for i in range(3)]
    assert all(1.95 < value < 2.15 for value in slopes)
    orientation = []
    for point in ((8, 0, 0), (6, 6, 0), (5, 5, 5)):
        radius = float(np.linalg.norm(point))
        value = 4*np.pi*radius**2*radial_force(final, point)
        orientation.append({'point': list(point), 'radius': round(radius, 9),
                            'force_amplitude_ratio': round(float(value), 10)})
        assert .95 < value < 1.06
    baseline = pair_energy(small, (3, 0, 0), 1)
    doubled = pair_energy(small, (3, 0, 0), 2)
    assert abs(doubled-4*baseline) < 1e-12
    # Rescaling the constitutive operator changes its Gauss response.
    rescaled_phi = phi/(1+.5)
    assert np.max(np.abs(l.divergence_adjoint(l.gradient(rescaled_phi))-rho)) > .3
    assert abs(4*rows[-1]['force_amplitude_ratio']-1) > 2
    return {'task': 'electrostatic_force', 'Gauss_max_error': round(gauss_error, 14),
            'field_Green_energy_error': round(abs(field_energy-green_energy), 14),
            'source_pair_energy_g1': round(baseline, 12),
            'source_pair_energy_g2': round(doubled, 12),
            'potential_force_rows': rows, 'force_log_slopes': [round(float(v), 10) for v in slopes],
            'orientation_controls': orientation,
            'checks': {'same_cochain_Gauss_minimizer': True, 'minimum_field_energy': True,
                       'potential_and_force_without_fit': True, 'near_inverse_square_force': True,
                       'declared_coupling_squared_scaling': True},
            'mutations': {'altered_Gauss_solution_rejected': True, 'hidden_coupling_change_rejected': True},
            'inputs': {'coupling_g': 1.0, 'lattice_spacing': 1.0,
                       'geometry': 'periodic cubic Gauss source with zero total charge'},
            'continuum_coefficient': 'g^2/(4 pi), conditional on the declared electric energy and charge normalization',
            'alpha_prediction': None,
            'scope': 'A force-level computation on the declared Gauss field architecture. The full interacting charged-particle state, observed e and alpha, and source-selected three-dimensional geometry are not identified by this calculation.'}
