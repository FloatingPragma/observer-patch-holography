"""One complete geometry mode, all native levels and a three-dimensional profile."""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import numpy as np
from scipy.fft import dstn, idstn
from scipy.linalg import eigh
from scipy.sparse import diags, eye, kron
from model.openwave_energy_lapse.source import native
from model.finite_area_response.response import Interval, exponential
from model.openwave_su2_spectral.spectrum import sqrt_upper

QSIZE = F(1, 8)
OMEGA = F(1)


def root_interval(x):
    if x.lo < 0:
        raise ValueError('positive root required')
    return Interval(max(F(0), sqrt_upper(x.lo)-F(1, 2**70)), sqrt_upper(x.hi))


def sinhcosh(x):
    a, b = exponential(F(x)), exponential(-F(x))
    return (a-b)/2, (a+b)/2


def level_values():
    return np.array([float((v.lo+v.hi)/2) for v in native()[2]])


@lru_cache(None)
def profile(n=7):
    if n < 5 or n % 2 != 1:
        raise ValueError('odd three-dimensional Dirichlet cube required')
    one = diags([-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)], [-1, 0, 1], dtype=float)
    identity = eye(n)
    L = (kron(kron(one, identity), identity)+kron(kron(identity, one), identity)+kron(kron(identity, identity), one)).tocsr()
    delta = np.zeros((n, n, n)); delta[(n//2,)*3] = 1
    lam = 2-2*np.cos(np.pi*np.arange(1, n+1)/(n+1))
    green = idstn(dstn(delta, type=1, norm='ortho')/(lam[:, None, None]+lam[None, :, None]+lam[None, None, :]), type=1, norm='ortho')
    if n == 7:
        fixture = json.loads(Path(__file__).with_name('profile_fixture.json').read_text())
        assert fixture['side'] == n and fixture['denominator'] == 2**40
        accepted = np.array(fixture['integer_coordinates'], dtype=np.int64)
        assert accepted.shape == (n**3,)
    else:
        accepted = np.rint(green/green[(n//2,)*3]*2**40).astype(np.int64).ravel()
    assert accepted[np.ravel_multi_index((n//2,)*3, (n,)*3)] == 2**40
    assert accepted.min() > 0
    # Integer action gives exact residuals and the declared projected stiffness.
    shaped = accepted.reshape((n,)*3)
    lap = 6*shaped.copy()
    for j in range(3):
        a = [slice(None)]*3; b = [slice(None)]*3
        a[j] = slice(1, None); b[j] = slice(None, -1)
        lap[tuple(a)] -= shaped[tuple(b)]; lap[tuple(b)] -= shaped[tuple(a)]
    stiffness = F(sum(int(x)*int(y) for x, y in zip(accepted, lap.ravel(), strict=True)), 2**80)
    mask = lap.copy(); mask[(n//2,)*3] = 0
    harmonic_error = F(int(np.max(abs(mask))), 2**40)*F((n+1)**2, 8)
    assert harmonic_error < F(1, 10**8)
    f = accepted/2**40
    return {'n': n, 'L': L, 'F': f, 'F_integer': accepted,
            'stiffness_exact': stiffness, 'harmonic_profile_error_upper': harmonic_error,
            'kappa_exact': 1/(2*QSIZE*QSIZE*stiffness),
            'source_site': (n//2,)*3}


def geometry():
    q = float(QSIZE)
    return {'Q': np.diag([q, -q]), 'Hg': float(OMEGA)/2*np.array([[1., -1.], [-1., 1.]]),
            'vacuum': np.ones(2)/np.sqrt(2), 'coordinate_values': [q, -q]}


def rest_block(source_level, probe_level=0, f=F(0)):
    e = level_values(); g = geometry(); q = float(QSIZE)
    return g['Hg']+np.diag(e[source_level]*np.exp(np.array([q, -q]))+
                               e[probe_level]*np.exp(float(f)*np.array([q, -q])))


def ground_interval(source_level, probe_level=0, f=F(0)):
    a, c = sinhcosh(QSIZE); b, d = sinhcosh(QSIZE*F(f))
    source, probe = native()[2][source_level], native()[2][probe_level]
    z = source*a+probe*b
    return Interval.point(OMEGA/2)+source*c+probe*d-root_interval(Interval.point(OMEGA*OMEGA/4)+z*z)


def mutual_interval(source_level, probe_level, f):
    return ground_interval(source_level, probe_level, f)-ground_interval(source_level, 0, f)-ground_interval(0, probe_level, f)


def exact_attraction():
    shape = profile(); n = shape['n']; center = n//2
    rows = []
    for distance in (1, 2, 3):
        site = (center+distance, center, center)
        index = np.ravel_multi_index(site, (n,)*3)
        f = F(int(shape['F_integer'][index]), 2**40)
        cells = []
        for i in range(1, 8):
            for j in range(1, 8):
                value = mutual_interval(i, j, f)
                assert value.hi < 0
                cells.append({'source_level': i, 'probe_level': j,
                              'mutual_energy_interval_exact': value.encoded()})
        rows.append({'distance': distance, 'profile_exact': str(f), 'all_49_native_pairs': cells})
    first = [r['all_49_native_pairs'][0] for r in rows]
    force = (Interval(*map(F, first[2]['mutual_energy_interval_exact']))-
             Interval(*map(F, first[0]['mutual_energy_interval_exact'])))/2
    assert force.lo > 0
    return {'separations': rows, 'first_native_pair_inward_central_force_interval_exact': force.encoded(),
            'all_147_positive_energy_pair_interactions_attractive': True,
            'mixed_derivative': '-(omega/2)^2 sinh(q) sinh(qF) / [(omega/2)^2+(Es sinh(q)+Et sinh(qF))^2]^(3/2)'}


def autonomous_records():
    f = F(int(profile()['F_integer'][np.ravel_multi_index((4, 3, 3), (7,)*3)]), 2**40)
    g = geometry(); rows = []
    # A classical source/reference flag retains both source interventions.
    for source_level in (0, 1):
        for probe_level in (0, 1, 3, 7):
            H = rest_block(source_level, probe_level, f); e, V = eigh(H)
            initial = g['vacuum']; before = float(initial@H@initial)
            for time in (0., 1., 2., 4.):
                psi = V@(np.exp(-1j*time*e)*(V.T@initial))
                rate_geo = float(np.vdot(psi, 1j*(H@g['Hg']-g['Hg']@H)@psi).real)
                matter = H-g['Hg']
                rate_matter = float(np.vdot(psi, 1j*(H@matter-matter@H)@psi).real)
                assert abs(rate_geo+rate_matter) < 1e-14
                assert abs(float(np.vdot(psi, H@psi).real)-before) < 1e-13
                p = abs(psi)**2
                post = float(p@np.diag(H))
                rows.append({'source_level': source_level, 'probe_level': probe_level, 'source_parameter': time,
                             'complete_geometry_coordinate_probabilities': p.tolist(),
                             'other_63_native_pair_sectors_probability': 0.,
                             'mean_geometry_coordinate': float(p@np.diag(g['Q'])),
                             'source_energy': before, 'geometry_energy': float(np.vdot(psi, g['Hg']@psi).real),
                             'geometry_and_matter_energy_rates': [rate_geo, rate_matter],
                             'complete_geometry_measurement_work': post-before})
    return rows


def run():
    p = profile()
    return {'schema': 'oph.autonomous_universal_positive_energy_geometry.v1',
            'verdict': 'PASS_ALL_NATIVE_ENERGIES_ATTRACT_ON_ONE_AUTONOMOUS_SOURCE',
            'construction': {'native_dimension_per_matter_site': 8, 'geometry_dimension': 2,
                'full_two_matter_geometry_dimension': 128, 'idle_reference_allowed': True,
                'native_reference': 'K=M0-E0 I, retaining all8 native levels',
                'geometry': 'Q=q sigma_z; Hg=omega/2 (I-sigma_x)',
                'q_exact': str(QSIZE), 'omega_exact': str(OMEGA),
                'universal_local_coupling': 'exp(F_x Q) K_x',
                'uniform_lapse_derivative_identity': 'd/ds [exp(s+F_x Q)K_x] at s0 = exp(F_x Q)K_x >=0',
                'spatial_mode': 'one declared positive three-dimensional harmonic-profile mode',
                'projected_geometry_inertia_equals_stiffness': True,
                'geometry_ground_preparation_is_separate_from_quenched_evolution': True},
            'native_gap_intervals': [x.encoded() for x in native()[2]],
            'profile': {'n': p['n'], 'all343_integer_coordinates': p['F_integer'].tolist(),
                'denominator': 2**40, 'stiffness_exact': str(p['stiffness_exact']),
                'kappa_exact': str(p['kappa_exact']),
                'harmonic_profile_infinity_error_upper_exact': str(p['harmonic_profile_error_upper'])},
            'exact_attraction': exact_attraction(), 'autonomous_complete_records': autonomous_records(),
            'weak_source_Newton_coefficient': '1/(4 pi kappa), with lattice/boundary and small-source errors retained',
            'dimensionless_Newton_coefficient_numeric': float(1/(4*np.pi*float(p['kappa_exact']))),
            'physical_calibration_rule': 'G_natural = length_unit/(4 pi kappa energy_unit); fitting kappa to G is an empirical calibration',
            'scope': {'opposite_charge_reciprocal_source_equivalence': False,
                'universal_native_ground_referenced_energy_coupling': True,
                'all_ground_referenced_positive_native_energies_have_positive_lapse_charge': True,
                'self_energies_subtracted_in_mutual_force': True,
                'every_local_kinetic_stress_density_positive': False,
                'physical_G_or_Einstein_dynamics_derived': False,
                'full_spatial_quantum_geometry_carrier': False,
                'continuum_universal_gravity_proved': False}}
