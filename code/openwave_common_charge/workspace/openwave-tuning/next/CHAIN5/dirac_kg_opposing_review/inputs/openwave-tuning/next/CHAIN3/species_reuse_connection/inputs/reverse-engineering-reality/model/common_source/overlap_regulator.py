"""Finite Euclidean overlap kernel on the retained M15 matter representation.

Neuberger's operator and Luescher's modified chiral projectors are applied
to explicit four-dimensional links. This is a finite regulator and native
query construction, not a chiral fermion measure or Lorentzian field theory.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from tools import strict_json
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('overlap_regulator_results.json')
MATTER = 'model/observer_particles/results/representation_identification_operators.json'
MATTER_OWNER = 'model/observer_particles/results/representation_identification.json'
SPLITS = 8192
COLOR_RATE = F(1, 10000)
PHASE_RATE = F(1, 20000)
GENERATOR_BOUND = 2*COLOR_RATE+24*PHASE_RATE
PLAQUETTE_BOUND = 16*GENERATOR_BOUND


def gamma_matrices():
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.]).astype(complex)
    gamma = tuple(np.kron(y, s) for s in (x, y, z))+(np.kron(x, np.eye(2)),)
    return gamma, np.kron(z, np.eye(2))


def _matrix(value, dimension, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(name+' cannot be masked')
    raw = np.asarray(value)
    if raw.shape != (dimension, dimension) or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all():
        raise ValueError(name+' must be a finite square matrix of the declared size')
    return np.asarray(raw, complex)


def _admit(value, ceiling, name):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= ceiling:
        raise ArithmeticError(name+' exceeds its numerical admission')
    return ceiling


@lru_cache(maxsize=8)
def _matter_basis(raw, owner_raw):
    owner = strict_json.loads(owner_raw.decode())
    if hashlib.sha256(raw).hexdigest() != owner['artifacts'][Path(MATTER).name]:
        raise ValueError('retained matter representation custody changed')
    document = strict_json.loads(raw.decode())
    basis = tuple(tuple(document['basis_subsets'][i]) for i in document['even_indices'])
    if len(basis) != 15 or {len(s) for s in basis} != {2, 4}:
        raise ValueError('the retained nonvacuum even exterior module is required')
    weights = (-2, -2, -2, 3, 3)
    charges = tuple(sum(weights[i] for i in s) for s in basis)
    q = document['even_identification_operators']['q']
    entries = {(i, j): int(v) for i, j, v in q['entries']}
    if q['shape'] != [15, 15] or entries != {(i, i): v for i, v in enumerate(charges)}:
        raise ValueError('overlap representation does not match the retained hypercharges')
    return basis


def matter_basis():
    # Hash validation is repeated even after warming the decoded-basis cache.
    return _matter_basis((ROOT/MATTER).read_bytes(), (ROOT/MATTER_OWNER).read_bytes())


def exterior_link(unitary):
    """The actual Lambda²(C5) plus Lambda⁴(C5), in retained basis order."""
    u = _matrix(unitary, 5, 'fundamental link')
    _admit(np.linalg.norm(u.conj().T@u-np.eye(5), 'fro'), 2e-11, 'fundamental unitarity')
    _admit(abs(np.linalg.det(u)-1), 2e-11, 'fundamental determinant')
    basis = matter_basis()
    return np.array([[np.linalg.det(u[np.ix_(a, b)]) if len(a) == len(b) else 0
                      for b in basis] for a in basis])


def _sites(side):
    if type(side) is not int or not 2 <= side <= 4:
        raise ValueError('a declared periodic side from two through four required')
    return tuple(product(range(side), repeat=4))


def sample_links(side=2):
    """A fixed smooth SM-subgroup background; no measured values enter."""
    sites = _sites(side)
    x = np.array([[0, 1], [1, 0]], complex)
    z = np.diag([1., -1.])
    t3 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], complex)
    s3 = np.diag([1., -1., 0.])
    links = []
    for point in sites:
        rows = []
        for mu in range(4):
            angle = float(COLOR_RATE)*(1+point[(mu+1) % 4] % 2)
            color = expm(1j*angle*(t3 if mu % 2 else s3))
            weak = expm(1j*angle*(x if mu % 2 else z))
            phase = np.exp(1j*float(PHASE_RATE)*(mu+1)*(1+point[mu] % 2))
            fundamental = np.zeros((5, 5), complex)
            fundamental[:3, :3] = phase**-2*color
            fundamental[3:, 3:] = phase**3*weak
            rows.append(exterior_link(fundamental))
        links.append(rows)
    return np.asarray(links)


def plaquette_defect(links, side=2):
    """Frobenius upper diagnostic for every represented plaquette norm."""
    sites = _sites(side)
    positions = {p: i for i, p in enumerate(sites)}
    colors = np.asarray(links).shape[2]
    maximum = 0.
    for i, site in enumerate(sites):
        for mu in range(4):
            next_mu = positions[tuple((v+1) % side if j == mu else v for j, v in enumerate(site))]
            for nu in range(mu+1, 4):
                next_nu = positions[tuple((v+1) % side if j == nu else v for j, v in enumerate(site))]
                loop = links[i, mu]@links[next_mu, nu]@links[next_nu, mu].conj().T@links[i, nu].conj().T
                maximum = max(maximum, float(np.linalg.norm(loop-np.eye(colors), 'fro')))
    return maximum


def wilson_operator(links, side=2):
    """Dimensionless a*D_W, with all four directions and both transports."""
    sites = _sites(side)
    raw = np.asarray(links)
    if (np.ma.isMaskedArray(links) or raw.ndim != 4 or raw.shape[:2] != (len(sites), 4)
            or raw.shape[2] != raw.shape[3] or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all()):
        raise ValueError('complete finite four-direction link matrices required')
    colors = raw.shape[2]
    for u in raw.reshape(-1, colors, colors):
        _admit(np.linalg.norm(u.conj().T@u-np.eye(colors), 'fro'), 2e-10, 'link unitarity')
    positions = {site: i for i, site in enumerate(sites)}
    gamma, g5 = gamma_matrices()
    block = 4*colors
    operator = 4*np.eye(len(sites)*block, dtype=complex)
    for index, site in enumerate(sites):
        left = slice(index*block, (index+1)*block)
        for mu in range(4):
            neighbor = tuple((v+1) % side if j == mu else v for j, v in enumerate(site))
            right_index = positions[neighbor]
            right = slice(right_index*block, (right_index+1)*block)
            forward = -.5*np.kron(np.eye(4)-gamma[mu], raw[index, mu])
            backward = -.5*np.kron(np.eye(4)+gamma[mu], raw[index, mu].conj().T)
            operator[left, right] += forward
            operator[right, left] += backward
    chirality = np.kron(np.eye(len(sites)), np.kron(g5, np.eye(colors)))
    return operator, chirality


def overlap(wilson, chirality):
    """a*D=I+gamma5 sign(gamma5(a*D_W-I)), rho=1."""
    n = len(wilson)
    w = _matrix(wilson, n, 'Wilson operator')
    g5 = _matrix(chirality, n, 'chirality')
    _admit(np.linalg.norm(g5-g5.conj().T, 'fro'), 2e-10, 'Hermitian chirality')
    _admit(np.linalg.norm(g5@g5-np.eye(n), 'fro'), 2e-10, 'involutive chirality')
    h = g5@(w-np.eye(n))
    _admit(np.linalg.norm(h-h.conj().T, 'fro'), 2e-10, 'Wilson Hermiticity')
    values, vectors = np.linalg.eigh(h)
    gap = float(np.min(abs(values)))
    if not math.isfinite(gap) or gap <= .1:
        raise ValueError('Wilson sign kernel has no admitted spectral gap')
    sign = (vectors*np.where(values > 0, 1., -1.))@vectors.conj().T
    d = np.eye(n)+g5@sign
    hat = -sign
    return d, hat, gap


def gauge_transform(links, gauges, side=2):
    sites = _sites(side)
    positions = {p: i for i, p in enumerate(sites)}
    raw_links = np.asarray(links)
    if (np.ma.isMaskedArray(links) or raw_links.ndim != 4
            or raw_links.shape[:2] != (len(sites), 4)
            or raw_links.shape[2] != raw_links.shape[3]
            or raw_links.dtype.kind not in 'iufc' or not np.isfinite(raw_links).all()):
        raise ValueError('complete finite four-direction link matrices required')
    if np.ma.isMaskedArray(gauges):
        raise ValueError('one finite gauge matrix at every site required')
    gauges = np.asarray(gauges)
    colors = raw_links.shape[2]
    if (np.ma.isMaskedArray(gauges) or gauges.shape != (len(sites), colors, colors)
            or gauges.dtype.kind not in 'iufc' or not np.isfinite(gauges).all()):
        raise ValueError('one finite gauge matrix at every site required')
    for gauge in gauges:
        _admit(np.linalg.norm(gauge.conj().T@gauge-np.eye(colors), 'fro'),
               2e-10, 'gauge unitarity')
    result = np.empty(raw_links.shape, dtype=complex)
    for i, site in enumerate(sites):
        for mu in range(4):
            neighbor = tuple((v+1) % side if j == mu else v for j, v in enumerate(site))
            result[i, mu] = gauges[i]@links[i, mu]@gauges[positions[neighbor]].conj().T
    whole = np.zeros((len(sites)*4*colors,)*2, complex)
    for i, gauge in enumerate(gauges):
        span = slice(i*4*colors, (i+1)*4*colors)
        whole[span, span] = np.kron(np.eye(4), gauge)
    return result, whole


def local_overlap_polynomial(wilson, chirality, degree=192, *, links, side=2):
    """Finite-range approximation with an exact ideal-arithmetic tail.

    The checked links must assemble exactly to the supplied nearest-neighbor
    Wilson operator and on-site chirality, with represented plaquettes below
    1/40. The fixed example has fundamental generator norm <=7/5000,
    exterior-link norm change <=7/1250 and plaquette norm <=14/625.
    Hernandez--Jansen--Luescher, eqs2.14--2.15, then gives
    H²>=1/4 at rho=1; their bound ||H||<=8 supplies H²<=64.
    The Legendre construction has t=15/17,kappa=4/17, hence
    ||D-D_degree||<=16*(15/17)**(degree+1), with range<=2*degree+1
    lattice links. No gauge-volume or refinement factor multiplies the tail.
    Numerical evaluation error is separate from this polynomial theorem.
    """
    if type(degree) is not int or not 0 <= degree <= 1024:
        raise ValueError('an integer degree from zero through1024 required')
    n = len(wilson)
    w = _matrix(wilson, n, 'Wilson operator')
    g5 = _matrix(chirality, n, 'chirality')
    local_w, local_g5 = wilson_operator(links, side)
    if local_w.shape != w.shape:
        raise ValueError('the supplied Wilson operator is not built from these local links')
    if not np.array_equal(w, local_w):
        raise ValueError('nearest-neighbor Wilson assembly must match exactly')
    if not np.array_equal(g5, local_g5):
        raise ValueError('on-site chirality must match exactly')
    _admit(plaquette_defect(links, side), 1/40, 'represented plaquette admissibility')
    _admit(np.linalg.norm(g5-g5.conj().T, 'fro'), 2e-10, 'Hermitian chirality')
    _admit(np.linalg.norm(g5@g5-np.eye(n), 'fro'), 2e-10, 'involutive chirality')
    h = g5@(w-np.eye(n))
    _admit(np.linalg.norm(h-h.conj().T, 'fro'), 2e-10, 'polynomial Hermiticity')
    values = np.linalg.eigvalsh(h)
    if min(values*values) < .25-1e-11 or max(values*values) > 64+1e-11:
        raise ValueError('operator lies outside the declared uniform spectral interval')
    z = (float(F(257, 4))*np.eye(n)-2*h@h)/float(F(255, 4))
    t, kappa = F(15, 17), F(4, 17)
    previous = np.eye(n, dtype=complex)
    inverse_root = float(kappa)*previous
    if degree:
        current = np.asarray(z, complex)
        inverse_root += float(kappa*t)*current
        for index in range(1, degree):
            following = ((2*index+1)*z@current-index*previous)/(index+1)
            inverse_root += float(kappa*t**(index+1))*following
            previous, current = current, following
    approximate_sign = h@inverse_root
    approximation = np.eye(n)+g5@approximate_sign
    tail = 16*t**(degree+1)
    return approximation, {'degree': degree, 'maximum_range_in_lattice_links': 2*degree+1,
                           'uniform_operator_tail_exact': str(tail),
                           'Ginsparg_Wilson_defect_bound_exact': str(tail*(2+tail)),
                           'floating_evaluator_error_included': False}


def free_symbol(momentum, spacing=1.):
    raw = np.asarray(momentum)
    if (np.ma.isMaskedArray(momentum) or raw.shape != (4,) or raw.dtype.kind not in 'iuf'
            or not np.isfinite(raw).all() or type(spacing) not in (int, float)
            or not math.isfinite(spacing) or spacing <= 0):
        raise ValueError('four finite Euclidean momenta and positive spacing required')
    gamma, _ = gamma_matrices()
    k = spacing*raw
    sine = np.sin(k)
    b = float(np.sum(1-np.cos(k))-1)
    radius = math.sqrt(b*b+float(sine@sine))
    v = (b*np.eye(4)+1j*sum((s*g for s, g in zip(sine, gamma, strict=True)), np.zeros((4, 4), complex)))/radius
    return (np.eye(4)+v)/spacing, v


def free_continuum_bound(momentum, spacing):
    """Analytic ||D(p)-i gamma.p|| <= a*|p|² for a*|p|<=1/4.

    With u_mu=1-cos(k_mu), r²=1+2 sum_(mu<nu)u_mu*u_nu>=1.
    Taylor remainder and r-1<=|k|^4/8 give error at most
    |k|²/2+|k|³/6+|k|^4/8+|k|^5/8 < |k|² on this band.
    This controls only the free Euclidean symbol, not an interacting limit.
    """
    d, _ = free_symbol(momentum, spacing)
    norm = float(np.linalg.norm(momentum))
    if spacing*norm > .25:
        raise ValueError('the declared continuum comparison band was exceeded')
    gamma, _ = gamma_matrices()
    target = 1j*sum((p*g for p, g in zip(momentum, gamma, strict=True)), np.zeros((4, 4), complex))
    bound = spacing*norm*norm
    _admit(np.linalg.norm(d-target, 2), bound+1e-13, 'free continuum comparison')
    return bound


def native_free_query(momentum, spacing=.125):
    """Complete native Hadamard query of aD/2, including its failure output.

    This is a four-spinor kernel query, not the full many-fermion field.
    The proper-code load/export and read services retain their existing
    accuracy assumptions. No success postselection changes the output.
    """
    d, v = free_symbol(momentum, spacing)
    gamma, _ = gamma_matrices()
    sine = np.sin(spacing*np.asarray(momentum))
    length = float(np.linalg.norm(sine))
    b = float(np.sum(1-np.cos(spacing*np.asarray(momentum)))-1)
    # Little-endian tensor Pauli strings for gamma1..gamma4.
    labels = ('XY', 'YY', 'ZY', 'IX')
    if length == 0:
        angle = 0. if b > 0 else math.pi
        step = compiler.compile_diagonal([0., 0., 0., 0., angle, angle, angle, angle])
        count = 1
    else:
        angle = math.atan2(length, b)
        coefficients = angle*sine/length
        step = []
        schedule = [(j, .5) for j in range(3)]+[(3, 1.)]+[(j, .5) for j in reversed(range(3))]
        for j, fraction in schedule:
            a = float(coefficients[j])*fraction/SPLITS
            step.extend(compiler.pauli_rotation(labels[j]+'I', a/2))
            step.extend(compiler.pauli_rotation(labels[j]+'Z', -a/2))
        count = SPLITS
    executed_step = compiler.execute(step, 3)
    middle = np.linalg.matrix_power(executed_step, count)
    hword = [compiler.pulse('H', math.pi/2, [2]), compiler.pulse('I', -math.pi/2, [2])]
    h = compiler.execute(hword, 3)
    actual = h@middle@h
    target_h = np.kron(np.array([[1, 1], [1, -1]])/math.sqrt(2), np.eye(4))
    target = target_h@np.block([[np.eye(4), np.zeros((4, 4))], [np.zeros((4, 4)), v]])@target_h
    ideal_step = np.eye(8, dtype=complex)
    if length:
        for j, fraction in schedule:
            factor = expm(1j*float(coefficients[j])*fraction/SPLITS*gamma[j])
            ideal_step = np.block([[np.eye(4), np.zeros((4, 4))], [np.zeros((4, 4)), factor]])@ideal_step
    else:
        ideal_step = np.block([[np.eye(4), np.zeros((4, 4))], [np.zeros((4, 4)), v]])
    # Local Taylor tails plus a telescoping product bound. pi<22/7 and
    # exp(2*pi/N)<2 give this rational bound for the ideal Strang word.
    truncation = F(0) if count == 1 else 2*F(44, 7)**3/(3*SPLITS**2)
    numerical = math.expm1(count*math.log1p(float(np.linalg.norm(executed_step-ideal_step, 2)))
                          +2*math.log1p(float(np.linalg.norm(h-target_h, 2))))
    _admit(np.linalg.norm(actual-target, 2), float(truncation)+numerical+1e-11, 'complete native query')
    _admit(np.linalg.norm(actual[:4, :4]-spacing*d/2, 2), float(truncation)+numerical+1e-11, 'native overlap block')
    _admit(np.linalg.norm(actual[:, :4].conj().T@actual[:, :4]-np.eye(4), 'fro'), 1e-7, 'complete native output mass')
    columns = [[[float(z.real), float(z.imag)] for z in row] for row in actual[:, :4]]
    return {'source_native_pulses': count*len(step)+2*len(hword),
            'ideal_trotter_operator_bound_exact': str(truncation),
            'source_numerical_operator_admission': _admit(numerical, 2e-7, 'source execution'),
            'both_quantum_outcomes_and_all_input_columns_checked': True,
            'complete_success_and_failure_reference_columns': columns,
            'kernel_query_is_a_full_chiral_field_simulation': False}, actual[:, :4]


def run():
    gamma, g5 = gamma_matrices()
    for a, left in enumerate(gamma):
        for b, right in enumerate(gamma):
            _admit(np.linalg.norm(left@right+right@left-(2*np.eye(4) if a == b else 0), 'fro'), 1e-12, 'Euclidean Clifford algebra')
    basis = matter_basis()
    weights = (-2, -2, -2, 3, 3)
    charges = [sum(weights[i] for i in s) for s in basis]
    if sum(charges) != 0 or sum(q**3 for q in charges) != 0:
        raise ArithmeticError('retained hypercharge anomaly traces changed')
    links = sample_links()
    if not PLAQUETTE_BOUND == F(14, 625) < F(1, 40):
        raise ArithmeticError('the supplied small-link domain changed')
    _admit(plaquette_defect(links), 1/40, 'small represented plaquettes')
    w, whole_g5 = wilson_operator(links)
    d, hat, gap = overlap(w, whole_g5)
    n = len(d)
    left = (np.eye(n)-hat)/2
    right = (np.eye(n)+whole_g5)/2
    admissions = {}
    for name, residual in {
        'Ginsparg_Wilson': whole_g5@d+d@whole_g5-d@whole_g5@d,
        'modified_chirality': hat@hat-np.eye(n),
        'Hermitian_chiral_projector': left-left.conj().T,
        'chiral_projector': left@left-left,
        'Weyl_map': d@left-right@d,
    }.items():
        admissions[name] = _admit(np.linalg.norm(residual, 'fro'), 2e-9, name)
    gauges = []
    for i in range(16):
        generators = np.zeros((5, 5), complex)
        generators[0, 1] = generators[1, 0] = .003*(i+1)
        generators[3, 4] = -.002j*(i+1)
        generators[4, 3] = .002j*(i+1)
        gauges.append(exterior_link(expm(1j*generators)))
    changed, gauge = gauge_transform(links, gauges)
    new_w, _ = wilson_operator(changed)
    new_d, new_hat, _ = overlap(new_w, whole_g5)
    admissions['Wilson_gauge_covariance'] = _admit(np.linalg.norm(new_w-gauge@w@gauge.conj().T, 'fro'), 2e-9, 'Wilson covariance')
    admissions['overlap_gauge_covariance'] = _admit(np.linalg.norm(new_d-gauge@d@gauge.conj().T, 'fro'), 2e-9, 'overlap covariance')
    admissions['chiral_projector_gauge_covariance'] = _admit(np.linalg.norm(new_hat-gauge@hat@gauge.conj().T, 'fro'), 2e-9, 'projector covariance')
    corners = []
    for bits in product((0, 1), repeat=4):
        symbol, _ = free_symbol(tuple(math.pi*x for x in bits))
        target = np.zeros((4, 4)) if not any(bits) else 2*np.eye(4)
        _admit(np.linalg.norm(symbol-target, 'fro'), 1e-12, 'all free corners')
        corners.append({'corner': list(bits), 'massless_Dirac_species': int(not any(bits)), 'dimensionless_D_at_corner': 0 if not any(bits) else 2})
    momenta = ((1/3, 0., 0., 0.), (1/3, 1/5, 1/7, 1/11))
    queries = [native_free_query(p)[0] for p in momenta]
    for spacing in (.125, .0625, .03125):
        for p in momenta:
            free_continuum_bound(p, spacing)
    # Independently evaluate the volume-uniform polynomial on the finite
    # four-direction free operator. The nonabelian background above is
    # checked separately, without representing its spectrum by free modes.
    scalar_links = np.ones((16, 4, 1, 1), complex)
    scalar_w, scalar_g5 = wilson_operator(scalar_links)
    polynomial, local = local_overlap_polynomial(scalar_w, scalar_g5, links=scalar_links)
    exact, _, _ = overlap(scalar_w, scalar_g5)
    _admit(np.linalg.norm(polynomial-exact, 2), float(F(local['uniform_operator_tail_exact']))+1e-11,
           'finite local overlap approximation')
    paths = (MATTER, MATTER_OWNER, 'model/common_source/overlap_regulator.py',
             'model/common_source/test_overlap_regulator.py', 'model/observer/source_scalar/compiler.py',
             *compiler.source_summary()['source_sha256'])
    return {'schema': 'oph.overlap_regulator.v1',
            'architecture': {'Euclidean_dimensions': 4, 'periodic_side': 2, 'sites': 16,
                             'internal_representation': 'Lambda2(C5)+Lambda4(C5)', 'internal_dimension': 15,
                             'Dirac_dimension': n, 'independent_families': 3,
                             'Weyl_components_per_site_with_families': 90,
                             'spacing_and_background_are_supplied': True},
            'representation': {'integer_hypercharges': charges, 'trace_q': sum(charges), 'trace_q_cubed': sum(q**3 for q in charges)},
            'finite_operator_admissions': admissions,
            'small_link_domain': {'fundamental_generator_norm_bound_exact': str(GENERATOR_BOUND),
                                 'represented_plaquette_norm_bound_exact': str(PLAQUETTE_BOUND),
                                 'admissibility_ceiling_exact': '1/40',
                                 'Wilson_H_squared_lower_bound_exact': '1/4',
                                 'Wilson_H_squared_upper_bound_exact': '64',
                                 'actual_background_plaquette_Frobenius_admission': 1/40,
                                 'domain_does_not_cover_all_gauge_fields': True},
            'finite_range_approximation': local,
            'Wilson_gap_admitted_above': .1,
            'left_projector_rank_per_family': int(round(float(np.trace(left).real))),
            'free_corners': corners,
            'free_continuum': {'parameter': 'Euclidean lattice spacing a at fixed model momentum',
                               'band': 'a*norm(p)<=1/4', 'analytic_operator_error_bound': 'a*norm(p)^2',
                               'interacting_continuum_limit': False},
            'native_free_kernel_queries': queries,
            'scope': {'finite_gauge_covariant_modified_chiral_operator_constructed': True,
                      'global_chiral_measure_phase_or_quantum_gauge_action_constructed': False,
                      'reflection_positivity_or_Lorentzian_continuation_established': False,
                      'full_many_fermion_field_source_executed': False,
                      'physical_masses_Yukawas_spacetime_or_hardware_accuracy_derived': False},
            'primary_sources': ['https://arxiv.org/abs/hep-lat/9707022', 'https://arxiv.org/abs/hep-lat/9802011',
                                'https://arxiv.org/html/hep-lat/9808010v2'],
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(set(paths))}}


def compare(retained, current, location='root'):
    """Exact contracts and custody; portable complete complex source outputs."""
    if type(retained) is not type(current):
        raise ValueError('overlap receipt type mismatch: '+location)
    if type(current) is dict:
        if retained.keys() != current.keys():
            raise ValueError('overlap receipt fields mismatch: '+location)
        for key in current:
            compare(retained[key], current[key], location+'.'+key)
    elif type(current) is list:
        if len(retained) != len(current):
            raise ValueError('overlap complete output length mismatch: '+location)
        for i, (left, right) in enumerate(zip(retained, current, strict=True)):
            compare(left, right, location+f'[{i}]')
    elif type(current) is float:
        if not math.isfinite(retained) or not math.isfinite(current):
            raise ValueError('nonfinite overlap evidence: '+location)
        if 'complete_success_and_failure_reference_columns' in location:
            if not math.isclose(retained, current, rel_tol=2e-9, abs_tol=2e-9):
                raise ValueError('complete complex native output changed: '+location)
        elif retained != current:
            raise ValueError('overlap numerical contract changed: '+location)
    elif retained != current:
        raise ValueError('overlap exact evidence changed: '+location)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    option = parser.add_mutually_exclusive_group(required=True)
    option.add_argument('--write', type=Path)
    option.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare(strict_json.load(args.check), result)
    else:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print('Finite overlap kernel, modified chirality, covariance and complete native queries passed')


if __name__ == '__main__':
    main()
