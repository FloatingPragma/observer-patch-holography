"""Continuum compatibility of the same native, gauge-linked staggered action.

Constant anti-Hermitian native currents define U_j(a)=exp(a A_j). The actual
nearest-neighbour source action has the complete eight-corner symbol
 m beta + c sum_j alpha_j sin(a K_j)/a, K_j=p_j I-i A_j.
No commutation of different K_j is assumed. Hermitian functional calculus
and unitary Duhamel give a uniform compact-band, finite-coordinate-time
bound. The full auxiliary register, two Dirac tastes and both chiralities
remain present. Independent site gauges act on the finite real-space law;
they generally leave the constant-background subclass.

This is a conditional background-field limit, not interacting continuum
QFT. Geometry, background/control paths, the coherent time parameter and
arbitrarily improving physical source accuracy are separate premises.
The finite matrix proof below is self-contained. For standard context see
Susskind, https://doi.org/10.1103/PhysRevD.16.3031; Catterall et al.,
https://arxiv.org/abs/2501.10862; Childs et al.,
https://doi.org/10.1103/PhysRevX.11.011020 (first-order product formula).
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from . import native_gauge_matter_transport as transport
from . import native_matter_attachment as matter
from . import native_staggered_field as field


ROOT = Path(__file__).resolve().parents[2]
CURRENT_INDICES = (0, 1, 2)
SPACINGS = (F(1, 4), F(1, 8), F(1, 16))
NUMERICAL_ADMISSION = 2e-9
SCHEMA = {'currents': (3, 32, 32), 'links': (3, 3, 32, 32),
          'site_gauge': (32, 32), 'site_fourier': (64, 64), 'family_mass': (3, 3)}


def _fraction(value, *, positive=False):
    if type(value) not in (int, F) or isinstance(value, bool):
        raise ValueError('an original exact integer or Fraction is required')
    value = F(value)
    if value < 0 or (positive and value == 0):
        raise ValueError('a nonnegative exact parameter is required')
    return value


@lru_cache(maxsize=1)
def architecture():
    r = matter.representation()
    currents = tuple(matter.exterior_derivative(matter.spin_lift(
        matter.canonical_currents()[j])) for j in CURRENT_INDICES)
    if any(not matter.zero(A+A.H) for A in currents):
        raise ArithmeticError('the actual lifted currents must be anti-Hermitian')
    e = r['encoder']
    if any(not matter.zero(A*e-e*(e.T*A*e)) for A in currents):
        raise ArithmeticError('selected matter/complement invariance failed')
    norms = []
    for A in currents:
        norm2 = s.expand(s.trace(A.H*A))
        upper = int(s.ceiling(s.sqrt(norm2)))
        if s.simplify(upper**2-norm2).is_nonnegative is not True:
            raise ArithmeticError('exact Frobenius norm majorant failed')
        norms.append(F(upper))
    commutators = tuple(matter.clean(currents[j]*currents[k]-currents[k]*currents[j])
                       for j in range(3) for k in range(j+1, 3))
    if all(matter.zero(x) for x in commutators):
        raise ArithmeticError('the declared non-Abelian control became commuting')
    return {'currents': currents, 'norm_bounds': tuple(norms), 'encoder': e,
            'commutators': commutators, 'clifford': field.clifford_matrices()}


def symbol(momentum, spacing, currents, *, mass=None, continuum=False):
    """Complete symbol; the input currents need not mutually commute."""
    a = float(_fraction(spacing, positive=True))
    A = np.asarray(currents, complex)
    if A.ndim != 3 or A.shape[0] != 3 or A.shape[1] != A.shape[2]:
        raise ValueError('three square current matrices required')
    if len(momentum) != 3 or not np.isfinite(momentum).all() or not np.isfinite(A).all():
        raise ValueError('finite momenta and currents required')
    if np.linalg.norm(A+A.conj().transpose(0, 2, 1)) > 1e-12:
        raise ValueError('anti-Hermitian currents required')
    d = A.shape[1]; p = field.native_parameters(); c = float(p['coordinate_speed'])
    m = float(p['mass'] if mass is None else mass)
    if not np.isfinite(m): raise ValueError('finite mass required')
    alpha = tuple(np.asarray(x, complex) for x in field.clifford_matrices())
    out = m*np.kron(alpha[3], np.eye(d))
    for j in range(3):
        K = float(momentum[j])*np.eye(d)-1j*A[j]
        if continuum:
            derivative = K
        else:
            U = expm(1j*a*K)
            derivative = (U-U.conj().T)/(2j*a)
        out += c*np.kron(alpha[j], derivative)
    return out


def source_symbol(momentum, spacing, links, *, mass=None):
    """Symbol directly from complete executed source links, without a log."""
    a = float(_fraction(spacing, positive=True)); links = np.asarray(links, complex)
    if links.ndim != 3 or links.shape[0] != 3 or links.shape[1] != links.shape[2]:
        raise ValueError('three complete square source links required')
    if len(momentum) != 3 or not np.isfinite(momentum).all() or not np.isfinite(links).all():
        raise ValueError('finite complete source links and momenta required')
    d = links.shape[1]; p = field.native_parameters()
    alpha = tuple(np.asarray(x, complex) for x in field.clifford_matrices())
    out = float(p['mass'] if mass is None else mass)*np.kron(alpha[3], np.eye(d))
    for j in range(3):
        forward = np.exp(1j*a*momentum[j])*links[j]
        out += float(p['coordinate_speed'])*np.kron(alpha[j],
            (forward-forward.conj().T)/(2j*a))
    return out


@lru_cache(maxsize=1)
def family_attachment():
    """Chosen common vectorlike attachment of the existing source response Y.

    This does not replace the chiral Higgs pairing or execute the joint
    field+finite-memory source action. Its three unequal gaps are the same
    target-free, empirically unsuccessful source hierarchy already retained.
    """
    from . import native_repair_mass as repair
    data = repair.architecture(); Y = data['Y']
    if not repair._zero(Y-Y.T) or s.trace(Y) != 1:
        raise ArithmeticError('the common source family response changed')
    if any(weight <= 0 for weight in data['stationary']):
        raise ArithmeticError('positive stationary source weights required')
    if any(not repair._zero(P*P-P) or s.trace(P) != 1 for P in data['projectors']):
        raise ArithmeticError('source family projector decomposition failed')
    return {'Y': Y, 'operator_norm_upper': F(1), 'source_record_probability': data['p']}


def family_symbol(momentum, spacing, currents, Y, *, continuum=False):
    """All eight corners x32 internal levels x3 families, no diagonalization."""
    Y = np.asarray(Y, complex)
    if Y.shape != (3,3) or not np.isfinite(Y).all() or np.linalg.norm(Y-Y.conj().T) > 1e-13:
        raise ValueError('a complete Hermitian three-family mass block required')
    d = np.asarray(currents).shape[-1]
    kinetic = symbol(momentum, spacing, currents, mass=0, continuum=continuum)
    beta = np.asarray(field.clifford_matrices()[3], complex)
    return np.kron(kinetic, np.eye(3))+np.kron(np.kron(beta, np.eye(d)), Y)


def refinement_contract(*, spacing, duration, momentum_cutoff, splitting_steps,
                        source_macro_count, source_macro_error, particles=1, mass_norm_upper=None):
    """Finite-particle, finite-volume/source contract, with no hidden precision limit.

    For seven one-body groups (onsite, six disjoint oriented bond matchings),
    Lambda=m+6J/a bounds the sum of norms. First-order splitting error is
    T^2 Lambda^2/(2N). A complete local Fock macro error epsilon contributes
    at most x/(1-x), x=G epsilon<1, by product telescoping; source accuracy is
    an independently supplied physical premise, not a floating-point test.
    """
    a = _fraction(spacing, positive=True); t = _fraction(duration)
    cutoff = _fraction(momentum_cutoff); epsilon = _fraction(source_macro_error)
    for value, name, low in ((splitting_steps, 'splitting steps', 1),
                             (source_macro_count, 'source macro count', 0),
                             (particles, 'fixed particle number', 0)):
        if type(value) is not int or value < low:
            raise ValueError('a finite nonnegative integer '+name+' is required')
    # pi>3 suffices for a deliberately strict, exact reduced-zone admission.
    if 2*a*cutoff >= 3:
        raise ValueError('the compact band must lie strictly inside the reduced zone')
    p = field.native_parameters(); c = p['coordinate_speed']
    bounds = architecture()['norm_bounds']
    delta = c*a*a*sum((cutoff+b)**3 for b in bounds)/6
    mass_upper = p['mass'] if mass_norm_upper is None else _fraction(mass_norm_upper)
    Lambda = mass_upper+6*p['exchange']/a
    splitting = t*t*Lambda*Lambda/(2*splitting_steps)
    x = source_macro_count*epsilon
    if x >= 1:
        raise ValueError('the accumulated source macro admission must be below one')
    source = x/(1-x)
    return {'spacing': str(a), 'source_coordinate_duration': str(t),
        'componentwise_momentum_cutoff': str(cutoff), 'fixed_particle_number': particles,
        'background_operator_norm_upper': list(map(str, bounds)),
        'declared_mass_block_operator_norm_upper': str(mass_upper),
        'one_body_discretization_error_upper': str(delta),
        'one_body_product_formula_error_upper': str(splitting),
        'complete_Fock_macro_count': source_macro_count,
        'complete_Fock_macro_error_assumption': str(epsilon),
        'complete_Fock_source_error_upper_conditional': str(source),
        'fixed_particle_propagator_error_upper_conditional': str(
            min(F(2), particles*(t*delta+splitting))+source),
        'finite_volume_domain': 'periodic cubic torus with an even site count >=4 on each axis; a=L/n, compatible fixed physical momenta inside the declared reduced-zone band',
        'background_domain': 'three constant anti-Hermitian native currents; noncommuting directions allowed; finite norm bounds uniform in a',
        'source_parameter': 'coherent source coordinate duration, not a derived physical clock',
        'precision_condition': 'G(a)*epsilon(a)->0 in addition to N(a)*a^2->infinity; fixed nonzero macro admission does not establish this limit',
        'sufficient_rate_example': 'at fixed physical volume and fixed word cost per site, N=O(a^-4), G=O(a^-7), epsilon=O(a^9) give O(a^2) error; actual routing/precision costs must be included in G',
        'source_precision_scaling_implemented': False,
        'arbitrarily_growing_source_register_implementation_supplied': False,
        'whole_growing_Fock_or_filled_sea_vacuum_limit_claimed': False}


def site_factors():
    coordinates = tuple(product(range(4), repeat=3)); index = {x: i for i, x in enumerate(coordinates)}
    forward = []
    for j in range(3):
        K = np.zeros((64, 64), complex)
        for i, x in enumerate(coordinates):
            y = list(x); y[j] = (y[j]+1) % 4
            K[i, index[tuple(y)]] = (-1)**sum(x[:j])
        forward.append(K)
    parity = np.diag([(-1)**sum(x) for x in coordinates])
    return tuple(forward), parity


def site_fourier():
    coordinates = np.asarray(tuple(product(range(4), repeat=3)))
    corners = np.asarray(tuple(product((0, 1), repeat=3)))
    return np.hstack([np.exp(1j*coordinates@(-np.pi/2*np.asarray(b)+np.pi*corners).T)/8
                     for b in product((0, 1), repeat=3)])


def apply_lattice(state, links, *, spacing=F(1), mass=None):
    """All real-space rows; last axes may carry any unknown reference."""
    a = float(_fraction(spacing, positive=True)); state = np.asarray(state, complex)
    links = np.asarray(links, complex); d = links.shape[-1]
    if links.shape not in ((3, d, d), (64, 3, d, d)) or state.shape[:2] != (64, d):
        raise ValueError('all64 sites and complete internal link/input columns required')
    if not np.isfinite(state).all() or not np.isfinite(links).all():
        raise ValueError('finite state and complete link arrays required')
    U = np.broadcast_to(links, (64, 3, d, d))
    forward, parity = site_factors(); p = field.native_parameters()
    out = float(p['mass'] if mass is None else mass)*np.einsum('xy,yi...->xi...', parity, state)
    J = float(p['exchange'])/a
    for axis, shift in enumerate(forward):
        for x, y in zip(*np.nonzero(shift), strict=True):
            z = -1j*J*shift[x, y]
            out[x] += z*np.einsum('ij,j...->i...', U[x, axis], state[y])
            out[y] += z.conjugate()*np.einsum('ij,j...->i...', U[x, axis].conj().T, state[x])
    return out


def transformed_links(links, gauges):
    links = np.asarray(links, complex); gauges = np.asarray(gauges, complex)
    d = links.shape[-1]
    if links.shape != (3, d, d) or gauges.shape != (64, d, d):
        raise ValueError('constant complete links and one complete gauge per site required')
    result = np.empty((64, 3, d, d), complex)
    for axis, shift in enumerate(site_factors()[0]):
        for x, y in zip(*np.nonzero(shift), strict=True):
            result[x, axis] = gauges[x]@links[axis]@gauges[y].conj().T
    return result


def execute():
    arrays = {'currents': np.asarray(architecture()['currents'], complex),
              'site_fourier': site_fourier(),
              'site_gauge': np.asarray(transport.current_path(((1, F(1, 3)),))['full32']),
              'family_mass': np.asarray(family_attachment()['Y'], complex)}
    arrays['links'] = np.asarray([[transport.current_path(((j, a),))['full32']
                                  for j in CURRENT_INDICES] for a in SPACINGS])
    validate_arrays(arrays)
    return arrays


def validate_arrays(arrays):
    if type(arrays) is not dict or set(arrays) != set(SCHEMA):
        raise ValueError('the complete continuum factor array menu is required')
    for key, shape in SCHEMA.items():
        value = arrays[key]
        if not isinstance(value, np.ndarray) or value.shape != shape or value.dtype != np.dtype(complex) or not np.isfinite(value).all():
            raise ValueError('invalid complete continuum factor: '+key)
    if np.linalg.norm(arrays['currents']-np.asarray(architecture()['currents'], complex)) > 1e-13:
        raise ValueError('retained currents differ from the exact native lift')
    if np.linalg.norm(arrays['site_fourier']-site_fourier()) > 1e-13:
        raise ValueError('retained full site/corner transform differs')
    if np.linalg.norm(arrays['family_mass']-np.asarray(family_attachment()['Y'], complex)) > 1e-13:
        raise ValueError('retained family mass differs from the fixed source response')
    return arrays


def report_from_arrays(arrays):
    validate_arrays(arrays)
    A = arrays['currents']; W = arrays['site_fourier']
    gauge_error = float(np.linalg.norm(arrays['site_gauge']-expm(A[1]/3), 2))
    if gauge_error > float(transport.CURRENT_ADMISSION):
        raise ValueError('retained complete site gauge fails its unchanged source admission')
    alpha = tuple(np.asarray(x, complex) for x in field.clifford_matrices())
    c = float(field.native_parameters()['coordinate_speed']); p = (0.2, -0.3, 0.1)
    rows = []
    for a, links in zip(SPACINGS, arrays['links'], strict=True):
        ideal = np.asarray([expm(float(a)*x) for x in A])
        link_error = max(float(np.linalg.norm(x-y, 2)) for x, y in zip(links, ideal, strict=True))
        if link_error > float(transport.CURRENT_ADMISSION):
            raise ValueError('retained complete source links fail their unchanged admission')
        h = source_symbol(p, a, links); hd = symbol(p, a, A, continuum=True)
        radius = float(float(a)**2*c*sum(np.linalg.norm(p[j]*np.eye(32)-1j*A[j], 2)**3 for j in range(3))/6)
        source_bound = 3*c*float(transport.CURRENT_ADMISSION)/float(a)
        error = float(np.linalg.norm(h-hd, 2))
        if error > radius+source_bound+NUMERICAL_ADMISSION:
            raise ArithmeticError('noncommuting background continuum bound failed')
        rows.append({'spacing': str(a), 'complete_link_error_diagnostic': link_error,
            'symbol_discretization_error_diagnostic': error,
            'sine_bound_diagnostic': radius, 'fixed_source_link_error_upper_conditional': source_bound})
    # Exact Fourier identities are tested on each site coefficient, before
    # tensoring arbitrary complete internal links and all their columns.
    forward, parity = site_factors(); fourier_error = np.linalg.norm(W.conj().T@parity@W-np.kron(np.eye(8), alpha[3]), 2)
    for j, shift in enumerate(forward):
        blocks = [np.exp(-1j*np.pi*b[j]/2)*alpha[j] for b in product((0, 1), repeat=3)]
        expected = np.zeros((64, 64), complex)
        for k, block in enumerate(blocks): expected[8*k:8*k+8, 8*k:8*k+8] = block
        fourier_error = max(fourier_error, np.linalg.norm(W.conj().T@shift@W-expected, 2))
    if fourier_error > 1e-12: raise ArithmeticError('the complete lattice/symbol map failed')
    # The continuum squared symbol retains the field-strength cross terms.
    K = [p[j]*np.eye(32)-1j*A[j] for j in range(3)]
    hd = symbol(p, SPACINGS[0], A, continuum=True)
    m = float(field.native_parameters()['mass'])
    scalar = m*m*np.eye(256)+c*c*np.kron(np.eye(8), sum(x@x for x in K))
    curvature = c*c*sum((np.kron(alpha[j]@alpha[k], K[j]@K[k]-K[k]@K[j])
                          for j in range(3) for k in range(j+1, 3)), np.zeros((256,256), complex))
    curvature_error = float(np.linalg.norm(hd@hd-scalar-curvature, 2))
    if curvature_error > 1e-10 or np.linalg.norm(curvature, 2) < 1e-4:
        raise ArithmeticError('the non-Abelian field-strength contribution was lost')
    # A fixed corner chirality is not an invariant real-space gauge subspace.
    chi = 1j*alpha[0]@alpha[1]@alpha[2]
    P = W@np.kron(np.eye(8), (np.eye(8)+chi)/2)@W.conj().T
    E = np.zeros((64, 64)); E[0, 0] = 1
    leak = float(np.linalg.norm((np.eye(64)-P)@E@P, 'fro')*
                 np.linalg.norm(arrays['site_gauge']-np.eye(32), 'fro'))
    if leak < 1e-4: raise ArithmeticError('the site-gauge chirality countercontrol vanished')
    gauges = np.repeat(np.eye(32, dtype=complex)[None], 64, axis=0)
    gauges[0] = arrays['site_gauge']; gauges[21] = arrays['site_gauge'].conj().T
    changed = transformed_links(arrays['links'][0], gauges)
    covariance_error = max(float(np.linalg.norm(G@G.conj().T-np.eye(32), 2)) for G in gauges)
    for j, shift in enumerate(forward):
        for x, y in zip(*np.nonzero(shift), strict=True):
            covariance_error = max(covariance_error, float(np.linalg.norm(
                changed[x,j]@gauges[y]-gauges[x]@arrays['links'][0,j], 2)))
    if covariance_error > NUMERICAL_ADMISSION:
        raise ArithmeticError('complete local-gauge covariance exceeds its numerical admission')
    return {'current_indices': list(CURRENT_INDICES), 'corner_dimension': 8,
        'scalar_mass_from_native_sector_zero': str(field.native_parameters()['mass']),
        'coordinate_speed_from_native_exchange': str(field.native_parameters()['coordinate_speed']),
        'full_auxiliary_internal_dimension': 32, 'selected_matter_dimension': 15,
        'retained_auxiliary_complement_dimension': 17, 'declared_family_factor': 3,
        'complete_finite_lattice_dimension_with_families': 64*32*3,
        'finite_lattice_operator_retained_as_complete_factors': True,
        'global_periodic_lattice_Fock_source_word_executed': False,
        'executed_source_scope': 'complete32-column links and site gauge from the same native backend; local physical-CAR transport factors are retained by the canonical transport receipt',
        'selected_Dirac_tastes_per_internal_component': 2,
        'selected_Dirac_fields_with_family_factor': 2*15*3,
        'source_link_refinements': rows,
        'source_link_paths': [[[[j,str(a)]] for j in CURRENT_INDICES] for a in SPACINGS],
        'complete_Fourier_factor_error_diagnostic': float(fourier_error),
        'complete_source_site_gauge_error_diagnostic': gauge_error,
        'maximum_unscaled_real_space_gauge_block_error_diagnostic': covariance_error,
        'local_covariance_domain': 'all192 oriented links, all64 onsite blocks and all32 internal columns; transformed links need not be constant',
        'field_strength_square_identity_error_diagnostic': curvature_error,
        'omitted_field_strength_operator_norm_diagnostic': float(np.linalg.norm(curvature, 2)),
        'fixed_chirality_site_gauge_leakage_Frobenius_diagnostic': leak,
        'chirality_control_scope': 'a fixed free corner projector cannot be kept while allowing arbitrary local gauge changes; a transformed projector describes a different subspace, not a selected chiral gauge theory',
        'continuum_bound': '||h_a(p)-h_D(p)|| <= c*a^2*sum_j ||p_j I-i A_j||^3/6; propagation <= |T| times this bound',
        'bound_derivation': [
            'anti-Hermitian A_j makes K_j=p_j I-i A_j Hermitian; the scalar integral Taylor bound |sin(x)-x|<=|x|^3/6 passes through its spectral calculus',
            'each alpha_j has norm one; summing three operator norms avoids any assumption that different K_j commute',
            'unitary Duhamel bounds the finite coherent-time propagator error by |T| times the Hamiltonian error; exterior products give the fixed-particle factor',
            'the six even/odd nearest-neighbour matchings and onsite group have summed one-body norms at most m+6J/a; the first-order commutator bound is T^2*(m+6J/a)^2/(2N)',
            'complete numerical source macros need their own telescope x/(1-x), x=G*epsilon<1, in addition to the ideal-unitary discretization and splitting errors'],
        'refinement_contract_example': refinement_contract(spacing=F(1,16), duration=F(1),
            momentum_cutoff=F(1), splitting_steps=16**4, source_macro_count=7*64**3*16**4,
            source_macro_error=F(1,16**2*2*7*64**3*16**4)),
        'noncommuting_backgrounds_or_source_leakage_dropped': False,
        'physical_chiral_SM_or_dynamical_gauge_continuum_claimed': False,
        'different_measured_species_masses_supplied_or_derived': False,
        'common_source_family_attachment': {
            'mass_block': 'beta tensor I32 tensor Y; gauge links U tensor I3',
            'Y': [[str(x) for x in row] for row in family_attachment()['Y'].tolist()],
            'operator_norm_upper': '1', 'native_event_probability': str(family_attachment()['source_record_probability']),
            'scalar_native_gap_control_retained': True,
            'discretization_and_finite_time_bounds_unchanged': True,
            'product_formula_mass_norm_upper': '1',
            'same_mass_block_for_all_internal_components': True,
            'chosen_attachment_not_the_chiral_Higgs_pairing': True,
            'joint_field_and_finite_memory_source_execution_claimed': False,
            'measured_particle_hierarchy_recovered': False},
        'coherent_precision_scaling_and_continuum_geometry_are_premises': True,
        'references': ['https://doi.org/10.1103/PhysRevD.16.3031',
                       'https://arxiv.org/abs/2501.10862',
                       'https://doi.org/10.1103/PhysRevX.11.011020']}


def source_bindings():
    from . import native_code_transport
    paths = set(matter.source_bindings()) | set(field.source_bindings()) | set(native_code_transport.source_bindings()) | {
        'model/common_source/native_repair_mass.py',
        'model/common_source/test_native_repair_mass.py',
        'model/common_source/native_gauge_matter_transport.py',
        'model/common_source/test_native_gauge_matter_transport.py',
        'model/common_source/native_gauge_continuum.py',
        'model/common_source/test_native_gauge_continuum.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def build():
    arrays = execute()
    return report_from_arrays(arrays), arrays


def _compare_report(expected, actual, path='continuum'):
    if type(expected) is not type(actual):
        raise ValueError('continuum report type differs at '+path)
    if isinstance(expected, dict):
        if set(expected) != set(actual): raise ValueError('continuum report fields differ at '+path)
        for key in expected: _compare_report(expected[key], actual[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(expected) != len(actual): raise ValueError('continuum report length differs at '+path)
        for i, (a,b) in enumerate(zip(expected,actual,strict=True)):
            _compare_report(a,b,path+'.'+str(i))
    elif isinstance(expected, float):
        # Only diagnostic leaves vary. Contract/admission values are exact.
        if '_diagnostic' not in path or not np.isfinite(actual) or abs(expected-actual) > 128*np.finfo(float).eps*max(1., abs(expected)):
            if expected != actual: raise ValueError('continuum report number differs at '+path)
    elif expected != actual:
        raise ValueError('continuum report value differs at '+path)


def verify_arrays(report, arrays):
    """Bind every reported number to the complete retained factor arrays."""
    _compare_report(report_from_arrays(arrays), report)
    return True
