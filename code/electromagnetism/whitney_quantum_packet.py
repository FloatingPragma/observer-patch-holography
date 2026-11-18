"""Full interacting Coulomb phase-space preparation and neutral Gaussian seeds.

Two classical samples supply preparation parameters, not a quantum history.
The kinetic reduction is the actual 68-coordinate Schur reduction. Every
56-dimensional seed direction is retained before exact circle projection.
Bessel evaluations and metric matrices are numerical diagnostics; the paper
proves exact initial cotangent/projection bounds and the operator domain.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import mpmath
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'runtime/whitney_quantum_packet_receipt.json'
SCHEMA = 'oph.whitney_quantum_packet.v1'
SCOPE = 'FULL_56D_NEUTRAL_PACKET_PREPARATION__NO_QUANTUM_PROPAGATION_CERTIFICATE'
PARENT_PATH = 'code/electromagnetism/runtime/whitney_charged_dynamics_receipt.json'
PIN_PATHS = (
    'Lean/Screen/SeamCurrentEdge30Moment.lean',
    'Lean/Screen/SeamCurrentCarrierQuotient.lean',
    'Lean/ObserverPatchHolography/CoreAxioms.lean',
    'code/electromagnetism/verify_cone_whitney_bridge.py',
    'code/electromagnetism/whitney_interacting_quantum.py',
    'code/electromagnetism/verify_whitney_quantum_state.py',
    'code/electromagnetism/verify_whitney_charged_dynamics.py',
    'code/electromagnetism/whitney_quantum_packet.py',
    'code/electromagnetism/verify_whitney_quantum_packet.py',
    'code/electromagnetism/test_whitney_quantum_packet.py',
    'code/electromagnetism/test_neutral_packet_observables.py',
    'code/electromagnetism/test_neutral_packet_projection.py',
    PARENT_PATH,
)


def sibling(name, canonical_name=None):
    """Load a pinned sibling independently of the caller's working directory."""
    path = HERE/(name+'.py')
    module_name = canonical_name or 'packet_preparation_'+name
    existing = sys.modules.get(module_name)
    if existing is not None:
        if Path(existing.__file__).resolve() != path.resolve():
            raise ImportError('unexpected sibling module origin')
        return existing
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError('missing sibling '+name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


sibling('verify_cone_whitney_bridge', 'verify_cone_whitney_bridge')
quantum = sibling('whitney_interacting_quantum')
parent = sibling('verify_whitney_charged_dynamics')


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode('ascii')


def real_vector(value, size, name):
    # These geometric/sampling callers are binary64 calculations. Validate
    # each original entry before coercion and refuse a changed input value.
    return np.array([_binary64(v, name) for v in _exact_vector(value, size)])


def positive(value, name):
    result = _binary64(_rational(value), name)
    if result <= 0:
        raise ValueError('positive finite '+name+' required')
    return result


def _binary64(value, name):
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(name+' must be exactly representable in binary64') from error
    if not np.isfinite(result) or Q(result) != value:
        raise ValueError(name+' must be exactly representable in binary64')
    return result


def phase_space(q68, velocity68, mesh=None, charge=.25):
    """True tangent rechart, then the full Schur cotangent at c=0.

    No Gauss condition is silently imposed. The difference between the
    transformed temporal potential and its Schur minimizer is returned.
    """
    mesh = quantum.geometry(4) if mesh is None else mesh
    q68 = real_vector(q68, 68, 'configuration')
    velocity68 = real_vector(velocity68, 68, 'velocity')
    charge = positive(charge, 'charge')
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            return _phase_space_values(q68, velocity68, mesh, charge)
    except FloatingPointError as error:
        raise ValueError('kinetic precision outside finite phase-space reporting range') from error


def _phase_space_values(q68, velocity68, mesh, charge):
    """Evaluate validated binary64 inputs within the public range guard."""
    a, psi = q68[:42], q68[42:55]+1j*q68[55:]
    av, pv = velocity68[:42], velocity68[42:55]+1j*velocity68[55:]
    b = mesh.mean_zero
    laplacian = b.T@mesh.d.T@mesh.mass@mesh.d@b
    xi = -b@np.linalg.solve(laplacian, b.T@mesh.d.T@mesh.mass@a)
    xidot = -b@np.linalg.solve(laplacian, b.T@mesh.d.T@mesh.mass@av)
    rotation = np.exp(1j*charge*xi)
    ac, pc = a+mesh.d@xi, rotation*psi
    gauge_edge_velocity = mesh.d@xidot
    gauge_scalar_velocity = 1j*charge*xidot*psi
    acv, pcv = av+gauge_edge_velocity, rotation*(pv+gauge_scalar_velocity)
    full_q, full_v = np.r_[ac, pc.real, pc.imag], np.r_[acv, pcv.real, pcv.imag]
    # A resolved metric cannot repair an already cancelled tangent. This
    # normwise rechart policy retains directly supplied tiny velocities but
    # refuses unresolved cancellation of large gauge/scalar summands. Use
    # max norms and compare before multiplying tiny scales by epsilon.
    # Use the same real-coordinate max norm as full_v. A complex modulus
    # may overflow although both supplied real components are representable.
    upstream = max(float(np.max(abs(component)))
                   for term in (av, pv, gauge_edge_velocity, gauge_scalar_velocity)
                   for component in (term.real, term.imag))
    reduced = float(np.max(abs(full_v)))
    allowance = quantum.KINETIC_RESOLUTION/(128*np.finfo(float).eps*np.linalg.cond(laplacian))
    if (not np.isfinite(upstream) or not np.isfinite(reduced)
            or not np.isfinite(allowance) or allowance <= 0
            or (upstream > 0 and (reduced == 0 or reduced/upstream < 1/allowance))):
        raise ValueError('kinetic precision insufficient for the recharted tangent')
    q, velocity = mesh.slice.T@full_q, mesh.slice.T@full_v
    reduction = quantum.reduced_kinetic(ac, pc, charge, .5, .25, mesh)
    # Componentwise cotangents need a resolved dense geometry, even when
    # square-root log densities can still be evaluated beyond that regime.
    gamma = reduction.dense_metric()
    momentum = reduction.factor.T@(reduction.factor@velocity)
    generator = np.r_[np.zeros(30), -pc.imag, pc.real]
    eta = mesh.mean_zero@reduction.eta_map@velocity
    result = {'q': q, 'velocity': velocity, 'momentum': momentum,
            'full_q': full_q, 'full_velocity': full_v, 'gauge_parameter': xi,
            'gauge_parameter_velocity': xidot, 'transformed_scalar_potential': -xidot,
            'schur_scalar_potential': eta, 'constant_moment_map': float(momentum@generator),
            'minimizer_defect': float(np.max(abs(eta+xidot))),
            'coulomb_defect': float(np.max(abs(mesh.d.T@mesh.mass@ac))),
            'cotangent_identity_defect': float(np.max(abs(momentum-gamma@velocity))),
            'log_rho': reduction.logdet()/2, 'potential': reduction.potential,
            'gamma_min_eigenvalue': float(reduction.singular_values[-1]**2)}
    if any(not np.isfinite(value).all() for value in result.values()):
        raise ValueError('kinetic precision outside finite phase-space reporting range')
    return result


def overlap_parameters(center, momentum, sigma, hbar=1):
    values = _scalar_observables(center, momentum, sigma, hbar)
    return {key: _reported(values[key], key) for key in ('A', 'B', 'norm_squared')}


def _rational(value):
    """Preserve each original real scalar before any array coercion."""
    if (np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, float, np.integer, np.floating, Q))):
        raise ValueError('finite real packet scalar required')
    try:
        # Fraction can retain NumPy integer components and thereby their
        # fixed-width overflow. Canonicalize before any rational operation.
        if isinstance(value, (int, np.integer)):
            return Q(int(value))
        if isinstance(value, Q):
            return Q(int(value.numerator), int(value.denominator))
        return Q(*value.as_integer_ratio())
    except (ValueError, OverflowError) as error:
        raise ValueError('finite real packet scalar required') from error


def _exact_vector(value, size=56):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked packet vector')
    raw = np.asarray(value, dtype=object)
    if raw.shape != (size,):
        raise ValueError(str(size)+' real packet coordinates required')
    return [_rational(v) for v in raw]


def _reported(value, name):
    """Require relative 1e-12 representability; a positive value is not zero."""
    result = float(value)
    if not np.isfinite(result) or (value != 0 and
            abs((result-value)/value) > 1e-12):
        raise ValueError(name+' outside reliable binary64 reporting range')
    return result


def _scalar_observables(center, momentum, sigma, hbar):
    """One original-input calculation for both neutral-state observables.

    Exact rational invariants decide A >= |B|, including equality. Scaled
    Bessel evaluation uses a private context with enough extra digits for
    z*(1-I1/I0). This is numerical evaluation, not an interval certificate.
    """
    x, p = _exact_vector(center)[30:], _exact_vector(momentum)[30:]
    s, h = _rational(sigma), _rational(hbar)
    if s <= 0 or h <= 0:
        raise ValueError('positive width and hbar required')
    xx = sum(v*v for v in x)/(4*s*s)
    pp = s*s*sum(v*v for v in p)/(h*h)
    b = sum(p[i+13]*x[i]-p[i]*x[i+13] for i in range(13))/h
    a, difference = xx+pp, xx-pp
    gram = 4*xx*pp-b*b
    if gram < 0:
        raise ValueError('invalid exact packet Gram determinant')
    mp = mpmath.mp.clone()
    # bit_length gives a conservative decimal digit budget for z <= A.
    mp.dps = 80+max(0, (a.numerator.bit_length()-a.denominator.bit_length()+3)//3)
    def real(q):
        return mp.mpf(q.numerator)/q.denominator
    aa, bb, ss, dd = map(real, (a, b, s*s, difference))
    z = mp.sqrt(real(difference*difference+gram))
    if z == 0:
        scaled_i0, deficit = mp.mpf(1), mp.mpf(0)
    else:
        i0 = mp.besseli(0, z)
        scaled_i0 = mp.exp(-z)*i0
        deficit = z*(1-mp.besseli(1, z)/i0)
    # Rationalize each dangerous subtraction rather than flooring it.
    gap = real(b*b)/(aa+z) if a else mp.mpf(0)
    positive_part = real(gram)/(z-dd) if difference < 0 else z+dd
    radius = 2*ss*(13+positive_part-deficit)
    if radius <= 0 or not mp.isfinite(radius):
        raise ValueError('unresolved positive packet radius')
    return {'A': aa, 'B': bb, 'norm_squared': mp.exp(-gap)*scaled_i0,
            'radius': radius}


def seed_log_half_density(point, center, momentum, sigma, hbar=1):
    """Seed logarithm with original-input displacement and unwrapped phase.

    Inputs retain the binary64 contract. Rational arithmetic prevents a
    cancelled dot product or a tiny squared displacement from being erased.
    The returned phase must also have absolute error at most 1e-12 radians;
    relative representability alone does not suffice for a phase.
    """
    point, center, momentum, sigma, hbar = _pointwise_inputs(
        point, center, momentum, sigma, hbar)
    delta = [y-x for y, x in zip(point, center)]
    decay = sum(v*v for v in delta)/(4*sigma*sigma)
    phase = sum(p*v for p, v in zip(momentum, delta))/hbar
    try:
        phase_float = float(phase)
    except OverflowError as error:
        raise ValueError('seed phase outside binary64 reporting range') from error
    if not np.isfinite(phase_float):
        raise ValueError('seed phase outside binary64 reporting range')
    phase_error = abs(Q(phase_float)-phase)
    if ((phase != 0 and phase_error > Q(1, 10**12)*abs(phase))
            or phase_error > Q(1, 10**12)):
        raise ValueError('seed phase outside absolute binary64 reporting precision')

    def evaluate(mp):
        return -14*(mp.log(2*mp.pi)+2*mp.log(_mp_rational(mp, sigma)))-_mp_rational(mp, decay)

    real = _pointwise_evaluation(evaluate, _pointwise_digits(decay),
                                 'seed logarithm', complex_output=False)
    return complex(real, phase_float)


def _pointwise_inputs(point, center, momentum, sigma, hbar):
    """Validate original values, then keep them exact for pointwise algebra."""
    vectors = []
    for value, name in ((point, 'packet evaluation point'),
                        (center, 'packet center'), (momentum, 'packet momentum')):
        exact = _exact_vector(value)
        for entry in exact:
            _binary64(entry, name)
        vectors.append(exact)
    sigma, hbar = _rational(sigma), _rational(hbar)
    positive(sigma, 'width')
    positive(hbar, 'hbar')
    return (*vectors, sigma, hbar)


def _mp_rational(mp, value):
    return mp.mpf(value.numerator)/value.denominator


def _pointwise_digits(*values):
    # Cover absolute errors in large phases and cancellation of large log
    # amplitudes. The binary64 input contract bounds this finite budget.
    sizes = [(v.numerator.bit_length()-v.denominator.bit_length()+3)//3 for v in values]
    return 80+max(0, *sizes)


def _pointwise_evaluation(evaluate, digits, name, *, complex_output=True):
    """Numerical precision agreement plus final relative reporting guard.

    For complex values the 1e-12 reporting tolerance uses the full modulus,
    rather than requiring relative precision in each Cartesian component.
    This is not an interval enclosure or a proof of an exact zero. A zero
    from the special-function calculation is unresolved, never a successful
    pointwise amplitude. Extra precision also checks cancellation near a
    Bessel zero. All arithmetic is isolated from the caller's mp context.
    """
    mp = mpmath.mp.clone()
    previous = None
    for extra in (0, 40, 100, 220, 460):
        mp.dps = digits+extra
        try:
            value = evaluate(mp)
        except (mp.NoConvergence, OverflowError, ZeroDivisionError) as error:
            raise ValueError(name+' unresolved at available numerical precision') from error
        if not mp.isfinite(value):
            raise ValueError(name+' outside finite numerical evaluation range')
        if (value != 0 and previous is not None and previous != 0
                and abs(value-previous) <= mp.mpf('1e-30')*abs(value)):
            try:
                if complex_output:
                    result = complex(value)
                    reported = mp.mpc(result.real, result.imag)
                else:
                    result = float(value)
                    reported = mp.mpf(result)
            except OverflowError as error:
                raise ValueError(name+' outside binary64 reporting range') from error
            if (not np.isfinite(result)
                    or abs(reported-value) > mp.mpf('1e-12')*abs(value)):
                raise ValueError(name+' outside reliable binary64 reporting range')
            return result
        previous = value
    raise ValueError(name+' unresolved at available numerical precision')


def rotate(vector, angle):
    vector = real_vector(vector, 56, 'circle vector')
    angle = _binary64(_rational(angle), 'angle')
    z = (vector[30:43]+1j*vector[43:])*np.exp(1j*angle)
    return np.r_[vector[:30], z.real, z.imag]


def projected_half_density(point, center, momentum, sigma, hbar=1, nodes=256):
    """Normalized neutral projection from the closed complex Bessel integral.

    ``nodes`` is a validated legacy argument and no longer affects the result.
    The full expression, including normalization, is evaluated before any
    binary64 conversion. Precision agreement is numerical, not a certified
    pointwise error bound. Unresolved zeros and reporting ranges are refused.
    """
    if type(nodes) is not int or nodes < 16:
        raise ValueError('at least sixteen integer circle nodes required')
    point, center, momentum, sigma, hbar = _pointwise_inputs(
        point, center, momentum, sigma, hbar)
    y, x, p = point[30:], center[30:], momentum[30:]
    s2 = sigma*sigma
    dot = lambda a, b: sum(v*w for v, w in zip(a, b))
    jx, jp = [-v for v in x[13:]]+x[:13], [-v for v in p[13:]]+p[:13]
    radiative_delta = [v-w for v, w in zip(point[:30], center[:30])]
    decay = -(dot(radiative_delta, radiative_delta)+dot(y, y)+dot(x, x))/(4*s2)
    phase = (dot(momentum[:30], radiative_delta)-dot(p, x))/hbar
    ur, ui = dot(y, x)/(2*s2), dot(p, y)/hbar
    vr, vi = dot(y, jx)/(2*s2), dot(jp, y)/hbar
    # W = U^2+V^2, with U and V the cosine and sine coefficients.
    wr, wi = ur*ur-ui*ui+vr*vr-vi*vi, 2*(ur*ui+vr*vi)
    a = dot(x, x)/(4*s2)+s2*dot(p, p)/(hbar*hbar)
    b = dot(p, jx)/hbar
    z2 = a*a-b*b
    if z2 < 0:
        raise ValueError('invalid exact packet Gram determinant')
    # The norm squared is exp(-A)*I0(sqrt(A^2-B^2)). Combine
    # decay+A/2 exactly: either separate exponential can be unreportable
    # while the normalized state is an ordinary centered Gaussian.
    normalized_decay = decay+a/2

    def evaluate(mp):
        real = lambda v: _mp_rational(mp, v)
        if wi == 0 and wr < 0:
            integral = mp.besselj(0, mp.sqrt(real(-wr)))
        elif wr == 0 and wi == 0:
            integral = mp.mpf(1)
        else:
            integral = mp.besseli(0, mp.sqrt(mp.mpc(real(wr), real(wi))))
        if integral == 0:
            return mp.mpc(0)
        norm_log = mp.log(mp.besseli(0, mp.sqrt(real(z2)))) if z2 else mp.mpf(0)
        log_size = (-14*(mp.log(2*mp.pi)+2*mp.log(real(sigma)))
                    +real(normalized_decay)+mp.log(abs(integral))-norm_log/2)
        return mp.exp(log_size)*mp.exp(mp.j*real(phase))*(integral/abs(integral))

    return _pointwise_evaluation(evaluate,
        _pointwise_digits(normalized_decay, a, phase, ur, ui, vr, vi), 'pointwise amplitude')


def scalar_radius_moment(center, momentum, sigma, hbar=1):
    """Exact nodal-radius formula, evaluated numerically; not a spatial L2 norm."""
    return _reported(_scalar_observables(center, momentum, sigma, hbar)['radius'], 'radius')


def exact_initial():
    momentum_norm = [Q(91, 6), Q(13, 2)]
    widths = []
    for sigma in (Q(1, 4), Q(1, 2), Q(1)):
        aa = [Q(13)/(4*sigma**2)+sigma**2*momentum_norm[0], sigma**2*momentum_norm[1]]
        # Vol<18, and 39 Vol^2/400 bounds the squared scalar momentum.
        upper = Q(13)/(4*sigma**2)+sigma**2*Q(39, 400)*18**2
        widths.append({'sigma': str(sigma), 'A_in_Qsqrt5': [str(x) for x in aa],
            'A_rational_upper': str(upper), 'B': '0', 'projection_norm_squared_lower': '1/64',
            'seed_position_variance': str(sigma**2), 'seed_momentum_variance_hbar1': str(1/(4*sigma**2))})
    return {'volume_in_Qsqrt5': ['10', '10/3'],
        'radial_velocity_in_Qsqrt5': ['7/80', '3/80'],
        'center_gauge_factor': '12/13', 'boundary_gauge_factor': '-1/13',
        'scalar_center_real': '1', 'scalar_boundary_real': '1',
        'momentum_center_imaginary_in_Qsqrt5': ['3', '1'],
        'momentum_boundary_imaginary_in_Qsqrt5': ['-1/4', '-1/12'],
        'scalar_momentum_norm_squared_in_Qsqrt5': [str(x) for x in momentum_norm],
        'radiative_center_and_momentum': '0', 'constant_moment_map': '0',
        'proof_labels': ['prop:whitney-packet-coulomb', 'prop:whitney-packet-neutral-projection'],
        'widths': widths}


def contract():
    return {'configuration_dimension': 56, 'radiative_dimension': 30, 'scalar_real_dimension': 26,
        'seed_covariance': 'sigma^2 I_56; full rank; no classical fixed-space restriction',
        'projected_covariance': 'not the seed covariance; circle averaging changes scalar moments',
        'seed_formula': '(2*pi*sigma^2)^(-14)*exp(-|x-q|^2/(4*sigma^2)+i*p.(x-q)/hbar)',
        'state_formula': 'F=rho^(-1/2)*P0(seed)/||P0(seed)||; rho=sqrt(det(gamma))',
        'quantization': 'same full interacting Laplace-Beltrami Hamiltonian and metric volume',
        'circle_projection': 'exact group integral; floating Bessel values are diagnostics',
        'physical_comparison': False, 'empirical_prediction': False, 'observer_preparation': False,
        'quantum_propagation_certified': False, 'covariance_evolution_computed': False,
        'full_56D_preparation': True, 'neutral_strong_domain_state': True,
        'classical_samples_are_quantum_means': False,
        'scalar_one_point_mean': '0 by residual U(1) invariance',
        'sample_role': 'two separate packet preparations from adjacent classical samples; not quantum propagation',
        'numeric_scope': 'float64 classical data, phase-space maps and exact-degree simplex integration; no interval trajectory enclosure',
        'requested_propagation_target': {'T': '1/40', 'hbar': '1', 'norm_error': '1/10',
            'status': 'NOT_CERTIFIED', 'residual_upper_bound': None, 'all_tail_bound': None}}


def build():
    raw = parent.load(ROOT/PARENT_PATH)
    parent.verify(raw)
    mesh = quantum.geometry(4)
    samples = []
    for index in (0, 1):
        source = raw['samples'][index]
        data = phase_space(source['q'], source['velocity'], mesh)
        row = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in data.items()}
        row.update({'parent_sample_index': index, 'model_time': source['t'], 'widths': []})
        for sigma in (.25, .5, 1.):
            values = overlap_parameters(data['q'], data['momentum'], sigma)
            values.update({'sigma': str(Q(sigma)), 'scalar_radius_numeric': scalar_radius_moment(data['q'], data['momentum'], sigma)})
            row['widths'].append(values)
        samples.append(row)
    return {'schema': SCHEMA, 'scope': SCOPE,
        'source_pins': {x: hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in PIN_PATHS},
        'parameters': {'e': '1/4', 'm_squared': '1/2', 'g': '1/4', 'hbar': '1'},
        'contract': contract(), 'initial_exact': exact_initial(),
        'orthonormal_coulomb_frame': mesh.slice[:42, :30].tolist(), 'samples': samples}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({'receipt': str(args.output), 'samples': 2, 'quantum_propagation_certified': False}))
