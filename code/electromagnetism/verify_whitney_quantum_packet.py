"""Independent Coulomb tangent, Schur cotangent and circle-projection replay.

Neither the packet producer nor interacting coefficient evaluator is imported.
The metric at the charged Coulomb representative a=0 is assembled from
simplex monomials; the gauge solve uses a grounded augmented Laplacian.
The parent trajectory is independently verified fresh after packet checks.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from math import factorial, prod
from pathlib import Path

import mpmath
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'runtime/whitney_quantum_packet_receipt.json'
PARENT_PATH = 'code/electromagnetism/runtime/whitney_charged_dynamics_receipt.json'
PIN_PATHS = {
    'Lean/Screen/SeamCurrentEdge30Moment.lean', 'Lean/Screen/SeamCurrentCarrierQuotient.lean',
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
}


def sibling(name):
    spec = importlib.util.spec_from_file_location('packet_independent_'+name, HERE/(name+'.py'))
    if spec is None or spec.loader is None:
        raise ImportError('missing independent source verifier')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


source = sibling('verify_whitney_quantum_state')
parent = sibling('verify_whitney_charged_dynamics')
require = source.require


def load(path=OUTPUT):
    return source.load(path)


def keys(value, expected, name):
    require(type(value) is dict and set(value) == set(expected), name+' key census')


def equal(actual, expected, name):
    try:
        require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
                json.dumps(expected, sort_keys=True, allow_nan=False), name)
    except (TypeError, ValueError) as error:
        raise ValueError(name) from error


def numbers(value, shape, name):
    def valid(item):
        return all(valid(v) for v in item) if type(item) is list else type(item) in (int, float)
    require(valid(value), name+' numeric types')
    result = np.asarray(value, dtype=float)
    require(result.shape == shape and np.isfinite(result).all(), name+' shape/finiteness')
    return result


def close(value, expected, name, tolerance=3e-10):
    actual = numbers(value, np.shape(expected), name)
    require(np.max(abs(actual-expected), initial=0) <= tolerance*(1+np.max(abs(expected), initial=0)), name)
    return actual


def simplex(powers):
    return Q(6*prod(factorial(x) for x in powers), factorial(3+sum(powers)))


def metric_and_potential(psi, xyz, edges, tetrahedra):
    """a=0 original-input moment assembly, reported in the legacy format."""
    parts = source.complex_components(psi)
    vertices = source.real_components(xyz, (13, 3))
    mp = source.replay_context([x for pair in parts for x in pair], vertices.flat)
    matter = [mp.mpc(source.mp_real(mp, re), source.mp_real(mp, im)) for re, im in parts]
    points = mp.matrix([[source.mp_real(mp, x) for x in row] for row in vertices])
    g, mass, potential = source.moment_system(mp, matter, points, edges, tetrahedra, potential=True)
    return (np.array(g.tolist(), dtype=float), np.array(mass.tolist(), dtype=float),
            source.reported(potential, 'potential'))


def replay_phase_space(q68, v68, xyz, edges, tetrahedra, frame):
    """Replay the original tangent before rounding a reduced Gram matrix.

    All field-dependent products, the gauge rechart, minimization and density
    remain in a private high-precision context. The producer instead uses its
    positive quadrature factor. Agreement is a finite numerical diagnostic.
    The a=0 moments require an exact original vertex-gradient edge field;
    small nonzero radiative fields are unsupported and explicitly refused.
    """
    configuration = source.real_components(q68, (68,))
    tangent = source.real_components(v68, (68,))
    vertices = source.real_components(xyz, (13, 3))
    coordinates = source.real_components(frame, (42, 30))
    original_potential = source.require_original_gradient(configuration[:42], edges)
    mp = source.replay_context(configuration, tangent, vertices.flat, coordinates.flat)
    real = lambda value: source.mp_real(mp, value)
    q68, v68 = mp.matrix([real(v) for v in configuration]), mp.matrix([real(v) for v in tangent])
    points = mp.matrix([[real(v) for v in row] for row in vertices])
    frame = mp.matrix([[real(v) for v in row] for row in coordinates])
    psi = [mp.mpc(q68[42+i], q68[55+i]) for i in range(13)]
    d = mp.matrix(42, 13)
    for row, (i, j) in enumerate(edges):
        d[row, i], d[row, j] = -1, 1
    _, mass, _ = source.moment_system(mp, [mp.mpf(0)]*13, points, edges, tetrahedra)
    laplacian = d.T*mass*d
    augmented = mp.matrix(14)
    augmented[:13, :13] = laplacian
    for i in range(13):
        augmented[i, 13] = augmented[13, i] = 1
    def gauge_parameter(edges_velocity):
        rhs = mp.matrix(14, 1)
        rhs[:13, :] = -d.T*mass*edges_velocity
        return mp.lu_solve(augmented, rhs)[:13, :]
    # Exact cycle closure already proves a=D*original_potential. Its Coulomb
    # shift is therefore the negative mean-zero potential and its edge field
    # is exactly zero, without a numerical solve whose residual could
    # underflow when reporting a large but physically pure-gauge input.
    xi = mp.matrix([-real(value) for value in original_potential])
    xidot = gauge_parameter(v68[:42, :])
    phase = [mp.exp(mp.j*xi[i]/4) for i in range(13)]
    scalar = [phase[i]*psi[i] for i in range(13)]
    scalar_velocity = [phase[i]*(mp.mpc(v68[42+i], v68[55+i])+mp.j*xidot[i]*psi[i]/4)
                       for i in range(13)]
    ac, av = mp.matrix(42, 1), v68[:42, :]+d*xidot
    full_q = mp.matrix(list(ac)+[mp.re(v) for v in scalar]+[mp.im(v) for v in scalar])
    full_v = mp.matrix(list(av)+[mp.re(v) for v in scalar_velocity]+[mp.im(v) for v in scalar_velocity])
    section = mp.matrix(68, 56)
    section[:42, :30], section[42:, 30:] = frame, mp.eye(26)
    require(max(abs(v) for v in frame.T*frame-mp.eye(30)) <= mp.mpf('4e-12'),
            'orthonormal Coulomb frame')
    require(max(abs(v) for v in d.T*mass*frame) <= mp.mpf('2e-12'), 'mass Coulomb frame')
    g, _, potential = source.moment_system(mp, scalar, points, edges, tetrahedra, potential=True)
    # A different, nonorthonormal gauge basis verifies basis independence.
    basis = mp.matrix([[-1]*12]+np.eye(12, dtype=int).tolist())
    vertical = mp.matrix(68, 13)
    vertical[:42, :] = d
    for i, value in enumerate(scalar):
        vertical[42+i, i], vertical[55+i, i] = -mp.im(value)/4, mp.re(value)/4
    vertical = vertical*basis
    gamma, eta_map, log_rho = source.reduced_moments(mp, g, vertical, section)
    q, velocity = section.T*full_q, section.T*full_v
    multiplier = eta_map*velocity
    eta = basis*multiplier
    # Obtain the cotangent from the minimized full action, and compare with
    # the separately formed reduced quadratic form in working precision.
    momentum = section.T*g*(section*velocity+vertical*multiplier)
    generator = mp.matrix([0]*30+[-mp.im(v) for v in scalar]+[mp.re(v) for v in scalar])
    scalar_values = {
        'constant_moment_map': (momentum.T*generator)[0],
        'minimizer_defect': max(abs(v) for v in eta+xidot),
        'coulomb_defect': max(abs(v) for v in d.T*mass*ac),
        'cotangent_identity_defect': max(abs(v) for v in momentum-gamma*velocity),
        'log_rho': log_rho, 'potential': potential,
        'gamma_min_eigenvalue': mp.eigsy(gamma, eigvals_only=True)[0],
    }
    vectors = {'q': q, 'velocity': velocity, 'momentum': momentum, 'full_q': full_q,
        'full_velocity': full_v, 'gauge_parameter': xi, 'gauge_parameter_velocity': xidot,
        'transformed_scalar_potential': -xidot, 'schur_scalar_potential': eta}
    return {**{key: np.array([source.reported(v, key) for v in value]) for key, value in vectors.items()},
            **{key: source.reported(value, key) for key, value in scalar_values.items()}}


def circle_observables(q, p, sigma, hbar=1):
    """Independent positive circle integrals, with no Bessel implementation.

    Analytic contour translation gives exp(z-A)*I0e(z). At large z the
    coordinate u=2*sqrt(z)*sin(theta/2) resolves the concentrated kernel.
    The cutoff u=16 leaves a Gaussian tail below the reporting tolerance;
    this high-precision quadrature is a diagnostic, not an interval proof.
    """
    x, p = [Q(v) for v in q[30:]], [Q(v) for v in p[30:]]
    s, h = Q(sigma), Q(hbar)
    xx, pp = sum(v*v for v in x), sum(v*v for v in p)
    b = sum(p[i]*(-x[i+13]) + p[i+13]*x[i] for i in range(13))/h
    a = xx/(4*s*s)+s*s*pp/(h*h)
    discriminant = a*a-b*b
    require(discriminant >= 0, 'nonnegative exact circle discriminant')
    mp = mpmath.mp.clone()
    mp.dps = 70+max(0, len(str(abs(a.numerator)))-len(str(a.denominator)))
    def real(v):
        return mp.mpf(v.numerator)/v.denominator
    aa, bb, ss = real(a), real(b), real(s*s)
    z = mp.sqrt(real(discriminant))
    if z > 128:
        def weight(u):
            return mp.exp(-u*u/2)/mp.sqrt(1-u*u/(4*z))
        intervals = [0, 1, 4, 8, 16]
        mass = mp.quad(weight, intervals)
        deficit = mp.quad(lambda u: u*u/2*weight(u), intervals)/mass
        norm = mp.exp(z-aa)*mass/(mp.pi*mp.sqrt(z))
    else:
        def weight(theta):
            return mp.exp(z*(mp.cos(theta)-1))
        intervals = [0, mp.pi/2, mp.pi]
        mass = mp.quad(weight, intervals)
        deficit = mp.quad(lambda t: z*(1-mp.cos(t))*weight(t), intervals)/mass
        norm = mp.exp(z-aa)*mass/mp.pi
    moment = 26*ss+real(xx)/2-2*ss*ss*real(pp/(h*h))+2*ss*(z-deficit)
    return {'A':float(aa),'B':float(bb),'norm_squared':float(norm),
            'scalar_radius_numeric':float(moment),'sigma':str(s)}


def close_observable(actual, expected, name):
    """Each scalar of the supplied packet has its own relative tolerance."""
    value = numbers(actual, (), name).item()
    require(value == 0 if expected == 0 else abs((value-expected)/expected) <= 2e-11, name)


def initial_algebra(volume):
    # Compute the two scalar mass rows from simplex moments and node incidence.
    c = 2*(3*simplex((2,0,0,0))-3*simplex((1,1,0,0)))
    boundary = 2*Q(1,4)*(3*simplex((1,1,0,0))-(Q(1,4)-simplex((1,1,0,0))))
    require(c == Q(3,10) and boundary == -Q(1,40), 'initial simplex cotangent')
    center = [c*x for x in volume]; edge = [boundary*x for x in volume]
    momentum_norm = [(c*c+12*boundary**2)*(volume[0]**2+5*volume[1]**2),
                     (c*c+12*boundary**2)*2*volume[0]*volume[1]]
    widths=[]
    for sigma in (Q(1,4),Q(1,2),Q(1)):
        upper=Q(13)/(4*sigma**2)+(c*c+12*boundary**2)*sigma**2*18**2
        require(upper < 64,'projection bound range')
        widths.append({'sigma':str(sigma),'A_in_Qsqrt5':[str(Q(13)/(4*sigma**2)+sigma**2*momentum_norm[0]),str(sigma**2*momentum_norm[1])],
            'A_rational_upper':str(upper),'B':'0','projection_norm_squared_lower':'1/64',
            'seed_position_variance':str(sigma**2),'seed_momentum_variance_hbar1':str(Q(1)/(4*sigma**2))})
    return {'volume_in_Qsqrt5':[str(x) for x in volume], 'radial_velocity_in_Qsqrt5':['7/80','3/80'],
        'center_gauge_factor':'12/13','boundary_gauge_factor':'-1/13',
        'scalar_center_real':'1','scalar_boundary_real':'1',
        'momentum_center_imaginary_in_Qsqrt5':[str(x) for x in center],
        'momentum_boundary_imaginary_in_Qsqrt5':[str(x) for x in edge],
        'scalar_momentum_norm_squared_in_Qsqrt5':[str(x) for x in momentum_norm],
        'radiative_center_and_momentum':'0','constant_moment_map':'0',
        'proof_labels':['prop:whitney-packet-coulomb','prop:whitney-packet-neutral-projection'], 'widths':widths}


def expected_contract():
    return {'configuration_dimension':56,'radiative_dimension':30,'scalar_real_dimension':26,
        'seed_covariance':'sigma^2 I_56; full rank; no classical fixed-space restriction',
        'projected_covariance':'not the seed covariance; circle averaging changes scalar moments',
        'seed_formula':'(2*pi*sigma^2)^(-14)*exp(-|x-q|^2/(4*sigma^2)+i*p.(x-q)/hbar)',
        'state_formula':'F=rho^(-1/2)*P0(seed)/||P0(seed)||; rho=sqrt(det(gamma))',
        'quantization':'same full interacting Laplace-Beltrami Hamiltonian and metric volume',
        'circle_projection':'exact group integral; floating Bessel values are diagnostics',
        'physical_comparison':False,'empirical_prediction':False,'observer_preparation':False,
        'quantum_propagation_certified':False,'covariance_evolution_computed':False,
        'full_56D_preparation':True,'neutral_strong_domain_state':True,
        'classical_samples_are_quantum_means':False,'scalar_one_point_mean':'0 by residual U(1) invariance',
        'sample_role':'two separate packet preparations from adjacent classical samples; not quantum propagation',
        'numeric_scope':'float64 classical data, phase-space maps and exact-degree simplex integration; no interval trajectory enclosure',
        'requested_propagation_target':{'T':'1/40','hbar':'1','norm_error':'1/10','status':'NOT_CERTIFIED','residual_upper_bound':None,'all_tail_bound':None}}


def verify(packet):
    keys(packet, {'schema','scope','source_pins','parameters','contract','initial_exact','orthonormal_coulomb_frame','samples'},'packet')
    equal(packet['schema'],'oph.whitney_quantum_packet.v1','schema')
    equal(packet['scope'],'FULL_56D_NEUTRAL_PACKET_PREPARATION__NO_QUANTUM_PROPAGATION_CERTIFICATE','scope')
    keys(packet['source_pins'],PIN_PATHS,'source pins')
    for path in PIN_PATHS:
        equal(packet['source_pins'][path],hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),'source pin '+path)
    equal(packet['parameters'],{'e':'1/4','m_squared':'1/2','g':'1/4','hbar':'1'},'parameters')
    equal(packet['contract'],expected_contract(),'preparation and non-propagation contract')
    volume,xyz,edges,tets=source.exact_geometry()
    equal(packet['initial_exact'],initial_algebra(volume),'exact initial algebra and projection bound')
    frame=numbers(packet['orthonormal_coulomb_frame'],(42,30),'Coulomb frame')
    require(type(packet['samples']) is list and len(packet['samples'])==2,'sample census')
    raw=parent.load(ROOT/PARENT_PATH)
    diagnostics=[]
    for index,row in enumerate(packet['samples']):
        require(type(row) is dict,'sample object')
        original=raw['samples'][index]
        expected=replay_phase_space(original['q'],original['velocity'],xyz,edges,tets,frame)
        keys(row,set(expected)|{'parent_sample_index','model_time','widths'},'sample')
        equal(row['parent_sample_index'],index,'parent index')
        equal(row['model_time'],index/40,'model sample time')
        for key,value in expected.items():
            close(row[key],value,'independent phase-space '+key)
        require(abs(expected['constant_moment_map'])<1e-9 and expected['minimizer_defect']<1e-9,'parent numerical Gauss compatibility')
        require(type(row['widths']) is list and len(row['widths'])==3,'width census')
        for width,sigma in zip(row['widths'],(.25,.5,1),strict=True):
            # The phase-space replay above authenticates the numerical state.
            # Evaluate that *reported* representative so a nearly cancelling
            # charge cannot hide behind an absolute tolerance on another solve.
            expected_width=circle_observables(row['q'],row['momentum'],sigma)
            keys(width,expected_width,'width')
            equal(width['sigma'],expected_width['sigma'],'width parameter')
            for key in ('A','B','norm_squared','scalar_radius_numeric'):
                close_observable(width[key],expected_width[key],'independent circle '+key)
        diagnostics.append({'model_time':row['model_time'],'gauss_abs':abs(expected['constant_moment_map']),
            'minimizer_defect':expected['minimizer_defect'],
            'projection_norm_squared_numeric':[x['norm_squared'] for x in row['widths']]})
    # Fresh full parent replay is deliberately last: malformed packets fail
    # before the longer numerical trajectory validation, never instead of it.
    parent_result=parent.verify(raw)
    require(parent_result['accepted'] is True,'independent parent trajectory')
    return {'accepted':True,'configuration_dimension':56,'radiative_dimension':30,
        'prepared_classical_samples':2,'widths':['1/4','1/2','1'],
        'initial_projection_norm_squared_lower':'1/64','full_56D_preparation':True,
        'neutral_strong_domain_state':True,'quantum_propagation_certified':False,
        'covariance_evolution_computed':False,'physical_comparison':False,
        'analytic_domain_proved_by_numeric_replay':False,'diagnostics':diagnostics}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,default=OUTPUT)
    args=parser.parse_args()
    print(json.dumps(verify(load(args.receipt)),sort_keys=True))
