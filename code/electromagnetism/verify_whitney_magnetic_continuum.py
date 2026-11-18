"""Independent replay of the fixed-background complex continuum packet.

No magnetic producer is imported. Meshes use the earlier independent lattice
chain construction; phases are direct affine line integrals, and Gauss--Jacobi
integration differs from the producer's Legendre rule. Continuum estimates
remain analytic paper proofs; the replay checks only finite numerical data.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from itertools import combinations
import hashlib
import json
from math import isfinite
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized, spsolve
from scipy.special import roots_jacobi

import verify_whitney_real_continuum as geometry

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent/'runtime/whitney_magnetic_continuum_receipt.json'
PINS = {
    'Lean/Screen/SeamCurrentEdge30Moment.lean',
    'Lean/ObserverPatchHolography/CoreAxioms.lean',
    'Lean/Screen/MagneticMatterStability.lean',
    'code/electromagnetism/whitney_real_continuum.py',
    'code/electromagnetism/verify_whitney_real_continuum.py',
    'code/electromagnetism/whitney_magnetic_continuum.py',
    'code/electromagnetism/verify_whitney_magnetic_continuum.py',
    'code/electromagnetism/test_whitney_magnetic_continuum.py',
}
SCOPE = 'FIXED_UNIFORM_MAGNETIC_BACKGROUND__CONDITIONAL_COMPLEX_MATTER_CONTINUUM__FINITE_NUMERICAL_CHECKS'
ANALYTIC = {
    'uniform_mesh_paper_proof': True,
    'conforming_magnetic_ritz_paper_proof': True,
    'conditional_complex_trajectory_bound': True,
    'continuum_reference_existence_assumed': True,
    'external_magnetic_background_fixed': True,
    'self_consistent_maxwell_backreaction': False,
    'numerical_error_interval_certified': False,
    'physical_source_or_clock_selected': False,
    'formalized_in_lean': False,
    'pointwise_cubic_bounds_formalized_in_lean': True,
}
PARAMETERS = {'e': '1/4', 'm2': '1/2', 'g': '1/4',
              'B': ['0', '0', '1'], 'c': ['0', '0', '0']}
LEAN_DECLARATIONS = ('norm_square_difference_le', 'norm_cubic_difference_le', 'norm_cubic_difference_on_ball')
LABELS = ('lem:whitney-magnetic-conforming-approximation', 'thm:whitney-magnetic-continuum-trajectory')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode('ascii')


def exact(actual, expected, message):
    require(canonical(actual) == canonical(expected), message)


def finite(value):
    require(type(value) in (float, int) and isfinite(value), 'finite numeric diagnostic required')
    return float(value)


def close(actual, expected, message, atol=2e-9, rtol=2e-7):
    require(abs(finite(actual)-expected) <= atol+rtol*abs(expected), message)


def load(path=OUTPUT):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def number(value):
        result = float(value)
        require(isfinite(result), 'nonfinite JSON number')
        return result
    def constant(_):
        raise ValueError('nonfinite JSON constant')
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs,
                      parse_float=number, parse_constant=constant)


def vector_potential(q, shift=None):
    result = np.stack((-q[..., 1]/2, q[..., 0]/2, np.zeros(q.shape[:-1])), axis=-1)
    return result if shift is None else result+shift


def rule(order):
    nodes, weights = [], []
    for power in (2, 1, 0):
        node, weight = roots_jacobi(order, 0, power)
        nodes.append((node+1)/2); weights.append(weight/2**(power+1))
    r, s, t = np.meshgrid(*nodes, indexing='ij')
    wr, ws, wt = np.meshgrid(*weights, indexing='ij')
    return (np.stack((1-r, r*(1-s), r*s*(1-t), r*s*t), axis=-1).reshape(-1, 4),
            (wr*ws*wt).ravel())


def direct_basis(x, grad, lam, shift=None, omit_phase_gradient=False):
    q = np.einsum('qi,tic->tqc', lam, x)
    ai = vector_potential(x, shift)
    theta = np.sum(ai[:, None]*(q[:, :, None]-x[:, None]), axis=-1)
    exponential = np.exp(.25j*theta)
    value = exponential*lam[None]
    ordinary = exponential[..., None]*grad[:, None]
    if not omit_phase_gradient:
        ordinary += .25j*value[..., None]*ai[:, None]
    derivative = ordinary-.25j*value[..., None]*vector_potential(q, shift)[:, :, None]
    return q, value, derivative, ordinary


def system(n, order=7, shift=None):
    vertices, cells, keys = geometry.mesh(n)
    x, grad, volume, _, _ = geometry.assemble(vertices, cells)
    lam, w = rule(order)
    q, basis, derivative, ordinary = direct_basis(x, grad, lam, shift)
    weights = volume[:, None]*6*w
    matrices = []
    for differentiated in (False, True):
        rows, cols, entries = [], [], []
        for i in range(4):
            for j in range(4):
                if differentiated:
                    term = np.sum(derivative[:, :, i].conj()*derivative[:, :, j], axis=-1)
                else:
                    term = basis[:, :, i].conj()*basis[:, :, j]
                entries.extend(np.sum(term*weights, axis=1))
                rows.extend(cells[:, i]); cols.extend(cells[:, j])
        matrices.append(coo_matrix((entries, (rows, cols)), shape=(len(vertices), len(vertices))).tocsr())
    return dict(vertices=vertices, cells=cells, keys=keys, x=x, grad=grad, lam=lam,
                q=q, value=basis, covariant=derivative, derivative=ordinary,
                weights=weights, mass=matrices[0], stiff=matrices[1])


def scatter(local, data):
    result = np.zeros(len(data['vertices']), dtype=complex)
    for i in range(4):
        np.add.at(result, data['cells'][:, i], local[:, i])
    return result


def reference(q):
    phase = .3*q[..., 0]-.2*q[..., 1]+.4*q[..., 2]
    exponential = np.cos(phase)+1j*np.sin(phase)
    amplitude = 1+np.sum(q*q, axis=-1)/5
    return (amplitude*exponential,
            exponential[..., None]*(2*q/5+1j*amplitude[..., None]*np.array([.3, -.2, .4])))


def fields(coefficients, data, derivative='covariant'):
    local = coefficients[data['cells']]
    return (np.sum(local[:, None]*data['value'], axis=-1),
            np.sum(local[:, None, :, None]*data[derivative], axis=2))


def ritz_check(data):
    ref, derivative = reference(data['q'])
    covariant = derivative-.25j*vector_potential(data['q'])*ref[..., None]
    local = np.sum(data['covariant'].conj()*covariant[:, :, None], axis=-1)
    local += .5*data['value'].conj()*ref[:, :, None]
    load = scatter(np.sum(local*data['weights'][..., None], axis=1), data)
    matrix = data['stiff']+.5*data['mass']
    projected = spsolve(matrix, load)
    result = {'ritz_weak_residual': float(np.max(abs(matrix@projected-load)))}
    for name, coefficients in [('nodal', reference(data['vertices'])[0]), ('ritz', projected)]:
        value, ordinary = fields(coefficients, data, 'derivative')
        error = value-ref
        grad_error = ordinary-derivative
        cov_error = grad_error-.25j*vector_potential(data['q'])*error[..., None]
        l2 = np.sum(abs(error)**2*data['weights'])
        result[name+'_l2_error'] = float(np.sqrt(l2))
        result[name+'_h1_error'] = float(np.sqrt(l2+np.sum(np.sum(abs(grad_error)**2, axis=-1)*data['weights'])))
        result[name+'_magnetic_error'] = float(np.sqrt(.5*l2+np.sum(np.sum(abs(cov_error)**2, axis=-1)*data['weights'])))
    return result


def nonlinear_load(u, data):
    values = np.sum(u[data['cells']][:, None]*data['value'], axis=-1)
    local = data['value'].conj()*(values.real**2+values.imag**2)[..., None]*values[..., None]
    return scatter(np.sum(local*data['weights'][..., None], axis=1), data)


def energy(u, velocity, data):
    values = np.sum(u[data['cells']][:, None]*data['value'], axis=-1)
    kinetic = np.real(np.vdot(velocity, data['mass']@velocity))/2
    quadratic = np.real(np.vdot(u, (data['stiff']+.5*data['mass'])@u))/2
    return float(kinetic+quadratic+np.sum((values.real**2+values.imag**2)**2*data['weights'])/16)


def initial_state(vertices):
    profile = np.maximum(1-(vertices*vertices).sum(axis=1), 0)**4/5
    angle = 3*vertices[:, 0]/10+vertices[:, 1]/5
    u = profile*(np.cos(angle)+1j*np.sin(angle))
    return np.concatenate((u, 2j*u/5))


def independent_trajectory(data):
    size = len(data['vertices'])
    solve = factorized(data['mass'].tocsc())
    matrix = data['stiff']+.5*data['mass']
    def rhs(_, state):
        return np.r_[state[size:], -solve(matrix@state[:size]+nonlinear_load(state[:size], data)/4)]
    initial = initial_state(data['vertices'])
    solution = solve_ivp(rhs, (0., .125), initial, method='DOP853', rtol=2e-12, atol=2e-13)
    require(solution.success, 'independent finite trajectory integration')
    return initial, solution.y[:, -1]


def controls_check():
    data = system(1)
    shift = np.array([1/5, -1/7, 1/9])
    q, basis, covariant, _ = direct_basis(data['x'], data['grad'], data['lam'], shift)
    coefficients = reference(data['vertices'])[0]
    values, derivative = fields(coefficients, data)
    shifted = coefficients*np.exp(.25j*data['vertices']@shift)
    local = shifted[data['cells']]
    transformed_values = np.sum(local[:, None]*basis, axis=-1)
    transformed_derivative = np.sum(local[:, None, :, None]*covariant, axis=2)
    phase = np.exp(.25j*q@shift)
    wrong = direct_basis(data['x'], data['grad'], data['lam'], omit_phase_gradient=True)[2]
    wrong_field = np.sum(coefficients[data['cells']][:, None, :, None]*wrong, axis=2)
    faces = defaultdict(list)
    for t, cell in enumerate(data['cells']):
        for face in combinations(cell, 3):
            faces[tuple(sorted(face))].append(t)
    differences = []
    for face, cells in faces.items():
        if len(cells) == 2:
            values_at_face = []
            for cell in cells:
                lam = np.array([[float(i in face)/3 for i in data['cells'][cell]]])
                value = direct_basis(data['x'][cell:cell+1], data['grad'][cell:cell+1], lam)[1]
                values_at_face.append(np.sum(value[0, 0]*coefficients[data['cells'][cell]]))
            differences.append(abs(values_at_face[0]-values_at_face[1]))
    return {'shared_faces_checked': len(differences), 'shared_face_trace_error': float(max(differences)),
            'gauge_value_error': float(np.max(abs(transformed_values-phase*values))),
            'gauge_covariant_error': float(np.max(abs(transformed_derivative-phase[..., None]*derivative))),
            'omitted_phase_gradient_gap': float(np.sqrt(np.sum(np.sum(abs(wrong_field-derivative)**2,
                                                                      axis=-1)*data['weights'])))}


def verify(packet):
    require(type(packet) is dict and set(packet) == {'schema', 'scope', 'source_pins', 'analytic_scope',
                                                   'parameters', 'mesh_checks', 'controls', 'trajectories'}, 'packet schema')
    exact(packet['schema'], 'oph.whitney_magnetic_continuum.v1', 'schema version')
    exact(packet['scope'], SCOPE, 'scientific scope')
    exact(packet['analytic_scope'], ANALYTIC, 'analytic versus finite evidence scope')
    exact(packet['parameters'], PARAMETERS, 'fixed external field and coupling parameters')
    require(type(packet['source_pins']) is dict and set(packet['source_pins']) == PINS, 'source pin census')
    for path in PINS:
        exact(packet['source_pins'][path], hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), 'source pin: '+path)
    lean = (ROOT/'Lean/Screen/MagneticMatterStability.lean').read_text(encoding='utf-8')
    for declaration in LEAN_DECLARATIONS:
        require('theorem '+declaration+' ' in lean, 'pointwise Lean declaration missing')
    require(type(packet['mesh_checks']) is list and len(packet['mesh_checks']) == 3, 'mesh census')
    mesh_keys = {'n', 'vertices', 'tetrahedra', 'mesh_sha256', 'whitney_reproduction_error', 'ritz_weak_residual',
                 'nodal_l2_error', 'nodal_h1_error', 'nodal_magnetic_error',
                 'ritz_l2_error', 'ritz_h1_error', 'ritz_magnetic_error'}
    for row, n, size in zip(packet['mesh_checks'], (1, 2, 4), (13, 55, 309), strict=True):
        require(type(row) is dict and set(row) == mesh_keys, 'mesh record schema')
        for key, value in [('n', n), ('vertices', size), ('tetrahedra', 20*n**3)]:
            exact(row[key], value, 'exact mesh count: '+key)
        for key in mesh_keys-{'n', 'vertices', 'tetrahedra', 'mesh_sha256'}:
            require(finite(row[key]) >= 0, 'nonnegative error diagnostic')
    control_keys = {'shared_faces_checked', 'shared_face_trace_error', 'gauge_value_error',
                    'gauge_covariant_error', 'omitted_phase_gradient_gap'}
    require(type(packet['controls']) is dict and set(packet['controls']) == control_keys, 'control schema')
    exact(packet['controls']['shared_faces_checked'], 30, 'shared face census')
    for key in control_keys-{'shared_faces_checked'}:
        require(finite(packet['controls'][key]) >= 0, 'nonnegative control diagnostic')
    require(type(packet['trajectories']) is list and len(packet['trajectories']) == 2, 'trajectory census')
    trajectory_keys = {'n', 'steps', 'time_interval', 'initialization', 'final_state_real', 'final_state_imag',
                       'initial_energy', 'final_energy', 'initial_charge_pairing', 'final_charge_pairing'}
    for row, steps in zip(packet['trajectories'], (32, 64), strict=True):
        require(type(row) is dict and set(row) == trajectory_keys, 'trajectory schema')
        for key, value in [('n', 2), ('steps', steps), ('time_interval', ['0', '1/8']),
                           ('initialization', 'nodal compact complex pulse; not Ritz data for a supplied continuum solution')]:
            exact(row[key], value, 'trajectory metadata: '+key)
        for key in ('final_state_real', 'final_state_imag'):
            require(type(row[key]) is list and len(row[key]) == 110, 'complex state census')
            for value in row[key]:
                finite(value)
    diagnostics, cached = [], None
    for row in packet['mesh_checks']:
        data = system(row['n'])
        exact(row['mesh_sha256'], hashlib.sha256(canonical({'keys': data['keys'], 'cells': data['cells'].tolist()})).hexdigest(),
              'independent conforming mesh hash')
        # Direct Whitney reconstruction independently checks the declared uniform B.
        reconstructed = np.zeros_like(data['q'])
        for i, j in combinations(range(4), 2):
            edge = np.sum(vector_potential(data['x'][:, i])*(data['x'][:, j]-data['x'][:, i]), axis=1)
            reconstructed += edge[:, None, None]*(data['lam'][None, :, i, None]*data['grad'][:, None, j]
                                                 -data['lam'][None, :, j, None]*data['grad'][:, None, i])
        require(np.max(abs(reconstructed-vector_potential(data['q']))) < 1e-12, 'Whitney reproduction')
        require(abs(finite(row['whitney_reproduction_error'])) < 1e-12, 'reported reproduction defect')
        expected = ritz_check(data)
        for key, value in expected.items():
            close(row[key], value, 'independent magnetic Ritz replay: '+key)
        require(row['ritz_magnetic_error'] <= row['nodal_magnetic_error']+1e-9, 'Ritz best approximation')
        require(row['ritz_weak_residual'] < 1e-10, 'Ritz weak residual')
        diagnostics.append({'n': row['n'], **expected})
        if row['n'] == 2:
            cached = data
    for coarse, fine in zip(packet['mesh_checks'], packet['mesh_checks'][1:]):
        require(fine['ritz_h1_error'] < .8*coarse['ritz_h1_error'], 'finite complex Ritz refinement control')
    controls = controls_check()
    for key, value in controls.items():
        if key != 'shared_faces_checked':
            close(packet['controls'][key], value, 'independent gauge/derivative/trace control: '+key)
    for key in ('shared_face_trace_error', 'gauge_value_error', 'gauge_covariant_error'):
        require(abs(finite(packet['controls'][key])) < 1e-10 and controls[key] < 1e-10, 'conforming gauge covariance')
    require(controls['omitted_phase_gradient_gap'] > .01, 'omitted phase derivative falsifier')
    initial, endpoint = independent_trajectory(cached)
    size = len(cached['vertices'])
    initial_energy = energy(initial[:size], initial[size:], cached)
    initial_charge = float(np.imag(np.vdot(initial[:size], cached['mass']@initial[size:])))
    require(initial_charge > 1e-5, 'nonzero external-field charge witness')
    endpoint_errors = []
    for row in packet['trajectories']:
        state = np.array(row['final_state_real'])+1j*np.array(row['final_state_imag'])
        error = float(np.max(abs(state-endpoint)))
        require(error < 2e-7, 'independent DOP853 endpoint replay')
        close(row['initial_energy'], initial_energy, 'initial energy')
        final_energy = energy(state[:size], state[size:], cached)
        close(row['final_energy'], final_energy, 'final energy')
        close(row['initial_charge_pairing'], initial_charge, 'initial charge pairing')
        final_charge = float(np.imag(np.vdot(state[:size], cached['mass']@state[size:])))
        close(row['final_charge_pairing'], final_charge, 'final charge pairing')
        require(abs(final_energy-initial_energy) < 2e-8, 'finite energy conservation control')
        require(abs(final_charge-initial_charge) < 2e-9, 'finite charge conservation control')
        endpoint_errors.append(error)
    require(endpoint_errors[1] < endpoint_errors[0]/8, 'RK4 refinement toward independent flow')
    return {'accepted': True, 'scope': SCOPE, **ANALYTIC,
            'analytic_proof_labels': list(LABELS), 'lean_declarations': list(LEAN_DECLARATIONS),
            'lean_theorem_count': 3, 'refinements': [1, 2, 4], 'tetrahedra': [20, 160, 1280],
            'vertices': [13, 55, 309], 'trajectory_count': 2, 'complex_state_length': 110,
            'nonzero_magnetic_field': True, 'nonzero_charge_witness': True,
            'numeric_diagnostics': {'ritz': diagnostics, 'controls': controls, 'endpoint_errors': endpoint_errors}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True, allow_nan=False))
