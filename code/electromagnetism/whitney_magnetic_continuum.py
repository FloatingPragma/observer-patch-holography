"""Finite checks of conforming complex matter in a fixed uniform magnetic field.

The paper proves conditional spatial trajectory convergence. This packet uses
numerical spatial quadrature and time integration, not interval enclosures.
The external potential is held fixed; Maxwell backreaction is not solved.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized, spsolve

import whitney_real_continuum as real

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / 'runtime/whitney_magnetic_continuum_receipt.json'
SCOPE = 'FIXED_UNIFORM_MAGNETIC_BACKGROUND__CONDITIONAL_COMPLEX_MATTER_CONTINUUM__FINITE_NUMERICAL_CHECKS'
PINS = (
    'Lean/Screen/SeamCurrentEdge30Moment.lean',
    'Lean/ObserverPatchHolography/CoreAxioms.lean',
    'Lean/Screen/MagneticMatterStability.lean',
    'code/electromagnetism/whitney_real_continuum.py',
    'code/electromagnetism/verify_whitney_real_continuum.py',
    'code/electromagnetism/whitney_magnetic_continuum.py',
    'code/electromagnetism/verify_whitney_magnetic_continuum.py',
    'code/electromagnetism/test_whitney_magnetic_continuum.py',
)
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


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode('ascii')


def potential(points, shift=None):
    return .5*np.cross(np.array([0., 0., 1.]), points)+(np.zeros(3) if shift is None else shift)


def basis(x, grad, lam, *, shift=None, omit_phase_gradient=False):
    """Full Whitney-edge phase derivative; arrays cell/point/vertex/component."""
    q = np.einsum('qi,tic->tqc', lam, x)
    delta = x[:, None, :, :]-x[:, :, None, :]
    edge = np.sum(potential((x[:, None]+x[:, :, None])/2, shift)*delta, axis=-1)
    theta = np.einsum('tij,qj->tqi', edge, lam)
    dtheta = np.einsum('tij,tjc->tic', edge, grad)
    phase = np.exp(.25j*theta)
    value = phase*lam[None]
    derivative = phase[..., None]*grad[:, None]
    if not omit_phase_gradient:
        derivative += .25j*value[..., None]*dtheta[:, None]
    covariant = derivative-.25j*value[..., None]*potential(q, shift)[:, :, None]
    reconstructed = np.zeros_like(q)
    for i, j in combinations(range(4), 2):
        reconstructed += edge[:, i, j, None, None]*(lam[None, :, i, None]*grad[:, None, j]
                                                    -lam[None, :, j, None]*grad[:, None, i])
    return q, value, covariant, derivative, reconstructed


def system(n, order=6, shift=None):
    vertices, cells, keys = real.mesh(n)
    x, grad, volume = real.geometry(vertices, cells)
    lam, weights = real.quadrature(order)
    q, value, covariant, derivative, reconstructed = basis(x, grad, lam, shift=shift)
    weights = 6*volume[:, None]*weights[None]
    mass_local = np.einsum('tqi,tqj,tq->tij', value.conj(), value, weights)
    stiff_local = np.einsum('tqic,tqjc,tq->tij', covariant.conj(), covariant, weights)
    rows = np.broadcast_to(cells[:, :, None], (len(cells), 4, 4)).ravel()
    cols = np.broadcast_to(cells[:, None, :], (len(cells), 4, 4)).ravel()
    shape = (len(vertices), len(vertices))
    mass = coo_matrix((mass_local.ravel(), (rows, cols)), shape=shape).tocsr()
    stiff = coo_matrix((stiff_local.ravel(), (rows, cols)), shape=shape).tocsr()
    return dict(vertices=vertices, cells=cells, keys=keys, x=x, grad=grad, lam=lam,
                q=q, weights=weights, value=value, covariant=covariant,
                derivative=derivative, mass=mass, stiff=stiff,
                whitney_error=float(np.max(abs(reconstructed-potential(q, shift)))))


def scatter(local, data):
    out = np.zeros(len(data['vertices']), dtype=complex)
    np.add.at(out, data['cells'].ravel(), local.ravel())
    return out


def reference(q):
    wave = np.array([.3, -.2, .4])
    phase = np.exp(1j*q@wave)
    amplitude = 1+.2*np.sum(q*q, axis=-1)
    value = amplitude*phase
    derivative = phase[..., None]*(.4*q+1j*amplitude[..., None]*wave)
    return value, derivative


def field(coefficients, data, derivative='covariant'):
    return (np.einsum('ti,tqi->tq', coefficients[data['cells']], data['value']),
            np.einsum('ti,tqic->tqc', coefficients[data['cells']], data[derivative]))


def ritz_record(n):
    data = system(n)
    ref, derivative = reference(data['q'])
    covariant = derivative-.25j*potential(data['q'])*ref[..., None]
    local = np.einsum('tqic,tqc,tq->ti', data['covariant'].conj(), covariant, data['weights'])
    local += .5*np.einsum('tqi,tq,tq->ti', data['value'].conj(), ref, data['weights'])
    load = scatter(local, data)
    matrix = data['stiff']+.5*data['mass']
    projected = spsolve(matrix, load)
    record = {'n': n, 'vertices': len(data['vertices']), 'tetrahedra': len(data['cells']),
              'mesh_sha256': hashlib.sha256(canonical({'keys': data['keys'], 'cells': data['cells'].tolist()})).hexdigest(),
              'whitney_reproduction_error': data['whitney_error'],
              'ritz_weak_residual': float(np.max(abs(matrix@projected-load)))}
    for name, coefficients in [('nodal', reference(data['vertices'])[0]), ('ritz', projected)]:
        value, gradient = field(coefficients, data, 'derivative')
        l2 = np.sum(abs(value-ref)**2*data['weights'])
        h1 = l2+np.sum(np.sum(abs(gradient-derivative)**2, axis=-1)*data['weights'])
        magnetic = .5*l2+np.sum(np.sum(abs(gradient-derivative-.25j*potential(data['q'])*(value-ref)[..., None])**2,
                                               axis=-1)*data['weights'])
        record[name+'_l2_error'] = float(np.sqrt(l2))
        record[name+'_h1_error'] = float(np.sqrt(h1))
        record[name+'_magnetic_error'] = float(np.sqrt(magnetic))
    return record


def nonlinear_load(u, data):
    values, _ = field(u, data)
    local = np.einsum('tqi,tq,tq->ti', data['value'].conj(), abs(values)**2*values, data['weights'])
    return scatter(local, data)


def energy(u, velocity, data):
    values, _ = field(u, data)
    return float(.5*np.real(np.vdot(velocity, data['mass']@velocity)
                            +np.vdot(u, (data['stiff']+.5*data['mass'])@u))
                 +.0625*np.sum(abs(values)**4*data['weights']))


def initial_state(vertices):
    bump = .2*np.maximum(1-np.sum(vertices**2, axis=1), 0)**4
    u = bump*np.exp(1j*(.3*vertices[:, 0]+.2*vertices[:, 1]))
    return np.r_[u, .4j*u]


def trajectory(steps):
    data = system(2)
    size = len(data['vertices'])
    solve = factorized(data['mass'].tocsc())
    matrix = data['stiff']+.5*data['mass']
    def rhs(state):
        u, velocity = state[:size], state[size:]
        return np.r_[velocity, -solve(matrix@u+.25*nonlinear_load(u, data))]
    state = initial_state(data['vertices'])
    initial_energy = energy(state[:size], state[size:], data)
    initial_charge = float(np.imag(np.vdot(state[:size], data['mass']@state[size:])))
    step = .125/steps
    for _ in range(steps):
        k1 = rhs(state); k2 = rhs(state+step*k1/2)
        k3 = rhs(state+step*k2/2); k4 = rhs(state+step*k3)
        state += step*(k1+2*k2+2*k3+k4)/6
    return {'n': 2, 'steps': steps, 'time_interval': ['0', '1/8'],
            'initialization': 'nodal compact complex pulse; not Ritz data for a supplied continuum solution',
            'final_state_real': state.real.tolist(), 'final_state_imag': state.imag.tolist(),
            'initial_energy': initial_energy, 'final_energy': energy(state[:size], state[size:], data),
            'initial_charge_pairing': initial_charge,
            'final_charge_pairing': float(np.imag(np.vdot(state[:size], data['mass']@state[size:])))}


def controls():
    data = system(1)
    shift = np.array([.2, -1/7, 1/9])
    shifted = system(1, shift=shift)
    coefficients = reference(data['vertices'])[0]
    value, covariant = field(coefficients, data)
    shifted_value, shifted_covariant = field(coefficients*np.exp(.25j*data['vertices']@shift), shifted)
    phase = np.exp(.25j*data['q']@shift)
    _, _, wrong, _, _ = basis(data['x'], data['grad'], data['lam'], omit_phase_gradient=True)
    wrong_gradient = np.einsum('ti,tqic->tqc', coefficients[data['cells']], wrong)
    face_cells = defaultdict(list)
    for cell_index, cell in enumerate(data['cells']):
        for face in combinations(cell, 3):
            face_cells[tuple(sorted(face))].append(cell_index)
    trace_errors = []
    for face, adjacent in face_cells.items():
        if len(adjacent) != 2:
            continue
        outputs = []
        for cell_index in adjacent:
            cell = data['cells'][cell_index]
            lam = np.array([[1/3 if i in face else 0. for i in cell]])
            result = basis(data['x'][cell_index:cell_index+1], data['grad'][cell_index:cell_index+1], lam)[1]
            outputs.append(np.sum(result[0, 0]*coefficients[cell]))
        trace_errors.append(abs(outputs[0]-outputs[1]))
    return {'shared_faces_checked': len(trace_errors), 'shared_face_trace_error': float(max(trace_errors)),
            'gauge_value_error': float(np.max(abs(shifted_value-phase*value))),
            'gauge_covariant_error': float(np.max(abs(shifted_covariant-phase[..., None]*covariant))),
            'omitted_phase_gradient_gap': float(np.sqrt(np.sum(np.sum(abs(wrong_gradient-covariant)**2,
                                                                      axis=-1)*data['weights'])))}


def build():
    return {'schema': 'oph.whitney_magnetic_continuum.v1', 'scope': SCOPE,
            'source_pins': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PINS},
            'analytic_scope': ANALYTIC, 'parameters': PARAMETERS,
            'mesh_checks': [ritz_record(n) for n in (1, 2, 4)],
            'controls': controls(), 'trajectories': [trajectory(steps) for steps in (32, 64)]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({'receipt': args.output.as_posix(), 'bytes': args.output.stat().st_size}))
