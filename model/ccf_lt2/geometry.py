#!/usr/bin/env python3
"""Per-vertex reciprocal volume/capacity/memory matrix primitive.

Tensor order: volume, active-capacity label, memory, source. The builder is
sparse and accepts the complete source matrix, not its fitted spectrum.
The standalone witness is a complete three-CAR single-vertex diagnostic.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse

TOL = 1e-10
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


@dataclass(frozen=True)
class GeometryConfig:
    volumes: tuple[float, ...] = (1.0, 1.5)
    capacity_max: int = 2
    lapse_coupling: float = .25
    volume_hopping: float = .5
    capacity_hopping: float = .125
    vacuum_strength: float = .03125
    volume_penalty: float = .0625
    capacity_penalty: float = .0625
    source_energy_reference: float = 0.0


@dataclass
class GeometryExtension:
    hamiltonian: sparse.csr_matrix
    source_hamiltonian: sparse.csr_matrix
    auxiliary_hamiltonian: sparse.csr_matrix
    lapse: sparse.csr_matrix
    volume: sparse.csr_matrix
    active_capacity: sparse.csr_matrix
    dimensions: tuple[int, int, int, int]
    source_gauss: tuple[sparse.csr_matrix, ...]
    auxiliary_outcomes: tuple[tuple[tuple[int, int, str], np.ndarray], ...]
    config: GeometryConfig
    source_piece_names: tuple[str, ...]

    @property
    def auxiliary_dimension(self):
        return int(np.prod(self.dimensions[:-1]))


def _norm(a):
    return float(sparse.linalg.norm(a)) if sparse.issparse(a) else float(np.linalg.norm(a))


def _matrix(a, name):
    a = sparse.csr_matrix(a, dtype=complex)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or a.shape[0] == 0:
        raise ValueError(f'{name} must be a nonempty square matrix')
    if not np.isfinite(a.data).all() or _norm(a-a.getH()) > TOL:
        raise ValueError(f'{name} must be finite and Hermitian')
    return a


def _path_laplacian(n):
    a = np.zeros((n, n))
    for i in range(n-1):
        a[i, i] += 1
        a[i+1, i+1] += 1
        a[i, i+1] = a[i+1, i] = -1
    return sparse.csr_matrix(a)


def build_geometry_extension(source_hamiltonian, gauss_operators=(), config=None,
                             source_pieces=None):
    """Couple every supplied source term to one positive finite lapse.

    `source_pieces` is an optional integrity dictionary. If supplied, its
    exact finite matrix sum must agree with source_hamiltonian numerically.
    No input state, source spectrum or target response is used by the builder.
    """
    cfg = config or GeometryConfig()
    h = _matrix(source_hamiltonian, 'source Hamiltonian')
    d = h.shape[0]
    v = np.asarray(cfg.volumes, dtype=float)
    numbers = [cfg.lapse_coupling, cfg.volume_hopping, cfg.capacity_hopping,
               cfg.vacuum_strength, cfg.volume_penalty, cfg.capacity_penalty]
    if (len(v) < 2 or not np.isfinite(v).all() or min(v) <= 0
            or not np.all(np.diff(v) > 0) or not all(np.isfinite(numbers))
            or not 0 <= cfg.lapse_coupling < 1 or min(numbers[1:]) < 0
            or not np.isfinite(cfg.source_energy_reference)
            or not isinstance(cfg.capacity_max, int) or cfg.capacity_max < 1):
        raise ValueError('positive ordered volumes, positive lapse and finite nonnegative couplings required')
    names = ()
    if source_pieces is not None:
        if not source_pieces:
            raise ValueError('the supplied source-piece dictionary is empty')
        pieces = [_matrix(p, f'source piece {name}') for name, p in source_pieces.items()]
        if any(p.shape != h.shape for p in pieces):
            raise ValueError('source pieces have different dimensions')
        if _norm(sum(pieces, sparse.csr_matrix(h.shape))-h) > TOL:
            raise ValueError('source pieces omit or alter part of the full source Hamiltonian')
        names = tuple(source_pieces)
    gauss = tuple(_matrix(g, 'Gauss operator') for g in gauss_operators)
    if any(g.shape != h.shape or _norm(g@h-h@g) > TOL for g in gauss):
        raise ValueError('the input source fails the declared Gauss identity')
    nv, nk, nm = len(v), cfg.capacity_max, 2**cfg.capacity_max
    k = np.arange(1, nk+1, dtype=float)
    iv, ik, im = sparse.eye(nv), sparse.eye(nk), sparse.eye(nm)
    x = 2*(v-v[0])/(v[-1]-v[0])-1
    lapse = sparse.kron(sparse.diags(1+cfg.lapse_coupling*x), sparse.kron(ik, im), format='csr')
    volume = sparse.kron(sparse.diags(v), sparse.kron(ik, im), format='csr')
    capacity = sparse.kron(iv, sparse.kron(sparse.diags(k), im), format='csr')
    potential = (cfg.vacuum_strength*v[:, None]/k[None, :]
                 + cfg.volume_penalty*(v[:, None]-v[0])**2
                 + cfg.capacity_penalty*(k[None, :]-1)**2)
    aux = (cfg.volume_hopping*sparse.kron(_path_laplacian(nv), sparse.kron(ik, im))
           + cfg.capacity_hopping*sparse.kron(iv, sparse.kron(_path_laplacian(nk), im))
           + sparse.kron(sparse.diags(potential.reshape(-1)), im)).tocsr()
    reference = cfg.source_energy_reference
    stress_energy = h-reference*sparse.eye(d)
    joint = (sparse.kron(aux, sparse.eye(d))+sparse.kron(lapse, stress_energy)
             +reference*sparse.eye(aux.shape[0]*d)).tocsr()
    outcomes = []
    for vi in range(nv):
        for ki in range(nk):
            bits = ki+1
            for word in range(2**bits):
                mask = np.zeros(nv*nk*nm)
                for memory in range(nm):
                    if memory % (2**bits) == word:
                        mask[(vi*nk+ki)*nm+memory] = 1
                outcomes.append(((vi, bits, format(word, f'0{bits}b')), mask))
    if not np.array_equal(sum(mask for _, mask in outcomes), np.ones(nv*nk*nm)):
        raise AssertionError('capacity instrument is incomplete')
    return GeometryExtension(joint, h, aux, lapse, volume, capacity,
                             (nv, nk, nm, d), gauss, tuple(outcomes), cfg, names)


def evolve(extension, state, time):
    """Actual full finite propagator; no projection or renormalization."""
    psi = np.asarray(state, dtype=complex).reshape(-1)
    if len(psi) != extension.hamiltonian.shape[0] or abs(np.vdot(psi, psi)-1) > TOL:
        raise ValueError('one normalized state on the entire joint carrier required')
    return sparse.linalg.expm_multiply((-1j*float(time))*extension.hamiltonian, psi)


def complete_source_records(extension, state, source_projectors):
    """Joint complete source/volume/capacity/active-memory instrument.

    Includes zero-probability outcomes and measurement energy work. The
    inactive memory bits remain inside each projector; they are not erased.
    """
    d = extension.dimensions[-1]
    projectors = [(str(label), _matrix(p, f'projector {label}')) for label, p in source_projectors]
    if not projectors or any(p.shape != (d, d) or _norm(p@p-p) > TOL for _, p in projectors):
        raise ValueError('complete orthogonal source projectors required')
    if _norm(sum((p for _, p in projectors), sparse.csr_matrix((d, d)))-sparse.eye(d)) > TOL:
        raise ValueError('source instrument omits outcomes')
    for i, (_, p) in enumerate(projectors):
        if any(_norm(p@q) > TOL for _, q in projectors[i+1:]):
            raise ValueError('source projectors overlap')
    psi = np.asarray(state, dtype=complex).reshape(-1)
    if len(psi) != extension.hamiltonian.shape[0] or abs(np.vdot(psi, psi)-1) > TOL:
        raise ValueError('normalized full joint state required')
    records, after_energy = [], 0.0
    source_marginal = {name: 0.0 for name, _ in projectors}
    volume_marginal = [0.0]*extension.dimensions[0]
    capacity_marginal = [0.0]*extension.dimensions[1]
    for (vi, bits, word), mask in extension.auxiliary_outcomes:
        for label, p in projectors:
            projector = sparse.kron(sparse.diags(mask), p, format='csr')
            branch = projector@psi
            probability = float(np.vdot(branch, branch).real)
            records.append({'volume_index': vi, 'active_bits': bits,
                            'memory_word': word, 'source_outcome': label,
                            'probability': probability})
            source_marginal[label] += probability
            volume_marginal[vi] += probability
            capacity_marginal[bits-1] += probability
            after_energy += float(np.vdot(branch, extension.hamiltonian@branch).real)
    before = float(np.vdot(psi, extension.hamiltonian@psi).real)
    ia = sparse.eye(extension.auxiliary_dimension)
    bare = float(np.vdot(psi, sparse.kron(ia, extension.source_hamiltonian)@psi).real)
    stress_energy = extension.source_hamiltonian-extension.config.source_energy_reference*sparse.eye(d)
    weighted = float(np.vdot(psi, sparse.kron(extension.lapse, stress_energy)@psi).real)
    auxiliary = float(np.vdot(psi, sparse.kron(extension.auxiliary_hamiltonian, sparse.eye(d))@psi).real)
    return {'records': records, 'source_marginal': source_marginal,
            'volume_marginal': volume_marginal, 'capacity_marginal': capacity_marginal,
            'probability_sum': sum(r['probability'] for r in records),
            'bare_source_energy': bare, 'lapse_weighted_source_relative_energy': weighted,
            'source_energy_reference': extension.config.source_energy_reference,
            'auxiliary_energy': auxiliary, 'total_energy': before,
            'terminal_measurement_work': after_energy-before}


def _projector(v):
    return np.outer(v, v.conj())


def diagnostic_source(neutral_mass=2.0):
    """Complete e,p,n CAR occupation space; not the master CCF-2 source."""
    n_e = np.diag([(i >> 0) & 1 for i in range(8)])
    n_p = np.diag([(i >> 1) & 1 for i in range(8)])
    n_n = np.diag([(i >> 2) & 1 for i in range(8)])
    pieces = {'charged_rest': n_e+n_p, 'neutral_rest': neutral_mass*n_n}
    h = sum(pieces.values())
    charge = -n_e+n_p
    projectors = []
    for n in range(2):
        for sign, label in ((1, 'plus'), (-1, 'minus')):
            state = np.zeros(8, complex)
            state[4*n] = 1/np.sqrt(2)
            state[4*n+3] = sign/np.sqrt(2)
            projectors.append((f'n{n}_{label}', _projector(state)))
        for offset, label in ((1, 'electron'), (2, 'positron')):
            state = np.eye(8)[4*n+offset]
            projectors.append((f'n{n}_{label}', _projector(state)))
    return h, pieces, charge, projectors


def diagnostic_initial(extension, neutral):
    source = np.zeros(8, complex)
    source[4*neutral] = source[4*neutral+3] = 1/np.sqrt(2)
    auxiliary = np.zeros(extension.auxiliary_dimension)
    auxiliary[0] = 1
    return np.kron(auxiliary, source)


def witness():
    """Frozen finite construction checks; all branches and controls retained."""
    rows, diagnostics = [], {}
    for control, config, neutral_mass in (
            ('coupled', GeometryConfig(), 2.0),
            ('uncoupled', replace(GeometryConfig(), lapse_coupling=0), 2.0),
            ('neutral_mass_changed', GeometryConfig(), 3.0)):
        h, pieces, charge, projectors = diagnostic_source(neutral_mass)
        extension = build_geometry_extension(h, [charge], config, pieces)
        lifted_charge = sparse.kron(sparse.eye(extension.auxiliary_dimension), charge)
        diagnostics[control] = {
            'dimension': extension.hamiltonian.shape[0],
            'hermiticity_residual': _norm(extension.hamiltonian-extension.hamiltonian.getH()),
            'gauss_commutator_residual': _norm(lifted_charge@extension.hamiltonian-extension.hamiltonian@lifted_charge),
            'min_lapse': float(extension.lapse.diagonal().real.min()),
            'configuration': {**asdict(config), 'volumes': list(config.volumes)},
            'neutral_mass': neutral_mass}
        for neutral in (0, 1):
            initial = diagnostic_initial(extension, neutral)
            initial_energy = float(np.vdot(initial, extension.hamiltonian@initial).real)
            for time in (0.0, .5, 1.0, 2.0):
                evolved = evolve(extension, initial, time)
                result = complete_source_records(extension, evolved, projectors)
                result.update(control=control, neutral_preparation=neutral, time=time,
                              energy_residual=result['total_energy']-initial_energy)
                result['charged_plus_probability'] = sum(
                    value for name, value in result['source_marginal'].items() if name.endswith('_plus'))
                rows.append(result)
    lookup = {(r['control'], r['neutral_preparation'], r['time']): r for r in rows}
    a, b = [lookup['coupled', n, 2.0] for n in (0, 1)]
    response = abs(a['charged_plus_probability']-b['charged_plus_probability'])
    geometry_response = max(abs(x-y) for x, y in zip(a['volume_marginal'], b['volume_marginal']))
    null = max(abs(lookup['uncoupled', 0, t]['charged_plus_probability']
                   - lookup['uncoupled', 1, t]['charged_plus_probability']) for t in (0., .5, 1., 2.))
    clock = max(abs(lookup['uncoupled', n, t]['charged_plus_probability']-(1+np.cos(2*t))/2)
                for n in (0, 1) for t in (0., .5, 1., 2.))
    changed = abs(lookup['neutral_mass_changed', 1, 2.]['charged_plus_probability']-b['charged_plus_probability'])
    h, pieces, charge, projectors = diagnostic_source()
    try:
        build_geometry_extension(h, [charge], source_pieces={'charged_rest': pieces['charged_rest']})
        omitted_rejected = False
    except ValueError:
        omitted_rejected = True
    # An isolated energy eigenstate has no observable global phase.
    state = np.eye(8, dtype=complex)[1]
    stationary_residual = 0.0
    for offset in (0., 3.):
        energies, frame = np.linalg.eigh(h+offset*np.eye(8))
        evolved = frame@(np.exp(-2j*energies)*(frame.conj().T@state))
        stationary_residual = max(stationary_residual, float(np.linalg.norm(_projector(evolved)-_projector(state))))
    cfg = GeometryConfig()
    eps, v0, v1, k0, k1 = cfg.vacuum_strength, *cfg.volumes, 1, 2
    cycle = [eps*(v1-v0)/k0, eps*v1*(1/k1-1/k0),
             eps*(v0-v1)/k1, eps*v0*(1/k0-1/k1)]
    slope = (eps*v1/k0-eps*v0/k0)/(v1-v0)
    # A conventional common shift changes only the full scalar. Holding the
    # source reference fixed instead changes the actual geometry potential.
    base = build_geometry_extension(h, [charge], cfg, pieces)
    offset = 3.0
    shifted_pieces = dict(pieces, scalar_offset=offset*np.eye(8))
    shifted = build_geometry_extension(h+offset*np.eye(8), [charge],
        replace(cfg, source_energy_reference=offset), shifted_pieces)
    fixed_ref = build_geometry_extension(h+offset*np.eye(8), [charge], cfg, shifted_pieces)
    reference_shift_residual = _norm(shifted.hamiltonian-base.hamiltonian-offset*sparse.eye(base.hamiltonian.shape[0]))
    fixed_ref_change = _norm(fixed_ref.hamiltonian-base.hamiltonian-offset*sparse.eye(base.hamiltonian.shape[0]))
    initial = diagnostic_initial(base, 0)
    original_record = complete_source_records(base, evolve(base, initial, 2.), projectors)
    shifted_record = complete_source_records(shifted, evolve(shifted, initial, 2.), projectors)
    fixed_record = complete_source_records(fixed_ref, evolve(fixed_ref, initial, 2.), projectors)
    shift_record_residual = max(abs(a['probability']-b['probability']) for a,b in zip(original_record['records'],shifted_record['records']))
    fixed_record_change = max(abs(a['probability']-b['probability']) for a,b in zip(original_record['records'],fixed_record['records']))
    gates = {
        'hermitian': max(x['hermiticity_residual'] for x in diagnostics.values()) <= TOL,
        'gauss_preserved': max(x['gauss_commutator_residual'] for x in diagnostics.values()) <= TOL,
        'probabilities_complete': max(abs(r['probability_sum']-1) for r in rows) <= TOL,
        'total_energy_conserved': max(abs(r['energy_residual']) for r in rows) <= TOL,
        'neutral_changes_charged_clock': response > 1e-6,
        'neutral_changes_geometry': geometry_response > 1e-6,
        'coupling_deleted_clock_response_zero': null <= TOL,
        'actual_pair_gap_clock': clock <= TOL,
        'neutral_mass_change_detected': changed > 1e-6,
        'source_piece_omission_rejected': omitted_rejected,
        'common_source_reference_shift_global_only': reference_shift_residual <= TOL and shift_record_residual <= TOL,
        'fixed_reference_vacuum_shift_is_physical': fixed_ref_change > 1e-6 and fixed_record_change > 1e-6,
        'isolated_global_phase_invisible': stationary_residual <= TOL,
        'volume_and_capacity_dynamical': a['volume_marginal'][1] > 1e-6 and a['capacity_marginal'][1] > 1e-6,
        'capacity_volume_cycle_closed': abs(sum(cycle)) <= TOL,
        'declared_fixed_capacity_pressure_identity': abs(slope-eps/k0) <= TOL,
        'all_outcomes_retained': all(len(r['records']) == 96 for r in rows)}
    gates = {name: bool(value) for name,value in gates.items()}
    return {'component': 'CCF-LT-2 finite reciprocal geometry',
            'scope': 'Complete three-CAR diagnostic source; generic builder accepts actual master source. Not Einstein/Newton/Lambda or calibrated physical time.',
            'tensor_order': ['volume', 'active_capacity', 'memory', 'source'],
            'diagnostics': diagnostics, 'gates': gates, 'all_passed': all(gates.values()),
            'neutral_clock_probability_difference': response,
            'neutral_geometry_probability_difference': geometry_response,
            'uncoupled_clock_difference': null, 'pair_gap_clock_residual': clock,
            'changed_neutral_mass_clock_difference': changed,
            'stationary_global_phase_density_residual': stationary_residual,
            'common_source_reference_shift_operator_residual': reference_shift_residual,
            'common_source_reference_shift_record_residual': shift_record_residual,
            'fixed_reference_vacuum_shift_record_change': fixed_record_change,
            'reference_shift_complete_records': {'base': original_record, 'common_shift': shifted_record, 'fixed_reference_shift': fixed_record},
            'vacuum_four_stroke_work': cycle, 'runs': rows,
            'remaining_obligations': ['spatial local metric and stress', 'inverse-square response and unfitted G',
                'RER joint-refinement observational maps', 'capacity selection',
                'horizon and vacuum-stress dictionary', 'unfitted Lambda and cosmic acceleration',
                'clock and energy calibration', 'actual master source replay']}

