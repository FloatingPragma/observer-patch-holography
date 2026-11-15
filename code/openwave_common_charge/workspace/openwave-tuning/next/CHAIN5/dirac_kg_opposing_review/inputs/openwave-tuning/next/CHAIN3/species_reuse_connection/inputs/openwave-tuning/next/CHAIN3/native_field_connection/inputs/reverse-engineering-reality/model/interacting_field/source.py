"""Completed native words, exact CAR routing, and complete pointer instruments.

The local primitive matrices come from source pulse contraction. Global
routing uses the proved FSWAP conjugation identity; its source residual and
all requested routing pulses are counted, not replayed as a giant qubit array.
"""
from __future__ import annotations
from functools import lru_cache
import math
import copy
import numpy as np
from scipy.linalg import expm
from ..observer.source_scalar import compiler
from ..common_source.chiral_candidate import gaussian_source as gaussian
from .lattice import real

COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'proper_code_load_events', 'proper_code_export_events', 'proper_code_service_events')


@lru_cache(maxsize=64)
def _contact(duration, interaction):
    word = compiler.compile_diagonal(np.array([0., 0., 0., -duration*interaction]))
    raw = compiler.execute(word, 2)
    target = np.diag(np.exp(-1j*duration*interaction*np.array([0., 0., 0., 1.])))
    residual = float(np.linalg.norm(raw-target, 2))
    if residual > 2e-10:
        raise ArithmeticError('contact word failed its complete Fock-space target')
    raw = np.frombuffer(raw.tobytes(), dtype=raw.dtype).reshape(raw.shape)
    return raw, {'resources': compiler.resources(word), 'source_target_residual': residual,
                 'contracted_Fock_columns': 4, 'word': word}


def contact(duration, interaction=1.):
    raw, info = _contact(real(duration, 'model duration'), real(interaction, 'interaction'))
    return raw, copy.deepcopy(info)


@lru_cache(maxsize=1)
def _occupation():
    word = compiler.cnot(0, 1)
    raw = compiler.execute(word, 2)
    kraus = (raw[:2, :2], raw[2:, :2])
    residual = max(float(np.linalg.norm(k-np.diag([1-b, b]), 2)) for b, k in enumerate(kraus))
    if residual > 2e-10:
        raise ArithmeticError('complete number pointer failed')
    kraus = tuple(np.frombuffer(k.tobytes(), dtype=k.dtype).reshape(k.shape) for k in kraus)
    return kraus, {'resources': compiler.resources(word), 'source_target_residual': residual,
                   'completeness_residual': float(np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(2), 2))}


def occupation():
    kraus, info = _occupation()
    return kraus, copy.deepcopy(info)


def read_number(space, state, mode):
    if type(mode) is not int or not 0 <= mode < space.modes:
        raise ValueError('mode outside declared space')
    state = space.state(state)
    occupied = (space.first == mode) | (space.second == mode)
    kraus, _ = occupation()
    # Off-diagonal source roundoff is outside exact sector compression; keep
    # its full local residual rather than treating the output as exact pulses.
    factors = [np.where(occupied, k[1, 1], k[0, 0]) for k in kraus]
    return tuple(state*(f if state.ndim == 1 else f[:, None]) for f in factors)


def hop(space, state, p, q, coefficient, duration):
    gate = gaussian.local_gate('hopping', coefficient, duration)
    return space.apply_pair(state, p, q, gate.fock_unitary)


def read_ramsey(space, state, p, q):
    reverse = gaussian.local_gate('hopping', -1., -math.pi/4)
    forward = gaussian.local_gate('hopping', -1., math.pi/4)
    before = space.apply_pair(state, p, q, reverse.fock_unitary)
    return tuple(space.apply_pair(branch, p, q, forward.fock_unitary)
                 for branch in read_number(space, before, p))


def phase_mode(space, state, mode, angle):
    # Pick a local spectator; the complete source onsite word is identity on it.
    if mode+1 < space.modes:
        p, q = mode, mode+1
        gate = gaussian.local_gate('onsite', -float(angle), 1.).fock_unitary
    else:
        p, q = mode-1, mode
        raw = gaussian.local_gate('onsite', -float(angle), 1.).fock_unitary
        swap = gaussian.local_gate('fswap').fock_unitary
        gate = swap@raw@swap
    return space.apply_pair(state, p, q, gate)


def evolve(lattice, state, duration, steps):
    duration = real(duration, 'model time')
    if type(steps) is not int or steps < 1:
        raise ValueError('positive split count required')
    out = lattice.state(state).copy()
    dt = duration/(2*steps)
    gates = {}
    for group in lattice.groups:
        for _, _, z in group:
            if z not in gates:
                raw = gaussian.local_gate('hopping', z, dt).fock_unitary
                # Factor the full scalar algebraically; never normalize a state.
                vacuum = raw[0, 0]
                gate = raw/vacuum; gate[0, 0] = 1
                gates[z] = (gate, vacuum)
    onsite, _ = contact(dt, lattice.interaction)
    diag = np.diag(onsite)
    # Execute every onsite block simultaneously by its diagonal identity.
    up = (lattice.first % 2 == 0).astype(int)+(lattice.second % 2 == 0)
    down = 2-up
    factors = diag[0]**lattice.volume*(diag[1]/diag[0])**up*(diag[2]/diag[0])**down
    factors[lattice.contact_indices] = diag[0]**(lattice.volume-1)*diag[3]
    scalar = 1.+0j
    for _ in range(steps):
        for group in lattice.groups:
            for p, q, z in group:
                gate, vacuum = gates[z]
                lattice.apply_pair_inplace(out, p, q, gate)
                scalar *= vacuum
        out *= factors if out.ndim == 1 else factors[:, None]
        out *= factors if out.ndim == 1 else factors[:, None]
        for group in reversed(lattice.groups):
            for p, q, z in reversed(group):
                gate, vacuum = gates[z]
                lattice.apply_pair_inplace(out, p, q, gate)
                scalar *= vacuum
    return scalar*out


def local_error(gate):
    """Floating diagnostic against an independently assembled Fock generator."""
    h=np.zeros((4,4),complex)
    if gate.kind=='hopping':
        h[1,2]=gate.coefficient;h[2,1]=gate.coefficient.conjugate()
        target=expm(-1j*gate.duration*h)
    elif gate.kind=='onsite':
        target=np.diag(np.exp(-1j*gate.duration*gate.coefficient.real*np.array([0,1,0,1])))
    elif gate.kind=='fswap':
        target=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]],complex)
    else:
        raise ValueError('unsupported diagnostic block')
    return float(np.linalg.norm(gate.fock_unitary-target,2))


def resources(lattice, duration, steps, segments=2):
    """Count a positive integer number of complete evolution segments.

    Resource admission uses the same duration and split-count domain as
    execution. Negative or fractional counts cannot represent a source word.
    """
    duration = real(duration, 'model time')
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer split count required')
    if type(segments) is not int or segments < 1:
        raise ValueError('positive integer segment count required')
    dt = duration/(2*steps)
    counts = {key: 0 for key in COUNT_KEYS}
    swap = gaussian.local_gate('fswap')
    routing = 0
    residual = 0.
    residual_sum = 0.
    swap_error = float(np.linalg.norm(swap.fock_unitary-np.array(
        [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]], complex), 2))
    for group in lattice.groups:
        for p, q, z in group:
            g = gaussian.local_gate('hopping', z, dt)
            routes = 2*(q-p-1)
            routing += routes*2*steps*segments
            for key in COUNT_KEYS:
                counts[key] += (g.resources[key]+routes*swap.resources[key])*2*steps*segments
            # Local exact target comparison independent of covariance extraction.
            h = np.zeros((4, 4), complex); h[1, 2] = z; h[2, 1] = z.conjugate()
            local_error = float(np.linalg.norm(g.fock_unitary-expm(-1j*dt*h), 2))
            residual = max(residual, local_error)
            residual_sum += (local_error+routes*swap_error)*2*steps*segments
    _, ci = contact(dt, lattice.interaction)
    for key in COUNT_KEYS:
        counts[key] += ci['resources'][key]*2*steps*segments*lattice.volume
    residual_sum += ci['source_target_residual']*2*steps*segments*lattice.volume
    return {**counts, 'routing_FSWAP_calls': routing,
            'floating_local_residual_sum_with_routing': residual_sum,
            'floating_residual_telescoping_estimate': math.expm1(residual_sum),
            'floating_residual_estimate_outward_certified': False,
            'number_sector_compression_residual_included_in_local_target_residuals': True,
            'completed_hopping_blocks': sum(map(len, lattice.groups))*2*steps*segments,
            'completed_contact_blocks': lattice.volume*2*steps*segments,
            'maximum_local_hopping_target_residual': residual,
            'contact_target_residual': ci['source_target_residual'],
            'FSWAP_target_residual': float(np.linalg.norm(swap.fock_unitary-np.array(
                [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]], complex), 2)),
            'routing_evaluation': 'exact CAR conjugation identity; local FSWAP source word contracted on all four inputs',
            'all_routing_pulses_numerically_replayed': False,
            'scope': 'worst-branch evolution only; separately retained apparatus counts must be added',
            'physical_hardware_error_included': False}
