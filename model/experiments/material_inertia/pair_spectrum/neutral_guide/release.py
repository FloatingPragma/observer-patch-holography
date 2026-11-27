"""Open every native seam and reference without changing the mounted material state."""

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from ... import neutral_defect
from .. import measure as spectrum
from .guide import Flight


def release_report(ring, rng):
    board = federation(0)
    excess, vacancy, adjacent, _, h, cost = neutral_defect.generators(board)
    size = len(ring)
    v_ports, offsets = np.indices((size, size - 1))
    exc_ports = (v_ports + offsets + 1) % size
    ids = neutral_defect.ordered_index(np.asarray(ring)[exc_ports], np.asarray(ring)[v_ports], board.n)
    columns, definitions = [], []
    for s in (1, -1):
        for j in (0, 1):
            character = size // 2 + j
            flight = Flight(size, character, 8, s)
            phase = np.exp(2j * np.pi * character * np.arange(size) / size) / np.sqrt(size)
            full = np.zeros(len(excess), complex)
            full[ids] = phase[:, None] * flight.initial
            columns.append(full)
            definitions.append((s, j, phase, flight))
    initial = np.column_stack(columns)
    reduction = spectrum.PairReduction(h, cost)
    count = reduction.inertia(spectrum.CUT)['upper_lines_below']
    if count != 24:
        raise ArithmeticError('unrestricted cluster has a different rank')
    energies, basis, residual, _ = reduction.lowest(h, count)
    projector_error = residual / min(spectrum.CUT - energies[-1], energies[0] - 30)
    admissions = (abs(basis.T @ initial)**2).sum(axis=0)
    on_ring = np.isin(excess, ring) & np.isin(vacancy, ring)
    codes = 4 * adjacent.astype(int) + 2 * on_ring.astype(int) + (excess // 12 == vacancy // 12)
    out = [dict(orientation=s, offset=j, full_band_weight=float(admissions[i]), exposures=[])
           for i, (s, j, _, _) in enumerate(definitions)]
    for tick in (0, 0.5, 2):
        if tick == 0:
            states = initial
        else:
            states = expm_multiply(-1j * tick * h, initial, traceA=-1j * tick * h.diagonal().sum())
        for i, (s, j, phase, flight) in enumerate(definitions):
            state = states[:, i]
            weights = abs(state)**2
            np.testing.assert_allclose(sum(weights), 1, atol=2e-12, rtol=0)
            bins = np.bincount(codes, weights=weights, minlength=8)
            character_weights = (abs(np.fft.fft(state[ids], axis=0))**2).sum(axis=1) / size
            inside = sum(character_weights)
            probs = np.maximum(np.r_[character_weights, 1 - inside], 0)
            ret = np.vdot(initial[:, i], state)
            guide_ret = np.vdot(flight.initial, flight.at(tick))
            out[i]['exposures'].append(dict(time=tick, joint_counts=rng.multinomial(2**20, bins / sum(bins)).tolist(),
                ring_character_and_outside_counts=rng.multinomial(2**20, probs / sum(probs)).tolist(),
                ring_character_probability=float(character_weights[size // 2 + j]),
                complete_ring_exit_probability=float(1 - inside),
                graph_adjacency_probability=float(sum(bins[4:])), incoming_return_probability=float(abs(ret)**2),
                guide_return_probability=float(abs(guide_ret)**2),
                return_amplitude=[float(ret.real), float(ret.imag)],
                guide_amplitude=[float(guide_ret.real), float(guide_ret.imag)]))
    phase_checks = []
    for ri in (0, 2):
        for a, b in zip(out[ri]['exposures'], out[ri + 1]['exposures']):
            def phase_ratio(field):
                return complex(*b[field]) * complex(*a[field]).conjugate()
            phase_checks.append(dict(orientation=out[ri]['orientation'], time=a['time'],
                    full_minus_guide_relative_return_phase=float(np.angle(phase_ratio('return_amplitude')
                                                                          / phase_ratio('guide_amplitude')))))
    return dict(full_states=len(excess), band_count=count, band_projector_norm_error=float(projector_error),
                joint_bin_code='4*graph_adjacent + 2*both_on_ring + common_owner', sources=out,
                relative_phase=phase_checks)
