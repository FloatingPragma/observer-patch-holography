"""Independent host-reference formation clock and native directed bridge-current read."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ... import neutral_defect
from .. import measure as spectrum
from ..sum_rule_inertia import measure as meter
from .free_reader import free_sources
from .guide import physical_rings


def with_refusal(weights, admission, rng):
    p = np.maximum(weights, 0)
    np.testing.assert_allclose(sum(p), 1, atol=2e-12, rtol=0)
    return meter.sample(np.r_[admission * p / sum(p), 1 - admission], rng)['counts']


def current_probabilities(left, right, outside=None):
    plus = float(np.linalg.norm((1j * left + right) / np.sqrt(2))**2)
    minus = float(np.linalg.norm((left + 1j * right) / np.sqrt(2))**2)
    return [plus, minus, max(0, 1 - plus - minus) if outside is None else outside]


def bridge_indices(a, b, n):
    other = np.setdiff1d(np.arange(n), [a, b])
    def at(p):
        return np.r_[neutral_defect.ordered_index(p, other, n),
                     neutral_defect.ordered_index(other, p, n)]
    return at(a), at(b)


def report():
    board = federation(0)
    _, rings, _ = physical_rings()
    _, _, _, _, h, cost = neutral_defect.generators(board)
    reduced = spectrum.PairReduction(h, cost)
    e, basis, residual, _ = reduced.lowest(h, 24)
    cut = (e[1] + e[2]) / 2
    if reduced.inertia(cut)['upper_lines_below'] != 2:
        raise ArithmeticError('incorrect carrier rest admission')
    rest_error = residual / (min(cut - e[1], spectrum.CUT - e[-1]) - residual) + 1e-12
    basis = basis[:, :2]
    full_energies = 96 + e[:2]
    # This independent row-sum bound includes every off-diagonal and the formation offset.
    abs_bound = float(96 + np.max(abs(h).sum(axis=1)))
    if abs_bound / 128 >= np.pi:
        raise ArithmeticError('unresolved shortest-wait formation alias')
    indices = [bridge_indices(a, b, board.n) for a, b in board.external]
    operator_bounds = []
    for left, right in indices:
        cross = basis[left].T @ basis[right]
        compressed = 1j * (cross.T - cross) / 4
        operator_bounds.append(float(np.linalg.norm(compressed, 2) + rest_error / 2))
    rng = np.random.default_rng(190029)
    native_clock = meter.clock_standard(rng)
    controls = {}
    for label, lft, rgt in [('right', 1 / np.sqrt(2), 1j / np.sqrt(2)),
                            ('left', 1 / np.sqrt(2), -1j / np.sqrt(2)), ('host', 0., 0.)]:
        controls[label] = with_refusal(current_probabilities(lft, rgt, 1. if label == 'host' else 0.), 1., rng)
    panels = []
    for sign, off, source in free_sources(rings[20]['ports'], board.n):
        g = basis.T @ source
        admission = float(np.vdot(g, g).real)
        g /= np.sqrt(admission)
        energies = []
        for tau in (1/128, 1/64, 1/32):
            amplitude = sum(abs(g)**2 * np.exp(-1j * full_energies * tau))
            counts = [with_refusal([(1 - v) / 2, (1 + v) / 2], admission, rng)
                      for v in (amplitude.real, amplitude.imag)]
            seen = complex(*[(p - m) / (p + m) for m, p, _ in counts])
            angle = -np.angle(seen)
            if energies:
                angle += 2 * np.pi * round((energies[0]['measured_cost'] * tau - angle) / (2 * np.pi))
            energies.append(dict(tau=tau, quadratures=counts, measured_cost=float(angle / tau)))
        currents = []
        for tau in (0, 10**6):
            arrived = basis @ (g * np.exp(-1j * (e[:2] - e[0]) * tau))
            counts = [with_refusal(current_probabilities(arrived[left], arrived[right]), admission, rng)
                      for left, right in indices]
            largest = max(abs(p - m) / (4 * (p + m + zero)) for p, m, zero, _ in counts)
            currents.append(dict(tau=tau, bridge_counts=counts, largest_abs_measured_current=largest))
        panels.append(dict(orientation=sign, offset=off, admission=admission, energies=energies, currents=currents))
    return dict(seed=190029, shots=meter.SHOTS, receiver_bridges=board.external, independent_clock=native_clock,
                bridge_calibration_counts=controls, current_exits=['toward second', 'toward first', 'neutral', 'refusal'],
                energy_exits=['minus', 'plus', 'refusal'], absolute_energy_upper_bound=abs_bound,
                rest_projector_norm_error=rest_error, all_time_doublet_bridge_current_bounds=operator_bounds,
                sources=panels)


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--output', type=Path, required=True)
    args = cli.parse_args()
    args.output.write_text(json.dumps(report(), indent=1, allow_nan=False) + '\n')
