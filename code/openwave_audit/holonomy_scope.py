"""Countercheck the regulator and observable domain of holonomy reduction."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from native_scope import triangle


def wall_triangle(cutoff=1, mass=.23):
    states = []
    for z in range(6):
        for r in range(-cutoff, cutoff + 1):
            f0 = r
            f1 = int(z == 3) - f0
            f2 = -int(z == 1) - f1
            f3 = int(z == 4) - f2
            f4 = -int(z == 2) - f3
            f5 = 1 - int(z == 0) - f0
            records = (f0, f1, f2, f3, f4, f5)
            if all(abs(f) <= cutoff for f in records):
                states.append((z, records))
    index = {s: i for i, s in enumerate(states)}
    h = np.diag([mass if z < 3 else -mass for z, f in states])
    for col, (z, records) in enumerate(states):
        for e, endpoints in enumerate(((0, 1), (1, 2), (2, 0))):
            for side, port in enumerate(endpoints):
                if z != port:
                    continue
                shifted = list(records)
                shifted[2 * e + side] += 1
                target = (3 + e, tuple(shifted))
                if target in index:
                    row = index[target]
                    h[row, col] = h[col, row] = 1 if side == 0 else -1
    fibre_sizes = np.bincount([z for z, f in states], minlength=6)
    v = np.zeros((len(states), 6))
    w = np.zeros((len(states), len(states)), dtype=int)
    for j, (z, f) in enumerate(states):
        v[j, z] = 1 / np.sqrt(fibre_sizes[z])
        shifted = tuple(x + (1 if k % 2 == 0 else -1) for k, x in enumerate(f))
        if (z, shifted) in index:
            w[index[z, shifted], j] = 1
    return states, h, v, w, fibre_sizes


def run():
    states, index, h, v, incidence = triangle()
    bare = v.T @ h @ v
    w = np.zeros_like(h, dtype=int)
    for j, (z, f) in enumerate(states):
        shifted = tuple((x + (1 if k % 2 == 0 else -1)) % 6 for k, x in enumerate(f))
        w[index[z, shifted], j] = 1
    observable = (w + w.T) / 2
    plus = v[:, 0]
    minus = np.array([((-1) ** f[0]) / np.sqrt(6) if z == 0 else 0 for z, f in states])
    plus_law = np.bincount([z for z, f in states], weights=abs(plus) ** 2, minlength=6)
    minus_law = np.bincount([z for z, f in states], weights=abs(minus) ** 2, minlength=6)
    observable_values = [float(s @ observable @ s) for s in (plus, minus)]
    wall_states, wall_h, wall_v, wall_w, sizes = wall_triangle()
    # A nonzero entry is an exact integer hop obstruction: diagonal mass
    # commutes with a cycle shift, so it cancels from this commutator.
    wall_commutator = wall_h @ wall_w - wall_w @ wall_h
    errors = {'modular_cycle_unitarity': float(np.max(abs(w.T @ w - np.eye(len(states))))),
              'modular_cycle_commutator': float(np.max(abs(h @ w - w @ h))),
              'observable_holonomy_commutator': float(np.max(abs(observable @ w - w @ observable))),
              'same_position_law': float(np.max(abs(plus_law - minus_law))),
              'wall_cycle_unitarity_failure': float(np.max(abs(wall_w.T @ wall_w - np.eye(len(wall_states))))),
              'wall_cycle_commutator': float(np.max(abs(wall_commutator))),
              'wall_uniform_intertwiner_failure': float(np.max(abs(wall_h @ wall_v - wall_v @ bare)))}
    def probabilities(h, vector, positions):
        energies, vectors = np.linalg.eigh(h)
        out = vectors @ (np.exp(-.37j * energies) * (vectors.T @ vector))
        return np.bincount(positions, weights=abs(out) ** 2, minlength=6)
    wall_prob = probabilities(wall_h, wall_v[:, 0], [z for z, f in wall_states])
    bare_prob = probabilities(bare, np.eye(6)[:, 0], list(range(6)))
    checks = {'finite_modular_reduction_premises_hold': all(errors[k] < 1e-12 for k in ('modular_cycle_unitarity', 'modular_cycle_commutator')),
              'commuting_observable_not_a_position_function': errors['observable_holonomy_commutator'] == 0 and errors['same_position_law'] == 0 and abs(observable_values[0] - observable_values[1]) > 1.99,
              'hard_wall_cycle_group_transfer_rejected': errors['wall_cycle_unitarity_failure'] == 1 and errors['wall_cycle_commutator'] == 1,
              'hard_wall_uniform_intertwining_rejected': errors['wall_uniform_intertwiner_failure'] > .1,
              'hard_wall_position_dynamics_differ': float(np.max(abs(wall_prob - bare_prob))) > 1e-3,
              'modulus_one_charge_lift_not_injective': 0 % 1 == 1 % 1}
    return {'schema': 'oph.audit.holonomy_scope.v1', 'checks': checks, 'errors': errors,
            'modular_observable': {'definition': '(W+W^T)/2', 'expectations_same_position_different_character': observable_values,
                                   'scope': 'Counterexample to the general commutant assertion. It is not among the narrower listed divF readouts.'},
            'hard_wall': {'cutoff': 1, 'states': len(wall_states), 'fibre_sizes': sizes.tolist(),
                          'bare_position_probabilities': bare_prob.tolist(), 'wall_position_probabilities': wall_prob.tolist(),
                          'maximum_position_difference': float(np.max(abs(wall_prob - bare_prob)))},
            'scope': 'Complete finite modular fibres support the character reduction for the preserved apparatus. Hard-wall truncation lacks that shift-group hypothesis. A general holonomy-commuting operator can depend on character, while the listed Gauss/divergence readouts reduce to position within a fixed sector. Modulus-one coarse evolution does not implement an injective modular charge decoder. No claim is made about new physical particles or platform status.',
            'reviewer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'errors', 'hard_wall')}, indent=2))
    raise SystemExit(0 if all(result['checks'].values()) else 1)
