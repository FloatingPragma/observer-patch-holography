"""Opposite native quadratures and alignment standards on known physical ports."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from model.experiments.material_inertia.pair_spectrum.native_source_probe import circuit
from model.experiments.interval_echo_consensus.signed_effects import (
    fixed_pointer as pointer,
)
from .material import axes, vertices
from .probe import Quadrupole, draw

HERE = Path(__file__).resolve().parent


def phase_standards():
    basis = circuit.apply(pointer.neutral_basis(1, 2), np.eye(27))
    word = pointer.diagonal_word(0, 2, Fraction(8, 31))
    native = circuit.apply(word, np.eye(27))
    response = []
    regional_phases = []
    for q, pair in zip((-1, 0, 1), circuit.CODE):
        h = basis[np.ix_(pair, pair)]
        d = native[np.ix_(pair, pair)]
        # 02 forward and 20 inverse eliminate the common local contact phase.
        relative = d[1, 1] * d[0, 0].conjugate()
        regional_phases.append(relative)
        for sign, ticks in enumerate((1, 3)):
            y = circuit.apply(
                [pointer.pulse("reference", 2, 2, Fraction(ticks, 64))], np.eye(27)
            )
            exits = abs(h.conj().T @ y[np.ix_(pair, pair)] @ d @ h[:, 0]) ** 2
            if sign == 0:
                response.append([])
            response[-1].append(exits[0] - exits[1])
        assert abs(relative - np.exp(2j * np.pi * 130 * q / 31)) < 1e-12
    response = np.array(response)
    assert np.max(abs(response[:, 0] + response[:, 1])) < 3e-12
    return response, np.array(regional_phases)


def run():
    rng = np.random.default_rng(161403)
    response, _ = phase_standards()
    direction = np.sign(response[-1, 0])
    response *= direction
    coordinates = vertices()
    frame = axes(coordinates[0])
    proj = abs(coordinates @ frame[0])
    marks = np.zeros((3, 12), int)
    for i, tip in enumerate((np.argmax(proj), np.argmin(proj))):
        opposite = np.argmin(coordinates @ coordinates[tip])
        marks[i, [tip, opposite]] = 1
    marks[2] = marks[0]
    words = marks + 1
    words[2] = 1
    results = []
    for j in (0, 1):
        reader = Quadrupole(words, marks, frame[:1], 0, response[-1, j])
        bank = {
            label: reader.bank(row, rng)["centroid_0"]
            for label, row in zip(("axial", "transverse", "absent"), np.eye(3))
        }
        results.append(
            {
                "y_reference_ticks": (1, 3)[j],
                "calibration": [
                    draw(rng, [(1 + r) / 2, (1 - r) / 2]) for r in response[:, j]
                ],
                "banks": bank,
            }
        )
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = run()
    if args.output:
        args.output.write_text(json.dumps(rows, indent=1) + "\n")
    else:
        assert rows == json.loads((HERE / "standards.json").read_text())
        print("Replayed both native read pulses and all three alignment standards.")
