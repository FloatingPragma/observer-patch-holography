"""Held clicks of a graded P1 module carried between two Interval Echo owners."""

import argparse
from functools import cache
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ..scalar_grammar import measure as pulse
from ..scalar_grammar.shell_line.apparatus import transport
from .bracket import ballot, sample

HERE = Path(__file__).resolve().parent
SOURCES = [
    ("upper",),
    ("lower",),
    ("ref",),
    *(
        ("pair", i, j, angle)
        for i, j in combinations(range(3), 2)
        for angle in (0.0, np.pi / 2)
    ),
]
READS = pulse.READS
SHOTS = 8192
SEED = 120421
IDENTITY = np.eye(3, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)


def vector(matrix):
    return np.array([matrix[2, 1], matrix[0, 2], matrix[1, 0]])


@cache
def rotations():
    tensor, q, matrices, _ = sample()
    first = vector(matrices[0])
    second = next(
        i
        for i, a in enumerate(matrices)
        if np.linalg.norm(np.cross(first, vector(a))) > 1e-8
    )
    return {
        "zero": [],
        "center": [("sum", 1.0)],
        "turn": [(0, 2 * np.pi)],
        "ordered": [(0, 0.73), (second, 1.09)],
        "reversed": [(second, 1.09), (0, 0.73)],
    }


def motion(appointments):
    tensor, q, matrices, _ = sample()
    word = []
    expected = IDENTITY.copy()
    for i, angle in appointments:
        if i == "sum":
            adj = matrices.sum(axis=0)
            duration = angle
            physical = tensor.sum(axis=0)
        else:
            adj = matrices[i]
            duration = angle / np.linalg.norm(vector(adj))
            physical = tensor[i]
        w = vector(adj)
        # Pauli words in the triangle are calibrated with the opposite i sign.
        word += pulse.field_word([*(-w / pulse.G), 0], duration)
        # Independent representation forecast. The adjoint-to-spinor
        # sign is fixed by continuous lifting from exposure zero, including 2π.
        half_axis = w @ PAULI.reshape(3, 4)
        step = expm(-0.5j * duration * half_axis.reshape(2, 2))
        full_adjoint = q.T @ expm(duration * physical) @ q
        for j in range(3):
            acted = step @ PAULI[j] @ step.conj().T
            np.testing.assert_allclose(
                acted, np.einsum("i,ijk->jk", full_adjoint[:, j], PAULI), atol=1e-12
            )
        expected[:2, :2] = step @ expected[:2, :2]
    return word, expected


def frame_program():
    # Predicted chart is independently given as ordinary three-coordinate rotations.
    word = pulse.ry(0, 2, 0.67) + pulse.rz(0, 1, 0.41) + pulse.ry(1, 2, -0.59)

    def real_turn(i, j, theta):
        t = theta / 2
        r = IDENTITY.copy()
        r[np.ix_([i, j], [i, j])] = [[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]]
        return r

    return word, real_turn(1, 2, -0.59) @ np.diag(
        np.exp(1j * np.array([0.41, -0.41, 0]) / 2)
    ) @ real_turn(0, 2, 0.67)


def with_refusal(prob):
    return np.r_[prob, 0.0]


def run():
    source = np.eye(3)[:, 0]
    site = transport()
    deliver = IDENTITY[:, site["connected_columns"]]
    turns, frame = frame_program()
    rng = np.random.default_rng(SEED)
    all_banks = []
    totals = {}
    # Include native, changed-law and cut acquisitions in the familywise count.
    bins = 4 * len(SOURCES) * len(READS) * (len(rotations()) * 2 + 2)
    width = np.sqrt(np.log(2 * bins / 0.01) / (2 * SHOTS)) + 3e-12
    for name, appointments in rotations().items():
        field, predicted = motion(appointments)
        _, wrong_sign = motion([(i, -angle) for i, angle in appointments])
        settings = ["native", "remounted"]
        if name == "ordered":
            settings += ["half_contrast", "cut_bridge"]
        for setting in settings:
            changed = setting == "half_contrast"
            remount = setting == "remounted"
            kappa = 4 if changed else pulse.candidate.KAPPA
            compiled = pulse.execute(field, kappa)
            # Physical route on the native law is a sector permutation. At κ=4
            # the same π swaps still have unit phase on all integer contrasts.
            incoming = deliver @ compiled
            if remount:
                incoming = pulse.execute(turns) @ incoming
            receiver_prefix = pulse.inverse(turns) if remount else []
            if remount:
                np.testing.assert_allclose(pulse.execute(turns), frame, atol=1e-12)
            series = []
            passed = np.zeros(4, dtype=int)
            for spec in SOURCES:
                src = pulse.execute(pulse.prep(spec), kappa) @ source
                ideal = pulse.standard_source(spec)
                for rd in READS:
                    analyzer = pulse.execute(receiver_prefix + pulse.reader(rd), kappa)
                    ideal_analyzer = pulse.reading_vector(rd)
                    p = with_refusal(abs(analyzer @ incoming @ src) ** 2)
                    if setting == "cut_bridge":
                        p = np.array([0, 0, 0, 1.0])
                    forecast = with_refusal(
                        abs(ideal_analyzer @ predicted @ ideal) ** 2
                    )
                    scalar_sham = with_refusal(abs(ideal_analyzer @ ideal) ** 2)
                    old_labels = with_refusal(
                        abs(ideal_analyzer @ frame @ predicted @ ideal) ** 2
                    )
                    sign_sham = with_refusal(
                        abs(ideal_analyzer @ wrong_sign @ ideal) ** 2
                    )
                    obs = pulse.grid_counts(rng, SHOTS, p)
                    passed += [
                        max(abs(obs / SHOTS - h)) <= width
                        for h in (forecast, scalar_sham, old_labels, sign_sham)
                    ]
                    series.append(obs.tolist())
            totals[f"{name}/{setting}"] = passed.tolist()
            all_banks.append(dict(setting=setting, word=name, counts=series))
    return dict(
        grammar=ballot(),
        bracket_frame_pivots=sample()[-1],
        stations=site["owners"],
        route_contacts=site["path_contacts"],
        words={
            k: [[i, round(float(a), 10)] for i, a in w] for k, w in rotations().items()
        },
        sources=[list(s) for s in SOURCES],
        analyzers=[list(r) for r in READS],
        raw_clicks=[r[1] for r in site["address_clicks"]]
        + [site["address_clicks"][0][2]],
        seed=SEED,
        starts_per_bank=SHOTS,
        simultaneous_radius=round(float(width), 10),
        forecast_order=[
            "module",
            "scalar",
            "unconverted_receiving_frame",
            "sign_reversal",
        ],
        covered=totals,
        banks=all_banks,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "records.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1) + "\n")
    elif json.loads(path.read_text()) != result:
        raise SystemExit("Matter class records differ")
    print(result["covered"])
