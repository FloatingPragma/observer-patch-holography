"""Prepare one load, apply four material generators, and compare occupied rails."""

import argparse
from functools import cache
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame.field import (
    transfer,
)
from model.experiments.interval_echo_consensus.signed_effects.sequential import (
    grid_counts,
)
from model.experiments.interval_echo_gauge.arrival_ladder.measure import read_loads

HERE = Path(__file__).resolve().parent
PORTS = (0, 5, 11)
G, GY = 1.5, 1.25
SHOTS = 2**15
SEED = 641533
DWELL = 0.83
HOLD = (0.6, -0.8, 0.45, 0.7)
THETA = (0.0, np.pi / 2)
SOURCES = [
    ("upper",),
    ("lower",),
    ("ref",),
    ("pair", 0, 1, 0.0),
    ("pair", 0, 1, np.pi / 2),
    ("pair", 0, 2, 0.0),
    ("pair", 1, 2, np.pi / 2),
]
READS = [
    ("ports",),
    *(("pair", i, j, p) for i, j in combinations(range(3), 2) for p in THETA),
]


def phase(j, angle):
    # Positive private dwell produces exp(+i angle) at the named occupied site.
    return ("private", j, j, (-angle) % (2 * np.pi) / candidate.KAPPA)


def rx(i, j, angle):
    # Blank 01/10 has H=(κ+1/2)I-X/2; 00 costs zero. Cancel its common
    # sector phase against 00 explicitly, since the third rail can observe it.
    t = angle % (4 * np.pi)
    c = (candidate.KAPPA + 0.5) * t
    return [("contact", i, j, t), phase(i, c), phase(j, c)]


def rz(i, j, angle):
    return [phase(i, angle / 2), phase(j, -angle / 2)]


def ry(i, j, angle):
    return rz(i, j, np.pi / 2) + rx(i, j, angle) + rz(i, j, -np.pi / 2)


def inverse(word):
    # Native contrast is 8. The contact eigenvalues are 0,8,9 and
    # the private eigenvalues are 0,8. Each complete turn is identity.
    return [(kind, a, b, (-t) % (2 * np.pi)) for kind, a, b, t in reversed(word)]


def field_word(direction, time=DWELL):
    x, y, z, b = direction
    radius = np.linalg.norm([x, y, z])
    turns = []
    if radius:
        polar = np.arctan2(np.hypot(x, y), z)
        azimuth = np.arctan2(y, x)
        axis = ry(0, 1, -polar) + rz(0, 1, -azimuth)
        turns = inverse(axis) + rz(0, 1, time * G * radius) + axis
    return turns + [phase(i, time * GY * b / 2) for i in (0, 1)]


@cache
def native(kappa):
    f = federation(0)
    assert all(f.degree[p] == 5 for p in PORTS)
    assert all(PORTS[j] in f.neighbours[PORTS[i]] for i, j in combinations(range(3), 2))
    law = candidate.native_law
    assert law.exact_checks()["checks"][
        "reference_pair_gives_n_squared_at_reference_zero"
    ]
    ex, con = law.local_term_parts()
    full = ex + kappa * con
    bases = np.eye(law.DIM)
    out = {}
    for i, j in combinations(range(3), 2):
        sub = bases[:, [law.SRC.index((int(a == i), int(a == j))) for a in range(3)]]
        h = sub.T @ full @ sub
        np.testing.assert_allclose(full @ sub, sub @ h, atol=1e-13)
        out["contact", i, j] = h
    for i in range(3):
        out["private", i, i] = np.diag([kappa * int(a == i) for a in range(3)])
    return out


@cache
def gate(kind, i, j, time, kappa):
    return expm(-1j * time * native(kappa)[kind, i, j])


def execute(word, kappa=candidate.KAPPA):
    u = np.eye(3, dtype=complex)
    for step in word:
        u = gate(*step, kappa) @ u
    return u


def prep(spec):
    if spec[0] in ("ports", "upper"):
        return []
    if spec[0] in ("lower", "ref"):
        return rx(0, 1 if spec[0] == "lower" else 2, np.pi)
    _, i, j, angle = spec
    park = [] if i == 0 else rx(0, i, np.pi)
    return park + ry(i, j, -np.pi / 2) + [phase(j, angle)]


def reader(spec):
    if spec[0] == "ports":
        return []
    _, i, j, angle = spec
    return inverse(ry(i, j, -np.pi / 2) + [phase(j, angle)])


def reading_vector(spec):
    """An ideal pair analyzer independently specified by its columns."""
    u = np.eye(3, dtype=complex)
    if spec[0] == "ports":
        return u
    _, i, j, angle = spec
    u[np.ix_([i, j], [i, j])] = np.array(
        [[1, -1], [np.exp(1j * angle), np.exp(1j * angle)]]
    ) / np.sqrt(2)
    return u.conj().T


def standard_source(spec):
    if spec[0] != "pair":
        return np.eye(3)[:, ("upper", "lower", "ref").index(spec[0])]
    _, i, j, p = spec
    v = np.zeros(3, complex)
    v[i], v[j] = 1 / np.sqrt(2), np.exp(1j * p) / np.sqrt(2)
    return v


def exterior(direction, time):
    # Independently use the original five-mode CAR action: weak coordinates 3,4.
    u = transfer(3, 4)
    x, y, z = u + u.T, -1j * (u - u.T), transfer(3, 3) - transfer(4, 4)
    hyper6 = -2 * sum(transfer(i, i) for i in range(3)) + 3 * (
        transfer(3, 3) + transfer(4, 4)
    )
    h = (
        G * (direction[0] * x + direction[1] * y + direction[2] * z) / 2
        + GY * direction[3] * hyper6 / 6
    )
    basis = np.eye(32)[:, [8, 16, 0]]
    h3 = basis.T @ h @ basis
    np.testing.assert_allclose(h @ basis, basis @ h3, atol=1e-14)
    return expm(1j * time * h3)


def report():
    assert read_loads() == [0, 1, 2]
    directions = [*np.eye(4).tolist(), HOLD]
    rng = np.random.default_rng(SEED)
    out = []
    source = np.eye(3)[:, 0]
    # For all bins, three outcomes and a union bound over all planned settings.
    delta = np.sqrt(
        np.log(2 * 3 * len(directions) * len(SOURCES) * len(READS) / 0.01) / (2 * SHOTS)
    )
    wrong_k = []
    for index, axis in enumerate(directions):
        word = field_word(axis)
        u, comparison = execute(word), exterior(axis, DWELL)
        np.testing.assert_allclose(u, comparison, atol=5e-13, rtol=0)
        flipped = exterior([axis[0], -axis[1], axis[2], axis[3]], DWELL)
        double_y = exterior([*axis[:3], 2 * axis[3]], DWELL)
        rows = []
        for src in SOURCES:
            s = execute(prep(src)) @ source
            ideal_s = standard_source(src)
            for rd in READS:
                read = execute(reader(rd))
                v = read @ u @ s
                p = abs(v) ** 2
                counts = grid_counts(rng, SHOTS, p)
                measured = counts / SHOTS
                predictions = [
                    abs(reading_vector(rd) @ m @ ideal_s) ** 2
                    for m in (comparison, flipped, double_y)
                ]
                matches = [
                    bool(max(abs(a - measured)) <= delta + 1e-11) for a in predictions
                ]
                if index == len(directions) - 1:
                    mutant = (
                        execute(reader(rd), 4)
                        @ execute(word, 4)
                        @ execute(prep(src), 4)
                        @ source
                    )
                    wrong_k.append(float(max(abs(abs(mutant) ** 2 - measured))))
                rows.append(
                    dict(
                        source=SOURCES.index(src),
                        read=READS.index(rd),
                        counts=counts.tolist(),
                        matched=matches,
                    )
                )
        out.append(dict(direction=list(axis), banks=rows))
    return dict(
        ports=list(PORTS),
        couplings=[G, GY],
        dwell=DWELL,
        shots=SHOTS,
        seed=SEED,
        simultaneous_radius=float(delta + 1e-11),
        sources=[list(s) for s in SOURCES],
        readers=[list(s) for s in READS],
        panels=out,
        matched_totals=[
            sum(b["matched"][j] for p in out for b in p["banks"]) for j in range(3)
        ],
        held_changed_law_excluded=int(sum(v > delta + 1e-11 for v in wrong_k)),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    result = report()
    dest = HERE / "results.json"
    if ap.parse_args().write:
        dest.write_text(json.dumps(result, indent=2) + "\n")
    elif json.loads(dest.read_text()) != result:
        raise SystemExit("Scalar grammar replay differs")
    print(result["matched_totals"], result["held_changed_law_excluded"])
