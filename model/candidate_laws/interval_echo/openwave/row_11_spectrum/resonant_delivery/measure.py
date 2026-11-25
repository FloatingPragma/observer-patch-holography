"""Resonant energy transfer from a charged vacancy dimer into a compensated dimer."""

import argparse
from functools import cache
from itertools import product
import json
from math import acos, atan, cos, pi, sin, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_01_charge.measure import (
    sharp_contrast,
)
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from .enclosure import certified_transfer

HERE = Path(__file__).resolve().parent
SHOTS = 2**22
CAL_SHOTS = 2**17
TICKS = 10**12
IDENTITY = np.eye(2)
P = np.diag([0.0, 1.0])


@cache
def fixture():
    b = federation(0)
    ga, gb = min(edge for edge in b.external if edge[0] // 12 == 0)

    def regular(g):
        return [u for u in b.neighbours[g] if u // 12 == g // 12 and b.degree[u] == 5]

    fa, pa = regular(ga)[:2]
    fb, pb = regular(gb)[:2]
    parked = next(
        u
        for u in range(12 * (gb // 12), 12 * (gb // 12 + 1))
        if b.degree[u] == 5 and u not in (fb, pb)
    )
    rails = [[fa, ga], [fb, gb]]
    return b, rails, [pa, pb], parked


@cache
def local_tables(strength):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + strength * contrast
    embed = np.eye(99)[:, :9]
    blank = embed.T @ h @ embed
    np.testing.assert_allclose(h @ embed, embed @ blank, atol=1e-13)
    stored = np.kron(label(2, 1)[:, None], np.eye(9))
    marked = stored.T @ h @ stored
    np.testing.assert_allclose(h @ stored, stored @ marked, atol=1e-13)
    expected = np.diag([4 * strength * int(w != (1, 1)) for w in law.SRC])
    np.testing.assert_allclose(marked, expected, atol=1e-13)
    silent = np.kron(label(0, 1)[:, None], np.eye(9))
    np.testing.assert_allclose(h @ silent, 0, atol=1e-13)
    blocks = [np.rint(2 * t).astype(int) for t in (blank, marked)]
    for block, t in zip(blocks, (blank, marked), strict=True):
        np.testing.assert_allclose(block / 2, t, atol=1e-13)
    return blocks


def reduced(which, setting, strength, *, compensated=True):
    """Rebuild two native energy matrices from the physical graph and all local outputs."""
    board, rails, posts, park = fixture()
    blank, marked = local_tables(strength)
    active = [tuple(sorted(r)) for r in rails]
    gateways = [arm[1] for arm in rails]
    crossing = tuple(sorted(gateways))
    fixed = [tuple(sorted((g, p))) for g, p in zip(gateways, posts, strict=True)]
    written = fixed if setting == "matched_cut" else [crossing]
    if setting == "detuned":
        written = [*written, fixed[1]]
    variable = {*rails[0], *rails[1], park}
    touched = {edge for edge in board.seams if variable.intersection(edge)}
    assert set(active + written) <= touched
    # Untouched seams see unit loads. Every other incident contact carries
    # an attainable inert zero-contrast mark; the fixed posts stay filled.
    silent = touched - set(active + written)
    sites = sorted(variable | set(posts))
    ids = {u: j for j, u in enumerate(sites)}
    inputs = []
    for positions in product(range(2), repeat=len(which)):
        word = [1] * len(sites)
        for arm, pos in zip(which, positions, strict=True):
            word[ids[rails[arm][pos]]] = 0
        if 1 in which and compensated:
            word[ids[park]] = 2
        inputs.append(tuple(word))
    indices = {w: i for i, w in enumerate(inputs)}
    matrix = np.zeros((len(inputs), len(inputs)), dtype=int)
    for j, w in enumerate(inputs):
        matrix[j, j] = sum(
            2 * strength * (6 - int(board.degree[u])) * (w[ids[u]] ** 2 - 1)
            for u in sites
        )
        for edge, term in [(e, blank) for e in active] + [(e, marked) for e in written]:
            a, b = (ids[u] for u in edge)
            for k, element in enumerate(term[:, 3 * w[a] + w[b]]):
                if element:
                    z = list(w)
                    z[a], z[b] = divmod(k, 3)
                    assert tuple(z) in indices
                    matrix[indices[tuple(z)], j] += element
    return matrix / 2, len(silent)


def optical_frame(h):
    """Implement a real eigenbasis rotation with the material dimer and a private phase."""
    d = h[1, 1] - h[0, 0]
    theta = atan(1 / d) / 2
    expected = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
    blank, _ = local_tables(int(candidate.KAPPA))
    ix = [3, 1]  # the two orders of a vacancy and a unit
    beam = expm(-1j * (2 * theta) * blank[np.ix_(ix, ix)] / 2)

    # Phase the far (degree-five) vacancy with a nonnegative native dwell
    # before and after the balanced exchange. A common phase drops out.
    def private(angle):
        dwell = (angle % (2 * pi)) / candidate.KAPPA
        return np.diag([np.exp(1j * candidate.KAPPA * dwell), 1])

    frame = private(pi / 2) @ beam @ private(-pi / 2)
    frame = frame / frame[0, 0] * abs(frame[0, 0])
    np.testing.assert_allclose(frame, expected, atol=2e-13)
    energies = np.diag(frame.conj().T @ h @ frame).real
    np.testing.assert_allclose(h @ frame, frame * energies, atol=2e-13)
    return energies, frame


def counts(probabilities, shots, rng):
    p = np.maximum(0, np.rint(np.asarray(probabilities) * TICKS).astype(np.int64))
    p[np.argmax(p)] += TICKS - int(sum(p))
    return rng.multinomial(shots, p / TICKS).tolist()


def clock(h, rng):
    plus = np.ones(2) / sqrt(2)
    banks = []
    for j in range(32):
        bright = abs(np.vdot(plus, expm(-1j * j * h / 32) @ plus)) ** 2
        banks.append(counts([bright, 1 - bright], CAL_SHOTS, rng))
    signal = np.diff([(c - d) / CAL_SHOTS for c, d in banks])
    cos_step = np.dot(signal[1:-1], signal[:-2] + signal[2:]) / (
        2 * np.dot(signal[1:-1], signal[1:-1])
    )
    return {"bright_other": banks, "gap_from_counts": round(32 * acos(cos_step), 8)}


def setting(name, strength, rng):
    reg = name if name in ("matched_cut", "detuned") else "native"
    ha, _ = reduced((0,), reg, strength)
    hb, _ = reduced((1,), reg, strength)
    h, silent = reduced((0, 1), reg, strength)
    ea, a = optical_frame(ha)
    eb, b = optical_frame(hb)
    frame = np.kron(a, b)
    h0 = np.kron(ha, IDENTITY) + np.kron(IDENTITY, hb)
    v = h - h0
    bridge = 0 if reg == "matched_cut" else -4 * strength
    np.testing.assert_allclose(v, bridge * np.kron(P, P), atol=2e-13)
    s, d = frame[:, 2], frame[:, 1]
    eig, vec = eigh(h)
    weight_a = abs(vec.conj().T @ s) ** 2
    weight_b = abs(vec.conj().T @ d) ** 2
    spectral = []
    # Sum tied eigenspaces before testing swap covariance.
    for i, energy in enumerate(eig):
        if i and abs(energy - eig[i - 1]) < 1e-9:
            spectral[-1][1] += weight_a[i]
            spectral[-1][2] += weight_b[i]
        else:
            spectral.append([energy, weight_a[i], weight_b[i]])
    if reg != "detuned":
        np.testing.assert_allclose(
            [r[1] for r in spectral], [r[2] for r in spectral], atol=1e-12
        )
    a_gap, b_gap = ea[1] - ea[0], eb[1] - eb[0]
    hj = np.kron(IDENTITY, hb)
    current = 1j * (h @ hj - hj @ h)
    settings = {
        "silent_incident_seams": silent,
        "local_charged_levels": ea.round(9).tolist(),
        "local_neutral_levels": eb.round(9).tolist(),
        "independent_parent_and_daughter_totals": [
            round(ea[1] + eb[0], 9),
            round(ea[0] + eb[1], 9),
        ],
        "shell_coupling": round(float(np.vdot(d, v @ s).real), 12),
        "shell_neutral_current_modulus": round(abs(np.vdot(d, current @ s)), 12),
        "full_energy_source_daughter_weights": np.round(spectral, 9).tolist(),
        "charged_clock": clock(ha, rng),
        "neutral_clock": clock(hb, rng),
        "joint_times": [],
    }
    if name == "native":
        settings["certified_delivery"] = certified_transfer((2 * h).astype(int), 64)
    for initial, vec_in in (("charged_excited", s), ("neutral_excited", d)):
        energy0 = float(np.vdot(vec_in, h @ vec_in).real)
        for t in (0, 32, 64, 128):
            y = expm(-1j * t * h) @ vec_in
            p = abs(frame.conj().T @ y) ** 2
            iq = abs(y[-1]) ** 2
            shots = counts(p, SHOTS, rng)
            contact_counts = counts([iq, 1 - iq], SHOTS, rng)
            # Separate banks measure bare local energies and the occupied
            # interface. Their sum estimates the complete conserved energy.
            empirical = (
                ea[0]
                + eb[0]
                + a_gap * (shots[2] + shots[3]) / SHOTS
                + b_gap * (shots[1] + shots[3]) / SHOTS
                + bridge * contact_counts[0] / SHOTS
            )
            e_after = (
                ea[0]
                + eb[0]
                + a_gap * (p[2] + p[3])
                + b_gap * (p[1] + p[3])
                + bridge * iq
            )
            np.testing.assert_allclose(e_after, energy0, atol=6e-10, rtol=0)
            settings["joint_times"].append(
                {
                    "initial": initial,
                    "exposure": t,
                    "energy_probability_gg_ge_eg_ee": p.round(10).tolist(),
                    "energy_counts_gg_ge_eg_ee": shots,
                    "both_gateway_vacant_other": contact_counts,
                    "neutral_excitation_from_counts": round(
                        b_gap * (shots[1] + shots[3]) / SHOTS, 8
                    ),
                    "conserved_total_energy": round(energy0, 9),
                    "total_energy_from_both_banks": round(empirical, 8),
                }
            )
    return settings


def run():
    b, rails, posts, park = fixture()
    rng = np.random.default_rng(84029)
    assert [sharp_contrast(i, 0) for i in range(3)] == [0, 1, 2]

    def charges(compensated):
        readings = set()
        for i, j in product(range(2), repeat=2):
            load = np.ones(b.n, dtype=int)
            load[[rails[0][i], rails[1][j]]] = 0
            load[park] = 2 if compensated else 1
            readings.add(
                tuple(
                    sum(
                        sharp_contrast(int(load[p]), 0) - 1
                        for p in range(12 * (r[1] // 12), 12 * (r[1] // 12 + 1))
                    )
                    for r in rails
                )
            )
        assert len(readings) == 1
        return list(readings.pop())

    nominal, _ = reduced((0, 1), "native", 8)
    uncompensated, _ = reduced((0, 1), "native", 8, compensated=False)
    np.testing.assert_array_equal(nominal - uncompensated, 24 * np.eye(4))
    empty, _ = reduced((0,), "native", 8)
    light_levels, light_basis = optical_frame(empty)
    np.testing.assert_allclose(
        empty @ light_basis[:, 1], light_levels[1] * light_basis[:, 1], atol=2e-13
    )
    return {
        "cells": b.cells,
        "rails_far_near": rails,
        "cut_posts": posts,
        "neutral_compensating_excess": park,
        "local_charge_records": charges(True),
        "uncompensated_receiver": {
            "local_charge_records": charges(False),
            "uniform_energy_subtracted": 24,
            "maximum_difference_after_shift": 0,
        },
        "initially_filled_receiver": {
            "invariant_charged_source_dimension": 2,
            "receiver_absorption_probability": 0,
        },
        "calibration_shots": CAL_SHOTS,
        "shots_per_bank": SHOTS,
        "seed": 84029,
        "settings": {
            name: setting(name, k, rng)
            for name, k in (
                ("native", 8),
                ("matched_cut", 8),
                ("detuned", 8),
                ("half_contrast", 4),
            )
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "results.json"
    data = run()
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("Resonant delivery receipt differs")
    for name, row in data["settings"].items():
        print(name, row["joint_times"][2])


if __name__ == "__main__":
    main()
