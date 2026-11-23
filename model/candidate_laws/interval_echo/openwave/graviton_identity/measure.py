"""Test the gravitational information in an Einstein decoder and native material records."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.linalg import eigh, expm
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "code"))
from matched_einstein_limit import frame
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field import dark_contact
from model.executable_model import cell, federation, rotations
from model.experiments.interval_echo_einstein_station import backreaction

LAW = candidate.native_law
KAPPA = candidate.KAPPA
TIMES = (1, 2, 4, 8)
SHOTS = 2**18


def metric_read(name, k, h):
    """Unit second phase derivative of h*f(k.x), signature (-,+,+,+)."""
    signature = (-1, 1, 1, 1)
    riemann = {
        (a, b, c, d): F(
            int(
                h[a, d] * k[b] * k[c]
                + h[b, c] * k[a] * k[d]
                - h[a, c] * k[b] * k[d]
                - h[b, d] * k[a] * k[c]
            ),
            2,
        )
        for a, b, c, d in product(range(4), repeat=4)
    }
    ricci = {
        (b, d): sum(signature[a] * riemann[a, b, a, d] for a in range(4))
        for b, d in product(range(4), repeat=2)
    }
    trace = sum(signature[a] * ricci[a, a] for a in range(4))
    tensor = tuple(
        ricci[b, d] - (trace * signature[b] / 2 if b == d else 0)
        for b, d in frame.PAIRS
    )
    charges = frame.apply(frame.design(), tensor)
    assert frame.decode(charges) == tensor
    # Relative acceleration in a geodesic arm, with twice its component to
    # normalize the plus and cross responses to one.
    arms = [(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, -1, 0), (0, 0, 1)]
    delays = [
        -2
        * sum(
            arm[i] * arm[j] * riemann[0, i + 1, 0, j + 1]
            for i, j in product(range(3), repeat=2)
        )
        / sum(v * v for v in arm)
        for arm in arms
    ]
    return {
        "profile": name,
        "wave_covector": list(k),
        "curvature_components_nonzero": sum(bool(v) for v in riemann.values()),
        "ten_Einstein_contractions": list(map(str, charges)),
        "decoded_components": list(map(str, tensor)),
        "twice_tidal_arm_acceleration_x_y_dplus_dminus_z": list(map(str, delays)),
    }


def tensor_records():
    xy = np.zeros((4, 4), dtype=int)
    xy[1, 2] = xy[2, 1] = 1
    plus = np.diag([0, 1, -1, 0])
    null = (1, 0, 0, -1)
    gauge = np.outer(null, (0, 1, 0, 0)) + np.outer((0, 1, 0, 0), null)
    rows = [
        metric_read(name, k, h)
        for name, k, h in (
            ("flat", null, np.zeros((4, 4), dtype=int)),
            ("vacuum_plus", null, plus),
            ("vacuum_cross", null, xy),
            ("coordinate_change", null, gauge),
            ("breathing_control", null, np.diag([0, 1, 1, 0])),
            ("off_null_plus_control", (1, 0, 0, -2), plus),
        )
    ]
    for row in rows[:4]:
        assert row["ten_Einstein_contractions"] == ["0"] * 10
    for row in rows[4:]:
        assert any(x != "0" for x in row["ten_Einstein_contractions"])
    assert rows[1]["twice_tidal_arm_acceleration_x_y_dplus_dminus_z"] == [
        "1",
        "-1",
        "0",
        "0",
        "0",
    ]
    assert rows[2]["twice_tidal_arm_acceleration_x_y_dplus_dminus_z"] == [
        "0",
        "0",
        "1",
        "-1",
        "0",
    ]
    assert rows[3]["curvature_components_nonzero"] == 0
    return {
        "ten_frame_inverse_verified": frame.inverse(frame.design())
        == frame.decoder_matrix(),
        "Einstein_read_rank_on_vacuum_polarizations": 0,
        "tidal_read_rank": 2,
        "profiles": rows,
    }


def source_calibration():
    terms = LAW.local_term_parts()
    h = terms[0] + KAPPA * terms[1]
    vac = np.eye(LAW.N_SRC)[:, LAW.SRC.index((0, 0))]
    prepared, gaps = [], []
    for d in (0, 2):
        columns = dark_contact.lift(d, 1, -1)
        np.testing.assert_array_equal(LAW._write_pattern().T @ columns[LAW.N_SRC :], 0)
        np.testing.assert_array_equal(terms[0] @ columns, columns * 0)
        np.testing.assert_array_equal(terms[1] @ columns, columns * (d * d))
        prepared.append(np.kron(backreaction.collar_vector(d, 1, -1), vac))
        gaps.append(KAPPA * d * d)
    zero, excited = prepared
    bright, dark = (zero + excited) / np.sqrt(2), (zero - excited) / np.sqrt(2)
    time = np.pi / 32
    control = {}
    for name, generator in (("registered", h), ("exchange_only", terms[0])):
        pulse = expm(-1j * time * generator)
        control[name] = [
            float(abs(dark.conj() @ pulse @ w) ** 2) for w in (zero, bright)
        ]
    np.testing.assert_allclose(control["registered"], [0.5, 1], atol=1e-12)
    np.testing.assert_allclose(control["exchange_only"], [0.5, 0], atol=1e-12)
    return {
        "occupied_local_energy_eigenvalues": gaps,
        "coherent_source_mean_energy": sum(gaps) / 2,
        "complete_load_columns_per_seam_word": LAW.N_SRC,
        "local_read_exposure": time,
        "Ramsey_yes_zero_word_then_superposed_word": control,
    }


def ring():
    board = cell()
    apex = 0
    start = min(board.neighbours[apex])
    g = min(
        g
        for g in rotations()
        if g[apex] == apex
        and g[start] in set(board.neighbours[apex]) & set(board.neighbours[start])
    )
    ports = [start]
    while g[ports[-1]] != start:
        ports.append(g[ports[-1]])
    assert len(ports) == 5
    return ports, [tuple(sorted((a, b))) for a, b in zip(ports, ports[1:] + ports[:1])]


def receiver_records(h, launches, groups):
    levels, modes = eigh(h)
    t0 = modes[launches].T
    result = []
    for time in TIMES:
        wave = modes @ (np.exp(-1j * time * levels)[:, None] * t0)
        np.testing.assert_allclose(
            wave.conj().T @ wave, np.eye(len(launches)), atol=2e-12
        )
        read = []
        for column, group in enumerate(groups):
            hits = abs(wave[group, column]) ** 2
            read.append(
                np.r_[hits, np.linalg.norm(np.delete(wave[:, column], group)) ** 2]
            )
        result.append(read)
    return np.array(result)


def check_full_seam(board, seams, contact, launches, groups, frozen, blank):
    """Keep all eleven states of the first seam in the independently built flight operator."""
    full = backreaction.lift_one_load(board, seams[0])
    memories = 1 + LAW.N_LAB
    # Only the other four seams are restricted to their invariant zero-energy word.
    other_contacts = np.zeros((board.n, board.n))
    for edge in seams[1:]:
        other_contacts[np.ix_(edge, edge)] += contact
    full -= sparse.kron(
        sparse.csr_matrix(other_contacts), sparse.eye(memories), format="csr"
    )
    zero = backreaction.collar_vector(0, 1, -1)
    excited = backreaction.collar_vector(2, 1, -1)
    words = np.column_stack(
        (zero, (zero + excited) / np.sqrt(2), np.eye(memories)[:, 0])
    )
    initial = np.zeros((board.n, memories, 3))
    initial[launches[0]] = words
    time = TIMES[1]
    centered = full - 6 * KAPPA * sparse.eye(full.shape[0], format="csr")
    evolved = expm_multiply(-1j * time * centered, initial.reshape(-1, 3))
    load = (abs(evolved.reshape(board.n, memories, 3)) ** 2).sum(axis=1)
    hits = load[groups[0]].T
    outcomes = np.column_stack((hits, np.delete(load, groups[0], axis=0).sum(axis=0)))
    expected = np.stack((frozen[1, 0], frozen[1, 0], blank[0, 1, 0]))
    np.testing.assert_allclose(outcomes, expected, atol=1e-11, rtol=0)
    return {
        "dimension": full.shape[0],
        "retained_memory_states": memories,
        "source_seam": seams[0],
        "launch": launches[0],
        "exposure": time,
        "word_order": ["zero", "coherent_source", "blank"],
        "outcomes": outcomes.tolist(),
        "max_absolute_difference_from_reduced_flight": float(
            np.max(abs(outcomes - expected))
        ),
    }


def source_and_tracer():
    board = federation(0)
    ports, seams = ring()
    matrix = candidate.generator(board).toarray()
    s = [LAW.SRC.index((1, 0)), LAW.SRC.index((0, 1))]
    contact = LAW.seam_term(KAPPA)[np.ix_(s, s)]
    # A zero occupied collar acts by zero while a blank collar has the material contact.
    for edge in seams:
        matrix[np.ix_(edge, edge)] -= contact
    launches = [b for a, b in board.external if a < 12]
    groups = [list(range((p // 12) * 12, (p // 12 + 1) * 12)) for p in launches]
    frozen = receiver_records(matrix, launches, groups)
    blank = []
    for edge in seams:
        gate = matrix.copy()
        gate[np.ix_(edge, edge)] += contact
        blank.append(receiver_records(gate, launches, groups))
    blank = np.stack(blank)
    full_seam_check = check_full_seam(
        board, seams, contact, launches, groups, frozen, blank
    )
    theta = 2 * np.pi * np.arange(5) / 5
    shapes = np.column_stack((np.cos(2 * theta), np.sin(2 * theta)))
    rng = np.random.default_rng(2026112301)
    mixtures = {"isotropic": np.ones(5) / 5}
    for j, shape in enumerate(("plus", "cross")):
        for sign, name in ((1, "positive"), (-1, "negative")):
            mixtures[f"{shape}_{name}"] = (1 + sign * shapes[:, j]) / 5
    rows, predicted = {}, {}
    for name, p in mixtures.items():
        np.testing.assert_allclose(p.sum(), 1, atol=1e-14)
        ramsey = (1 + p) / 2
        q = p @ shapes
        blank_read = np.tensordot(p, blank, axes=(0, 0))
        predicted[name] = blank_read
        source_counts = rng.binomial(SHOTS, ramsey)
        cohorts = {}
        for arm, reading in (
            ("energy_pattern", frozen),
            ("one_blank_seam", blank_read),
        ):
            counts = [
                [rng.multinomial(SHOTS, prob / prob.sum()).tolist() for prob in span]
                for span in reading
            ]
            cohorts[arm] = {"outcome_probabilities": reading.tolist(), "counts": counts}
        rows[name] = {
            "mixture_weights": p.tolist(),
            "quadrupole_coordinates": q.tolist(),
            "source_ramsey_probabilities": ramsey.tolist(),
            "source_ramsey_yes_counts": source_counts.tolist(),
            "inferred_coordinates_from_source_counts": (
                (2 * source_counts / SHOTS - 1) @ shapes
            ).tolist(),
            "remote_tracer": cohorts,
        }
    contrast_tv = {}
    for shape in ("plus", "cross"):
        tv = 0.5 * np.abs(
            predicted[f"{shape}_positive"] - predicted[f"{shape}_negative"]
        ).sum(axis=-1)
        contrast_tv[shape] = tv.tolist()
        assert float(tv.max()) > 2e-4
    return {
        "native_KAPPA": KAPPA,
        "board_cells": board.cells,
        "pentagon_ports": ports,
        "source_seams": seams,
        "launches": launches,
        "remote_receiver_ports": groups,
        "quiet_outcome_index": 12,
        "native_exposures": TIMES,
        "shots_per_setting": SHOTS,
        "simulated_count_seed": 2026112301,
        "source_calibration": source_calibration(),
        "eleven_state_crosscheck": full_seam_check,
        "all_time_response_rank_on_frozen_energy_quadrupoles": 0,
        "exact_remote_tv_either_energy_quadrupole": 0,
        "blank_control_tv_by_time_and_launch": contrast_tv,
        "source_arms": rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group()
    group.add_argument(
        "--output", type=Path, help="write a fresh receipt to a new file"
    )
    group.add_argument(
        "--check",
        type=Path,
        default=HERE / "results.json",
        help="compare to a saved receipt",
    )
    args = ap.parse_args()
    out = {
        "einstein_and_tidal_comparison": tensor_records(),
        "registered_contact_run": source_and_tracer(),
        "sha256": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), HERE / "DESIGN.md")
        },
    }
    serialized = json.dumps(out, indent=2, allow_nan=False)
    if args.output:
        with args.output.open("x") as f:
            f.write(serialized + "\n")
    else:
        backreaction.compare(json.loads(serialized), json.loads(args.check.read_text()))
        print(
            "Replayed six continuum profiles and five native angular sources with their controls."
        )


if __name__ == "__main__":
    main()
