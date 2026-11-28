"""Replay oriented face-transfer records on the same tagged local arrivals."""

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import pairwise
import json
from math import inf, nextafter, pi, sqrt
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.stats import beta

from ...fixture import FanPath, archived_fixture
from ...spectral import verified_inertia


HERE = Path(__file__).resolve().parent
READS = HERE.parent / "word_reads.json"
COUNTS = HERE.parent / "witness.json"
LABEL = [1045, 0, 334, 365]
Y = np.array([[0, -1j], [1j, 0]])
TOL = 2e-12


def local_face(path, label):
    _, index, a, b = label
    words = path.blocks[3, 0].words[3]
    bits = path.field.bits
    faces = []
    for ring in path.field.rings:
        route = list(pairwise((*ring, ring[0])))
        mask = sum(bits[tuple(sorted(edge))] for edge in route)
        clockwise = sum(bits[u, v] for u, v in route if u < v)
        h = np.zeros((len(words), len(words)))
        for i, w in enumerate(words):
            if w & mask in (clockwise, clockwise ^ mask):
                h[i, i] = 0.5
                h[words.index(w ^ mask), i] = -0.5
        faces.append(h)
    np.testing.assert_array_equal(sum(faces), path.field.ring_laplacian(words))
    pair = [words.index(a), words.index(b)]
    other = [i for i in range(len(words)) if i not in pair]
    np.testing.assert_array_equal(faces[index][np.ix_(pair, other)], 0)
    h = faces[index][np.ix_(pair, pair)]
    np.testing.assert_array_equal(h, [[0.5, -0.5], [-0.5, 0.5]])
    return h


def initial_stocks(path):
    """Build the original bypassed Bell input and the delayed write after a band pass."""
    source = path.source()
    rng = np.random.default_rng(923)
    energy, modes = eigsh(
        path.h, k=4, which="SA", v0=rng.standard_normal(len(path.basis)), tol=1e-13
    )
    assert np.count_nonzero(energy < 4.25) == 2
    v = modes[:, :2]
    stock = v @ (v.T @ source)
    admitted = float(np.vdot(stock, stock).real / 2)
    source_block = path.blocks[0, 0]
    spectator_words = defaultdict(list)
    for i in range(source_block.rows.start, source_block.rows.stop):
        _, words = path.basis[i]
        spectator_words[words[1:]].append((i, source_block.words[0].index(words[0])))
    fresh = np.zeros((len(source), 8))
    for locations in spectator_words.values():
        old = sum(
            np.outer(path.code[word], stock[i]) for i, word in locations
        ).ravel() / sqrt(2)
        for i, word in locations:
            fresh[i] = np.outer(old, path.code[word]).ravel()
    written = float(np.vdot(fresh, fresh).real / 2)
    # The two separately energy-resolved reference branches, mixed over labels.
    stationary = np.column_stack(
        [np.outer(v[:, j], v[:, j] @ source) for j in range(2)]
    )
    return {
        "clock_bypassed": (source, 1.0, 1.0),
        "binary_clock": (fresh, admitted, written),
    }, (stationary, admitted, admitted)


def pair_amplitudes(path, wave, label):
    port, _, a, b = label
    paired = defaultdict(dict)
    cell_success = 0
    hub_success = 0
    for i, (x, words) in enumerate(path.basis):
        if path.owner[i] == 3:
            weight = np.vdot(wave[i], wave[i]).real / 2
            cell_success += weight
            if x == port:
                hub_success += weight
                if words[-1] in (a, b):
                    paired[words[:-1]][words[-1]] = wave[i].reshape(-1, 2)
    assert paired and all(set(data) == {a, b} for data in paired.values())
    amps = np.concatenate(
        [np.stack([group[a], group[b]], axis=1) for group in paired.values()]
    )
    return amps, float(cell_success), float(hub_success)


def head_meters(amplitudes, h, refusals, dephase_words=False):
    """Propagate the isolated two-word contact and read both binary records."""
    cohorts = {}
    banks = []
    for pulse in (1, 3):
        unitary = expm(-1j * pulse * pi / 2 * h)
        reference_counts = []
        for q in (1, -1):
            reference = np.array([1, -q * 1j]) / sqrt(2)
            projected = amplitudes @ reference
            if dephase_words:
                out = (abs(unitary) ** 2 @ (abs(projected) ** 2).T).sum(axis=1) / 2
            else:
                out = (abs(projected @ unitary.T) ** 2).sum(axis=0) / 2
            reference_counts.append(out.tolist())
        raw = np.ravel(reference_counts).tolist()
        np.testing.assert_allclose(sum(raw) + sum(refusals), 1, atol=TOL, rtol=0)
        banks.append(
            {"pulse_pi_over_two": pulse, "joint_records": raw, "refusals": refusals}
        )
    openings = np.array([row["joint_records"] for row in banks]).reshape(2, 2, 2)
    np.testing.assert_allclose(
        openings[0].sum(axis=1), openings[1].sum(axis=1), atol=TOL
    )
    rate = 0.5
    for i, name in enumerate(
        ("reference_y_plus", "reference_y_minus", "reference_traced")
    ):
        values = openings[:, i, :] if i < 2 else openings.sum(axis=1)
        denominator = values[0].sum()
        current = rate * (values[0, 1] - values[1, 1])
        cohorts[name] = {
            "face_admission_per_original_start": float(denominator),
            "current_per_original_start": float(current),
            "current_given_face_and_reference": float(current / denominator),
        }
    return {"banks": banks, "current": cohorts}


def run_flight(path, h_face, label, flight, initial):
    wave, admitted, written = initial
    h = path.parked + (flight["material_bridge_rate"] / 4) * path.bridge_term
    elapsed = flight["appointment"]
    wave = expm_multiply(-1j * elapsed * h, wave)
    amps, cell, hub = pair_amplitudes(path, wave, label)
    accepted = float(np.vdot(amps, amps).real / 2)
    refusals = [
        1 - admitted,
        admitted - written,
        written - cell,
        cell - hub,
        hub - accepted,
    ]
    assert min(refusals) >= -TOL
    observed = head_meters(amps, h_face, refusals)
    off = head_meters(amps, np.zeros_like(h_face), refusals)
    word_dephased = head_meters(amps, h_face, refusals, dephase_words=True)
    for control in (off, word_dephased):
        assert (
            max(
                abs(c["current_given_face_and_reference"])
                for c in control["current"].values()
            )
            < TOL
        )
    pair = next(
        row for row in flight["pair_readings"] if row["port_face_and_words"] == LABEL
    )
    pauli = pair["receiver_reference_pauli_moments"]
    pad = flight["observable_absolute_error_upper_per_pair"]
    for direction, bank in zip((1, -1), observed["banks"], strict=True):
        expected = [
            (
                pauli[0][0]
                + ref * pauli[0][2]
                + meter * direction * (pauli[2][0] + ref * pauli[2][2])
            )
            / 4
            for ref in (1, -1)
            for meter in (-1, 1)
        ]
        np.testing.assert_allclose(bank["joint_records"], expected, atol=TOL, rtol=0)
    for name, flag in zip(observed["current"], (1, -1, 0), strict=True):
        d = pauli[0][0] + flag * pauli[0][2]
        j = (pauli[2][0] + flag * pauli[2][2]) / 2
        predicted = j / d
        actual = observed["current"][name]["current_given_face_and_reference"]
        assert abs(actual - predicted) < 3e-8
        # Each moment is enclosed individually; nonzero flags use two moments.
        nd = (1 + abs(flag)) * (pad + TOL)
        assert d > nd
        error = (0.5 + abs(predicted)) * nd / (d - nd)
        observed["current"][name]["moment_prediction_interval"] = [
            nextafter(predicted - error, -inf),
            nextafter(predicted + error, inf),
        ]
    return {
        "flight_bridge_rate": flight["material_bridge_rate"],
        "appointment": elapsed,
        "isolated_face": observed,
        "face_turned_off": off,
        "word_directions_dephased": word_dephased,
    }


def count_intervals(data, flights):
    """Condition each already-saved multinomial bank on the two Y subcounts."""
    assert data["counter_order"] == ["++", "+-", "-+", "--", "refused"]
    alpha = 0.01 / (2 * 2 * len(data["cohorts"]))
    output = []
    for cohort in data["cohorts"]:
        row = next(b for b in cohort["settings"] if b["axes"] == "YY")
        counts = row["counts"]
        assert sum(counts) == data["preparations_per_setting"]
        source = next(
            f
            for f in flights
            if (f["material_bridge_rate"], f["appointment"])
            == (cohort["rate"], cohort["wait"])
        )
        pad = source["observable_absolute_error_upper_per_pair"]
        moments = next(
            p["receiver_reference_pauli_moments"]
            for p in source["pair_readings"]
            if p["port_face_and_words"] == LABEL
        )
        result = []
        for j, flag in enumerate((1, -1)):
            wins = counts[j]
            n = wins + counts[j + 2]
            assert 0 < wins < n
            lo = float(beta.ppf(alpha, wins, n - wins + 1)) - 0.5
            hi = float(beta.isf(alpha, wins + 1, n - wins)) - 0.5
            accepted_y = (
                moments[0][0]
                + (0 if cohort["reference_dephased"] else flag * moments[0][2])
            ) / 2
            assert accepted_y > pad
            systematic = 2 * pad / (accepted_y - pad)
            result.append(
                {
                    "reference_y": flag,
                    "outcome_count": n,
                    "sample_mean_current": wins / n - 0.5,
                    "simultaneous_99pct_count_interval": [lo, hi],
                    "including_inherited_model_probability_error": [
                        lo - systematic,
                        hi + systematic,
                    ],
                }
            )
        output.append(
            {
                "flight_bridge_rate": cohort["rate"],
                "appointment": cohort["wait"],
                "reference_dephased": cohort["reference_dephased"],
                "raw_counts": counts,
                "current_by_reference": result,
            }
        )
    return output


def run():
    receipt = json.loads(READS.read_text())
    sample = json.loads(COUNTS.read_text())
    assert sample["source_receipt_sha256"] == sha256(READS.read_bytes()).hexdigest()
    assert (
        sample["pre_count_protocol_sha256"]
        == sha256((HERE.parent / "WITNESS.md").read_bytes()).hexdigest()
    )
    case = next(row for row in receipt["cases"] if row["fan_count"] == 4)
    with archived_fixture() as a:
        positive = FanPath(a, 4, 1, asymmetric=True)
        negative = FanPath(a, 4, -1, asymmetric=True)
        for cut in (Fraction(17, 4), Fraction(145, 32)):
            assert verified_inertia(positive.h, cut)["count_below"] == 2
        mask = sum(positive.field.bits.values())
        order = [
            negative.ids[x, tuple(w ^ mask for w in words)]
            for x, words in positive.basis
        ]
        assert (positive.h != negative.h[order][:, order]).nnz == 0
        np.testing.assert_array_equal(negative.source()[order], -positive.source())
        flights = []
        for path in (positive, negative):
            label = (
                LABEL
                if path.sign == 1
                else [*LABEL[:2], LABEL[2] ^ mask, LABEL[3] ^ mask]
            )
            assert label[0] == path.cells[-1][0]
            h_face = local_face(path, label)
            stock_arms, stable = initial_stocks(path)
            amps, cell, hub = pair_amplitudes(path, stable[0], label)
            face = float(np.vdot(amps, amps).real / 2)
            stationary_refusals = [
                1 - stable[1],
                0,
                stable[1] - cell,
                cell - hub,
                hub - face,
            ]
            stationary = head_meters(amps, h_face, stationary_refusals)
            assert (
                max(
                    abs(c["current_given_face_and_reference"])
                    for c in stationary["current"].values()
                )
                < TOL
            )
            preparations = []
            for record in case["stock_arms"]:
                name = record["preparation"]
                np.testing.assert_allclose(
                    stock_arms[name][1:],
                    [
                        record["original_clock_admission"],
                        record["fresh_write_per_original"],
                    ],
                    atol=TOL,
                    rtol=0,
                )
                preparations.append(
                    {
                        "preparation": name,
                        "flights": [
                            run_flight(path, h_face, label, flight, stock_arms[name])
                            for flight in record["flights"]
                            if flight["appointment"] > 0
                        ],
                    }
                )
            flights.append(
                {
                    "source_sign": path.sign,
                    "oriented_label": label,
                    "stationary_real_separately_filtered_lines": stationary,
                    "preparations": preparations,
                }
            )
    bypass = next(
        s["flights"] for s in case["stock_arms"] if s["preparation"] == "clock_bypassed"
    )
    return {
        "original_face_receipt_sha256": sha256(READS.read_bytes()).hexdigest(),
        "original_counts_sha256": sha256(COUNTS.read_bytes()).hexdigest(),
        "joint_record_order": [
            "reference Y+, word a",
            "reference Y+, word b",
            "reference Y-, word a",
            "reference Y-, word b",
        ],
        "refusal_order": [
            "original clock",
            "source write",
            "outside receiver",
            "off hub",
            "different head words",
        ],
        "source_signs": flights,
        "saved_yy_subcohorts": count_intervals(sample, bypass),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as stream:
            stream.write(
                json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
            )
    else:
        compare = run_path(str(HERE.parents[2] / "replay.py"))["compare"]
        compare(result, json.loads((HERE / "results.json").read_text()), "face_current")
        print("Verified the isolated face current, controls and saved-count intervals.")


if __name__ == "__main__":
    main()
