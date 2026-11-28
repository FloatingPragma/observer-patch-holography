"""Measure motion and component energies on the face-0 arrival at hub 1045."""

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from math import asin, atan2, cos, hypot, log, pi, sin, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix

from ...fan_path.fixture import FanPath, archived_fixture
from ...fan_path.local_launch.measure import PAULI
from ...fan_path.measure import conjugacy
from ...field_charge_records.asymmetric_clock import (
    STEP,
    UNDERFLOW_PAD,
    step_error,
    taylor_step,
)
from ...field_charge_records.electric_collar.fresh_arrival.measure import (
    ARRIVAL,
    LABEL,
    source_receipt,
)

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[1]
EVENTS = (
    "null",
    "minus_1",
    "plus_1",
    "minus_2",
    "plus_2",
    "minus_3",
    "plus_3",
    "sharp_voltage",
    "face_phase",
)
TIMES = (1 / 128, 1 / 64, 1 / 8)
SHOTS = 2**30
LOG = log(2 * 8192 / 0.01)
SEED = 1045308125


def interval(k, n, pad=0):
    """Invert bounded-sample Bernstein tails with the observed proportion."""
    assert 0 <= k <= n and n
    p = k / n
    h = LOG / n
    radius = sqrt(2 * p * h) + 4 * h + pad
    return [max(0, p - radius), min(1, p + radius)]


def difference(a, b):
    return [a[0] - b[1], a[1] - b[0]]


def divide(a, b):
    assert b[0] * b[1] > 0
    outcomes = [x / y for x, y in product(a, b)]
    return [min(outcomes), max(outcomes)]


def mass(values):
    return float(np.vdot(values, values).real / 2)


def draw(rng, values):
    assert min(values) >= -2e-12 and abs(sum(values) - 1) <= 2e-12
    probabilities = np.round(np.clip(values, 0, 1), 13)
    return rng.multinomial(SHOTS, probabilities / probabilities.sum()).tolist()


def hit(rng, probability, gate):
    return draw(rng, [probability, gate - probability, 1 - gate])


def probability(counts, pad):
    return interval(counts[0], counts[0] + counts[1], pad)


def phase(x, y, pad):
    measurements = []
    errors = []
    for cells in [x, y]:
        n = sum(cells[:2])
        p = cells[0] / n
        band = interval(cells[0], n, pad)
        measurements.append(2 * p - 1)
        errors.append(2 * max(p - band[0], band[1] - p))
    center = atan2(measurements[1], measurements[0])
    error = asin(hypot(*errors) / hypot(*measurements))
    return [center - error, center + error]


def evolve(h, vector, time):
    """Use dyadic steps no larger than the independently enclosed carrier step."""
    pieces = 2 if time > STEP else 1
    scale = time / (pieces * float(STEP))
    assert scale in (1 / 8, 1 / 4, 1)
    scaled_h = h * scale
    result = vector
    for _ in range(pieces):
        result = taylor_step(scaled_h, result)
    np.testing.assert_allclose(mass(result), mass(vector), atol=1e-13, rtol=0)
    return result


def arrival(path, gauge):
    frame = (path.source() @ gauge).astype(complex)
    assert int(ARRIVAL / STEP) * STEP == ARRIVAL
    _, _, truncation, roundoff = step_error(path.h)
    ticks = int(ARRIVAL / STEP)
    assert ticks * (truncation + roundoff + UNDERFLOW_PAD) < 1 / 2
    for _ in range(ticks):
        frame = taylor_step(path.h, frame)
    gate = np.zeros_like(frame)
    hub = path.blocks[path.size - 1, 0].rows
    mask = sum(path.field.bits.values())
    words = [w ^ mask for w in LABEL[-2:]]
    raw = frame[hub].reshape(-1, len(path.words_at[0]), 2)
    axes = [path.words_at[0].index(w) for w in words]
    detected = raw[:, axes, :].reshape(-1, 4)
    rho = detected.T @ detected.conj() / 2
    moments = [[np.trace(rho @ np.kron(a, b)).real for b in PAULI] for a in PAULI]
    np.testing.assert_allclose(moments, source_receipt(), rtol=0, atol=3e-12)
    body = gate[hub].reshape(raw.shape)
    body[:, axes, :] = raw[:, axes, :]
    assert mass(gate) > 0
    wave_error = float(4 * ticks * (truncation + roundoff + UNDERFLOW_PAD))
    pad = (2 * sqrt(mass(gate)) + wave_error) * wave_error + 1e-10
    return gate, words, float(pad)


def bridge(path):
    right_rows = np.arange(
        path.blocks[path.size - 1, 0].rows.start,
        path.blocks[path.size - 1, 0].rows.stop,
    )
    left_port, right_port = path.bridges[-1]
    assert [left_port, right_port] == [1010, 1045]
    left_rows = np.array([path.ids[left_port, path.basis[j][1]] for j in right_rows])
    assert len(set(left_rows)) == len(left_rows)
    for i, j in zip(left_rows, right_rows, strict=True):
        row = path.bridge_term.getrow(j)
        assert row.nnz == 2 and row[0, j] == 4 and row[0, i] == -4
        assert path.bridge_term[i, j] == -4
    h = coo_matrix(
        (
            np.concatenate([np.full(len(left_rows), s) for s in (4, 4, -4, -4)]),
            (
                np.r_[left_rows, right_rows, left_rows, right_rows],
                np.r_[left_rows, right_rows, right_rows, left_rows],
            ),
        ),
        shape=path.h.shape,
    ).tocsr()
    assert (path.h - path.parked - path.bridge_term).nnz == 0
    # A current across the terminal owner has exactly one boundary term.
    is_terminal = path.owner == path.size - 1
    assert (path.h - h)[is_terminal][:, ~is_terminal].nnz == 0
    assert [path.board.degree[i] for i in (left_port, right_port)] == [6, 6]
    raw = expm(-1j * pi / 16 * np.array([[4, -4], [-4, 4]]))
    np.testing.assert_allclose(
        raw,
        np.exp(-1j * pi / 4) * np.array([[1, 1j], [1j, 1]]) / sqrt(2),
        rtol=0,
        atol=1e-14,
    )
    return left_rows, right_rows, h


def components(path):
    """Recover every addressed material pair and circulating face from H."""
    moves = defaultdict(list)
    rings = defaultdict(list)
    off = path.h.tocoo()
    for i, j, rate in zip(off.row, off.col, off.data, strict=True):
        if i >= j:
            continue
        a, wa = path.basis[i]
        b, wb = path.basis[j]
        if a != b:
            assert rate == -4
            moves[tuple(sorted((a, b)))].append((i, j))
        else:
            assert rate == -0.5
            (c,) = [k for k in range(path.size) if wa[k] != wb[k]]
            face = path.field.masks.index(wa[c] ^ wb[c])
            rings[c, face].append((i, j))
    full_diagonal = np.array(
        [8 * (6 - int(path.board.degree[p])) for p, _ in path.basis], float
    )
    for group, cost in [(moves, 4), (rings, 0.5)]:
        for rows in group.values():
            assert len(np.unique(rows)) == 2 * len(rows)
            np.add.at(full_diagonal, np.asarray(rows).ravel(), cost)
    np.testing.assert_array_equal(full_diagonal, path.h.diagonal())
    return moves, rings


def neutral_check(path, rings, captured, words):
    """Check that coherent erasure leaves exactly the same ring read."""
    block = path.blocks[path.size - 1, 0]
    amplitudes = captured[block.rows].reshape(*block.shape, 2)
    model_ring_sum = 0.0
    for k, sector in enumerate(block.words):
        ring = path.field.ring_laplacian(sector)
        view = np.moveaxis(amplitudes, k, 0).reshape(len(sector), -1)
        from_field = np.vdot(view, ring @ view).real / 2
        model_ring_sum += from_field
        direct = 0.0
        for (c, _), rows in rings.items():
            if c == k:
                i, j = np.asarray(rows).T
                direct += mass(captured[i] - captured[j]) / 2
        np.testing.assert_allclose(direct, from_field, atol=2e-16, rtol=0)
    assert words[0] ^ words[-1] == path.field.masks[0]
    return float(model_ring_sum)


def make_packet(captured, left_rows, right_rows, name, word_mask):
    # The preceding phase pulse is on a sharp position, hence global.
    prep = captured * (1j if name == "sharp_voltage" else 1)
    packet = np.zeros_like(captured)
    packet[left_rows] = prep[right_rows] / sqrt(2)
    packet[right_rows] = prep[right_rows] / sqrt(2)
    if name == "face_phase":
        packet[left_rows[word_mask]] *= 1j
        packet[right_rows[word_mask]] *= 1j
    elif name not in ("null", "sharp_voltage"):
        angle = int(name[-1]) * (1 if name.startswith("plus") else -1) * pi / 4
        packet[right_rows] *= np.exp(1j * angle)
    return packet


def meters(rng, packet, left_rows, right_rows, gate, pad):
    def quarter(phase_shift):
        left, right = packet[left_rows], packet[right_rows] * np.exp(1j * phase_shift)
        # Read gateway side, entry side, or original gate failure.
        return draw(
            rng,
            [
                mass((right + 1j * left) / sqrt(2)),
                mass((left + 1j * right) / sqrt(2)),
                1 - gate,
            ],
        )

    x, y = quarter(pi / 2), quarter(0)
    dy = [2 * p - 1 for p in probability(y, pad)]
    return {
        "x_counts": x,
        "y_counts": y,
        "material_bridge_phase_interval": phase(x, y, pad),
        "signed_quarter_current_pointer_interval": dy,
    }


def energy(rng, packet, moves, rings, gate, pad, law):
    bank = []
    total = [0.0, 0.0]
    for name, keyed, weight in [("matter", moves, 8), ("face", rings, 1)]:
        for loc, rows in keyed.items():
            i, j = np.array(rows).T
            # The anti-symmetric mode of the laplacian gives the minus flag
            # in an independently addressed controlled full exchange.
            odd = mass(packet[i] - packet[j]) / 2
            counts = hit(rng, odd, gate)
            bounds = probability(counts, pad)
            scale = (
                weight / 2 if law == "half_bridge" and loc == (1010, 1045) else weight
            )
            total = [a + scale * b for a, b in zip(total, bounds, strict=True)]
            bank.append(
                {"kind": name, "address": loc, "odd_energy": scale, "counts": counts}
            )
    return {"local_odds": bank, "energy_interval": total}


def calibrate(path, rng):
    ring = path.field.ring_laplacian(path.field.words[path.heads])
    first = (path.vacuum + path.code[:, 0]) / sqrt(2)
    final = expm(-1j * ring) @ first
    z = 2 * np.conj(path.vacuum @ final) * (path.code[:, 0] @ final)
    x = draw(rng, [(1 + z.real) / 2, (1 - z.real) / 2, 0])
    y = draw(rng, [(1 - z.imag) / 2, (1 + z.imag) / 2, 0])
    clock = phase(x, y, 1e-12)
    result = {
        "neutral_ring_x_counts": x,
        "neutral_ring_positive_phase_y_counts": y,
        "neutral_clock_frequency_interval": clock,
    }
    for law, g in (("native", 4), ("half_bridge", 2)):
        # A distinct sharp-port cohort measures the first increasing half turn.
        sub = np.array([[g, -g], [-g, g]], float)
        prob = abs(expm((-1j / 8) * sub)[0, 1]) ** 2
        counts = hit(rng, prob, 1)
        rate = [8 * asin(sqrt(p)) for p in probability(counts, 1e-12)]
        result[law] = {
            "bare_bridge_transfer_counts": counts,
            "bare_bridge_rate_interval": rate,
        }
    return result


def speed_squared(e0, work, momentum):
    assert min(e0) > 0 and min(work) > 0
    a, b = work
    return divide(
        [2 * e0[0] * a + a * a, 2 * e0[1] * b + b * b],
        [momentum[0] ** 2, momentum[1] ** 2],
    )


def motion(
    path,
    captured,
    left_rows,
    right_rows,
    target,
    cal,
    absolute_pad,
    terms,
    neutral,
    rng,
):
    gate = mass(captured)
    pad = 4 * absolute_pad / (gate - absolute_pad)
    moves, rings = terms
    cell = path.owner == path.size - 1
    word_mask = np.array([path.basis[i][1][-1] == target for i in right_rows])
    records = {}
    for law, h in [("native", path.h), ("half_bridge", path.h - bridge(path)[2] / 2)]:
        outcomes = {}
        _, _, truncate, arithmetic = step_error(h)
        assert 2 * (truncate + arithmetic + UNDERFLOW_PAD) < 1e-11
        rate = cal[law]["bare_bridge_rate_interval"]
        clock = cal["neutral_clock_frequency_interval"]
        for name in EVENTS:
            packet = make_packet(captured, left_rows, right_rows, name, word_mask)
            np.testing.assert_allclose(mass(packet), gate, rtol=0, atol=2e-16)
            quadratures = meters(rng, packet, left_rows, right_rows, gate, pad)
            displacements = []
            for t in TIMES:
                flown = evolve(h, packet, t)
                counts = hit(rng, mass(flown[cell]), gate)
                displacements.append(
                    {
                        "dwell": t,
                        "cell_counts": counts,
                        "cell_interval_given_gate": probability(counts, pad),
                    }
                )
            energy_read = energy(rng, packet, moves, rings, gate, pad, law)
            outcomes[name] = {
                "quadratures": quadratures,
                "owner_displacements": displacements,
                "component_energy_read": energy_read,
            }
        original = outcomes["null"]["component_energy_read"]["energy_interval"]
        formation = divide(difference(original, neutral), clock)
        comparisons = []
        for kick in (1, 2, 3):
            plus, minus = [outcomes[f"{label}_{kick}"] for label in ("plus", "minus")]
            yp = plus["quadratures"]["signed_quarter_current_pointer_interval"]
            ym = minus["quadratures"]["signed_quarter_current_pointer_interval"]
            current = [
                bridge_rate * j / 2
                for bridge_rate, j in product(rate, difference(yp, ym))
            ]
            assert min(current) > 0
            impulses = []
            for a, b in zip(
                plus["owner_displacements"], minus["owner_displacements"], strict=True
            ):
                dt = a["dwell"]
                slope = [
                    p / (2 * dt)
                    for p in difference(
                        a["cell_interval_given_gate"], b["cell_interval_given_gate"]
                    )
                ]
                # Taylor derivative error from ||[H,P]||=rate and ||H|| by row sums.
                taylor_error = float(abs(h).sum(axis=1).max()) * rate[1] * dt
                assert slope[0] <= max(current) + taylor_error
                assert slope[1] >= min(current) - taylor_error
                impulses.append(
                    {
                        "wait": dt,
                        "conditional_opposite_kick_slope": slope,
                        "continuity_remainder_ceiling": taylor_error,
                    }
                )
            assert impulses[0]["conditional_opposite_kick_slope"][0] > 0
            energies = [
                arm["component_energy_read"]["energy_interval"] for arm in (plus, minus)
            ]
            symmetric_energy = [(energies[0][j] + energies[1][j]) / 2 for j in range(2)]
            work = divide(difference(symmetric_energy, original), clock)
            angle = plus["quadratures"]["material_bridge_phase_interval"]
            speeds = speed_squared(formation, work, angle)
            lattice_sq = [4 * sin(x / 2) ** 2 for x in angle]
            measured_mass = divide(lattice_sq, [2 * w for w in work])
            expected_mass = divide(clock, rate)
            assert max(measured_mass[0], expected_mass[0]) <= min(
                measured_mass[1], expected_mass[1]
            )
            clocked_rate = divide(rate, clock)
            expected_work = [clocked_rate[j] * (1 - cos(angle[j])) for j in range(2)]
            assert max(expected_work[0], work[0]) <= min(expected_work[1], work[1])
            comparisons.append(
                {
                    "kick_index": kick,
                    "positive_phase_impulse_interval": angle,
                    "instantaneous_opposite_kick_current_in_clock_units": divide(
                        [min(current), max(current)], clock
                    ),
                    "independent_owner_difference": impulses,
                    "energy_increment_in_clock_units": work,
                    "independent_clock_and_phase_expected_work": expected_work,
                    "independent_bridge_expected_local_mass": expected_mass,
                    "local_quadratic_mass_using_lattice_phase": measured_mass,
                    "implied_relativistic_speed_squared_in_one_link_units": speeds,
                }
            )
        s1, s3 = [
            comparisons[k]["implied_relativistic_speed_squared_in_one_link_units"]
            for k in (0, 2)
        ]
        assert max(s1[0], s3[0]) > min(s1[1], s3[1])
        for name in ("null", "sharp_voltage", "face_phase"):
            lo, hi = outcomes[name]["quadratures"][
                "signed_quarter_current_pointer_interval"
            ]
            assert lo < 0 < hi
        records[law] = {
            "formation_mean_in_clock_units": formation,
            "full_setting_reads": outcomes,
            "opposite_kick_comparison": comparisons,
        }
    a, b = [
        records[law]["formation_mean_in_clock_units"]
        for law in ("native", "half_bridge")
    ]
    assert max(a[0], b[0]) < min(a[1], b[1])
    a, b = [
        records[law]["opposite_kick_comparison"][1][
            "instantaneous_opposite_kick_current_in_clock_units"
        ]
        for law in ("native", "half_bridge")
    ]
    assert a[0] > b[1]
    return records


def run():
    rng = np.random.default_rng(SEED)
    with archived_fixture() as archive:
        plus = FanPath(archive, 4, 1, asymmetric=True)
        path = FanPath(archive, 4, -1, asymmetric=True)
        source_gauge, _ = conjugacy(plus, path)
        captured, words, pad = arrival(path, np.array(source_gauge))
        moves, faces = components(path)
        left_rows, right_rows, h = bridge(path)
        packet = make_packet(captured, left_rows, right_rows, "null", None)
        np.testing.assert_array_equal(h @ packet, 0)
        assert all(
            path.board.degree[path.basis[i][0]] == 6
            for i in np.flatnonzero(np.any(packet != 0, axis=1))
        )
        expected_neutral = neutral_check(path, faces, captured, words)
        gate = mass(captured)
        action = path.h @ packet
        level_mean = np.vdot(packet, action).real / (2 * gate)
        variance = mass(action) / gate - level_mean**2
        assert variance > 0
        neutrals = energy(
            rng, captured, {}, faces, gate, 4 * pad / (gate - pad), "neutral"
        )
        assert (
            neutrals["energy_interval"][0]
            <= expected_neutral / gate
            <= neutrals["energy_interval"][1]
        )
        clock = calibrate(path, rng)
        laws = motion(
            path,
            captured,
            left_rows,
            right_rows,
            words[-1],
            clock,
            pad,
            (moves, faces),
            neutrals["energy_interval"],
            rng,
        )
    return {
        "native_fixture_formation_baseline": "filled same head words after compensated charge removal",
        "label_in_positive_chart": LABEL,
        "armed_vacancy_words": words,
        "source_and_receiving_gate_hubs": [98, 1045],
        "arrival_wait": ARRIVAL,
        "kick_link": [1010, 1045],
        "original_gate_probability": gate,
        "original_clock_and_source_write_refusal": 0,
        "original_face_refusal": 1 - gate,
        "counts_order": ["hit_or_plus", "miss_or_minus", "original_capture_refused"],
        "original_launches_per_bank": SHOTS,
        "register_confidence_at_least": 0.99,
        "seed": SEED,
        "probability_error_upper_per_original": pad,
        "source_identity_sha256": sha256(
            (FIELD / "fan_path" / "local_launch" / "witness.json").read_bytes()
        ).hexdigest(),
        "native_balanced_packet_energy_variance_audit": float(variance),
        "neutral_energy_reference": neutrals,
        "independent_time_and_bridge_rate": clock,
        "laws": laws,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    content = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as handle:
            handle.write(content)
    else:
        from ...replay import compare

        compare(
            json.loads(content),
            json.loads((HERE / "results.json").read_text()),
            "fresh_kick",
        )
        print("Reproduced the same-arrival kick, formation, motion and rival counts.")


if __name__ == "__main__":
    main()
