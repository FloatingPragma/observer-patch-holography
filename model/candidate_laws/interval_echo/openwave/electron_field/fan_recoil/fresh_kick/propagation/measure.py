"""Independent spatial, temporal and arrival counts on the commissioned route."""

import argparse
from hashlib import sha256
import json
from math import asin, atan2, cos, factorial, hypot, pi, sin, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags

from ..measure import LOG, arrival, divide, interval, mass, phase
from ....fan_path.fixture import FanPath, archived_fixture
from ....fan_path.measure import conjugacy
from ....field_charge_records.asymmetric_clock import (
    STEP,
    UNDERFLOW_PAD,
    step_error,
    taylor_step,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**42
CELLS = 48
RATE = 4
COST = 8
TIMES = (1 / 64, 2 / 64, 3 / 64)
COEFFICIENTS = (1.5, -0.3, 1 / 30)
WAITS = (0, 4)
SEED = 104510521047
SETTINGS = (
    ("rest", 0),
    ("plus_low", pi / 3),
    ("minus_low", -pi / 3),
    ("plus_held", 5 * pi / 6),
    ("minus_held", -5 * pi / 6),
)


def counts(rng, p, gate=1):
    values = np.r_[np.asarray(p) * gate, 1 - gate]
    assert min(values) >= -1e-12 and abs(sum(values) - 1) <= 2e-12
    values = np.round(np.clip(values, 0, 1), 13)
    return rng.multinomial(SHOTS, values / values.sum()).tolist()


def signed(bank, pad):
    bounds = interval(bank[0], sum(bank[:-1]), pad)
    return [2 * p - 1 for p in bounds]


def symbol(k, r, v):
    b = np.array(
        [
            [2 * r, -r, -r * np.exp(-1j * k)],
            [-r, 2 * r + v, -r],
            [-r * np.exp(1j * k), -r, 2 * r],
        ]
    )
    eigenvalues, eigenvectors = np.linalg.eigh(b)
    u = eigenvectors[:, 0]
    return (
        float(eigenvalues[0]),
        u,
        float(2 * r * np.imag(u[2].conj() * u[0] * np.exp(1j * k))),
    )


def trace_route(path, captured, words):
    board = path.board
    mates = {v: u for a, b in board.external for u, v in ((a, b), (b, a))}
    first = (1045, 1052, 1047)
    forbidden = {c[0] // 12 for c in path.cells}

    def extend(cells, seen):
        if len(cells) == CELLS:
            return cells
        port = mates[cells[-1][-1]]
        owner = port // 12
        if owner in seen:
            return None
        for m in sorted(board.neighbours[port]):
            if m // 12 != owner or board.degree[m] != 5:
                continue
            for q in sorted(board.neighbours[m]):
                if q == port or q // 12 != owner or board.degree[q] != 6:
                    continue
                if q not in mates or mates[q] // 12 in (seen | {owner}):
                    continue
                out = extend(cells + [(port, m, q)], seen | {owner})
                if out is not None:
                    return out
        return None

    cells = extend([first], forbidden)
    ports = [int(p) for cell in cells for p in cell]
    degrees = np.asarray(board.degree)[ports]
    np.testing.assert_array_equal(degrees, np.tile([6, 5, 6], CELLS))
    assert len(set(p // 12 for p in ports[::3])) == CELLS
    assert all(
        b in board.neighbours[a] for a, b in zip(ports[:-1], ports[1:], strict=True)
    )
    first_bit = path.field.bits[(0, path.cells[-1].index(1052))]
    assert (
        int(ports[0]) == path.cells[-1][0]
        and words[0] ^ words[1] == path.field.masks[0]
    )
    old = {
        frozenset((cell[a], cell[b]))
        for cell in path.cells
        for a, b in path.field.edges
    } | set(map(frozenset, path.bridges))
    route_edges = list(zip(ports[:-1], ports[1:], strict=True))
    assert old.intersection(map(frozenset, route_edges)) == {frozenset(route_edges[0])}
    support = np.flatnonzero(np.any(captured != 0, axis=1))
    assert {path.basis[i][0] for i in support} == {1045}
    assert {path.basis[i][1][-1] for i in support} == set(words)
    assert all(path.basis[i][1][-1] & first_bit == 0 for i in support)
    # Reverse the prefix as the vacancy moves. Every hop is toward its head.
    baseline = np.bincount(np.arange(1, len(ports)), minlength=len(ports))
    masks = []
    for j in range(len(ports)):
        heads = [i if i < j else i + 1 for i in range(len(ports) - 1)]
        got = np.bincount(heads, minlength=len(ports)) - baseline
        want = (np.arange(len(ports)) == 0).astype(int) - (np.arange(len(ports)) == j)
        np.testing.assert_array_equal(got, want)
        assert all(heads[i] != j for i in (j - 1, j) if 0 <= i < len(heads))
        masks.append((1 << j) - 1)
    assert all((masks[j] ^ masks[j + 1]) == 1 << j for j in range(len(masks) - 1))
    return ports, degrees, support.size


def generator(degrees, removed=False):
    links = np.full(len(degrees) - 1, RATE)
    if removed:
        links[2::3] = 0
    onsite = 8 * (6 - degrees)
    diagonal = onsite + np.r_[0, links] + np.r_[links, 0]
    return diags([-links, diagonal, -links], [-1, 0, 1], format="csr", dtype=float)


def packet(k):
    _, u, _ = symbol(k, RATE, COST)
    env = sqrt(2 / 32) * np.sin(np.arange(1, 32) * pi / 32)
    prepared = np.zeros((CELLS, 3), complex)
    prepared[8:39] = env[:, None] * np.exp(1j * k * np.arange(8, 39))[:, None] * u
    return prepare(prepared.ravel())


def prepare(target):
    stock = np.zeros_like(target)
    stock[0] = 1
    for j in range(len(stock) - 1):
        remaining = np.linalg.norm(target[j:])
        if not remaining:
            break
        theta = atan2(np.linalg.norm(target[j + 1 :]), abs(target[j]))
        # A native adjacent exchange for theta/4, then two addressed phases.
        a, b = stock[j : j + 2]
        pulse = np.exp(-1j * theta) * np.array(
            [
                cos(theta) * a + 1j * np.sin(theta) * b,
                1j * np.sin(theta) * a + cos(theta) * b,
            ]
        )
        stock[j : j + 2] = pulse * np.exp(
            1j * np.array([theta + np.angle(target[j]), theta - pi / 2])
        )
    stock[-1] *= np.exp(1j * np.angle(target[-1]))
    np.testing.assert_allclose(stock, target, rtol=0, atol=2e-12)
    return stock


def translation(rng, u, gate, pad):
    data = {}
    for axis in ("x", "y"):
        banks = []
        for parity in (0, 1):
            a = (3 * np.arange(parity, CELLS - 1, 2)[:, None] + np.arange(3)).ravel()
            b = a + 3
            assert len(set(a) | set(b)) == 2 * len(a)
            left, right = u[a], u[b] * (1j if axis == "x" else 1)
            out_right = (right + 1j * left) / sqrt(2)
            out_left = (left + 1j * right) / sqrt(2)
            pr, pl = (float(np.vdot(x, x).real) for x in (out_right, out_left))
            banks.append(counts(rng, [pr, pl, 1 - pr - pl], gate))
        bands = [difference_of_pointer(bank, pad) for bank in banks]
        data[axis] = {
            "counts": banks,
            "quadrature_interval": [sum(t) / 2 for t in zip(*bands)],
        }
    xy = [data[q]["quadrature_interval"] for q in ("x", "y")]
    x, y = (sum(q) / 2 for q in xy)
    center = atan2(y, x)
    radius = hypot(*(q[1] - q[0] for q in xy)) / 2
    error = asin(radius / hypot(x, y))
    data["wave_number_interval"] = [center - error, center + error]
    data["overlap_modulus_interval"] = [
        hypot(x, y) - radius,
        min(1, hypot(x, y) + radius),
    ]
    return data


def difference_of_pointer(bank, pad):
    a, b, other, _ = bank
    n = a + b + other
    x, y = interval(a, n, pad), interval(b, n, pad)
    return [x[0] - y[1], x[1] - y[0]]


def evolve(h, u, t):
    steps = max(1, int(t / float(STEP)))
    scale = t / (steps * float(STEP))
    assert 0 <= scale <= 1
    scaled = h * scale
    _, _, truncation, roundoff = step_error(scaled)
    assert steps * float(truncation + roundoff + UNDERFLOW_PAD) < 1e-10
    state = u
    for _ in range(steps):
        state = taylor_step(scaled, state)
    return state


def read_phase(rng, z, gate=1):
    return [counts(rng, [(1 + a) / 2, (1 - a) / 2], gate) for a in (z.real, -z.imag)]


def calibration(path, rng):
    clock_ring = path.field.ring_laplacian(path.field.words[path.heads])
    amp = expm(-1j * clock_ring) @ ((path.vacuum + path.code[:, 0]) / sqrt(2))
    c = 2 * np.conj(path.vacuum @ amp) * (path.code[:, 0] @ amp)
    clock_counts = read_phase(rng, c)
    clock_scale = phase(*clock_counts, 1e-10)
    bridge = np.array([[4, -4], [-4, 4]])
    p = abs(expm(-1j * bridge / 8)[1, 0]) ** 2
    hop = counts(rng, [p, 1 - p])
    rate_interval = divide(
        [8 * asin(sqrt(q)) for q in interval(hop[0], sum(hop[:-1]), 1e-10)], clock_scale
    )
    stored = expm(-1j * np.diag([0, COST]) / 16) @ (np.ones(2) / sqrt(2))
    onsite = 2 * stored[0].conj() * stored[1]
    onsite_counts = read_phase(rng, onsite)
    onsite_interval = divide(
        [16 * q for q in phase(*onsite_counts, 1e-10)], clock_scale
    )
    return {
        "neutral_ring_counts": clock_counts,
        "clock_scale_interval": clock_scale,
        "sharp_bridge_counts": hop,
        "transfer_rate_interval": rate_interval,
        "parked_degree_phase_counts": onsite_counts,
        "degree_cost_interval": onsite_interval,
    }


def formation(rng, h, u, gate, pad, clock):
    norms = float(np.max(np.asarray(abs(h).sum(axis=1))))
    assert norms <= 24
    error = (
        norms**7
        * TIMES[0] ** 6
        * sum(abs(w) * j**7 for j, w in enumerate(COEFFICIENTS, 1))
        / factorial(7)
    )
    bank = [read_phase(rng, np.vdot(u, evolve(h, u, t)), gate) for t in TIMES]
    terms = [
        sorted(w * q / TIMES[0] for q in signed(rows[1], pad))
        for w, rows in zip(COEFFICIENTS, bank, strict=True)
    ]
    lo, hi = [sum(t) for t in zip(*terms)]
    spread = float(np.linalg.norm(h @ u - np.vdot(u, h @ u) * u))
    blank = read_phase(rng, 1 + 0j, gate)
    return {
        "ramsey_counts": bank,
        "parked_filled_reference_counts": blank,
        "energy_interval": divide([lo - error, hi + error], clock),
        "ideal_energy_standard_deviation": spread,
        "derivative_truncation_bound": error,
    }


def stationary(rng, h, gate, pad, clock):
    # The equal exterior slots make the end-contact Laplacian vanish on the trial.
    anchor = 10 - 2 * sqrt(17)
    trial = np.tile([1, 1 - anchor / 4, 1], CELLS)
    u = prepare((trial / np.linalg.norm(trial)).astype(complex))
    np.testing.assert_allclose(h @ u, anchor * u, rtol=0, atol=2e-13)
    waits = (1 / 4, 1 / 2)
    ramsey = [read_phase(rng, np.vdot(u, evolve(h, u, t)), gate) for t in waits]
    bands = [
        divide([q / t for q in phase(*bank, pad)], clock)
        for t, bank in zip(waits, ramsey, strict=True)
    ]
    line = [max(b[0] for b in bands), min(b[1] for b in bands)]
    assert line[0] < anchor < line[1]
    return {
        "positive_formation_line_interval": line,
        "ramsey_waits": list(waits),
        "ramsey_counts": ramsey,
        "translation": translation(rng, u, gate, pad),
        "motion": owner_motion(rng, h, u, gate, pad, clock),
        "ideal_sharp_capture_weight_in_line": float(abs(u[0]) ** 2),
    }


def owner_motion(rng, h, u, gate, pad, clock):
    output = []
    for t in WAITS:
        evolved = evolve(h, u, t)
        p = np.sum(abs(evolved.reshape(CELLS, 3)) ** 2, axis=1)
        bank = counts(rng, p, gate)
        n = sum(bank[:-1])
        center = np.arange(CELLS) @ np.asarray(bank[:-1]) / n
        error = (CELLS - 1) * (sqrt(LOG / (2 * n)) + pad)
        output.append(
            {
                "wait": t,
                "owner_counts_then_gate_refusal": bank,
                "centroid_interval": [center - error, center + error],
            }
        )
    a, b = (d["centroid_interval"] for d in output)
    return {
        "shots": output,
        "velocity_interval": divide(
            [b[0] - a[1], b[1] - a[0]], [WAITS[-1] * q for q in clock]
        ),
    }


def compare(cal, row):
    k, ek = (
        row["translation"]["wave_number_interval"],
        row["formation"]["energy_interval"],
    )
    ri, vi = cal["transfer_rate_interval"], cal["degree_cost_interval"]
    center = sum(k) / 2
    r = sum(ri) / 2
    omega, v = symbol(center, r, sum(vi) / 2)[::2]
    err = 2 * (ri[1] - ri[0]) + (vi[1] - vi[0]) / 2 + ri[1] * (k[1] - k[0]) / 2
    fringe = ri[1] * (1 - row["translation"]["overlap_modulus_interval"][0])
    expected = [omega - err - fringe, omega + err + fringe]
    ev = row["motion"]["velocity_interval"]
    assert max(ek[0], expected[0]) <= min(ek[1], expected[1])
    assert max(ev[0], v - 0.06) <= min(ev[1], v + 0.06)
    assert (
        row["contact_removed"]["velocity_interval"][0]
        < 0
        < row["contact_removed"]["velocity_interval"][1]
    )
    return {
        "held_energy_interval_from_calibrated_symbol": expected,
        "bulk_group_velocity_center": v,
        "finite_packet_velocity_allowance": 0.06,
    }


def run():
    rng = np.random.default_rng(SEED)
    with archived_fixture() as archive:
        path = FanPath(archive, 4, -1, True)
        positive = FanPath(archive, 4, 1, True)
    captured, words, absolute_pad = arrival(path, conjugacy(positive, path)[0])
    gate = mass(captured)
    pad = 4 * absolute_pad / (gate - absolute_pad) + 1e-9
    ports, degrees, n = trace_route(path, captured, words)
    h = generator(degrees)
    cut = generator(degrees, True)
    calibration_data = calibration(path, rng)
    timeunit = calibration_data["clock_scale_interval"]
    cut_entries = cut.tocoo()
    assert cut_entries.nnz and np.all(cut_entries.row // 3 == cut_entries.col // 3)
    records = {}
    for name, k in SETTINGS:
        u = packet(k)
        row = {
            "translation": translation(rng, u, gate, pad),
            "formation": formation(rng, h, u, gate, pad, timeunit),
            "motion": owner_motion(rng, h, u, gate, pad, timeunit),
            "contact_removed": owner_motion(rng, cut, u, gate, pad, timeunit),
        }
        row["comparison"] = compare(calibration_data, row)
        records[name] = row
    line = stationary(rng, h, gate, pad, timeunit)
    e0 = line["positive_formation_line_interval"]
    for name in ("plus_low", "plus_held"):
        row = records[name]
        e = row["formation"]["energy_interval"]
        k = row["translation"]["wave_number_interval"]
        vel = row["motion"]["velocity_interval"]
        e2diff = [e[0] ** 2 - e0[1] ** 2, e[1] ** 2 - e0[0] ** 2]
        ce = divide(e2diff, [q * q for q in k])
        cv = divide([e[0] * vel[0], e[1] * vel[1]], k)
        row["canonical_energy_squared_speed_interval"] = ce
        row["canonical_drift_squared_speed_interval"] = cv
        row["sine_lattice_energy_squared_speed_interval"] = divide(
            e2diff, [(2 * sin(q / 2)) ** 2 for q in k]
        )
    a, b = (
        records[k]["canonical_energy_squared_speed_interval"]
        for k in ("plus_low", "plus_held")
    )
    assert a[1] < b[0]
    held = records["plus_held"]
    assert held["canonical_drift_squared_speed_interval"][1] < b[0]
    a_lattice, b_lattice = (
        records[k]["sine_lattice_energy_squared_speed_interval"]
        for k in ("plus_low", "plus_held")
    )
    assert a_lattice[1] < b_lattice[0]
    k = held["translation"]["wave_number_interval"]
    pad_e = calibration_data["transfer_rate_interval"][1] * (
        1 - held["translation"]["overlap_modulus_interval"][0]
    )
    energy = held["formation"]["energy_interval"]
    velocity = held["motion"]["velocity_interval"]
    pad_v = held["comparison"]["finite_packet_velocity_allowance"]
    eb = [energy[0] - pad_e, energy[1] + pad_e]
    vb = [velocity[0] - pad_v, velocity[1] + pad_v]
    wide_energy = divide(
        [eb[0] ** 2 - e0[1] ** 2, eb[1] ** 2 - e0[0] ** 2], [q * q for q in k]
    )
    wide_drift = divide([eb[0] * vb[0], eb[1] * vb[1]], k)
    assert wide_drift[1] < wide_energy[0]
    held["widened_energy_squared_speed_interval"] = wide_energy
    held["widened_drift_squared_speed_interval"] = wide_drift
    for level in ("low", "held"):
        assert records["plus_" + level]["motion"]["velocity_interval"][0] > 0
        assert records["minus_" + level]["motion"]["velocity_interval"][1] < 0
    witness = HERE.parents[2] / "fan_path" / "local_launch" / "witness.json"
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "original_arrival_witness_sha256": sha256(witness.read_bytes()).hexdigest(),
        "original_start_count_per_bank": SHOTS,
        "count_rng_seed": SEED,
        "conditional_probability_absolute_allowance": pad,
        "gate_probability_interval": [gate - absolute_pad, gate + absolute_pad],
        "gate_words": words,
        "ideal_clock_refusals": 0,
        "ideal_source_write_refusals": 0,
        "spectator_coordinate_support": n,
        "ports": ports,
        "calibration": calibration_data,
        "runs": records,
        "stationary_compensated_anchor": line,
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    records = run()
    receipt = HERE / "results.json"
    if args.write:
        receipt.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
    else:
        saved = json.loads(receipt.read_text())

        def compare_record(a, b):
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for k in a:
                    compare_record(a[k], b[k])
            elif isinstance(a, list):
                assert len(a) == len(b)
                for x, y in zip(a, b, strict=True):
                    compare_record(x, y)
            elif isinstance(a, float):
                assert abs(a - b) <= 1e-10, (a, b)
            else:
                assert a == b, (a, b)

        compare_record(records, saved)
        print("Verified 1045 propagation counts, holds and controls.")


if __name__ == "__main__":
    main()
