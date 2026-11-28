"""Check native bond-current readings on the lowest three-arm winding lines."""

import argparse
from hashlib import sha256
from itertools import pairwise
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from ..apparatus import LAW, blocks, fixture
from ..certificate import isolate
from ..measure import BASE, CY, ONE
from ..recoil_search.guides import detour
from ..recoil_search.measure import clock
from ..recoil_search.threshold import DELAYS, DT


HERE = Path(__file__).resolve().parent
THRESHOLD = HERE.parent / "recoil_search" / "threshold.json"
TOL = 5e-11
FLUX = np.pi / 8
WAITS = (1, 4, 16)


def verify_local_meter():
    """Check the destination effect of both pulse lengths on every local block."""
    for kind in ("blank", (1, 1)):
        native = blocks()[kind] / 8
        local = 9 * np.arange(len(native) // 9)[:, None]
        source = (local + LAW.SRC.index((0, 1))).ravel()
        dest = (local + LAW.SRC.index((1, 0))).ravel()
        positions = np.r_[source, dest]
        outside = np.ones(len(native), dtype=bool)
        outside[positions] = False
        np.testing.assert_array_equal(native[outside][:, positions], 0)
        h = native[np.ix_(positions, positions)]
        dst = np.diag(np.isin(positions, dest))
        magnitudes = abs(h[: len(source), len(source) :])
        amplitude = magnitudes.max()
        assert amplitude > 0
        assert np.all((magnitudes == 0) | (magnitudes == amplitude))
        assert np.all((magnitudes > 0).sum(axis=0) <= 1)
        assert np.all((magnitudes > 0).sum(axis=1) <= 1)
        for i, j in zip(*np.nonzero(magnitudes), strict=True):
            assert h[i, i] == h[j + len(source), j + len(source)]
        effects = []
        for multiplier in (1, 3):
            pulse = expm(-1j * h * multiplier * np.pi / (4 * amplitude))
            effects.append(pulse.conj().T @ dst @ pulse)
        expected = 1j * (h @ dst - dst @ h) / amplitude
        np.testing.assert_allclose(effects[0] - effects[1], expected, rtol=0, atol=TOL)


def bonds(h, actions, ports, routes, marks):
    """Attach every off-diagonal word transition to its physical seam."""
    index = {port: i for i, port in enumerate(ports)}
    blank = [edge for route in routes for edge in pairwise(route)]
    ordered = list(dict.fromkeys([*blank, *marks]))
    ids = {tuple(sorted(edge)): j for j, edge in enumerate(ordered)}
    assert len(ids) == len(ordered)
    arm_of = {port: j for j, r in enumerate(routes) for port in r[1:-1]}
    hub, tip = routes[0][0], routes[0][-1]
    data = []
    allowed = 0
    for edge in ordered:
        u, v = edge
        arm = arm_of[v if u in (hub, tip) else u]
        others = [i for i in range(3) if i != arm]
        permutation = list(range(3))
        permutation[others[0]], permutation[others[1]] = reversed(others)
        reflection = actions[tuple(permutation)]
        rows_u = 64 * index[u] + np.arange(64)
        rows_v = 64 * index[v] + np.arange(64)
        part = h[rows_v][:, rows_u].tocoo()
        assert (part.nnz, part.data.min(), part.data.max()) == (
            (32, -4, -4) if edge in marks else (64, -0.5, -0.5)
        )
        origin, target = rows_u[part.col], rows_v[part.row]
        key = np.argsort(origin)
        origin, target = origin[key], target[key]
        pairs = set(zip(origin, target, strict=True))
        assert pairs == set(zip(reflection[origin], reflection[target], strict=True))
        if edge in marks:
            np.testing.assert_array_equal(
                (origin % 64) ^ (target % 64), 1 << marks.index(edge)
            )
            assert len(set(origin % 64 & (1 << marks.index(edge)))) == 1
        else:
            np.testing.assert_array_equal(origin % 64, target % 64)
        data.append((edge, origin, target, -float(part.data[0])))
        allowed += 2 * len(origin)
    assert allowed == h.nnz - np.count_nonzero(h.diagonal())
    return data, ids


def verify_phase_classes(ports, routes, marks, data, ids):
    """Reduce arbitrary physical edge phases in integer coefficient arithmetic."""
    index = {port: i for i, port in enumerate(ports)}
    coefficients = np.eye(len(data), dtype=np.int16)
    potential = np.zeros((len(ports), len(data)), dtype=np.int16)
    totals = []
    for route in routes:
        running = np.zeros(len(data), dtype=np.int16)
        for u, v in pairwise(route):
            running += coefficients[ids[tuple(sorted((u, v)))]]
            if v != route[-1]:
                potential[index[v]] = running
        totals.append(running)
    potential[index[routes[0][-1]]] = totals[0]
    memory = np.zeros((len(marks), len(data)), dtype=np.int16)
    residual = {}
    for j, mark in enumerate(marks):
        edge, origin, target, _ = data[ids[tuple(sorted(mark))]]
        assert edge == mark
        u, v = map(index.get, edge)
        delta = ((target[0] >> j) & 1) - ((origin[0] >> j) & 1)
        memory[j] = delta * (
            coefficients[ids[tuple(sorted(mark))]] - potential[v] + potential[u]
        )
    for i, route in enumerate(routes):
        residual[(route[-2], route[-1])] = totals[i] - totals[0]
    zero = np.zeros(len(data), dtype=np.int16)
    for edge, origin, target, _ in data:
        u, v = map(index.get, edge)
        phase = potential[v] - potential[u] + residual.get(edge, zero)
        if edge in marks:
            j = marks.index(edge)
            delta = ((target >> j) & 1) - ((origin >> j) & 1)
            actual = phase + delta[:, None] * memory[j]
            expected = coefficients[ids[tuple(sorted(edge))]]
            np.testing.assert_array_equal(
                actual, np.broadcast_to(expected, actual.shape)
            )
        else:
            np.testing.assert_array_equal(phase, coefficients[ids[tuple(sorted(edge))]])
    independent = [np.flatnonzero(total - totals[0]).tolist() for total in totals[1:]]
    assert independent[0] != independent[1] and all(independent)
    return independent


def head_dependent_control(h, actions, routes, data):
    """Change phases conditioned on the local written head and preserve arm permutations."""
    first_detours = [(r[1], r[2]) for r in routes]
    rows, cols, values = [], [], []
    for edge, source, dest, magnitude in data:
        if edge not in first_detours:
            continue
        bit = first_detours.index(edge)
        phase = FLUX * (1 - 2 * ((source >> bit) & 1))
        delta = -magnitude * np.expm1(1j * phase)
        rows.extend([*dest, *source])
        cols.extend([*source, *dest])
        values.extend([*delta, *delta.conj()])
    changed = h + coo_matrix((values, (rows, cols)), shape=h.shape).tocsr()
    assert (changed - changed.T.conj()).nnz == 0
    assert changed.nnz == h.nnz
    np.testing.assert_allclose((abs(changed) - abs(h)).data, 0, rtol=0, atol=TOL)
    for permutation in actions.values():
        assert (changed[permutation][:, permutation] - changed).nnz == 0
    return changed


def blank_loop_controls(routes):
    """Recombine two prescribed all-blank tours relative to a calibrated route."""
    local = blocks()["blank"] / 8
    rail = [LAW.SRC.index(w) for w in ((0, 1), (1, 0))]
    native = local[np.ix_(rail, rail)] - np.eye(2) * local[4, 4]
    np.testing.assert_array_equal(native, [[8.5, -0.5], [-0.5, 8.5]])
    light = []
    for theta in (0.0, FLUX):
        transfers = []
        for arm, route in enumerate(routes):
            amplitude = 1 + 0j
            for e, _ in enumerate(pairwise(route)):
                phase = theta if arm == 0 and e == 0 else 0.0
                h = native.astype(complex)
                h[1, 0] *= np.exp(1j * phase)
                h[0, 1] = h[1, 0].conjugate()
                amplitude *= expm(-1j * np.pi * h)[1, 0]
            transfers.append(amplitude)
        outcomes = []
        for other in (1, 2):
            phasor = transfers[0] * transfers[other].conjugate()
            np.testing.assert_allclose(phasor, np.exp(1j * theta), rtol=0, atol=TOL)
            angle = np.angle(phasor)
            p = (1 + np.array([np.cos(angle), np.sin(angle)])) / 2
            outcomes.append(np.column_stack((p, 1 - p)).tolist())
        light.append(
            {"arm_zero_phase": theta, "in_phase_and_quadrature_outcomes": outcomes}
        )
    return light


def meter(wave, data, port_indices, clock_probability):
    """Complete port probabilities in each quarter-pulse arm, including refusals."""
    table = []
    largest = 0.0
    port_weights = {
        p: float(np.vdot(wave[indices], wave[indices]).real)
        for p, indices in port_indices.items()
    }
    for (u, v), origin, target, _ in data:
        ps = wave[origin]
        pt = wave[target]
        other_u = np.setdiff1d(port_indices[u], origin, assume_unique=True)
        other_v = np.setdiff1d(port_indices[v], target, assume_unique=True)
        untouched_u = np.vdot(wave[other_u], wave[other_u]).real
        untouched_v = np.vdot(wave[other_v], wave[other_v]).real
        elsewhere = sum(weight for p, weight in port_weights.items() if p not in (u, v))
        reading = []
        for sign in (1, -1):
            pu = untouched_u + np.linalg.norm(ps + sign * 1j * pt) ** 2 / 2
            pv = untouched_v + np.linalg.norm(pt + sign * 1j * ps) ** 2 / 2
            result = np.array([pu, pv, elsewhere])
            assert abs(result.sum() - 1) < TOL and min(result) >= -TOL
            reading.append(result)
        difference = float(reading[0][1] - reading[1][1])
        largest = max(largest, abs(difference))
        full = [
            [*(r * clock_probability / 3), 2 / 3, (1 - clock_probability) / 3]
            for r in reading
        ]
        for row in full:
            np.testing.assert_allclose(sum(row), 1, rtol=0, atol=TOL)
        table.append(
            {
                "edge": [u, v],
                "original_start_laws": full,
                "conditional_current_score": difference,
            }
        )
    return {"maximum_conditional_score": largest, "edges": table}


def certified_clock_weight(h, enclosure, center, target_amplitude):
    """Bound the finite-filter purity without treating the solver eigenvalues as exact."""
    detuning = max(abs(x - center) for x in enclosure["interval"]) + TOL
    gap = enclosure["nearest_other_standard_distance"] - detuning
    source_lower = abs(target_amplitude) - enclosure["eigenvector_error"] - TOL
    assert gap > 0 and source_lower > 0
    # All frequency differences are in the rising half of sin(tick*|E-center|/2).
    max_energy = float(np.max(abs(h).sum(axis=1)))
    assert DT * (max_energy + abs(center)) < np.pi
    # Center the geometric sum at its midpoint and bound the real quadrature.
    angle = DT * (DELAYS - 1) * detuning / 2
    assert angle < np.pi / 2
    target_lower = (source_lower * np.cos(angle)) ** 2
    leakage_upper = min(1.0, (DELAYS * np.sin(DT * gap / 2)) ** -2)
    weight = target_lower / (target_lower + leakage_upper) - TOL
    return {
        "line_weight_lower": weight,
        "unconditional_filter_leakage_upper_within_source": leakage_upper,
        "normalized_current_allowance": 2 * np.sqrt(1 - weight),
    }


def measure(length_index, original):
    _, ports, legs, h, actions, lift, small = detour(length_index)
    stock = fixture()
    routes = [[stock["hub"], *leg, stock["tip"]] for leg in legs]
    data, ids = bonds(h, actions, ports, routes, stock["marked"])
    gauge = verify_phase_classes(ports, routes, stock["marked"], data, ids)
    changed = head_dependent_control(h, actions, routes, data)
    ev, vec = eigh(small, driver="evd")
    line = int(ev.searchsorted(16.0))
    enclosure = isolate(h, actions, lift, vec, ev, line)
    center = original["lowest_and_next_harmonic"][0]["energy"]
    assert abs(ev[line] - center) < 1e-11
    assert len(routes[0]) - 1 == original["length"]
    row = lift @ vec[:, line]
    c = actions[CY]
    partner = (row[c] - row[c[c]]) / np.sqrt(3)
    frame = np.column_stack([row, partner])
    np.testing.assert_allclose(frame.T @ frame, np.eye(2), rtol=0, atol=TOL)
    np.testing.assert_allclose(h @ frame, ev[line] * frame, rtol=0, atol=TOL)
    beta = BASE[:, 0] @ (lift[ONE] @ vec)
    prepared = beta * clock(ev, dict(center=center, tick=DT, number_delays=DELAYS))
    clock_probability = float(np.vdot(prepared, prepared).real)
    purity = float(abs(prepared[line]) ** 2 / clock_probability)
    allowance = 2 * np.sqrt(1 - purity) + TOL
    certified = certified_clock_weight(h, enclosure, center, float(beta[line]))
    assert 0 < certified["line_weight_lower"] <= purity
    initial = (lift @ (vec @ prepared)) / np.sqrt(clock_probability)
    if float(partner[ONE] @ BASE[:, 1]) * beta[line] < 0:
        partner, handedness = -partner, -1
    else:
        handedness = 1
    basis_indices = {p: np.arange(64 * i, 64 * (i + 1)) for i, p in enumerate(ports)}
    ideal = [
        meter(
            (row + sign * 1j * partner) / np.sqrt(2),
            data,
            basis_indices,
            clock_probability,
        )
        for sign in (1, -1)
    ]
    assert max(m["maximum_conditional_score"] for m in ideal) < TOL
    records = []
    for t in (0, *WAITS):
        native0 = lift @ (vec @ (prepared * np.exp(-1j * t * (ev - ev[line]))))
        native0 /= np.sqrt(clock_probability)
        native1 = handedness * (native0[c] - native0[c[c]]) / np.sqrt(3)
        reads = [
            meter(
                (native0 + sign * 1j * native1) / np.sqrt(2),
                data,
                basis_indices,
                clock_probability,
            )
            for sign in (1, -1)
        ]
        assert max(r["maximum_conditional_score"] for r in reads) <= allowance
        records.append({"native_exchange_exposure": t, "circular_readers": reads})
    initial1 = handedness * (initial[c] - initial[c[c]]) / np.sqrt(3)
    launch = np.column_stack(
        [(initial + s * 1j * initial1) / np.sqrt(2) for s in (1, -1)]
    )
    control = []
    for t in WAITS:
        later = expm_multiply(-1j * t * changed, launch)
        reads = [meter(col, data, basis_indices, clock_probability) for col in later.T]
        control.append({"native_exchange_exposure": t, "circular_readers": reads})
    difference = max(
        max(
            abs(a["conditional_current_score"] - b["conditional_current_score"])
            for a, b in zip(altered["edges"], native["edges"], strict=True)
        )
        for family, base in zip(control, records[1:], strict=True)
        for altered, native in zip(
            family["circular_readers"], base["circular_readers"], strict=True
        )
    )
    assert difference > 1e-8
    return {
        "arm_links": len(routes[0]) - 1,
        "physical_edges": len(data),
        "phase_edge_order": [edge for edge, *_ in data],
        "two_blank_phase_supports": gauge,
        "blank_loop_control": blank_loop_controls(routes),
        "line_enclosure": enclosure,
        "conditional_filter_weight": purity,
        "conditional_clock_acceptance": clock_probability,
        "computed_conditional_all_time_score_allowance": allowance,
        "certified_clock_weight": certified,
        "zero_intrinsic_current": ideal,
        "native": records,
        "changed_law": control,
        "maximum_changed_law_difference": difference,
    }


def run():
    verify_local_meter()
    original = json.loads(THRESHOLD.read_text(encoding="utf-8"))
    assert original["clock_tick"] == DT and original["clock_delays"] == DELAYS
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "threshold_receipt_sha256": sha256(THRESHOLD.read_bytes()).hexdigest(),
        "source_owner": fixture()["hub"] // 12,
        "receiving_owner": fixture()["tip"] // 12,
        "outcome_order": [
            "source port",
            "destination port",
            "elsewhere",
            "source refusal",
            "filter refusal",
        ],
        "cases": [measure(i, data) for i, data in enumerate(original["lengths"])],
    }


def compare(actual, expected):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for key in actual:
            compare(actual[key], expected[key])
    elif isinstance(actual, (list, tuple)):
        for a, b in zip(actual, expected, strict=True):
            compare(a, b)
    elif isinstance(actual, float):
        np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=5e-12)
    else:
        assert actual == expected


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", type=Path)
    args = cli.parse_args()
    data = run()
    if args.record:
        with args.record.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    else:
        compare(data, json.loads((HERE / "results.json").read_text(encoding="utf-8")))
    for case in data["cases"]:
        highest = max(
            r["maximum_conditional_score"]
            for bank in case["native"]
            for r in bank["circular_readers"]
        )
        print(
            f"{case['arm_links']} links, {case['physical_edges']} current meters, "
            f"native max {highest:.4g}, changed law {case['maximum_changed_law_difference']:.4g}"
        )


if __name__ == "__main__":
    main()
