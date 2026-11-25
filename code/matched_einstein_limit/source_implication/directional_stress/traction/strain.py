"""Two contact actuators subject to the same first-arrival strain."""

import argparse
import json

import numpy as np
from scipy.sparse import diags
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from .. import bridge
from .measure import HERE


def observables(ports, route):
    a, b = (ports.index(p) for p in route[1:3])
    size = len(ports)
    d = np.zeros((size, size), dtype=np.int64)
    d[a, a] = d[b, b] = 1
    x = np.zeros_like(d)
    x[a, b] = x[b, a] = 1
    ex = d - x
    out = {"full_seam": ex + 16 * d, "exchange_only": ex, "contact_occupation": 2 * d}
    assert np.array_equal(ex @ ex, 2 * ex)
    # Two times each read. Same eigenspaces: outside, symmetric, antisymmetric.
    return {k: np.kron(effect, np.eye(4, dtype=np.int64)) for k, effect in out.items()}


def clock_check(rng):
    """Independent arrival counts from an isolated sharp endpoint at finite strain."""
    trials = 2**26
    d = np.eye(2)
    ex = np.array([[1, -1], [-1, 1]]) / 2
    rows = {}
    for ratio in (1.0, 0.75):
        labels = {}
        for name, h in (
            ("full_seam", ratio * (ex + 8 * d)),
            ("exchange_only", ratio * ex + 8 * d),
        ):
            energy, basis = np.linalg.eigh(h)
            p = []
            for t in (1, 2, 4):
                u = (basis * np.exp(-1j * t * energy)) @ basis.T
                p.append(float(abs(u[1, 0]) ** 2))
            np.testing.assert_allclose(
                p, np.sin(ratio * np.array([1, 2, 4]) / 2) ** 2, atol=2e-14
            )
            counts = rng.binomial(trials, p)
            width = np.sqrt(np.log(2 * 12 / 0.01) / (2 * trials))
            assert np.max(abs(counts / trials - p)) < width
            labels[name] = dict(hit_counts=counts.tolist(), population=p)
        rows[str(ratio)] = labels
    return dict(
        times=[1, 2, 4],
        trials_per_setting=trials,
        coverage=0.99,
        simultaneous_halfwidth=width,
        scale_settings=rows,
    )


def material_standard(rng):
    """Measure the dwell phase relative to an outside rail and the exchange gap independently."""
    shots = 2**26
    width = np.sqrt(np.log(2 * 8 / 0.01) / (2 * shots))
    records = {}
    for s in (1.0, 0.75):
        for name, dwell in (("full_seam", 8 * s), ("exchange_only", 8.0)):
            intervals = {}
            counts = {}
            for label, freq, t in (
                ("dwell", dwell, np.pi / 16),
                ("exchange", s, np.pi / 2),
            ):
                # The two clocks start in balanced arms of the corresponding
                # energy eigenspaces and finish in their original superposition.
                p = (1 + np.cos(freq * t)) / 2
                sample = int(rng.binomial(shots, p))
                est = sample / shots
                interval = (
                    np.arccos(
                        np.clip(2 * np.array([est + width, est - width]) - 1, -1, 1)
                    )
                    / t
                )
                assert interval[0] < freq < interval[1]
                counts[label] = dict(plus=sample, exposure=t)
                intervals[label] = interval.tolist()
            a, b = intervals["dwell"], intervals["exchange"]
            index = [a[0] / b[1], a[1] / b[0]]
            records[f"{name}_{s}"] = dict(
                settings=counts, frequency_intervals=intervals, material_index=index
            )
    return dict(trials_per_setting=shots, coverage=0.99, records=records)


def fronts(board, route):
    """Check all destination ports at the first three graph orders for the four launches."""
    start = bridge.beam.indices(route)
    dist = shortest_path(
        board.adjacency(), indices=route, unweighted=True, directed=False
    ).T
    num = 0
    b = bridge.bridges(board, [route[1:3]])
    diag = np.zeros(board.n, dtype=np.int64)
    diag[list(route[1:3])] = 16
    d2 = diags(np.repeat(diag, 4), dtype=np.int64, format="csr")
    for q in (-1, 1):
        for odd in (False, True):
            h = bridge.native.generator(board, route, q, odd)
            # 8 times each strained full-board H at s=3/4.
            full = 4 * h - b
            ex = full + d2
            w = np.zeros((h.shape[0], len(start)), dtype=np.int64)
            w[start] = np.eye(len(start), dtype=np.int64)
            z = w.copy()
            assert (
                max(abs(full).sum(axis=1).max(), abs(ex).sum(axis=1).max()) ** 3 < 2**63
            )
            for order in (1, 2, 3):
                w = full @ w
                z = ex @ z
                mask = np.repeat(np.repeat(dist == order, 4, axis=0), 4, axis=1)
                np.testing.assert_array_equal(w[mask], z[mask])
                num += int(np.count_nonzero(mask))
    return dict(
        replayed_graph_orders=[1, 2, 3],
        launch_ports=route,
        compared_integer_entries=num,
    )


def neighboring_clocks(menus, state, rng):
    """Disjoint same-appointment clocks with every outside arrival counted."""
    population = np.array(
        [bridge.evaluate(panel, {"ref": state})["ref"] for panel in menus]
    )
    assert population.min() > 0 and population.sum() < 1
    shots = 2**26
    time = np.pi / 2
    results = {}
    for name, rates in (
        ("bridge_stretch", (1.0, 0.75, 1.0)),
        ("common_lapse", (0.75, 0.75, 0.75)),
    ):
        enrolled = rng.multinomial(shots, [*population, 1 - population.sum()])
        crossed = rng.binomial(enrolled[:3], np.sin(time * np.array(rates) / 2) ** 2)
        intervals = []
        for n, hits, freq in zip(enrolled[:3], crossed, rates):
            # Conditional on each observed nonzero enrollment count. A union
            # bound covers 3 dimers times 2 protocols times 2 carrier signs.
            width = np.sqrt(np.log(2 * 12 / 0.01) / (2 * n))
            p = hits / n
            lo, hi = (
                np.arccos(np.clip(1 - 2 * np.array([p - width, p + width]), -1, 1))
                / time
            )
            assert lo < freq < hi
            intervals.append([float(lo), float(hi)])
        results[name] = dict(
            enrolled=list(map(int, enrolled[:3])),
            outside=int(enrolled[3]),
            crossed=list(map(int, crossed)),
            frequency_intervals=intervals,
        )
    return dict(
        dimer_order=["left", "bridge", "right"],
        endpoint_flip_exposure=time,
        initial_occupation=population.tolist(),
        shots_per_protocol=shots,
        conditional_frequency_coverage=0.99,
        settings=results,
    )


def run():
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    ports = list(route[1:3])
    pos = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    effects = observables(ports, route)
    nominal = bridge.tangents()
    prepared = bridge.action.preparations()
    rng = np.random.default_rng(4192638)
    trials = 2**56
    rows = {}
    all_angles = json.loads((bridge.HERE.parent / "signed_transfers.json").read_text())[
        "tests"
    ]
    angles = {
        ("ordinary" if row["source"] == "ordinary" else "pairing"): row[
            "angle_per_unit"
        ]
        for row in all_angles
    }
    clock_rng = np.random.default_rng(34192638)
    dimers = (rings[0][1:3], route[1:3], rings[1][1:3])
    assert len(set(p for edge in dimers for p in edge)) == 6
    assert all(
        j in board.neighbours[i] and board.degree[i] == board.degree[j]
        for i, j in dimers
    )
    for q in (-1, 1):
        clock_menus = [[] for _ in dimers]
        meters = {name: [] for name in effects}
        for odd in (False, True):
            h = bridge.native.generator(board, route, q, odd)
            start = np.zeros((h.shape[0], len(pos)))
            start[pos] = np.eye(len(pos))
            ht = -2j * h
            full_wave = expm_multiply(ht, start, traceA=ht.diagonal().sum())
            wave = full_wave[bridge.beam.indices(ports)]
            for menu, edge in zip(clock_menus, dimers):
                segment = full_wave[bridge.beam.indices(edge)]
                menu.append(
                    np.mean(
                        [arm.conj().T @ arm for arm in np.hsplit(segment, 10)], axis=0
                    )
                )
            for name, h2 in effects.items():
                vals = [arm.conj().T @ (h2 / 2) @ arm for arm in np.hsplit(wave, 10)]
                meters[name].append(np.mean(vals, axis=0))
        values = {name: bridge.evaluate(head, nominal) for name, head in meters.items()}
        pminus = bridge.evaluate(meters["exchange_only"], prepared)
        pin = bridge.evaluate(meters["contact_occupation"], prepared)
        counts = {}
        for name in prepared:
            prob = [1 - pin[name], pin[name] - pminus[name], pminus[name]]
            counts[name] = rng.multinomial(trials, prob).tolist()
        reads = {}
        for name, scores in (
            ("full_seam", [0, 8, 9]),
            ("exchange_only", [0, 0, 1]),
            ("contact_occupation", [0, 1, 1]),
        ):
            source = {k: np.asarray(row) @ scores / trials for k, row in counts.items()}
            margin = max(scores) * np.sqrt(np.log(2 * 3 * 2 * 5 / 0.01) / (2 * trials))
            reads[name] = {}
            for k, t in angles.items():
                d = (
                    source[k]
                    - source["equilibrium"]
                    + t * (source["plus"] - source["minus"]) / 2
                )
                w = margin * (2 + abs(t))
                assert abs(d - values[name][k]) < w
                reads[name][k] = dict(estimate=d, simultaneous_halfwidth=w)
        rows[str(q)] = dict(
            means=values,
            outcome_counts=counts,
            count_slopes=reads,
            neighbor_clocks=neighboring_clocks(
                clock_menus, prepared["equilibrium"], clock_rng
            ),
        )
    return dict(
        time=4,
        radius_reference="isolated blank one-hop front",
        neighbor_clock_seed=34192638,
        disjoint_dimers=dimers,
        trials_per_preparation=trials,
        seed=4192638,
        coverage=0.99,
        outcome_order=["outside", "contact_symmetric", "contact_antisymmetric"],
        charges=rows,
        isolated_strain_clock=clock_check(rng),
        independent_material_standard=material_standard(rng),
        leading_flight_comparison=fronts(board, route),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "strain.json").write_text(data)
    else:
        print(data, end="")
