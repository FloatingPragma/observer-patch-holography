"""Charge de-excitation and neutral orbital excitation arriving in a filled owner."""

import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    blocks,
)
from ..conjugate_flight.measure import charge_reflection, closed_walk, seam_data
from ..emission_channels.apparatus import (
    frame_columns,
    read_local,
    rotate_frame,
    sector,
)
from ..emission_channels.measure import SHOTS, samples

HERE = Path(__file__).resolve().parent
PORTS = (3, 2, 4, 5, 0, 1, 82, 72, 73, 77)
SOURCE = ((1, 2), (4, 5))
DAUGHTER = ((1, 2), (6, 7))
TIMES = (0, 32, 128, 512, 1024, 4096, 16384, 17739, 32768)
ORDERS = ((0, 1), (0, 2))
K = int(candidate.KAPPA)


def geometry():
    board = federation(0)
    graph = nx.Graph(board.seams)
    assert len(set(PORTS)) == 10
    physical = {tuple(sorted(bond)) for bond in zip(PORTS, PORTS[1:])}
    assert physical <= set(board.seams)
    assert [p // 12 for p in PORTS] == [0] * 6 + [6] * 4
    boundary = sorted(
        e for e in board.seams if set(PORTS).intersection(e) and e not in physical
    )
    shutters = {PORTS[3]: 9, PORTS[8]: 79}
    assert all(board.degree[ref] == 5 and ref not in PORTS for ref in shutters.values())
    shutter_edges = {tuple(sorted(e)) for e in shutters.items()}
    assert shutter_edges <= set(boundary)
    outside = set(range(board.n)) - {p for e in boundary for p in e} - set(PORTS)
    store = sorted(outside)[:3]
    loads = np.ones(board.n, dtype=int)
    loads[store] = [0, 0, 2]

    def transpose(i, j):
        path = nx.shortest_path(graph, i, j)
        for a, b in closed_walk(path):
            loads[a], loads[b] = loads[b], loads[a]

    private = [p for p in PORTS if board.degree[p] == 5]
    for p in private:
        transpose(store[0], p)
        assert loads[p] == 0
        transpose(store[0], p)
    for e in boundary:
        if e not in shutter_edges:
            assert all(loads[p] == 1 for p in e)
            graph.remove_edge(*e)
    for s, p in zip(store, [PORTS[1], PORTS[4], PORTS[5]], strict=True):
        transpose(s, p)
    access = []
    for orb, buffer in [(1, 3), (4, 3), (6, 8)]:
        ref = shutters[PORTS[buffer]]
        local = graph.subgraph([p for p in PORTS if p // 12 == ref // 12] + [ref])
        access.append(nx.shortest_path(local, PORTS[orb], ref))
    for e in shutter_edges:
        assert all(loads[p] == 1 for p in e)
        graph.remove_edge(*e)
    assert set(nx.node_connected_component(graph, PORTS[0])) == set(PORTS)
    assert sorted(physical) == sorted(
        tuple(sorted(e)) for e in graph.subgraph(PORTS).edges
    )
    # The native full swap on blank contacts and silent writes for every load.
    seam_data(K)
    return dict(
        ordered_ports=PORTS,
        bridge=[1, 82],
        silent_seams=boundary,
        private_zero_writes=private,
        phase_accesses=access,
        initial_store=store,
    )


def regional_frame(order):
    """Measure the local generator compression with independent material inventory."""
    vacancies, excesses = (1, int(order == (0, 2)))
    words, h = sector(10, vacancies, excesses, K)
    positions = (1, 6) if excesses == 0 else (4, 6)
    matrices = []
    for p in positions:
        indices = []
        for pair in (order, order[::-1]):
            word = [1] * 10
            word[p : p + 2] = pair
            indices.append(words.index(tuple(word)))
        matrices.append(h[np.ix_(indices, indices)])
    np.testing.assert_array_equal(*matrices)
    return matrices[0], rotate_frame(matrices[0], K, order)


def filled_head_certificate():
    """Exact commutation on an orthogonal decomposition of the entire seam space."""
    filled = candidate.native_law.SRC.index((1, 1))
    levels = {}
    tested = 0
    for sector_name, h in blocks().items():
        heads = len(h) // 9
        b = h[filled::9, filled::9]
        host = np.kron(b, np.eye(9, dtype=int))
        p = np.tile(np.arange(9) == filled, heads)
        np.testing.assert_array_equal(h @ host, host @ h)
        assert not (h - host)[:, p].any()
        # All spectral values, certified by a square-free exact polynomial.
        roots = sorted(set(np.rint(np.linalg.eigvalsh(b) / 8).astype(int).tolist()))
        ident = np.eye(heads, dtype=int)
        poly = ident
        for root in roots:
            poly = poly @ (b - 8 * root * ident)
        assert not poly.any()
        levels[str(sector_name)] = roots
        tested += len(h)
    assert tested == 99
    return dict(columns=tested, filled_seam_levels=levels)


def project_channels(words, frames, target, y):
    reads = [
        (left, right, buffer, order, frame)
        for (left, right), buffer, order, frame in zip(
            target, (3, 8 if target == DAUGHTER else 3), ORDERS, frames, strict=True
        )
    ]
    # Swap rotations do not alter the filled status of either buffer.
    return read_local(y, words, reads)


def exact_transition(words, h, source, target):
    """Integer Taylor coefficients of the amplitude, in the physical parity frames."""
    twice_h = np.rint(2 * h).astype(int)
    np.testing.assert_allclose(2 * h, twice_h)
    a = np.rint(2 * source).astype(int)
    b = np.rint(2 * target).astype(int)
    np.testing.assert_allclose(a / 2, source, atol=1e-13)
    np.testing.assert_allclose(b / 2, target, atol=1e-13)
    # Python integers retain every cancellation in the products.
    values = a.astype(object)
    rows, cols = twice_h.nonzero()
    entries = [(int(i), int(j), int(twice_h[i, j])) for i, j in zip(rows, cols)]
    coefficients = []
    for order in range(len(words)):
        moment = int(b @ values)
        coefficients.append(moment)
        if moment:
            return coefficients
        next_values = np.zeros(len(words), dtype=object)
        for i, j, value in entries:
            next_values[i] += value * values[j]
        values = next_values
    return coefficients


def boundary_inventory(words, h):
    """Only the inter-owner seam contributes to either receiver inventory current."""
    _, cut = sector(10, 2, 1, K, 5)
    connected = h - cut
    counts = np.array([[w[6:].count(load) for w in words] for load in (0, 2)])
    for inventory in counts:
        change = inventory[None, :] - inventory[:, None]
        np.testing.assert_array_equal(h * change, connected * change)
    return counts


def moments(h, x):
    hx = h @ x
    mean = x @ hx
    return [
        round(float(mean), 10),
        round(float(np.linalg.norm(hx - mean * x) ** 2), 10),
    ]


def run_setting(frames, charged_level, rng, flux_rng, *, cut=False, reverse=False):
    words, h = sector(10, 2, 1, K, 5 if cut else None)
    src, end = (DAUGHTER, SOURCE) if reverse else (SOURCE, DAUGHTER)
    source, positions = frame_columns(words, src, frames)
    target, destinations = frame_columns(words, end, frames)
    for slots in (positions, destinations):
        assert all(sum(x != 1 for x in words[i]) == 3 for i in slots)
    initial, final = (1, 2) if reverse else (2 * charged_level, 1)
    s, d = source[:, initial], target[:, final]
    # All accepted regional tags exhaust the inventory. The third label includes
    # every missed detection, independently for each of the two owners.
    eigenvalues, v = eigh(h)
    incoming = v.T @ s
    inventories = boundary_inventory(words, h)
    joint_bins = 2 * inventories[0] + inventories[1]
    rows = []
    for t in TIMES:
        y = v @ (np.exp(-1j * eigenvalues * t) * incoming)
        p = project_channels(words, frames, end, y)
        wanted = divmod(final, 2)
        np.testing.assert_allclose(p[wanted], abs(d @ y) ** 2, atol=1e-12)
        counts = samples(p, rng)
        inventory_p = np.bincount(joint_bins, weights=abs(y) ** 2, minlength=6).reshape(
            3, 2
        )
        inventory_counts = samples(inventory_p, flux_rng)
        rows.append(
            dict(
                exposure=t,
                counts=counts.tolist(),
                receiver_inventory=inventory_counts.tolist(),
                accepted=int(counts[wanted]),
                probability=round(float(p[wanted]), 10),
            )
        )
    if not cut:
        expected = np.kron(regional_frame(ORDERS[0])[0], np.eye(2)) + np.kron(
            np.eye(2), regional_frame(ORDERS[1])[0]
        )
        for ids in (positions, destinations):
            np.testing.assert_array_equal(h[np.ix_(ids, ids)], expected)
    return dict(
        initial_levels=divmod(initial, 2),
        target_levels=divmod(final, 2),
        source_mean_variance=moments(h, s),
        target_mean_variance=moments(h, d),
        dimension=len(h),
        records=rows,
        integer_moments=exact_transition(words, h, s, d)
        if not (reverse or cut) and charged_level == 1
        else None,
    )


def pair_alone(frames):
    words, h = sector(10, 1, 1, K)
    index = {w: i for i, w in enumerate(words)}
    perm = [index[tuple(2 - x for x in w)] for w in words]
    assert charge_reflection()["native_dimension"] == 99
    np.testing.assert_array_equal(h[np.ix_(perm, perm)], h)
    src = frame_columns(words, [SOURCE[1]], [frames[1]])[0][:, 0]
    end = frame_columns(words, [DAUGHTER[1]], [frames[1]])[0][:, 1]
    np.testing.assert_allclose(src[perm], src, atol=1e-13)
    np.testing.assert_allclose(end[perm], -end, atol=1e-13)
    return dict(
        dimension=len(words),
        incoming_reflection=1,
        outgoing_reflection=-1,
        all_time_probability=0,
    )


def run():
    banks = [regional_frame(order) for order in ORDERS]
    frames = [f for _, f in banks]
    rng = np.random.default_rng(113077)
    flux_rng = np.random.default_rng(113078)
    calibration = []
    for i, (h, f) in enumerate(banks):
        records = []
        base = K * (2 if i == 0 else 6) + (0.5 if i == 0 else 1)
        for j in (0, 1):
            pminus = abs(f[0, j] - f[1, j]) ** 2 / 2
            counts = samples([1 - pminus, pminus], rng)
            records.append(
                dict(
                    parity_counts=counts.tolist(),
                    measured_energy=base + counts[1] / SHOTS,
                )
            )
        calibration.append(
            dict(twice_compression=np.rint(h * 2).astype(int).tolist(), records=records)
        )
    return dict(
        kappa=K,
        starts_per_setting=SHOTS,
        geometry=geometry(),
        calibration=calibration,
        raised=run_setting(frames, 1, rng, flux_rng),
        ground=run_setting(frames, 0, rng, flux_rng),
        sealed=run_setting(frames, 1, rng, flux_rng, cut=True),
        reverse=run_setting(frames, 1, rng, flux_rng, reverse=True),
        pair_alone=pair_alone(frames),
        head_only=filled_head_certificate(),
    )


if __name__ == "__main__":
    (HERE / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
