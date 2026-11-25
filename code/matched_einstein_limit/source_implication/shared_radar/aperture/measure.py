"""Equal-hop phase-sensitive thresholds on two blank coherent apertures."""

import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from ...common_front import measure as front
from ...full_algebra import station
from ...material_clocks.measure import blank_clock
from ... import signed_transfers
from .. import measure as survey

HERE = Path(__file__).resolve().parent
APERTURES = ((0, 5, 39), (39, 40, 0))
NORM = 64 * 36


def hopping(board, route, q, odd, strength):
    if strength == 8:
        return front.hop(board, route, q, odd)
    h = (2 * station.iso.ranging.generator(board, route, q, strength, odd)).tocoo()
    keep = (h.row // 4 != h.col // 4) & (h.data != 0)
    integer = (-h.data[keep]).astype(np.int64)
    np.testing.assert_array_equal(integer, -h.data[keep])
    return sparse.csr_matrix((integer, (h.row[keep], h.col[keep])), shape=h.shape)


def operators(board, route, strength=8):
    """Raw integer numerator, with 4 times the probability normalization."""
    graph = shortest_path(board.adjacency(), unweighted=True, directed=False)
    edges = {tuple(sorted(e)) for e in board.seams}
    for q in (-1, 1):
        gate = expm(-1j * np.pi / 2 * blank_clock(q, 8))
        np.testing.assert_allclose(
            gate / gate[0, 0], np.array([[1, 1j], [1j, 1]]), atol=4e-14, rtol=0
        )
    for x, y, z in APERTURES:
        assert graph[x, z] == graph[y, z] == 3
        assert tuple(sorted((x, y))) in edges - set((route[:2], route[2:]))
    terms = np.empty((2, 2, 4, 4), complex)
    for odd in (0, 1):
        hops = [hopping(board, route, q, odd, strength) for q in (-1, 1)]
        assert all(int(abs(h).sum(axis=1).max()) ** 3 < 2**62 for h in hops)
        powers = [h**3 for h in hops]
        for k, (x, y, z) in enumerate(APERTURES):
            pair = []
            for h in powers:
                u, v = [
                    h[4 * z : 4 * z + 4, 4 * w : 4 * w + 4].toarray() for w in (x, y)
                ]
                back = [
                    h[4 * w : 4 * w + 4, 4 * z : 4 * z + 4].toarray() for w in (x, y)
                ]
                pair.append((u, v))
                np.testing.assert_array_equal(back, [u.T, v.T])
            np.testing.assert_array_equal(pair[0], np.transpose(pair[1], (0, 2, 1)))
            u, v = pair[1]
            # All previously calculated matrices are exact integers with entries <= 72.
            assert max(np.max(abs(u)), np.max(abs(v))) <= 72
            terms[k, odd] = (u.T @ u + v.T @ v + u @ u.T + v @ v.T) + 1j * (
                u.T @ v - v.T @ u + u @ v.T - v @ u.T
            )
    return terms


def directions():
    _, _, y, z = station.phase.AXES
    result = []
    for a, b in ((y, z), (z, y)):
        d = np.zeros((2, 4, 4), complex)
        d[1] = np.kron(a, b) / 2
        result.append(d)
    return result


def source_and_changes():
    states = survey.preparations()
    chi = states["equilibrium"]
    row = json.loads(
        (Path(signed_transfers.__file__).with_name("signed_transfers.json")).read_text()
    )
    angle = {r["source"]: r["angle_per_unit"] for r in row["tests"]}
    rot = (states["plus"] - states["minus"]) / 2
    matched = {
        "ordinary": states["ordinary"] - chi + angle["ordinary"] * rot,
        "unpairing": states["pairing"] - chi + angle["collar"] * rot,
    }
    basis = directions()
    menu = station.radar_menus()[6]
    dV = station.range_gradient(station.pair(menu, chi))[1]
    components = np.array([dV @ station.pair(menu, d) for d in basis])
    normal = np.linalg.norm(components)
    weights = np.array([-components[1], components[0]]) / normal
    matched["phase"] = sum(w * d for w, d in zip(weights, basis))
    K = station.log_blocks(states["pairing"]) - station.log_blocks(chi)
    charges = {
        name: dict(
            collar=float(station.pair(K, d)),
            full=float(-station.pair(station.log_blocks(chi), d)),
            energy=float(np.trace(d[1]).real),
        )
        for name, d in matched.items()
    }
    return (
        chi,
        matched,
        charges,
        dict(volume_basis=components.tolist(), phase_coefficients=weights.tolist()),
    )


def population_read(board, route, aperture, q, time, states):
    """Free finite-time flight; backward analyzer is the complex conjugate of the launch."""
    x, y, z = aperture
    indices = np.array([4 * p + k for p in (x, y, z) for k in range(4)])
    effect = np.zeros((2, 4, 4), complex)
    for odd in (0, 1):
        h = front.native.generator(board, route, q, odd) / 2
        init = np.zeros((h.shape[0], 12), complex)
        init[indices] = np.eye(12)
        g = (-1j * time) * h
        matrix = expm_multiply(g, init, traceA=g.diagonal().sum())[indices]
        a = (matrix[8:12, :4] + 1j * matrix[8:12, 4:8]) / np.sqrt(2)
        b = (matrix[:4, 8:12] + 1j * matrix[4:8, 8:12]) / np.sqrt(2)
        effect[odd] = (a.conj().T @ a + b.conj().T @ b) / 2
    return [float(station.pair(effect, rho)) for rho in states]


def blank_link(board, route):
    graph = board.adjacency().astype(np.int64)
    d = shortest_path(graph, unweighted=True, directed=False)
    cube = graph**3
    x, y = next(
        (a, b)
        for a, b in zip(*np.where(np.triu(d == 3, 1)))
        if d[a, list(route)].min() > 3 and cube[a, b] == 1
    )
    for q in (-1, 1):
        for odd in (0, 1):
            cube = front.hop(board, route, q, odd) ** 3
            block = cube[4 * y : 4 * y + 4, 4 * x : 4 * x + 4].toarray()
            np.testing.assert_array_equal(block, np.eye(4))
    return [int(x), int(y)]


def count_replay(board, route, phase_weights):
    # Unconditional pure even and 50:50 parity preparations synthesize the odd-block
    # linear slope; the individual source settings are native, positive preparations.
    _, _, yy, zz = station.phase.AXES
    settings = []
    contrasts = np.zeros((2, 16))
    for col, (a, b) in enumerate(((yy, zz), (zz, yy))):
        for i in (1, -1):
            for j in (1, -1):
                rho = np.kron((np.eye(2) + i * a) / 2, (np.eye(2) + j * b) / 2)
                for blocks, weight in (
                    ([rho, 0 * rho], -i * j / 2),
                    ([rho / 2, rho / 2], i * j),
                ):
                    contrasts[col, len(settings)] = weight
                    settings.append(np.stack(blocks))
    np.testing.assert_array_equal(
        np.einsum("ij,jbkl->ibkl", contrasts, settings), directions()
    )
    rng = np.random.default_rng(6202758)
    n = 2**68
    chunk = 2**60
    appointments = (1 / 64, 1 / 128)
    log_level = np.log(
        2 * len(APERTURES) * 2 * len(appointments) * len(settings) / 0.01
    )
    weights = np.array(phase_weights) @ contrasts
    rows = []
    for aperture in APERTURES:
        for q in (-1, 1):
            for t in appointments:
                probabilities = population_read(board, route, aperture, q, t, settings)
                assert min(probabilities) >= 0 and max(probabilities) <= 1
                clicks = [
                    int(rng.binomial(chunk, p, n // chunk).sum()) for p in probabilities
                ]
                est = np.array(clicks, dtype=float) / n
                # Simultaneous two-sided binomial Bernstein/Chernoff bound on Np in
                # terms of the observed click count X, including the case X=0.
                err = (
                    np.sqrt(2 * np.array(clicks, dtype=float) * log_level)
                    + 2 * log_level
                ) / n
                rows.append(
                    dict(
                        aperture=aperture,
                        q=q,
                        time=t,
                        click=clicks,
                        miss=[n - h for h in clicks],
                        phase_slope_population=float(weights @ probabilities / t**6),
                        phase_slope_count=float(weights @ est / t**6),
                        simultaneous_99pct_halfwidth=float(abs(weights) @ err / t**6),
                    )
                )
    return dict(shots_per_setting=n, seed=6202758, rows=rows)


def experiment():
    b = front.optics.pairing.formation.federation(1)
    route = front.optics.pairing.formation.region(b)
    numerators = operators(b, route)
    control = operators(b, route, 0) / (4 * NORM)
    assert np.all(control.imag == 0)
    physical = numerators / (4 * NORM)
    basis = directions()
    transfer = np.array([station.pair(physical, d) for d in basis]).T
    np.testing.assert_array_equal(transfer, [[-1 / 8, 0], [0, -9 / 64]])
    chi, ds, charges, geom = source_and_changes()
    rows = []
    for idx, op in enumerate(physical):
        f = float(station.pair(op, chi))
        length = (NORM * f) ** (-1 / 6)
        df = {tag: float(station.pair(op, v)) for tag, v in ds.items()}
        dl = {tag: float(-length * val / (6 * f)) for tag, val in df.items()}
        fringe = {tag: float(station.pair(1j * op.imag, v)) for tag, v in ds.items()}
        calibrations = {}
        for charge in ("collar", "full"):
            gain = dl["ordinary"] / charges["ordinary"][charge]
            calibrations[charge] = dict(
                ordinary_gain=gain,
                held={
                    tag: dict(
                        predicted=gain * charges[tag][charge],
                        observed=dl[tag],
                        observed_over_predicted=dl[tag] / (gain * charges[tag][charge]),
                    )
                    for tag in ("unpairing", "phase")
                },
            )
        rows.append(
            dict(
                aperture=APERTURES[idx],
                probability_coefficient=f,
                normalized_threshold_length=length,
                coefficient_slopes=df,
                threshold_slopes=dl,
                signed_fringe_slopes=fringe,
                calibrations=calibrations,
            )
        )
    return dict(
        phase_arm_transfer=transfer.tolist(),
        exchange_only_phase_transfer=np.stack(
            [station.pair(control, d) for d in basis]
        ).T.tolist(),
        blank_three_hop_standard=blank_link(b, route),
        finite_count_replay=count_replay(b, route, geom["phase_coefficients"]),
        source_geometry=geom,
        parked_source_slopes=charges,
        channels=rows,
    )


if __name__ == "__main__":
    print(json.dumps(experiment(), indent=2))
