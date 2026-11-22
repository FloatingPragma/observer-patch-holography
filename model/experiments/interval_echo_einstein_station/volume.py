"""Read boundary histories and regional entropy for two freely evolving load modes."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from . import cap

ARCHIVE = Path(__file__).with_name("volume.json")
SEED = 26112243


def native_modes():
    """Integer certificate over Q(sqrt(5)) and compilation of its two preparations."""
    car = cell()
    adjacency = car.adjacency().toarray().astype(int)
    # Twice the desired vector is base + lambda * direction, lambda^2=5.
    base = np.array([0, 0, 0, 0, 1, 1, -1, -1, -2, 2, -2, 2])
    direction = np.array([0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0])
    np.testing.assert_array_equal(adjacency @ base, 5 * direction)
    np.testing.assert_array_equal(adjacency @ direction, base)
    np.testing.assert_array_equal(np.r_[base[:3], direction[:3]], 0)
    eigenvalues = np.array([np.sqrt(5), -np.sqrt(5)])
    modes = base[:, None] + direction[:, None] * eigenvalues
    modes /= np.linalg.norm(modes, axis=0)
    np.testing.assert_allclose(modes.T @ modes, np.eye(2), atol=2e-15)

    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    port_columns = np.eye(law.DIM)[
        :, [law.SRC.index(w) for w in [(1, 0), (0, 1), (0, 0)]]
    ]
    echo = law.echo()
    labels = np.repeat([-9, *(d for _, d in law.LABELS)], law.N_SRC)
    interior = nx.Graph(car.seams).subgraph(set(range(12)) - {0, 1, 2})
    tree = list(nx.bfs_edges(interior, 4, sort_neighbors=sorted))
    slow, high = modes.T
    targets = dict(
        slow=slow,
        high=high,
        in_phase=(slow + high) / np.sqrt(2),
        quadrature=(slow - 1j * high) / np.sqrt(2),
    )
    readers, pulses = {}, {}
    for name, target in targets.items():
        work = abs(target) if np.iscomplexobj(target) else target.copy()
        compilation = []
        for a, b in reversed(tree):
            if work[b] == 0:
                continue
            angle = 2 * np.arctan2(work[b], work[a])
            compilation.append((a, b, angle))
            work[a], work[b] = np.hypot(work[a], work[b]), 0
        preparation = np.eye(12, dtype=complex)
        for a, b, angle in reversed(compilation):
            time = angle % (4 * np.pi)
            common = (candidate.KAPPA + 0.5) * time
            injection = port_columns if a < b else port_columns[:, [1, 0, 2]]
            second = labels == (-1 if a < b else 1)
            # An echo, signed-label phase, inverse echo addresses either occupied
            # endpoint while the third injected column checks the vacant spectator.
            before = np.exp(0.5j * np.pi * second)
            after = np.exp(1j * (common * (abs(labels) == 1) - np.pi / 2 * second))
            columns = echo @ (before[:, None] * (echo @ injection))
            columns = expm(-1j * time * h) @ columns
            columns = echo @ (after[:, None] * (echo @ columns))
            rotation = injection.T @ columns
            np.testing.assert_allclose(
                columns, injection @ rotation, atol=2e-13, rtol=0
            )
            np.testing.assert_allclose(rotation[:, 2], [0, 0, 1], atol=2e-14)
            preparation[[a, b]] = rotation[:2, :2] @ preparation[[a, b]]
        phases = 0
        if np.iscomplexobj(target):
            for a in np.flatnonzero(target):
                b = min(interior.neighbors(a))
                tag = labels == (1 if a < b else -1)
                columns = echo @ (
                    np.exp(1j * np.angle(target[a]) * tag)[:, None]
                    * (echo @ port_columns)
                )
                expected = np.ones(3, dtype=complex)
                expected[0 if a < b else 1] = np.exp(1j * np.angle(target[a]))
                np.testing.assert_allclose(columns, port_columns * expected, atol=2e-14)
                preparation[a] *= (port_columns.T @ columns).diagonal()[
                    0 if a < b else 1
                ]
                phases += 1
        np.testing.assert_allclose(preparation[:, 4], target, atol=3e-12, rtol=0)
        inverse = preparation.conj().T
        np.testing.assert_allclose(
            inverse @ preparation, np.eye(12), atol=3e-12, rtol=0
        )
        readers[name] = inverse
        pulses[name] = dict(rotations=len(compilation), port_phases=phases)
    return (
        modes,
        readers,
        dict(
            integer_part=base.tolist(),
            sqrt5_part=direction.tolist(),
            native_interior_pulses=pulses,
        ),
    )


def census(modes):
    """Classify all gateway triples by their antipodal pairs; check both full-board lifts."""
    local = cell()
    adjacency = local.adjacency().toarray().astype(int)
    ones = np.ones((12, 12), dtype=int)
    opposite = ones - ((np.eye(12, dtype=int) + adjacency + adjacency @ adjacency) > 0)
    np.testing.assert_array_equal(opposite.sum(axis=1), 1)
    # Exact spectral projector onto either rank-three band is (M +/- N/sqrt5)/4.
    # Unit diagonals after rescaling; an antipodal pair has inner product -1;
    # all other distinct vectors have inner product +/-1/sqrt5.
    M = np.eye(12, dtype=int) - opposite
    N = adjacency @ M
    np.testing.assert_array_equal(adjacency @ N, 5 * M)
    np.testing.assert_array_equal(M @ M + N @ N // 5, 4 * M)
    np.testing.assert_array_equal(N @ N % 5, 0)
    np.testing.assert_array_equal(M @ N + N @ M, 4 * N)
    np.testing.assert_array_equal(abs(N), ones - np.eye(12, dtype=int) - opposite)
    # The remaining nonscalar eigenspace is antipodally even, with zero mean.
    even = 6 * (np.eye(12, dtype=int) + opposite) - ones
    np.testing.assert_array_equal(adjacency @ even, -even)
    np.testing.assert_array_equal(even @ even, 12 * even)
    np.testing.assert_equal(np.trace(even), 60)
    rows = []
    for level in (0, 1, 2):
        board = federation(level)
        triples = [np.flatnonzero(d == 6) for d in board.degree.reshape(-1, 12)]
        antipodal = sum(bool(opposite[np.ix_(t, t)].any()) for t in triples)
        home = next(
            (i for i, t in enumerate(triples) if np.array_equal(t, [0, 1, 2])), None
        )
        if home is not None:
            packet = np.zeros((board.n, 2))
            packet[12 * home : 12 * home + 12] = modes
            energies = (
                6 * candidate.KAPPA + (5 - np.array([np.sqrt(5), -np.sqrt(5)])) / 2
            )
            np.testing.assert_allclose(
                candidate.generator(board) @ packet,
                packet * energies,
                atol=2e-14,
                rtol=0,
            )
        rows.append(
            dict(
                cells=board.cells,
                antipodal_gateway_triples=antipodal,
                matching_home=home,
            )
        )
    return rows


def cap_charge(modes):
    """Apply the existing 64-sweep reference and independent material meter at pole 0."""
    weights, _ = cap.cap_geometry(0)
    indices = [cap.seam_indices(e) for e in cell().seams]
    native = cap.NativeFilter()
    zero, one, odd, *_ = cap.exact_stationary_source(1, cap.standard_count())[
        "mode_probabilities"
    ]
    odds, beta = np.array([one, odd]) / zero, np.log(one / odd)
    filters = [
        (idx, native.success(odds, w / 256))
        for idx, w in zip(indices, weights, strict=True)
    ]
    sigma, charge, acceptance, weak, weak_acceptance = cap.prepared_charge(filters, 64)
    states = np.zeros((2, 79, 79))
    for rho, vector in zip(states, modes.T, strict=True):
        rho[1:13, 1:13] = np.outer(vector, vector)
    delta = states[1] - states[0]
    exact = np.sum(charge * delta)
    entropy_contrasts = []
    for fraction in (2**-16, 2**-20, 2**-24):
        entropies = []
        for rho in states:
            # The vacuum and two-load blocks agree and cancel at each fraction.
            mixture = (1 - fraction) * sigma[1:13, 1:13] + fraction * rho[1:13, 1:13]
            eigen = np.linalg.eigvalsh(mixture)
            entropies.append(-np.dot(eigen, np.log(eigen)))
        entropy_contrasts.append(
            dict(
                fraction=fraction, nats=float((entropies[1] - entropies[0]) / fraction)
            )
        )
    if abs(entropy_contrasts[-1]["nats"] - exact) > 2e-4:
        raise ValueError(
            "Two physical mixtures failed their independent entropy comparison"
        )
    energy, word_meter, projectors = cap.native_energy(indices, weights)
    weak_result, fitted, *_ = cap.weak_tomography(
        weak, delta, 64, beta, weak_acceptance
    )
    rv = [cap.read_charge(rho, beta * fitted, SEED + j) for j, rho in enumerate(states)]
    ev = [
        cap.read_charge(rho, energy, SEED + 2 + j, (word_meter, weights, projectors))
        for j, rho in enumerate(states)
    ]
    measured = rv[1][0] - rv[0][0]
    se = np.sqrt(
        rv[0][1] + rv[1][1] + (beta * weak_result["estimated_standard_error"]) ** 2
    )
    if abs(measured - exact) > 6 * se or exact < 0.4:
        raise ValueError("The measured signed entropy slope failed")
    return dict(
        pole=0,
        native_source_acceptance=acceptance,
        common_reference_minimum=float(np.linalg.eigvalsh(sigma)[0]),
        high_minus_slow_nats=float(exact),
        finite_entropy_contrasts=entropy_contrasts,
        measured_nats=float(measured),
        standard_error_nats=float(se),
        weak_reference_shots=weak_result["total_accepted_preparations"],
        weak_reference_count_digest=weak_result["packed_count_sha256"],
        independent_weighted_energy=float(np.sum(energy * delta)),
        measured_weighted_energy=float(ev[1][0] - ev[0][0]),
        energy_standard_error=float(np.sqrt(ev[0][1] + ev[1][1])),
    )


def crossing_probe(modes, home):
    """Add one unit on the neighboring gateway and release both on the whole native sphere."""
    board = federation(1)
    gateway = 12 * home
    remote = next(
        b if a == gateway else a for a, b in board.external if gateway in (a, b)
    )
    sector = PairSector(board)
    h = sector.hamiltonian(candidate.KAPPA)
    packet = np.zeros((sector.P, 2), complex)
    packet[[sector.index(gateway + j, remote) for j in range(12)]] = modes
    patch = (sector.aa // 12 == remote // 12).astype(int) + (
        sector.bb // 12 == remote // 12
    )
    bins = np.where(
        patch == 2,
        3,
        np.where(
            patch == 0, 2, np.where((sector.aa == remote) | (sector.bb == remote), 0, 1)
        ),
    )
    result = []
    for time in (0.5, 1.0):
        packet = expm_multiply(-0.5j * h, packet, traceA=-0.5j * h.diagonal().sum())
        outcomes = np.array(
            [np.bincount(bins, weights=abs(v) ** 2, minlength=4) for v in packet.T]
        )
        np.testing.assert_allclose(outcomes.sum(axis=1), 1, rtol=0, atol=2e-12)
        trials = 2**22
        counts = np.random.default_rng([SEED, 4, int(2 * time)]).multinomial(
            trials, cap.grid_distribution(outcomes)
        )
        empty = counts[:, 2] / trials
        error = np.sqrt(np.sum(empty * (1 - empty)) / trials)
        delta = outcomes[1, 2] - outcomes[0, 2]
        if abs((empty[1] - empty[0]) - delta) > 6 * error:
            raise ValueError("The crossing signal failed its independent count read")
        result.append(
            dict(
                exposure=time,
                probabilities=outcomes.tolist(),
                counts=counts.tolist(),
                quiet_difference=float(delta),
                measured_difference=float(empty[1] - empty[0]),
                standard_error=float(error),
            )
        )
    return dict(
        home=home,
        gateway=gateway,
        receiver_port=remote,
        cells=board.cells,
        two_load_states=sector.P,
        receiver_outcomes=["one_at_entry", "one_elsewhere", "quiet", "two_in_receiver"],
        appointments=result,
    )


def interior_clock(readers, home):
    """Read the freely rotating in-phase and quarter-phase hands entirely on the interior patch."""
    sphere = federation(1)
    region = slice(12 * home, 12 * home + 12)
    h = candidate.generator(sphere)
    times, trials = (0.25, 0.75, 4.0), 2**20
    packet = np.zeros((sphere.n, 3), complex)
    packet[region] = np.column_stack(
        [readers[name].conj().T[:, 4] for name in ("in_phase", "slow", "high")]
    )
    rows = []
    for j, time in enumerate(times):
        reached = expm_multiply(
            -1j * time * h, packet, traceA=-1j * time * h.diagonal().sum()
        )
        inside = reached[region]
        np.testing.assert_allclose(
            np.linalg.norm(inside, axis=0), 1, atol=3e-12, rtol=0
        )
        np.testing.assert_allclose(inside[:3], 0, atol=2e-12, rtol=0)
        read = np.array(
            [
                abs((readers[name] @ inside)[4]) ** 2
                for name in ("in_phase", "quadrature", "slow")
            ]
        )
        coherent, dephased = read[:, 0], read[:, 1:].mean(axis=1)
        angle = np.sqrt(5) * time
        np.testing.assert_allclose(
            coherent,
            [(1 + np.cos(angle)) / 2, (1 + np.sin(angle)) / 2, 0.5],
            atol=3e-12,
            rtol=0,
        )
        np.testing.assert_allclose(dephased, 0.5, atol=3e-12, rtol=0)
        returns = np.column_stack((coherent, dephased))
        counts = np.random.default_rng([SEED, 5, j]).binomial(
            trials, np.rint(returns * 2**31) / 2**31
        )
        cx, cy = 2 * counts[:2, 0] / trials - 1
        phase = np.arctan2(cy, cx)
        se = np.sqrt((cx**2 * (1 - cy**2) + cy**2 * (1 - cx**2)) / trials) / (
            cx**2 + cy**2
        )
        rows.append(
            dict(
                exposure=time,
                coherent_probabilities=coherent.tolist(),
                dephased_probabilities=dephased.tolist(),
                return_counts=counts.tolist(),
                phase=float(phase),
                phase_standard_error=float(se),
            )
        )
    gap, gap_se = (rows[0][key] / times[0] for key in ("phase", "phase_standard_error"))
    for row in rows[1:]:
        predicted = gap * row["exposure"]
        residual = np.angle(np.exp(1j * (row["phase"] - predicted)))
        se = np.hypot(row["phase_standard_error"], gap_se * row["exposure"])
        if abs(residual) > 6 * se:
            raise ValueError(
                "The held patch clock disagrees with its calibrated frequency."
            )
        row.update(
            held_phase_residual=float(residual), held_phase_standard_error=float(se)
        )
    return dict(
        home=home,
        return_port=12 * home + 4,
        trials_per_arm_and_setting=trials,
        settings=["in_phase", "quadrature", "slow_population"],
        arms=["coherent", "dephased"],
        calibrated_gap=float(gap),
        calibrated_gap_standard_error=float(gap_se),
        appointments=rows,
    )


def run():
    modes, readers, certificate = native_modes()
    screens = census(modes)
    home = screens[1]["matching_home"]
    # The previous cap compiler routed all local seam filters through one
    # adjacent cell. Transport its verified two-cell gate station by incidence.
    old = federation(0)
    old_bridge = min(old.external)
    new = federation(1)
    new_gateway = 12 * home
    new_remote = next(
        b if a == new_gateway else a for a, b in new.external if new_gateway in (a, b)
    )
    graph = nx.Graph(cell().seams)
    permutations = list(
        nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter()
    )
    old_remote_home = old_bridge[1] // 12
    new_remote_home = new_remote // 12
    maps = [
        next(p for p in permutations if p[old_bridge[j] % 12] == r % 12)
        for j, r in enumerate((new_gateway, new_remote))
    ]
    mapping = {
        12 * source + j: 12 * dest + p[j]
        for source, dest, p in zip(
            [0, old_remote_home], [home, new_remote_home], maps, strict=True
        )
        for j in range(12)
    }
    station = cap.routes([old_bridge])
    mapped_edges = {
        tuple(sorted((mapping[i], mapping[j])))
        for i, j in old.seams
        if i in mapping and j in mapping
    }
    actual_edges = {
        e for e in new.seams if e[0] in mapping.values() and e[1] in mapping.values()
    }
    if (
        len(station["paths"]) != 30
        or len(mapped_edges) != 61
        or mapped_edges != actual_edges
    ):
        raise ValueError(
            "The compiled cap services do not fit the selected physical cells"
        )
    certificate["cap_station_transport"] = dict(
        old_bridge=old_bridge,
        actual_bridge=[new_gateway, new_remote],
        registered_edges=len(mapped_edges),
        seam_filter_paths=len(station["paths"]),
    )
    return dict(
        preparation=certificate,
        screens=screens,
        cap_charge=cap_charge(modes),
        interior_clock=interior_clock(readers, home),
        entering_probe=crossing_probe(modes, home),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    if args.write:
        ARCHIVE.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        cap.compare(data, json.loads(ARCHIVE.read_text()))
    print(data["cap_charge"])
    print(data["entering_probe"]["appointments"][-1])
