"""Read boundary histories and regional entropy for two freely evolving load modes."""

import argparse
import hashlib
import json
from math import lgamma, log
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.special import jv

from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    PairSector,
    pair_index,
)
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.experiments.contact import radar
from model.experiments.material_inertia.dark_collar import (
    common_pair_generator,
    label,
    native_readouts,
)
from .backreaction import ordered_collar_control
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
        opposite=(slow - high) / np.sqrt(2),
        quadrature=(slow - 1j * high) / np.sqrt(2),
        conjugate=(slow + 1j * high) / np.sqrt(2),
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


def marked_flights(base, starts, marks, gaps, times):
    """Simultaneous Chebyshev flights with the two native contact energies in each branch."""
    diagonal = base.diagonal()
    off_radius = np.asarray(abs(base).sum(axis=1)).ravel() - abs(diagonal)
    lower = diagonal - off_radius
    lower[marks] -= np.max(gaps, axis=1)
    lo, hi = lower.min(), (diagonal + off_radius).max()
    center, radius = (hi + lo) / 2, (hi - lo) / 2
    scaled = (base - center * sparse.eye(base.shape[0], format="csr")) / radius

    def apply(v):
        w = scaled @ v
        w[marks] -= gaps / radius * v[marks]
        return w

    # Poisson's integral bounds |J_k(x)| by (x/2)^k/k!. Bound the omitted
    # tail by a geometric series once its consecutive ratio is less than one.
    limits = []
    for t in times:
        x = radius * t
        k = max(2, int(x))
        while log(2) + (k + 1) * log(x / 2) - lgamma(k + 2) - log(
            1 - x / (2 * (k + 2))
        ) > log(1e-13):
            k += 1
        limits.append(k)
    coefficients = [
        jv(np.arange(k + 1), radius * t) for k, t in zip(limits, times, strict=True)
    ]
    waves = [(coef[0] * starts).astype(complex) for coef in coefficients]
    previous, current = starts, apply(starts)
    for k in range(1, max(limits) + 1):
        for wave, coef in zip(waves, coefficients, strict=True):
            if k < len(coef):
                component = wave.imag if k % 2 else wave.real
                component += (2 if k % 4 in (0, 3) else -2) * coef[k] * current
        if k < max(limits):
            previous, current = current, 2 * apply(current) - previous
    # The omitted common exp(-i center*t) has no effect on any register.
    for wave in waves:
        np.testing.assert_allclose(np.linalg.norm(wave, axis=0), 1, atol=2e-12, rtol=0)
    return waves, dict(
        gershgorin_bounds=[float(lo), float(hi)], orders=limits, norm_tail_bound=1e-13
    )


def clock_and_detector_measurement(sector, home, remote, readers):
    """Commuting local effects, with every occupancy retained under all four analyzer choices."""
    port = 12 * home
    outside = np.r_[np.arange(port), np.arange(port + 12, sector.n)]
    hi = np.arange(port, port + 12)[:, None]
    one = pair_index(sector.n, np.minimum(hi, outside), np.maximum(hi, outside))
    nh = (sector.aa // 12 == home).astype(int) + (sector.bb // 12 == home)
    nr = (sector.aa // 12 == remote // 12).astype(int) + (
        sector.bb // 12 == remote // 12
    )
    # Remote ordering: no unit, two units, one at entry, one elsewhere.
    r = np.where(
        nr < 1,
        0,
        np.where(
            nr > 1, 1, np.where((sector.aa == remote) | (sector.bb == remote), 2, 3)
        ),
    )
    outside_r = np.where(
        outside // 12 != remote // 12, 0, np.where(outside == remote, 2, 3)
    )
    rays = np.stack(
        [readers[key][4] for key in ("in_phase", "opposite", "quadrature", "conjugate")]
    )

    def read(wave):
        prob = abs(wave) ** 2
        block = wave[one]
        target = abs(rays @ block) ** 2
        table = np.zeros((4, 4, 4))
        for n in (0, 2):
            q = nh == n
            table[:, n // 2, :] = np.bincount(r[q], weights=prob[q], minlength=4)
        for j in (0, 2, 3):
            allowed = outside_r == j
            table[:, 2, j] = target[:, allowed].sum(axis=1)
            table[:, 3, j] = np.sum(abs(block[:, allowed]) ** 2) - table[:, 2, j]
        np.testing.assert_allclose(table.sum(axis=(1, 2)), 1, atol=3e-12, rtol=0)
        np.testing.assert_allclose(
            table.sum(axis=1),
            np.broadcast_to(table[0].sum(axis=0), (4, 4)),
            atol=3e-12,
            rtol=0,
        )
        return table

    return read


def hand_summary(table, clock, trials=None):
    clicks = table[:, 2, :]
    r = table.mean(axis=0).sum(axis=0)

    def hand(p):
        z = p[0] - p[1] + 1j * (p[2] - p[3])
        mass = p.sum() / 2
        result = dict(
            mass=float(mass),
            phase=float(np.angle(z)) if abs(z) > 1e-14 else None,
            visibility=float(abs(z) / mass) if mass > 1e-14 else None,
        )
        if result["phase"] is None:
            return result
        frequency = clock["calibrated_gap"]
        reading = (result["phase"] % (2 * np.pi)) / frequency
        result["native_exposure_mod_period"] = reading
        if trials is not None:
            variance = p * (1 - p) / trials
            noise = variance.sum()
            # Do not quote a delta-method phase error for a near-empty
            # conditional bin or a phasor whose length is compatible with zero.
            if trials * mass >= 100 and abs(z) ** 2 > 25 * noise:
                se = (
                    np.sqrt(
                        z.imag**2 * variance[:2].sum() + z.real**2 * variance[2:].sum()
                    )
                    / abs(z) ** 2
                )
                result["phase_se"] = float(se)
                result["native_exposure_se"] = float(
                    np.hypot(se, reading * clock["calibrated_gap_standard_error"])
                    / frequency
                )
        return result

    xplus = table[0]
    return dict(
        all_trials=hand(clicks.sum(axis=1)),
        receiver_probabilities=r.tolist(),
        by_receiver=[hand(p) for p in clicks.T],
        xplus_joint_vs_independent_tv=float(
            abs(xplus - np.outer(xplus.sum(axis=1), xplus.sum(axis=0))).sum() / 2
        ),
    )


def joint_clock_experiment(readers, home, clock):
    board = federation(1)
    gateway = 12 * home
    remote = next(
        b if a == gateway else a for a, b in board.external if gateway in (a, b)
    )
    regular = [
        a
        for a in board.neighbours[remote]
        if a // 12 == remote // 12 and board.degree[a] == 5
    ]
    a, c = next(
        (a, c) for a in regular for c in regular if a < c and c in board.neighbours[a]
    )
    edges = [(a, remote), (remote, c)]
    source = json.loads(Path(__file__).with_name("backreaction.json").read_text())[
        "correlated_contacts"
    ]["ordered_three_levels"]
    levels = np.rint(source["source_energies"]).astype(int)
    mixtures = (
        np.array([row["source_units_of_1_over_36"] for row in source["arms"]]) / 36
    )
    np.testing.assert_array_equal(levels, [0, 32, 33])
    np.testing.assert_allclose(
        mixtures + mixtures.transpose(0, 2, 1),
        np.broadcast_to(mixtures[1] + mixtures[1].T, mixtures.shape),
    )
    law = candidate.native_law
    ex, con = law.local_term_parts()
    native = ex + candidate.KAPPA * con
    rails = [law.SRC.index(w) for w in ((1, 0), (0, 1))]
    prepared = ordered_collar_control(
        board,
        a,
        remote,
        c,
        native,
        native_readouts()[0],
        label,
        expm(-0.5j * np.pi * native)[np.ix_(rails, rails)],
        expm(-1.5j * np.pi * native)[np.ix_(rails, rails)],
        source["shots_per_setting"],
    )
    np.testing.assert_array_equal(
        [row["source_units_of_1_over_36"] for row in prepared["arms"]],
        np.rint(36 * mixtures),
    )
    inputs = np.eye(law.N_SRC)[
        :, [law.SRC.index(w) for w in [(0, 0), (1, 0), (0, 1), (1, 1)]]
    ]
    for energy, rec in zip(
        levels, [label(0, 1), label(2, 1), label(2, -1)], strict=True
    ):
        columns = np.kron(rec[:, None], inputs)
        np.testing.assert_allclose(
            native @ columns,
            columns * (energy * np.array([1, 1, 1, 0])),
            atol=3e-13,
            rtol=0,
        )
    # A native signed echo phase on (4,5) gives the symmetry-breaking control.
    signs = np.repeat([-9, *(d for _, d in law.LABELS)], law.N_SRC)
    injection = np.eye(law.DIM, law.N_SRC) @ inputs[:, :3]
    phase = law.echo() @ ((1 - 2 * (signs == 1))[:, None] * (law.echo() @ injection))
    np.testing.assert_allclose(
        phase, injection * np.array([1, -1, 1]), rtol=0, atol=1e-14
    )
    ref = np.arange(board.n)
    ref[gateway : gateway + 12] = gateway + np.array(
        [0, 1, 2, 3, 6, 7, 4, 5, 9, 8, 11, 10]
    )
    if {tuple(sorted((int(ref[i]), int(ref[j])))) for i, j in board.seams} != set(
        board.seams
    ):
        raise ValueError("The proposed reflection is not a full-federation symmetry")
    psi = readers["in_phase"].conj().T[:, 4]
    np.testing.assert_allclose(psi.imag, 0, atol=2e-12, rtol=0)
    psi = psi.real
    np.testing.assert_allclose(
        psi[ref[gateway : gateway + 12] - gateway], -psi, rtol=0, atol=2e-12
    )
    sector = PairSector(board)
    free = sector.hamiltonian(candidate.KAPPA)
    base = free.copy()
    for edge in edges:
        base += common_pair_generator(sector, edge) - free
    del free
    marks = np.array([sector.index(*e) for e in edges])
    reflection = pair_index(
        sector.n,
        np.minimum(ref[sector.aa], ref[sector.bb]),
        np.maximum(ref[sector.aa], ref[sector.bb]),
    )
    np.testing.assert_array_equal(reflection[marks], marks)
    if (base[reflection][:, reflection] != base).nnz:
        raise ValueError(
            "The complete two-load base does not commute with the reflection"
        )
    starts = np.zeros((sector.P, 10))
    ids = [sector.index(gateway + k, remote) for k in range(12)]
    starts[ids] = psi[:, None]
    starts[ids[4], 1:] *= -1
    gaps = np.array([[0, 0], *[(int(i), int(j)) for i in levels for j in levels]]).T
    times = [0.25, 0.75, 4.0]
    waves, certificate = marked_flights(base, starts, marks, gaps, times)
    read = clock_and_detector_measurement(sector, home, remote, readers)
    trials = 2**22
    rows = []
    for it, (time, wave) in enumerate(zip(times, waves, strict=True)):
        np.testing.assert_allclose(wave[reflection, 0], -wave[:, 0], rtol=0, atol=2e-12)
        np.testing.assert_allclose(wave[marks, 0], 0, rtol=0, atol=2e-12)
        null = read(wave[:, 0])
        alternatives = np.stack([read(v) for v in wave[:, 1:].T])
        contacts = (abs(wave[marks, 1:]) ** 2).T
        rv = []
        for control, branches in enumerate(
            (np.broadcast_to(null, alternatives.shape), alternatives)
        ):
            arm_results = []
            for arm, weights in enumerate(mixtures):
                joint = np.einsum("w,wsij->sij", weights.ravel(), branches)
                hist = []
                word_hist = []
                digest = hashlib.sha256()
                for s in range(4):
                    enrolled = weights.ravel()[:, None, None] * branches[:, s]
                    prob = cap.grid_distribution(enrolled.reshape(-1))
                    result = (
                        np.random.default_rng([SEED, 6, it, control, arm, s])
                        .multinomial(trials, prob)
                        .reshape(enrolled.shape)
                    )
                    digest.update(result.astype("<u8").tobytes())
                    hist.append(result.sum(axis=0).tolist())
                    word_hist.append(result.sum(axis=(1, 2)).tolist())
                arm_results.append(
                    dict(
                        direction=arm - 1,
                        joint_counts=hist,
                        source_counts=word_hist,
                        source_and_joint_count_digest=digest.hexdigest(),
                        exact=hand_summary(joint, clock),
                        sampled=hand_summary(np.array(hist) / trials, clock, trials),
                    )
                )
            diff = np.einsum(
                "w,wsij->sij", (mixtures[2] - mixtures[0]).ravel(), branches
            )
            upper = 0.0
            if control:
                for i, j in ((0, 1), (0, 2), (1, 2)):
                    v, w = wave[:, 1 + 3 * i + j], wave[:, 1 + 3 * j + i]
                    shift = np.vdot(w, v)
                    # Pure-state trace distance <= the aligned vector distance;
                    # unchanged diagonal weights and opposite (i,j)/(j,i) terms cancel.
                    upper += abs(
                        mixtures[2, i, j] - mixtures[0, i, j]
                    ) * np.linalg.norm(v - w * shift / abs(shift))
            first, last = (
                np.asarray(arm_results[k]["joint_counts"])[0, 2, 2] / trials
                for k in (0, 2)
            )
            rv.append(
                dict(
                    arms=arm_results,
                    maximum_reader_source_tv=float(
                        abs(diff).sum(axis=(1, 2)).max() / 2
                    ),
                    any_load_measurement_source_tv_upper=float(
                        upper + 2 * certificate["norm_tail_bound"] if control else 0
                    ),
                    xplus_entry_change=float(diff[0, 2, 2]),
                    sampled_xplus_entry_change=float(last - first),
                    sampled_xplus_entry_standard_error=float(
                        np.sqrt((first * (1 - first) + last * (1 - last)) / trials)
                    ),
                )
            )
        rows.append(
            dict(
                exposure=time,
                symmetric=rv[0],
                local_phase_control=rv[1],
                local_phase_contact_probabilities=contacts.tolist(),
            )
        )
    return dict(
        home=home,
        receiver_port=remote,
        marks=edges,
        source_preparation_spares=prepared["spare_ports"],
        declared_signed_charge=[row["charge_change"] for row in prepared["arms"]],
        levels=levels.tolist(),
        source_units=np.rint(36 * mixtures).astype(int).tolist(),
        trials_per_arm_setting=trials,
        home_bins=["zero", "two", "ray_click", "other_one"],
        receiver_bins=["zero", "two", "one_at_entry", "one_elsewhere"],
        settings=["x_plus", "x_minus", "y_plus", "y_minus"],
        phase_control_port=gateway + 4,
        solver=certificate,
        appointments=rows,
    )


def radar_source_and_hand(board):
    """Lift both source collars and an independent physical hand through the native gate."""
    previous = json.loads(Path(__file__).with_name("backreaction.json").read_text())[
        "correlated_contacts"
    ]
    source = previous["ordered_three_levels"]
    a, b, c = previous["ports"]
    marks = [(a, b), (b, c)]
    hand = min(
        edge for edge in board.seams if edge[1] < 12 and not set(edge) & {a, b, c}
    )
    marks.append(hand)
    sites = tuple(sorted(set(range(12)) - {a, b, c, *hand})[:6])
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    native = exchange + candidate.KAPPA * contrast
    energies = np.rint(source["source_energies"]).astype(int)
    np.testing.assert_array_equal(energies, [0, 32, 33])
    words = [law.SRC.index(w) for w in [(0, 0), (1, 0), (0, 1), (1, 1)]]
    records = [label(0, 1), label(2, 1), label(2, -1)]
    # Each marked word is scalar even when the single load reaches either endpoint.
    for energy, memory in zip(energies, records, strict=True):
        cols = np.kron(memory[:, None], np.eye(law.N_SRC)[:, words[:3]])
        np.testing.assert_allclose(native @ cols, energy * cols, rtol=0, atol=2e-13)
    h = candidate.generator(board).toarray()
    blank = law.seam_term(candidate.KAPPA)[np.ix_(words[1:3], words[1:3])]
    for edge in marks:
        h[np.ix_(edge, edge)] -= blank
    mixtures = (
        np.array([row["source_units_of_1_over_36"] for row in source["arms"]]) / 36
    )
    tangent = (mixtures[2] - mixtures[0]) / 2
    np.testing.assert_allclose(mixtures[1] + tangent, mixtures[2], atol=1e-15)
    entropy_slope = -float(np.sum(tangent * np.log(mixtures[1])))
    np.testing.assert_allclose(entropy_slope, np.log(3) / 12, atol=2e-15)
    # Encode a sharp signed word by echo and remove the two metering units.
    # At a successful reunion restore (1,1), read echo, then the enrolled quadrature.
    both = np.eye(law.N_SRC)[:, words[3]]
    levels = np.column_stack([np.kron(r, both) for r in records[1:]])
    sharp = law.SRC.index((2, 0))
    encoded = law.echo()[:, sharp]
    initial = levels.T @ encoded
    np.testing.assert_allclose(levels @ initial, encoded, rtol=0, atol=2e-14)
    gap = energies[2] - energies[1]
    n = np.arange(2, 2 * radar.TICKS + 1)
    evolved = initial[:, None] * np.exp(
        -1j * (energies[1:] - energies[1])[:, None] * radar.DT * n
    )
    decoded = law.echo() @ levels @ evolved
    sine_read = expm(-1.5j * np.pi * native) @ decoded
    yes = abs(np.vstack([decoded[sharp], sine_read[sharp]])) ** 2
    np.testing.assert_allclose(
        2 * yes - 1,
        [np.cos(gap * radar.DT * n), np.sin(gap * radar.DT * n)],
        atol=2e-13,
    )
    # Supply one metering load at each private spare after a successful return.
    # Native swaps then feed the hand without touching the parked radar load.
    swap = expm(-1j * np.pi * native)
    rails = words[1:3]
    np.testing.assert_allclose(
        abs(swap[np.ix_(rails, rails)]) ** 2, [[0, 1], [1, 0]], atol=2e-13
    )
    spares = []
    for site in sites:
        eligible = set(range(12)) - {a, b, c, *hand, site}
        spares.append(
            min(
                (p, q)
                for p in board.neighbours[hand[0]]
                if p in eligible
                for q in board.neighbours[hand[1]]
                if q in eligible and p != q
            )
        )
    return (
        h,
        sites,
        yes,
        dict(
            source_ports=[a, b, c],
            hand_seam=hand,
            terminal_spares_by_origin=spares,
            source_word_energies=energies.tolist(),
            hand_gap=int(gap),
            source_units=np.rint(36 * mixtures).astype(int).tolist(),
            interior_source_entropy_tangent=entropy_slope,
            source_charge_tangent=source["arms"][2]["charge_change"],
            source_derivatives="per unit epsilon at the middle preparation",
        ),
    )


def five_site_rank_bound(lower, upper, relative_tolerance):
    """All ten squared distances of a five-point subset vary independently in this box."""
    subset = (0, 1, 2, 4, 5)
    i, j = np.triu_indices(5, 1)
    lo = ((1 - relative_tolerance) * lower[np.ix_(subset, subset)][i, j]) ** 2
    hi = ((1 + relative_tolerance) * upper[np.ix_(subset, subset)][i, j]) ** 2
    selectors = np.arange(1 << len(i))[:, None] & (1 << np.arange(len(i)))
    squares = np.zeros((len(selectors), 5, 5))
    squares[:, i, j] = squares[:, j, i] = np.where(selectors, hi, lo)
    radii = squares[:, 0, 1:]
    gram = (radii[:, :, None] + radii[:, None, :] - squares[:, 1:, 1:]) / 2
    values, modes = np.linalg.eigh(gram)
    rebuilt = (modes * values[:, None, :]) @ modes.transpose(0, 2, 1)
    residual = float(np.linalg.norm(gram - rebuilt, axis=(1, 2)).max())
    orthogonal = modes.transpose(0, 2, 1) @ modes - np.eye(4)
    orthogonality_error = float(np.linalg.norm(orthogonal, axis=(1, 2)).max())
    guard = 1e-9 * max(1, np.linalg.norm(gram, axis=(1, 2)).max())
    if residual + abs(values).max() * orthogonality_error > guard / 10:
        raise ValueError("Corner eigensolver exceeds its numerical guard")
    # The anchor Gram depends affinely on the squared distances. Every box
    # point is a convex combination of its vertices. A positive lower eigenvalue
    # at every vertex forces four independent edge vectors throughout the box.
    return dict(
        site_indices=list(subset),
        tolerance=relative_tolerance,
        vertices=len(gram),
        smallest_vertex_eigenvalue=float(values[:, 0].min()),
        worst_vertex=int(values[:, 0].argmin()),
        maximum_eigendecomposition_residual=residual,
        maximum_orthogonality_error=orthogonality_error,
        numerical_guard=float(guard),
        excludes_three_dimensions=bool(np.all(values[:, 0] > guard)),
    )


def approximate_radar_embedding(distances, precision_run):
    """Check a supplementary in-sample construction and its global lower bound."""
    # Fifteen-distance minimax fit by SLSQP. Reproduction needs only these
    # coordinates, the measured distances and the independent Gram lower bound.
    xyz = np.array(
        [
            [0.50788223, 0.98866613, -1.72882015],
            [1.00278738, -0.04574380, 1.27493344],
            [-0.78366558, -0.39789148, 1.46656876],
            [0.21597213, -1.65505773, 0.01993595],
            [-0.56376765, -0.46392804, -2.00488774],
            [-0.37920849, 1.57395501, 0.97226975],
        ]
    )
    fit = np.linalg.norm(xyz[:, None] - xyz, axis=2)
    i, j = np.triu_indices(6, 1)
    errors = abs(fit[i, j] / distances[i, j] - 1)
    fit_bound = five_site_rank_bound(distances, distances, 0.0442663)
    sample_bound = five_site_rank_bound(
        np.array(precision_run["distance_confidence_lower"]),
        np.array(precision_run["distance_confidence_upper"]),
        0.034,
    )
    if (
        not fit_bound["excludes_three_dimensions"]
        or not sample_bound["excludes_three_dimensions"]
    ):
        raise ValueError(
            "Five-site rank obstruction did not exceed its numerical guard"
        )
    if errors.max() > 0.0442665:
        raise ValueError("The coordinates exceed the claimed all-dyad relative error")
    return dict(
        sites_used_in_fit=6,
        distances_used_in_fit=15,
        coordinates=xyz.tolist(),
        relative_errors=errors.tolist(),
        maximum_relative_error=float(errors.max()),
        exact_ensemble_lower_bound=fit_bound,
        material_hand_interval_lower_bound=sample_bound,
    )


def radar_geometry():
    """One-load projective first-hit histories with an independently decoded material clock."""
    board = federation(0)
    h, sites, hand_yes, source = radar_source_and_hand(board)
    dt, ticks = radar.DT, radar.TICKS
    unitary = radar._unitary(h, dt)
    legs = radar.first_hit(unitary, sites, ticks)
    np.testing.assert_allclose(
        legs, radar.renewal_prediction(unitary, sites, ticks), rtol=0, atol=5e-13
    )
    joint = radar.round_trip_law(legs)
    exact, _ = radar.distance_statistics(joint, dt)
    geometry = radar.geometry_test(exact)
    x = np.arange(6) / 5
    xyz = np.column_stack((x, x**2, x**3))
    calibration = radar.geometry_test(np.linalg.norm(xyz[:, None] - xyz, axis=2))
    if not calibration["all_pairs_accepted"]:
        raise ValueError(
            "The radar calibration failed on a known three-dimensional geometry"
        )

    acquisitions = [
        sample_material_radar(source, hand_yes, legs, joint, shots, seed)
        for shots, seed in ((131072, 2026112214), (4194304, 2026112215))
    ]
    return dict(
        cells=board.n // 12,
        ports=board.n,
        sites=sites,
        source=source,
        interval=dt,
        attempts_per_leg=ticks,
        exact_distances=exact.tolist(),
        exact_held_geometry=geometry,
        euclidean_calibration_accepted=calibration["accepted_pairs"],
        exact_spectral_bound=radar.euclidean_error_bound(exact),
        acquisitions=acquisitions,
        approximate_embedding=approximate_radar_embedding(exact, acquisitions[1]),
    )


def sample_material_radar(source, hand_yes, legs, joint, shots, seed):
    """Decode returned material hands and bound the fifteen conditional distance means."""
    rng = np.random.default_rng(seed)
    dt, ticks = radar.DT, radar.TICKS
    n = np.arange(2, 2 * ticks + 1)
    handset = np.zeros((2, len(n), 2), dtype=np.int64)
    by_pair = np.zeros((6, 6, len(n)), dtype=np.int64)
    bins = np.add.outer(np.arange(ticks), np.arange(ticks)).ravel()
    failures = np.zeros((3, 2, 6, 6, 3), dtype=np.int64)
    digest = hashlib.sha256()
    successes = []
    for left in range(6):
        for right in range(6):
            if left == right:
                continue
            p = joint[left, right]
            returned = np.bincount(bins, weights=p[:-1, :-1].ravel(), minlength=len(n))
            successes.append(float(returned.sum()))
            for arm in range(3):
                for axis, yes in enumerate(hand_yes):
                    categories = np.r_[
                        p[-1, -1],
                        p[:-1, -1].sum(),
                        returned * yes,
                        returned * (1 - yes),
                    ]
                    if int(np.argmax(categories)) > 1:
                        raise ValueError(
                            "The rounding-remainder bin must lie outside the clock observations"
                        )
                    weights = radar.acquisition_weights(categories)
                    counts = rng.multinomial(shots, weights / radar.LATTICE)
                    digest.update(counts.astype("<i8").tobytes())
                    retained = counts[2:].reshape(2, len(n)).T
                    handset[axis] += retained
                    by_pair[left, right] += retained.sum(axis=1)
                    failures[arm, axis, left, right] = [
                        counts[0],
                        counts[1],
                        retained.sum(),
                    ]
    np.testing.assert_array_equal(
        failures.sum(axis=-1)[:, :, ~np.eye(6, dtype=bool)], shots
    )
    hand_counts = handset.sum(axis=-1)
    if np.any(hand_counts == 0):
        raise ValueError("A returned duration bin has no physical hand records")
    xy = 2 * handset[:, :, 0] / hand_counts - 1
    phases = np.unwrap(np.arctan2(xy[1], xy[0]))
    # 46 conditional binary means share half the family error budget. Rounding
    # the acquisition table changes each conditional binary probability by at
    # most three lattice units divided by its physical duration-bin mass.
    masses = [
        np.convolve(legs[i, j, :-1], legs[j, i, :-1])
        for i in range(6)
        for j in range(6)
        if i != j
    ]
    min_bin = np.min(masses, axis=0)
    if np.any(min_bin <= 3 / radar.LATTICE):
        raise ValueError(
            "Insufficient true probability to bound material-hand rounding bias"
        )
    alpha = radar.CONFIDENCE_ALPHA
    xy_radius = 2 * (
        np.sqrt(np.log(4 * hand_counts.size / alpha) / (2 * hand_counts))
        + 3 / (radar.LATTICE * min_bin - 3)
    )
    disk = np.linalg.norm(xy_radius, axis=0)
    if np.any(disk >= 1):
        raise ValueError("Quadratures fail to resolve the enrolled duration bins")
    phase_radius = np.arcsin(disk)
    if abs(phases[0]) + phase_radius[0] >= np.pi or np.any(
        abs(np.diff(phases)) + phase_radius[:-1] + phase_radius[1:] >= np.pi
    ):
        raise ValueError(
            "The adjacent material clock intervals do not identify a winding"
        )
    durations = phases / source["hand_gap"]
    sym = by_pair + by_pair.transpose(1, 0, 2)
    population = sym.sum(axis=2)
    weights = sym / np.maximum(1, population)[:, :, None]
    measured = weights @ (durations / 2)
    # Half-trip times are in [dt,ticks*dt]. The existing two-preparation bound
    # covers 30 dyads; using 15 consumes at most the other half of alpha.
    calendar = weights @ (dt * n / 2)
    # Categories rounded to the acquisition lattice have <= 48/LATTICE total L1 error.
    # A returned conditional mean has bias <= range*error/(P(return)-error).
    error = (2 + 2 * len(n)) / radar.LATTICE
    rounding = (ticks - 1) * dt * error / (min(successes) - error)
    calendar_intervals, _, _ = radar.confidence_envelopes(
        calendar, population, dt, ticks
    )
    uncertainty = (
        np.array(calendar_intervals["statistical_half_width"])
        + weights @ (phase_radius / (2 * source["hand_gap"]))
        + rounding * (1 - np.eye(6))
    )
    lower, upper = np.maximum(0, measured - uncertainty), measured + uncertainty
    sampled_bound = radar.euclidean_error_bound(measured, lower, upper)
    return dict(
        shots_per_arm_direction_and_hand_setting=shots,
        seed=seed,
        clock_settings=["cosine", "sine"],
        acquisition_digest=digest.hexdigest(),
        count_columns=["outbound_timeout", "return_timeout", "both_hits"],
        counts_by_arm_and_clock=failures.sum(axis=(2, 3)).tolist(),
        counts_by_direction=failures.sum(axis=(0, 1)).tolist(),
        duration_bins=n.tolist(),
        hand_counts_setting_duration_yes_no=handset.tolist(),
        measured_hand_phase=phases.tolist(),
        phase_confidence_radii=phase_radius.tolist(),
        successes_by_unordered_pair_duration=[
            sym[i, j].tolist() for i in range(6) for j in range(i + 1, 6)
        ],
        exact_return_probability_range=[min(successes), max(successes)],
        distance_units="half round-trip material hand radians / native hand gap",
        sampled_hand_distances=measured.tolist(),
        distance_confidence_lower=lower.tolist(),
        distance_confidence_upper=upper.tolist(),
        simultaneous_failure_probability=alpha,
        sampled_hand_spectral_bound=sampled_bound,
    )


def clocked_radar_cohorts(
    names, joint_laws, pairs, hand_yes, gap, batch_size, batches, seed, alpha
):
    """Decode material-clock half trips, retaining both no-hit records on every start."""
    shots = batch_size * batches
    if batch_size > 2**42 or 4 * len(pairs) * shots >= np.iinfo(np.int64).max:
        raise ValueError(
            "Multinomial batch or accumulated enrollment exceeds its integer budget"
        )
    rng = np.random.default_rng(seed)
    n = np.arange(2, 2 * radar.TICKS + 1)
    bins = np.add.outer(np.arange(radar.TICKS), np.arange(radar.TICKS)).ravel()
    # One simultaneous bound covers all hand means and dyads in this acquisition.
    log_factor = np.log(2 * len(names) * (2 * len(n) + len(pairs)) / alpha)
    lattice_error = (2 + 2 * len(n)) / radar.LATTICE
    archive = []
    for name, joint in zip(names, joint_laws, strict=True):
        hand_counts = np.zeros((2, len(n), 2), dtype=np.int64)
        conditional_bias_sum = np.zeros((2, len(n)))
        retained = np.zeros((len(pairs), len(n)), dtype=np.int64)
        totals = np.zeros((len(pairs), 2, 2, 3), dtype=np.int64)
        success_floor = np.ones(len(pairs))
        digest = hashlib.sha256()
        for k, (i, j) in enumerate(pairs):
            # A persistent source word gives the convolution of its two directional
            # first hits. Reversing the order preserves every successful duration.
            directions = [joint[i, j], joint[j, i]]
            masses = [
                np.bincount(bins, weights=p[:-1, :-1].ravel(), minlength=len(n))
                for p in directions
            ]
            np.testing.assert_allclose(masses[0], masses[1], rtol=1e-12, atol=1e-17)
            for direction, (p, ok) in enumerate(zip(directions, masses, strict=True)):
                success_floor[k] = min(success_floor[k], ok.sum())
                for axis, yes in enumerate(hand_yes):
                    probabilities = np.r_[
                        p[-1, -1], p[:-1, -1].sum(), ok * yes, ok * (1 - yes)
                    ]
                    if int(np.argmax(probabilities)) > 1:
                        raise ValueError(
                            "The grid correction must fall in a timeout bin"
                        )
                    weights = radar.acquisition_weights(probabilities)
                    probabilities = weights / radar.LATTICE
                    if batches == 1:
                        counts = rng.multinomial(batch_size, probabilities)
                    else:
                        draws = rng.multinomial(batch_size, probabilities, size=batches)
                        counts = draws.sum(axis=0, dtype=np.int64)
                    digest.update(counts.astype("<i8").tobytes())
                    read = counts[2:].reshape(2, len(n)).T
                    hand_counts[axis] += read
                    duration_counts = read.sum(axis=1)
                    retained[k] += duration_counts
                    rounded = weights[2:].reshape(2, len(n))
                    rounded_yes = rounded[0] / np.maximum(1, rounded.sum(axis=0))
                    conditional_bias_sum[axis] += duration_counts * abs(
                        rounded_yes - yes
                    )
                    totals[k, direction, axis] = [*counts[:2], read.sum()]
        np.testing.assert_array_equal(totals.sum(axis=-1), shots)
        exposure = hand_counts.sum(axis=-1)
        if np.any(exposure == 0):
            raise ValueError("A duration bin has no observed material-hand quadrature")
        xy = 2 * hand_counts[:, :, 0] / exposure - 1
        phase = np.unwrap(np.arctan2(xy[1], xy[0]))
        # Given all realized bin counts, the retained hand bits are independent.
        # Weight each rounded conditional bias by its realized sample population
        # to keep the joint confidence guarantee even on extremely rare paths.
        xy_error = (
            np.sqrt(2 * log_factor / exposure) + 2 * conditional_bias_sum / exposure
        )
        disk_radius = np.linalg.norm(xy_error, axis=0)
        if np.any(disk_radius >= 1):
            raise ValueError("Insufficient independent clock reads in a retained bin")
        phase_bound = np.arcsin(disk_radius)
        if abs(phase[0]) + phase_bound[0] >= np.pi or np.any(
            abs(np.diff(phase)) + phase_bound[1:] + phase_bound[:-1] >= np.pi
        ):
            raise ValueError("The material clock does not identify its winding")
        populations = retained.sum(axis=1)
        frac = retained / np.maximum(1, populations)[:, None]
        half_time = frac @ phase / (2 * gap)
        width = (radar.TICKS - 1) * radar.DT
        # Sparse dyads can have no returns, and a quantity smaller than the
        # acquisition grid cannot bound a sampled conditional time nontrivially.
        round_bias = width * np.minimum(
            1, lattice_error / np.maximum(success_floor - lattice_error, lattice_error)
        )
        count_radius = width * np.sqrt(log_factor / (2 * np.maximum(1, populations)))
        hand_radius = frac @ phase_bound / (2 * gap)
        radius = count_radius + round_bias + hand_radius
        exact, _ = radar.distance_statistics(joint, radar.DT)
        ideal = np.array([exact[i, j] for i, j in pairs])

        def retained_values(values):
            return [
                float(v) if k > 0 else None
                for v, k in zip(values, populations, strict=True)
            ]

        archive.append(
            dict(
                name=name,
                ideal_distance=ideal.tolist(),
                material_hand_distance=retained_values(half_time),
                simultaneous_radius=retained_values(radius),
                min_max_return_probability=[min(success_floor), max(success_floor)],
                hand_counts_by_axis_and_duration=hand_counts.tolist(),
                returned_by_pair_duration=retained.tolist(),
                all_terminal_counts=(
                    totals.astype(str).tolist() if shots > 2**53 else totals.tolist()
                ),
                acquisition_sha256=digest.hexdigest(),
            )
        )
        if batches > 1:
            archive[-1]["radius_terms"] = dict(
                duration_sampling=count_radius.tolist(),
                conditional_duration_lattice=round_bias.tolist(),
                material_phase=hand_radius.tolist(),
            )

    return archive


def ordinary_load_comparison(archive, precision, pairs, held_pairs, sites):
    """Compare complete terminal probabilities in the same two independently acquired dyads."""
    # Four source arms, two dyads, two directions, three mutually exclusive
    # outcomes. Both clock axes share these endpoint probabilities. Pool them.
    cases = []
    family = 4 * len(held_pairs) * 2 * 3
    rounded_categories = 2 + 2 * (2 * radar.TICKS - 1)
    for held, pair in enumerate(held_pairs):
        for direction in range(2):
            means, errors = [], []
            for arm, index in zip(
                [*archive[:2], *precision],
                [pairs.index(pair)] * 2 + [held] * 2,
                strict=True,
            ):
                counts = np.asarray(arm["all_terminal_counts"][index][direction], dtype=np.int64).sum(axis=0)
                n = int(counts.sum())
                means.append(counts / n)
                errors.append(np.sqrt(np.log(2 * family / 0.01) / (2 * n)) + rounded_categories / radar.LATTICE)
            addresses = [sites[k] for k in (pair if direction == 0 else reversed(pair))]
            cases.append(dict(sites=addresses, ordinary_delta=(means[1] - means[0]).tolist(),
                              ordinary_radius=float(errors[1] + errors[0]),
                              signed_delta=(means[3] - means[2]).tolist(),
                              signed_radius=float(errors[3] + errors[2])))
    return dict(terminal_outcomes=["outbound timeout", "return timeout after outward hit", "round trip completed"],
                confidence=0.99, method="Simultaneous two-sided Hoeffding bound for all 48 arm probabilities, including acquisition grid bias.",
                cases=cases)


def occupied_source_radar():
    """Let a travelling unit range past two written contacts sharing a parked load."""
    board = federation(0)
    baseline, _, hand_yes, source = radar_source_and_hand(board)
    a, b, c = source["source_ports"]
    hand = tuple(source["hand_seam"])
    archived = json.loads(Path(__file__).with_name("backreaction.json").read_text())
    tag_record = archived["correlated_contacts"]["tagged_fluctuations"]
    tag = tuple(tag_record["native_tag_seam"])
    fixed = {b, *tag}
    marks = {(a, b), (b, c), hand}
    idle = [e for e in board.seams if fixed.intersection(e) and e not in marks]
    law = candidate.native_law
    words = [law.SRC.index(w) for w in ((0, 0), (1, 0), (0, 1), (1, 1))]
    full = sum(
        f * m for f, m in zip((1, candidate.KAPPA), law.local_term_parts(), strict=True)
    )
    blank = full[np.ix_(words[1:3], words[1:3])]
    for edge in idle:
        baseline[np.ix_(edge, edge)] -= blank
    free = sorted(set(range(board.n)) - fixed)
    np.testing.assert_allclose(baseline[np.ix_(free, sorted(fixed))], 0, atol=1e-14)
    gaps = np.array(source["source_word_energies"])
    for gap, record in zip(gaps, (label(0, 1), label(2, 1), label(2, -1)), strict=True):
        lifts = np.kron(record[:, None], np.eye(law.N_SRC)[:, words])
        np.testing.assert_allclose(
            full @ lifts, lifts * [gap, gap, gap, 0], rtol=0, atol=3e-13
        )
    (bridge,) = [e for e in board.external if e[0] < 12 and e not in idle]
    remote = bridge[1]
    sites = [
        4,
        7,
        8,
        remote,
        *sorted(p for p in board.neighbours[remote] if p // 12 == remote // 12)[:2],
    ]
    protected = {a, b, c, *hand, *tag}
    spares = [
        min(
            (p, q)
            for p in board.neighbours[hand[0]]
            for q in board.neighbours[hand[1]]
            if p != q
            and not {p, q}.intersection(protected | {origin})
            and max(p, q) < 12
        )
        for origin in sites
    ]
    pairs = [(i, j) for i in range(len(sites)) for j in range(i + 1, len(sites))]
    mixtures = np.array(source["source_units"]).reshape(3, 9) / 36
    np.testing.assert_array_equal(
        2 * np.array(source["source_units"])[1],
        np.array(source["source_units"])[[0, 2]].sum(axis=0),
    )

    def trip(energy_a, energy_c):
        h = baseline.copy()
        h[a, a] -= energy_a
        h[c, c] -= energy_c
        unitary = radar._unitary(h, radar.DT)
        legs = radar.first_hit(unitary, sites, radar.TICKS)
        # Independent renewal amplitudes check complete hits and both timeouts.
        alternate = radar.renewal_prediction(
            expm(-1j * radar.DT * h), sites, radar.TICKS
        )
        np.testing.assert_allclose(legs, alternate, rtol=0, atol=6e-13)
        return radar.round_trip_law(legs), radar.round_trip_law(alternate)

    by_word, alternate_by_word = (
        np.array(values)
        for values in zip(*(trip(x, y) for x in gaps for y in gaps), strict=True)
    )
    occupied = np.einsum("sw,wijab->sijab", mixtures, by_word)
    # The source word persists through outward transit, reflection and reunion.
    np.testing.assert_allclose(
        (occupied[0] + occupied[2]) / 2, occupied[1], rtol=0, atol=5e-15
    )
    joint_laws = [by_word[0], occupied[1], occupied[0], occupied[2]]
    names = [
        "vacant_reference",
        "occupied_reference",
        "occupied_minus",
        "occupied_plus",
    ]
    shots, seed, alpha = 2**36, 2026112240, 0.01
    archive = clocked_radar_cohorts(
        names, joint_laws, pairs, hand_yes, source["hand_gap"], shots, 1, seed, alpha
    )
    # An independent precision acquisition on the ordinary calibration dyad and a
    # remote control. The same supplied probe, appointments and material clock apply.
    held_pairs = [
        (sites.index(left), sites.index(right)) for left, right in ((7, 8), (48, 53))
    ]
    batch_size, batches, precision_seed = 2**42, 2**12, 2026112243
    precision = clocked_radar_cohorts(
        names[2:],
        joint_laws[2:],
        held_pairs,
        hand_yes,
        source["hand_gap"],
        batch_size,
        batches,
        precision_seed,
        alpha,
    )
    contrasts = []
    alternate_mixtures = np.einsum(
        "sw,wijab->sijab", mixtures[[0, 2]], alternate_by_word
    )
    alternate_ranges = [
        radar.distance_statistics(joint, radar.DT)[0] for joint in alternate_mixtures
    ]
    for k, (i, j) in enumerate(held_pairs):
        lower, upper = precision
        change = upper["material_hand_distance"][k] - lower["material_hand_distance"][k]
        radius = upper["simultaneous_radius"][k] + lower["simultaneous_radius"][k]
        contrasts.append(
            dict(
                sites=[sites[i], sites[j]],
                measured_change=change,
                native_change=upper["ideal_distance"][k] - lower["ideal_distance"][k],
                matrix_exponential_renewal_change=alternate_ranges[1][i, j]
                - alternate_ranges[0][i, j],
                simultaneous_interval=[change - radius, change + radius],
            )
        )
    # A fresh measurement after inspecting the first counts: repeat the outbound
    # detector comparison on all 30 directions. Each +/- pair comprises two
    # independent complete released preparations. Given m discordances, under
    # equality their signs are independent fair bits, regardless of the common
    # arrival rate. The conditional two-sided tail bound includes all quiet trials.
    repeat_shots, repeat_seed = 2**39, 2026112241
    repeat_rng = np.random.default_rng(repeat_seed)
    replications = []
    for i, j in pairs:
        for start, end in ((i, j), (j, i)):
            minus, plus = [
                float(1 - arm[start, end, -1, -1]) for arm in occupied[[0, 2]]
            ]
            qm, qp = np.rint(np.array([minus, plus]) * radar.LATTICE) / radar.LATTICE
            thresholds = np.array(
                [(1 - qm) * qp, qm * (1 - qp), (1 - qm) * (1 - qp), qm * qp]
            )
            # Round each arm separately before multiplying independent probabilities.
            # Equal native rates round identically; the contrast bias is <= 1/LATTICE.
            counts = repeat_rng.multinomial(repeat_shots, thresholds)
            discord = int(counts[0] + counts[1])
            difference = int(counts[0] - counts[1])
            bound = min(1.0, 2 * np.exp(-(difference**2) / (2 * max(discord, 1))))
            replications.append(
                dict(
                    origin=sites[start],
                    target=sites[end],
                    ordered_counts=counts.tolist(),
                    native_difference=plus - minus,
                    sampled_difference=difference / repeat_shots,
                    equal_rates_p_upper=float(bound),
                    rejects_equal_rates=bool(bound < alpha / 30),
                )
            )
    return dict(
        cells=board.n // 12,
        ports=board.n,
        source_ports=[a, b, c],
        shared_ordinary_port=b,
        protected_tag_seam=tag,
        hand_seam=hand,
        spares_by_origin=spares,
        remaining_bridge=bridge,
        sites=sites,
        dyads=[[sites[i], sites[j]] for i, j in pairs],
        idle_blank_seams=idle,
        source_units=source["source_units"],
        signed_log_score_shifts_nats=[-np.log(3) / 12, 0, np.log(3) / 12],
        clock_gap=source["hand_gap"],
        dt=radar.DT,
        attempts_per_leg=radar.TICKS,
        shots_per_direction_and_quadrature=shots,
        seed=seed,
        simultaneous_failure_probability=alpha,
        terminal_counts_axes=[
            "dyad",
            "direction",
            "clock quadrature",
            "outbound timeout / return timeout / both hits",
        ],
        acquisitions=archive,
        precision_signed_ranges=dict(
            seed=precision_seed,
            enrollments_per_batch=batch_size,
            batches_per_menu=batches,
            effective_enrollments_per_arm_direction_and_quadrature=str(
                batch_size * batches
            ),
            terminal_counts_format="Decimal integer strings retain exact totals above the JSON safe-integer range.",
            simultaneous_failure_probability=alpha,
            acquisition_lattice=radar.LATTICE,
            dyads=[[sites[i], sites[j]] for i, j in held_pairs],
            contrasts=contrasts,
            ordinary_load_probability_comparison=ordinary_load_comparison(
                archive, precision, pairs, held_pairs, sites
            ),
            acquisitions=precision,
        ),
        largest_tag_blind_record_difference=float(
            abs((occupied[0] + occupied[2]) / 2 - occupied[1]).max()
        ),
        independent_outward_replication=dict(
            paired_enrollments_per_direction=repeat_shots,
            seed=repeat_seed,
            family_false_positive_bound=alpha,
            directional_tests=30,
            maximum_contrast_rounding_bias=1 / radar.LATTICE,
            categories=[
                "only plus arrived",
                "only minus arrived",
                "both timed out",
                "both arrived",
            ],
            directions=replications,
        ),
        method="Full first-detection law per persistent ordered word. In the occupied arms the native two-load columns subtract the gap at each contact. Terminal material-hand winding uses only retained attempt ordering. Tag is a protected spectator during all flight.",
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
    clock = interior_clock(readers, home)
    return dict(
        preparation=certificate,
        screens=screens,
        cap_charge=cap_charge(modes),
        interior_clock=clock,
        entering_probe=crossing_probe(modes, home),
        joint_clock=joint_clock_experiment(readers, home, clock),
        radar_geometry=radar_geometry(),
        occupied_source_radar=occupied_source_radar(),
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
