"""Free propagation and native two-wing pointer records on the registered screen."""

import argparse
import json
from itertools import product
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_04_compact.experiment import free
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import free_seam

K = candidate.KAPPA
LOCAL = ((1, 0), (0, 1), (0, 0), (1, 1))
WORDS = [a + b for a, b in product(LOCAL, repeat=2)]
SCORES = np.array([1, -1, 1])  # positive, negative, quiet
COARSE = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1]])
CHSH = np.array([[1, 1], [1, -1]])
N = 2**28
HERE = Path(__file__).resolve().parent


def local_gate(time):
    law = candidate.native_law
    ids = [law.SRC.index(x) for x in LOCAL]
    return expm(-1j * time * free_seam(K))[np.ix_(ids, ids)]


def verify_read():
    law = candidate.native_law
    blank = np.eye(law.DIM)[:, [law.SRC.index(x) for x in LOCAL]]
    echo = law.echo()
    for bins in ((1,), (-1,), (0,), (-1, 1)):
        labels = np.repeat([False, *(d in bins for _, d in law.LABELS)], law.N_SRC)
        back = echo.T @ (labels[:, None] * (echo @ blank))
        np.testing.assert_allclose(
            back, blank * [a - b in bins for a, b in LOCAL], atol=2e-14
        )
    np.testing.assert_allclose(local_gate(np.pi), np.eye(4)[[1, 0, 2, 3]], atol=2e-14)


def sample(n, prob, rng):
    p = np.asarray(prob, float)
    if p.min() < -1e-10 or not np.isclose(p.sum(), 1, atol=2e-10, rtol=0):
        raise ValueError("Incomplete pointer distribution")
    p = np.maximum(p, 0).round(13)
    return rng.multinomial(n, p.ravel() / p.sum()).reshape(p.shape)


def interference(board, rng):
    # In the one-load sector the degree reference completes the contrast to 6K I.
    h = board.laplacian().toarray() / 2
    u = expm(-1j * h)
    # Decode the final position through local signed seam records. On this
    # one-load sector the oriented sign at each port's first seam singles it out.
    ids = np.arange(board.n)
    indicators = []
    for site, neighbours in enumerate(board.neighbours):
        left, right = sorted((site, neighbours[0]))
        sign = (ids == left).astype(int) - (ids == right).astype(int)
        indicators.append(sign == (1 if site == left else -1))
    np.testing.assert_array_equal(indicators, np.eye(board.n, dtype=bool))
    source = 0
    a, b = [
        i for i in board.neighbours[source] if board.degree[i] == 5 and i // 12 == 0
    ][:2]
    midway = u[:, source]
    first, second = u[:, a] * midway[a], u[:, b] * midway[b]
    panels = {}
    shots = 2**20

    def measure(name, probability):
        p = np.r_[probability, max(0, 1 - probability.sum())]
        panels[name] = sample(shots, p, rng).tolist()

    measure("first_only", abs(first) ** 2)
    measure("second_only", abs(second) ** 2)
    measure("which_aperture", abs(first) ** 2 + abs(second) ** 2)
    for j in range(8):
        # Waiting j*pi/(4*K) on the first degree-five reference gives this phase.
        measure(f"phase_{j}", abs(first * np.exp(-1j * np.pi * j / 4) + second) ** 2)
    measure("free", abs(u @ midway) ** 2)
    measure("free_unread_midpoint", abs(u) ** 2 @ abs(midway) ** 2)

    def ret(name):
        return panels[name][source] / shots

    return {
        "source": source,
        "apertures": [a, b],
        "shots_per_setting": shots,
        "screen_counts_then_absorbed": panels,
        "zero_phase_return": ret("phase_0"),
        "sum_single_returns": ret("first_only") + ret("second_only"),
        "unread_aperture_return": ret("which_aperture"),
        "free_return": ret("free"),
        "free_unread_midpoint_return": ret("free_unread_midpoint"),
    }


def orient(board, edge):
    return tuple(sorted(edge, key=lambda i: (board.degree[i], i)))


def free_source(board):
    edges = [e for e in board.external if e[0] // 12 == 0][:2]
    cells = [e[1] // 12 for e in edges]
    banks = [
        [
            orient(board, e)
            for e in board.seams
            if e[0] // 12 == e[1] // 12 == c and min(board.degree[list(e)]) == 5
        ]
        for c in cells
    ]
    placements = list(product(*banks))
    sector = PairSector(board)
    ports = np.array(
        [
            [[sector.index(a, b) for b in bob] for a in alice]
            for alice, bob in placements
        ]
    )
    start = np.zeros(sector.P)
    start[sector.index(edges[0][0], edges[1][0])] = 1
    h = sector.hamiltonian(K)
    history = []
    best = None
    for time in (1, 2, 4, 8):
        state = free(h, start, time)
        joint = state[ports]
        yields = (abs(joint) ** 2).sum(axis=(1, 2))
        concurrence = np.divide(
            2 * abs(np.linalg.det(joint)),
            yields,
            out=np.zeros_like(yields),
            where=yields > 0,
        )
        expected = 2 * np.sqrt(1 + concurrence**2)
        merit = yields * (expected - 2) ** 2
        i = int(np.argmax(merit))
        history.append(
            {
                "time": time,
                "eligible_placements": len(placements),
                "best_rails": placements[i],
                "joint_yield": float(yields[i]),
                "ideal_selected_CHSH": float(expected[i]),
            }
        )
        if best is None or merit[i] > best[0]:
            best = merit[i], time, placements[i], state
    _, time, rails, state = best
    changed = free(sector.hamiltonian(0), start, time)
    return (
        sector,
        rails,
        state,
        changed,
        {
            "source_ports": [e[0] for e in edges],
            "cells": cells,
            "survey": history,
            "selected_time": time,
            "rails": rails,
        },
    )


def reduce_wings(sector, rails, state):
    # Trace the actual exterior, retaining both single captures and all quiet sectors.
    sites = sum(rails, ())
    bits = [(site == sector.aa) | (site == sector.bb) for site in sites]
    code = sum(bit.astype(int) << i for i, bit in enumerate(bits))

    def word_code(word):
        return sum(x << i for i, x in enumerate(word))

    codes = [word_code(w) for w in WORDS]
    local = np.isin(sector.aa, sites).astype(int) + np.isin(sector.bb, sites)
    block = np.zeros((16, sector.n + 1), complex)
    for i, c in enumerate(codes):
        inside = (code == c) & (local == 2)
        block[i, -1] = state[inside].sum()
        one = (code == c) & (local == 1)
        outside = np.where(
            np.isin(sector.aa[one], sites), sector.bb[one], sector.aa[one]
        )
        block[i, outside] = state[one]
    rho = block @ block.conj().T
    rho[WORDS.index((0, 0, 0, 0)), WORDS.index((0, 0, 0, 0))] += (
        np.linalg.norm(state[local == 0]) ** 2
    )
    np.testing.assert_allclose(rho.trace(), 1, atol=5e-12)
    return rho


def axis_gate(axis):
    # On these seams the first rail has degree five. Apply its native degree
    # phase, then the local exchange rotation, then the signed echo.
    x, y, z = axis
    theta = np.arccos(np.clip(z, -1, 1))
    phi = (-np.arctan2(y, x) - np.pi / 2) % (2 * np.pi)
    g = local_gate(theta) @ np.diag(np.exp(-1j * phi * np.array([1, 0, 0, 1])))
    read = g[:2, :2].conj().T @ np.diag([1, -1]) @ g[:2, :2]
    np.testing.assert_allclose(read, [[z, x - 1j * y], [x + 1j * y, -z]], atol=2e-12)
    return g


def pointer(rho, alice, bob):
    # Each observer's effect acts only on its own two ports after isolation.
    joint = np.kron(alice, bob)
    sharp = np.diag(joint @ rho @ joint.conj().T).real.reshape(4, 4)
    return COARSE @ sharp @ COARSE.T


def axes_from_records(rho, rng, shots=2**26):
    gates = [axis_gate(v) for v in np.eye(3)]
    rows, tensor = [], np.empty((3, 3))
    for i, j in product(range(3), repeat=2):
        counts = sample(shots, pointer(rho, gates[i], gates[j]), rng)
        sub = counts[:2, :2]
        tensor[i, j] = (np.array([1, -1]) @ sub @ np.array([1, -1])) / sub.sum()
        rows.append(counts.tolist())
    u, s, vh = np.linalg.svd(tensor)
    # Singular vectors only select the settings. The test uses new shots.
    # Fix a reproducible paired sign convention for the two laboratory axes.
    for i in range(2):
        if u[np.argmax(abs(u[:, i])), i] < 0:
            u[:, i] *= -1
            vh[i] *= -1
    c, d = s[:2] / np.linalg.norm(s[:2])
    alice = [u[:, 0], u[:, 1]]
    bob = [c * vh[0] + d * vh[1], c * vh[0] - d * vh[1]]
    return (
        [axis_gate(v) for v in alice],
        [axis_gate(v) for v in bob],
        {
            "counts_XYZ_order": rows,
            "shots_per_axis_pair": shots,
            "alice_axes": np.array(alice).tolist(),
            "bob_axes": np.array(bob).tolist(),
        },
    )


def capture_axes_bound(training):
    """All projective two-rail axes on the same event-ready population."""
    sub = np.array(training["counts_XYZ_order"]).reshape(3, 3, 3, 3)[:, :, :2, :2]
    n = sub.sum(axis=(2, 3))
    signs = np.array([1, -1])
    tensor = np.einsum("i,xyij,j->xy", signs, sub, signs) / n
    singular = np.linalg.svd(tensor, compute_uv=False)
    estimate = 2 * np.linalg.norm(singular[:2])
    # The change in the leading-two singular-value norm is bounded by the
    # Frobenius norm of the error. Nine two-sided Hoeffding intervals suffice.
    error = 2 * np.sqrt((2 * np.log(18 / 0.01) / n).sum())
    return {
        "estimated_maximum": float(estimate),
        "simultaneous_99_upper": float(estimate + error),
    }


def all_axes_bound(training):
    """Bound any choice of the four axes with fixed quiet output +1."""
    shots = training["shots_per_axis_pair"]
    p = np.array(training["counts_XYZ_order"]).reshape(3, 3, 3, 3) / shots
    neither = float(p[:, :, 2, 2].mean())
    both = float(p[:, :, :2, :2].sum(axis=(2, 3)).mean())
    alice = (p[:, :, 0, 2] - p[:, :, 1, 2]).mean(axis=1)
    bob = (p[:, :, 2, 0] - p[:, :, 2, 1]).mean(axis=0)
    bound = 2 * (
        neither + np.linalg.norm(alice) + np.linalg.norm(bob) + np.sqrt(2) * both
    )
    # Eight quantities: two event probabilities, six expectations in [-1,1].
    log_tail = np.log(16 / 0.01)
    p_error = np.sqrt(log_tail / (2 * 9 * shots))
    component_error = np.sqrt(2 * log_tail / (3 * shots))
    error = 2 * ((1 + np.sqrt(2)) * p_error + 2 * np.sqrt(3) * component_error)
    return {
        "neither": neither,
        "both": both,
        "alice_only_vector": alice.tolist(),
        "bob_only_vector": bob.tolist(),
        "estimated_upper_bound": float(bound),
        "simultaneous_99_upper_bound": float(bound + error),
    }


def scores(tables):
    stats = {}
    for name, values in (("all_launches", SCORES), ("joint_capture", SCORES[:2])):
        width = len(values)
        selected = tables[:, :, :width, :width]
        counts = selected.sum(axis=(2, 3))
        if np.any(counts == 0):
            continue
        e = np.einsum("i,xyij,j->xy", values, selected, values) / counts
        # Conditional on the enrolled setting counts, the products are in [-1,1].
        radius = np.sqrt(2 * np.log(8 / 0.01) / counts).sum()
        stats[name] = {
            "S": float(np.sum(CHSH * e)),
            "simultaneous_99_radius": float(radius),
            "counts": counts.tolist(),
            "correlations": e.tolist(),
        }
    return stats


def bell(rho, ga, gb, rng, shots=N):
    p_heads = abs(local_gate(np.pi / 2)[0, 0]) ** 2
    p = np.array([[pointer(rho, a, b) for b in gb] for a in ga])
    err = max(
        np.max(abs(p[:, 0].sum(axis=-1) - p[:, 1].sum(axis=-1))),
        np.max(abs(p[0].sum(axis=-2) - p[1].sum(axis=-2))),
    )
    np.testing.assert_allclose(err, 0, atol=1e-12)
    for x, y in product(range(2), repeat=2):
        left, right = np.kron(ga[x], np.eye(4)), np.kron(np.eye(4), gb[y])
        np.testing.assert_allclose(left @ right, right @ left, atol=1e-13)
    # Both parity flags are written before settings. Keep unsuccessful flags
    # even when only one wing receives a load. The late analyzers commute with
    # the local flag, so erasing the flag gives the unconditioned sign law.
    yes = np.diag([1.0, 1.0, 0.0, 0.0])
    flags = [np.kron(a, b) for a, b in product((np.eye(4) - yes, yes), repeat=2)]
    branches = [f @ rho @ f for f in flags]
    weights = np.array([b.trace().real for b in branches])
    cohorts = sample(shots, weights, rng)
    tables = np.zeros((2, 2, 3, 3), dtype=int)
    predicted = np.zeros_like(p)
    records = []
    for k, (n, weight, branch) in enumerate(
        zip(cohorts, weights, branches, strict=True)
    ):
        n_a = sample(n, [p_heads, 1 - p_heads], rng)
        settings = [sample(m, [p_heads, 1 - p_heads], rng) for m in n_a]
        observed = np.zeros_like(tables)
        for x, y in product(range(2), repeat=2):
            probs = pointer(branch, ga[x], gb[y])
            predicted[x, y] += probs
            if weight > 1e-14:
                observed[x, y] = sample(settings[x][y], probs / weight, rng)
            elif settings[x][y]:
                raise ValueError("Enrolled empty capture branch")
        tables += observed
        records.append(
            {
                "alice_capture": bool(k // 2),
                "bob_capture": bool(k % 2),
                "enrolled_before_settings": int(n),
                "counts": observed.tolist(),
            }
        )
    np.testing.assert_allclose(p, predicted, atol=2e-12)
    np.testing.assert_array_equal(
        np.array(records[-1]["counts"])[:, :, :2, :2], tables[:, :, :2, :2]
    )
    return {
        "counts_xy_positive_negative_quiet": tables.tolist(),
        "scores": scores(tables),
        "capture_then_settings_then_local_records": records,
        "joint_capture_fraction": float(cohorts[-1] / shots),
        "max_exact_remote_marginal_change": float(err),
    }


def routed_source(board):
    graph = nx.Graph(board.seams)
    source_edges = [e for e in board.seams if e[0] // 12 == e[1] // 12 == 0]
    e = source_edges[0]
    f = next(f for f in source_edges if set(e).isdisjoint(f))
    initial = [e[0], f[0], e[1], f[1]]
    cell_graph = nx.Graph((u // 12, v // 12) for u, v in board.external)
    distances = nx.single_source_shortest_path_length(cell_graph, 0)
    far = max(sorted(distances), key=distances.get)
    target = list(e) + list(
        next(s for s in board.seams if s[0] // 12 == s[1] // 12 == far)
    )
    current = initial.copy()
    schedule = []
    for i, t in enumerate(target):
        path = nx.shortest_path(graph, current[i], t)
        edges = list(zip(path, path[1:]))
        for a, b in edges + edges[-2::-1]:
            schedule.append([a, b])
            current = [b if x == a else a if x == b else x for x in current]
    assert current == target
    state = np.eye(16, dtype=complex)[:, WORDS.index((1, 1, 0, 0))]
    for first, second in ((0, 2), (1, 3)):
        out = np.zeros_like(state)
        g = local_gate(np.pi / 2)
        for j, w in enumerate(WORDS):
            for i, local in enumerate(LOCAL):
                dest = list(w)
                dest[first], dest[second] = local
                out[WORDS.index(tuple(dest))] += (
                    g[i, LOCAL.index((w[first], w[second]))] * state[j]
                )
        state = out
    return np.outer(state, state.conj()), {
        "initial_rails": initial,
        "terminal_rails": target,
        "cell_distance": distances[far],
        "swap_schedule": schedule,
    }


def midpoint_publications(sector, rails, start, state, gates, rng):
    h = sector.hamiltonian(K)
    middle = free(h, start, 4)
    released = free(h, middle, 4)
    np.testing.assert_allclose(abs(np.vdot(released, state)), 1, atol=5e-12)
    gates_at_home = [e for e in sector.car.external if e[0] // 12 == 0][:2]
    tags = []
    for left, right in gates_at_home:
        a = (sector.aa == left) | (sector.bb == left)
        b = (sector.aa == right) | (sector.bb == right)
        tags.append(a.astype(int) - b.astype(int))
    rho = np.zeros((16, 16), complex)
    records = []
    components = []
    for a, b in product((-1, 0, 1), repeat=2):
        branch = np.where((tags[0] == a) & (tags[1] == b), middle, 0)
        weight = np.vdot(branch, branch).real
        if weight > 0:
            after = free(h, branch / np.sqrt(weight), 4)
            local = reduce_wings(sector, rails, after)
        else:
            local = np.zeros_like(rho)
        rho += weight * local
        components.append(local)
        records.append({"gateway_signs": [a, b], "probability": float(weight)})
    np.testing.assert_allclose(np.trace(rho), 1, atol=4e-12)
    weights = [r["probability"] for r in records]
    launches = sample(N, weights, rng)
    # Keep the signs attached to each final record table before pooling. Random
    # setting draws within each branch remain independent of its earlier marks.
    all_counts = np.zeros((2, 2, 3, 3), int)
    all_flags = np.zeros((4, 2, 2, 3, 3), int)
    for n, local, record in zip(launches, components, records, strict=True):
        record["launches"] = int(n)
        if n:
            run = bell(local, *gates, rng, shots=int(n))
            flags = np.array(
                [r["counts"] for r in run["capture_then_settings_then_local_records"]]
            )
            record["flag_xy_local_counts"] = flags.tolist()
            all_flags += flags
            all_counts += np.array(run["counts_xy_positive_negative_quiet"])
    followup_rng = np.random.default_rng(178044)
    a, b, training = axes_from_records(rho, followup_rng, shots=N)
    # The moment both signs are summed, rho is exactly the state of the
    # retained material. No midpoint outcome is used to pick a setting.
    optimized = bell(rho, a, b, followup_rng)
    return {
        "seams": gates_at_home,
        "pooled_recalibration": {
            "seed": 178044,
            "training": training,
            "bound": capture_axes_bound(training),
            "held": optimized,
        },
        "publications": records,
        "scores": scores(all_counts),
        "signed_counts": all_counts.tolist(),
        "joint_capture_fraction": float(all_flags[3].sum() / N),
    }


def run():
    verify_read()
    board = federation(0)
    rng = np.random.default_rng(178041)
    interferometer = interference(board, rng)
    sector, rails, state, changed, fixture = free_source(board)
    rho = reduce_wings(sector, rails, state)
    ga, gb, training = axes_from_records(rho, rng)
    autonomous = {
        "apparatus": fixture,
        "axis_training": training,
        "held": bell(rho, ga, gb, rng),
        "unread_rails": bell(np.diag(rho.diagonal()), ga, gb, rng),
        "exchange_only": bell(reduce_wings(sector, rails, changed), ga, gb, rng),
    }
    routed, routing = routed_source(board)
    ga = [local_gate(x) for x in (0, np.pi / 2)]
    gb = [local_gate(x) for x in (3 * np.pi / 4, 5 * np.pi / 4)]
    autonomous["all_axes_fixed_quiet_bound"] = all_axes_bound(training)
    nearby = []
    repeat_rng = np.random.default_rng(178042)
    start = np.zeros(sector.P)
    start[sector.index(*fixture["source_ports"])] = 1
    fixed_a = [axis_gate(v) for v in training["alice_axes"]]
    fixed_b = [axis_gate(v) for v in training["bob_axes"]]
    for time in (7.5, 8.5):
        wave = free(sector.hamiltonian(K), start, time)
        measurement = bell(
            reduce_wings(sector, rails, wave), fixed_a, fixed_b, repeat_rng
        )
        nearby.append({"time": time, "measurement": measurement})
    autonomous["midpoint_publication"] = midpoint_publications(
        sector, rails, start, state, (fixed_a, fixed_b), np.random.default_rng(178043)
    )
    autonomous["fixed_axes_nearby_times"] = nearby
    return {
        "seed": 178041,
        "held_launches_per_panel": N,
        "nearby_times_seed": 178042,
        "midpoint_seed": 178043,
        "interference": interferometer,
        "free_pair": autonomous,
        "routed": {
            "apparatus": routing,
            "held": bell(routed, ga, gb, rng),
            "unread_rails": bell(np.diag(routed.diagonal()), ga, gb, rng),
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(record, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(
            json.loads(json.dumps(record)),
            json.loads((HERE / "results.json").read_text()),
            "quantum",
        )
        print("Reproduced both interference controls and both Bell experiments.")
