"""Compare signed whole-cell head records with the exact finite Coulomb field."""

import argparse
from collections import deque
from fractions import Fraction
from hashlib import sha256
from itertools import pairwise
import json
from math import exp, factorial, log, sqrt
from pathlib import Path
from runpy import run_path

import networkx as nx
import numpy as np
from scipy import sparse

from ........executable_model import cell
from ....fan_path.fixture import FanPath, archived_fixture
from ..radial_probe.measure import LABEL, READS, late_arrival


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[7]
FRAME = HERE.parents[2] / "dressed_exchange/intact_path/face_turn/cell_frame.json"
TRIALS = 2**34
WAITS = (0, 1 / 8, 1, 4)
PULSES = 4096
UNIT = 2**-53
# Selected as nonzero support of a convex feasibility control, in reached-word
# order for the displaced and fixed arms. The solve and checks below are exact.
CONTROL_SUPPORTS = (
    (
        (93, 172, 298, 383, 1294, 1978, 2206, 2983, 4109, 4446, 4833),
        (127, 1692, 2135, 2665, 4089, 4570, 5021, 5028, 5235, 5848),
    ),
    (
        (172, 250, 428, 591, 714, 1294, 1752, 2545, 2651, 3210, 3882, 4349, 4446, 4900),
        (140, 1138, 1692, 3541, 3972, 5982, 6330),
    ),
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def rational_solve(columns, target):
    """Solve an exact full-column-rank rectangular linear system."""
    rows = [
        [Fraction(int(x)) for x in row] + [Fraction(y)]
        for row, y in zip(columns, target, strict=True)
    ]
    width = len(rows[0]) - 1
    for j in range(width):
        i = next(i for i in range(j, len(rows)) if rows[i][j])
        rows[j], rows[i] = rows[i], rows[j]
        normal = rows[j][j]
        rows[j] = [x / normal for x in rows[j]]
        for i, row in enumerate(rows):
            coeff = row[j]
            if i != j and coeff:
                rows[i] = [x - coeff * y for x, y in zip(row, rows[j], strict=True)]
    assert all(not any(row) for row in rows[width:])
    return [row[-1] for row in rows[:width]]


def is_positive(rows):
    """Exact symmetric elimination by positive pivots."""
    data = [list(map(Fraction, row)) for row in rows]
    for j in range(len(data)):
        diagonal = data[j][j]
        assert diagonal > 0
        for i in range(j + 1, len(data)):
            weight = data[i][j] / diagonal
            for k in range(j + 1, len(data)):
                data[i][k] -= weight * data[j][k]


def cell_geometry(shape):
    """Fix orientations and prove a norm bound for the complete cycle meter."""
    edges = shape.seams
    edge_index = {edge: j for j, edge in enumerate(edges)}
    inc = np.zeros((12, 30), dtype=int)
    for i, (a, b) in enumerate(edges):
        inc[a, i], inc[b, i] = -1, 1
    faces = shape.faces()
    triangle = np.zeros((20, 30), dtype=int)
    for i, f in enumerate(faces):
        for a, b in pairwise((*f, f[0])):
            triangle[i, edge_index[tuple(sorted((a, b)))]] = 1 if a < b else -1
    face_sign = {0: 1}
    next_faces = [0]
    for i in next_faces:
        for e in np.flatnonzero(triangle[i]):
            j = next(j for j in np.flatnonzero(triangle[:, e]) if j != i)
            sign = -face_sign[i] * triangle[i, e] * triangle[j, e]
            if j not in face_sign:
                face_sign[j] = sign
                next_faces.append(j)
            assert face_sign[j] == sign
    assert len(face_sign) == 20
    triangle *= [[face_sign[i]] for i in range(20)]
    np.testing.assert_array_equal(triangle.sum(axis=0), 0)
    np.testing.assert_array_equal(inc @ triangle.T, 0)
    # 20*(C C^T + 11^T/20 - 3I/4) is exact integer SPD. The rows of C
    # therefore span the 19 cycle directions and control their norm.
    is_positive(
        (20 * (triangle @ triangle.T) + 1 - 15 * np.eye(20, dtype=int)).tolist()
    )
    return inc, triangle, edge_index


def green_table(shape, inc):
    """Identify the physical seam graph with the exact pinned theorem table."""
    api = run_path(str(REPO / "code/electromagnetism/discrete_coulomb_green.py"))
    carrier, theorem = [api["load_pinned_source"](x) for x in api["SOURCE_CONTRACT"]]
    left, right = api["seam_table"](carrier)
    green = api["parse_int_matrix"](theorem, "greenNumZ", "theorem incidence_product")
    physical = nx.Graph(shape.seams)
    formal = nx.Graph(zip(left, right, strict=True))
    permutations = [
        tuple(m[i] for i in range(12))
        for m in nx.algorithms.isomorphism.GraphMatcher(
            physical, formal
        ).isomorphisms_iter()
    ]
    assert len(permutations) == 120
    mapping = min(permutations)
    g = np.array(green, dtype=int)[mapping, :][:, mapping]
    np.testing.assert_array_equal((inc @ inc.T) @ g, 180 * np.eye(12, dtype=int) - 15)
    np.testing.assert_array_equal(g, g.T)
    np.testing.assert_array_equal(g.sum(axis=0), 0)
    return g, mapping


class Component:
    def __init__(self, start, triangle):
        flips = [
            (
                sum(1 << int(j) for j in np.flatnonzero(c)),
                sum(1 << int(j) for j in np.flatnonzero(c == 1)),
            )
            for c in triangle
        ]
        self.flips = flips
        reached = {start[0]: None}
        pending = deque([start[0]])
        graph = []
        while pending:
            w = pending.popleft()
            for i, (mask, arrows) in enumerate(flips):
                if w & mask in (arrows, arrows ^ mask):
                    target = w ^ mask
                    graph.append((w, target))
                    if target not in reached:
                        reached[target] = (w, i)
                        pending.append(target)
        self.parents = reached
        self.words = sorted(reached)
        index = {w: j for j, w in enumerate(self.words)}
        self.start = [index[w] for w in start]
        left, right = np.array([(index[w], index[v]) for w, v in graph]).T
        self.heads = 1 - ((np.array(self.words)[:, None] >> np.arange(30)) & 1)
        off = sparse.coo_matrix(
            (np.ones(len(left)), (left, right)), shape=(len(index),) * 2
        ).tocsr()
        self.degree = degree = np.asarray(off.sum(axis=1)).ravel()
        assert (off != off.T).nnz == 0 and int(off.max()) == 1
        self.h = (sparse.diags(degree, format="csr") - off) / 2
        self.mix = sparse.eye(len(index), format="csr") - self.h / 16
        assert self.mix.min() == 0

    def path(self, i):
        word = self.words[i]
        route = []
        while (preceding := self.parents[word]) is not None:
            word, face = preceding
            route.append(face)
        return route[::-1]

    def apply(self, word, route):
        for face in route:
            mask, arrows = self.flips[face]
            assert word & mask in (arrows, arrows ^ mask)
            word ^= mask
        return word

    def probability_banks(self, density):
        n = len(self.words)
        waves = np.zeros((n, 2), dtype=complex)
        waves[self.start] = np.eye(2)
        banks = [np.einsum("xi,ij,xj->x", waves, density, waves.conj()).real]
        scaled = -1j * self.h / 32
        for a, b in pairwise(WAITS):
            for _ in range(int(32 * (b - a))):
                term = waves
                total = waves.copy()
                for k in range(1, 19):
                    term = (scaled @ term) / k
                    total += term
                waves = total
            banks.append(np.einsum("xi,ij,xj->x", waves, density, waves.conj()).real)
        # Each Taylor step has norm <= exp(b). Bound all arithmetic along
        # any summand by 4096u and the omitted terms by the norm tail.
        # Errors are composed with the exact unitary of norm one.
        b = self.degree.max() / 32
        assert b <= 20 / 32
        guard = 4096 * UNIT
        step = (guard + b**19 / factorial(19)) * exp(b)
        wave_error = exp(128 * step) * 128 * step
        assert wave_error < 2**-30
        banks.append(banks[0].copy())
        mixed = banks[0]
        for _ in range(PULSES):
            mixed = self.mix @ mixed
        assert exp(PULSES * 64 * UNIT) * (PULSES * 64 * UNIT) < 2**-30
        banks.append(mixed)
        return [np.maximum(x, 0) for x in banks]


def specified_mixture(reference, moved, support, canonical):
    j, k = support
    columns = np.vstack(
        [
            np.concatenate(
                [moved.heads[list(j)].T, -reference.heads[list(k)].T], axis=1
            ),
            [1] * len(j) + [0] * len(k),
            [0] * len(j) + [1] * len(k),
        ]
    )
    coefficients = rational_solve(columns, [*canonical, 1, 1])
    assert all(x > 0 for x in coefficients)
    profiles = []
    for indices, weights, part in (
        (j, coefficients[: len(j)], moved),
        (k, coefficients[len(j) :], reference),
    ):
        assert sum(weights) == 1
        reset = part.path(part.start[1])[::-1]
        root = part.words[part.start[0]]
        assert part.apply(part.words[part.start[1]], reset) == root
        paths = [
            [part.words[i], str(w), part.path(i)]
            for i, w in zip(indices, weights, strict=True)
        ]
        assert all(
            part.apply(root, path) == destination for destination, _, path in paths
        )
        population = np.zeros(len(part.words))
        population[list(indices)] = [float(x) for x in weights]
        profiles.append(
            (
                population,
                {
                    "second_admitted_word_to_first": reset,
                    "target_word_probability_and_route": paths,
                },
            )
        )
    return profiles


def sample(rng, probabilities, component, acceptance, sign):
    assert abs(float(probabilities.sum()) - 1) <= 2e-9
    sample_mass = acceptance * probabilities / probabilities.sum()
    counts = rng.multinomial(TRIALS, [*sample_mass, 1 - acceptance])
    n = int(counts[:-1].sum())
    assert n + counts[-1] == TRIALS and n > 0
    heads = component.heads if sign == 1 else 1 - component.heads
    total_heads = counts[:-1] @ heads
    complemented = (1 << 30) - 1
    words = (
        component.words if sign == 1 else [w ^ complemented for w in component.words]
    )
    return {
        "words_and_counts": [
            [w, int(k)] for w, k in zip(words, counts[:-1], strict=True) if k > 0
        ],
        "seam_head_count_toward_larger_port": total_heads.tolist(),
        "face_admission": n,
        "face_refusal": int(counts[-1]),
    }


def fit_difference(moved, reference, cycles, margin, tested_banks):
    vectors = [
        np.array(x["seam_head_count_toward_larger_port"]) / x["face_admission"]
        for x in (moved, reference)
    ]
    field = vectors[0] - vectors[1]
    width = (
        sum(
            sqrt(log(2 * tested_banks / 0.01) / (2 * x["face_admission"]))
            for x in (moved, reference)
        )
        + margin
    )
    return {
        "observed_mean_seam_difference": field.tolist(),
        "observed_face_circulation_difference": (cycles @ field).tolist(),
        "simultaneous_99pct_seam_half_width": width,
        "simultaneous_99pct_circulation_half_width": 3 * width,
    }


def arrival_input(positive, negative, frame):
    data = json.loads(READS.read_text())
    flight = next(
        f
        for case in data["cases"]
        if case["fan_count"] == 4
        for s in case["stock_arms"]
        if s["preparation"] == "clock_bypassed"
        for f in s["flights"]
        if f["material_bridge_rate"] == 4 and f["appointment"] == 30.8125
    )
    receipt = next(
        w for w in flight["pair_readings"] if w["port_face_and_words"] == LABEL
    )
    raw = []
    for i, path in enumerate((positive, negative)):
        assert path.cells[-1][0] == frame["fixed_charge_port"] == 1045
        words = frame["proper_pose_checkpoints"][i][
            "head_words_at_five_successive_axial_poses"
        ][0]
        fly_edges = [
            (path.cells[-1][a] % 12, path.cells[-1][b] % 12)
            for a, b in path.field.edges
        ]
        chart = []
        for w in words:
            chart.append(
                sum(
                    1 << j
                    for j, (a, b) in enumerate(fly_edges)
                    if (a < b)
                    == bool(w & (1 << cell().seams.index(tuple(sorted((a, b))))))
                )
            )
        assert chart == (LABEL[-2:] if i == 0 else [177, 146])
        states, admissions, *_ = late_arrival(path, [*LABEL[:2], *chart], receipt)
        acceptance = float(sum(admissions))
        density = states[0].sum(axis=0) / acceptance
        raw.append((words, acceptance, density))
    a, b = raw
    assert [w ^ ((1 << 30) - 1) for w in a[0]] == b[0]
    np.testing.assert_allclose(a[1], b[1], atol=2e-12, rtol=0)
    np.testing.assert_allclose(a[2], b[2], atol=2e-11, rtol=0)
    source_error = 4 * (flight["observable_absolute_error_upper_per_pair"] + 2e-12)
    assert a[1] - source_error > 0.015
    return a, source_error


def run():
    rng = np.random.default_rng(10451052)
    frame = json.loads(FRAME.read_text())
    shape = cell()
    inc, cycles, idx = cell_geometry(shape)
    g, mapping = green_table(shape, inc)
    with archived_fixture() as archive:
        pos, neg = [FanPath(archive, 4, sign, asymmetric=True) for sign in (1, -1)]
        (words, acceptance, density), source_error = arrival_input(pos, neg, frame)
        offset = 1044
        assert all((offset + a, offset + b) in pos.board.seams for a, b in shape.seams)
    origin = 1045 - offset
    base = Component(words, cycles)
    base_bank = base.probability_banks(density)
    moves = []
    for (a, b), j in idx.items():
        if origin not in (a, b):
            continue
        if all(bool(w & (1 << j)) == (origin == a) for w in words):
            target = b if origin == a else a
            moves.append((target, j, 1 if origin == a else -1))
    assert [offset + t for t, *_ in moves] == [1051, 1052]
    records = []
    false_means = float("inf")
    protocol_names = ["simultaneous_free_" + str(t) for t in WAITS] + [
        "all_faces_shuttered_4",
        "equal_face_pulses_4096",
        "programmed_Green_control",
    ]
    # Dependence across different face sums from the same multihead read is
    # permitted by the union bound. Each protocol has two separate banks.
    tests = 2 * len(moves) * len(protocol_names) * 2 * 50
    margin = 4 * source_error / (acceptance - source_error) + 2**-26
    for (target, edge_id, direction), support in zip(
        moves, CONTROL_SUPPORTS, strict=True
    ):
        shifted_words = [w ^ (1 << edge_id) for w in words]
        moved = Component(shifted_words, cycles)
        load = direction * inc[:, edge_id]
        np.testing.assert_array_equal(base.heads @ inc.T - base.heads[0] @ inc.T, 0)
        np.testing.assert_array_equal(
            moved.heads @ inc.T - base.heads[0] @ inc.T,
            np.broadcast_to(load, (len(moved.words), 12)),
        )
        green_numerator = inc.T @ g @ load
        np.testing.assert_array_equal(cycles @ green_numerator, 0)
        np.testing.assert_array_equal(inc @ green_numerator, 180 * load)
        canonical = [Fraction(int(x), 180) for x in green_numerator]
        assert canonical[edge_id] == direction * Fraction(11, 30)
        uniform = [
            Fraction(int(a), len(moved.words)) - Fraction(int(b), len(base.words))
            for a, b in zip(
                moved.heads.sum(axis=0), base.heads.sum(axis=0), strict=True
            )
        ]
        uniform_circulations = [
            sum(int(k) * x for k, x in zip(row, uniform, strict=True)) for row in cycles
        ]
        assert max(map(abs, uniform_circulations)) > Fraction(29, 1000)
        (tailored_moved, mp), (tailored_ref, rp) = specified_mixture(
            base, moved, support, canonical
        )
        moved_banks = moved.probability_banks(density) + [tailored_moved]
        base_banks = base_bank + [tailored_ref]
        comparisons = []
        for sign in (1, -1):
            banks = []
            for name, moved_prob, ref_prob in zip(
                protocol_names, moved_banks, base_banks, strict=True
            ):
                r0 = sample(rng, ref_prob, base, acceptance, sign)
                r1 = sample(rng, moved_prob, moved, acceptance, sign)
                for r, p in ((r0, origin), (r1, target)):
                    r["port_and_local_signed_load_counts"] = [
                        [offset + v, sign if v == p else 0, r["face_admission"]]
                        for v in range(12)
                    ]
                comparison = fit_difference(r1, r0, cycles, margin, tests)
                np.testing.assert_allclose(
                    inc @ comparison["observed_mean_seam_difference"],
                    sign * load,
                    atol=1e-14,
                    rtol=0,
                )
                m = max(map(abs, comparison["observed_face_circulation_difference"]))
                pad = comparison["simultaneous_99pct_circulation_half_width"]
                if name == "programmed_Green_control":
                    assert m <= pad
                else:
                    assert m > pad
                    false_means = min(false_means, m - pad)
                banks.append(
                    {"arm": name, "fixed": r0, "moved": r1, "comparison": comparison}
                )
            comparisons.append({"signed_load": sign, "banks": banks})
        records.append(
            {
                "moved_occupied_port": offset + target,
                "fixed_and_moved_ring_component_sizes": [
                    len(base.words),
                    len(moved.words),
                ],
                "positive_arrival_port_load_difference": load.tolist(),
                "positive_Green_seam_target": [str(x) for x in canonical],
                "positive_exact_uniform_seam_difference": [str(x) for x in uniform],
                "positive_exact_uniform_face_difference": [
                    str(x) for x in uniform_circulations
                ],
                "maximum_uniform_seam_error": str(
                    max(abs(x - y) for x, y in zip(uniform, canonical, strict=True))
                ),
                "specified_Green_control_moved_and_reference_word_probability_and_path": [
                    mp,
                    rp,
                ],
                "comparisons": comparisons,
            }
        )
    return {
        "protocol_sha256": digest(HERE / "PROTOCOL.md"),
        "parent_word_reads_sha256": digest(READS),
        "full_cell_frame_sha256": digest(FRAME),
        "array_port_order": list(range(offset, offset + 12)),
        "Lean_port_index_by_physical_port": mapping,
        "seams_in_physical_orientation": [
            [a + offset, b + offset] for a, b in shape.seams
        ],
        "triangles_in_physical_order": [[v + offset for v in f] for f in shape.faces()],
        "oriented_triangle_matrix": cycles.tolist(),
        "rigorous_positive_cycle_eigenvalue_lower": 0.75,
        "acceptance": acceptance,
        "original_launches_per_bank": TRIALS,
        "number_of_simultaneous_mean_bounds": tests,
        "guaranteed_minimum_rejected_face_circulation_among_native_waits_and_scheduled_mixes": false_means,
        "two_route_comparisons": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(results, separators=(",", ":")) + "\n")
    else:
        compare = run_path(str(HERE.parents[2] / "replay.py"))["compare"]
        compare(results, json.loads((HERE / "results.json").read_text()), "cycle_flux")
        print(
            "Verified signed Gauss, triangle reads, and exact Green passing controls."
        )


if __name__ == "__main__":
    main()
