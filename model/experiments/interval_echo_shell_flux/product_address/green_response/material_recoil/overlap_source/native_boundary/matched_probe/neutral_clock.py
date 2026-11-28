"""An independent native near/far clock, calibrated with the material probe installed."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm

from . import measure as m


def clock_register():
    law = m.source.LAW

    def ids(entries):
        return np.array([0] + [entries.get((c, d), 0) for c, d in law.LABELS], int)

    # Exact, unnormalized columns for the existing native-emission Ramsey stock.
    zero = ids({(0, 0): 1, (1, 0): -1})
    even = ids({(c, d): (-1) ** c for c in (0, 1) for d in (-2, 2)})
    odd = ids({(c, d): (-1) ** c * (d // 2) for c in (0, 1) for d in (-2, 2)})
    clock = np.stack([zero, even, odd], axis=1)
    energy = {}
    for k in (0, 8):
        full = m.source.formation.emission.four_h(k)
        vals = [0, 4 * k, 4 * k + 1]
        embedded = np.kron(clock, np.eye(9, dtype=int))
        np.testing.assert_array_equal(
            full @ embedded, embedded * np.repeat(4 * np.array(vals), 9)
        )
        energy[str(k)] = vals
    # Even an arbitrary magnitude-two memory separates in the original
    # one-binary-unit sector: its writable endpoint pair (1,1) is unavailable.
    native = m.source.formation.emission.four_h(8)
    tags = [1 + i for i, (c, d) in enumerate(law.LABELS) if abs(d) == 2]
    pairs = [(0, 0), (1, 0), (0, 1)]
    chosen = [9 * i + law.SRC.index(w) for i in tags for w in pairs]
    bare_dirs = native[np.ix_(chosen[::3], chosen[::3])]
    np.testing.assert_array_equal(
        native[np.ix_(chosen, chosen)], np.kron(bare_dirs, np.eye(3, dtype=int))
    )
    assert not np.any(native[np.ix_(np.setdiff1d(np.arange(99), chosen), chosen)])
    normalized = clock / np.linalg.norm(clock, axis=0)
    ramsey = (normalized[:, 0] + normalized[:, 1]) / np.sqrt(2)
    sharp = (normalized[:, 1] + normalized[:, 2]) / np.sqrt(2)
    np.testing.assert_allclose(
        sharp[1:][[d != 2 for c, d in law.LABELS]], 0, atol=1e-15
    )
    reference = m.source.formation.emission.memory_vector
    lo = reference([(0, 0, 1), (1, 0, -1)])
    hi = reference([(0, 2, 1), (0, -2, 1), (1, 2, -1), (1, -2, -1)])
    np.testing.assert_allclose((lo + hi) / np.sqrt(2), ramsey, atol=1e-15)
    times = [0, np.pi / 32, np.pi / 16, 1, np.pi, 4, 2 * np.pi, 16, 64]
    records = {}
    for k, vals in energy.items():
        u = np.stack([expm(-1j * t * np.diag(vals)) for t in times])

        def pointer(basis):
            v = normalized.T @ basis
            return np.abs(np.einsum("i,tij,j->t", v, u, v)) ** 2

        fast = pointer(ramsey)
        slow = pointer(sharp)
        np.testing.assert_allclose(
            fast, np.cos(2 * int(k) * np.array(times)) ** 2, atol=2e-13
        )
        np.testing.assert_allclose(slow, np.cos(np.array(times) / 2) ** 2, atol=5e-13)
        records[k] = dict(ramsey_yes=fast.tolist(), sharp_sign_same=slow.tolist())
    return dict(
        register_levels=energy,
        verified_native_columns=27,
        all_magnitude_two_memories_in_single_source_columns=12,
        relative_near_far_rate=[1, 1],
        source_induced_frequency_change=[0, 0],
        times=times,
        hand_calibration=records,
    )


def clock_edges(board, seam, dist, origins):
    edges = []
    for gateway, far in zip(seam, origins, strict=True):
        near = min(edge for edge in board.seams if gateway in edge and edge != seam)
        remote = min(
            edge
            for edge in board.seams
            if far in edge and min(dist[list(edge)]) == dist[far] - 1
        )
        edges.extend([near, remote])
    assert len(set(edges)) == 4
    return edges


def material_h(board, seam, total, edges):
    h, words, bits = m.fp.generator(board, seam, total, 8, True)
    # The clock term factorizes. Subtract the previously blank ordinary bond
    # on its physical address; degree completion belongs to the original ports.
    ix = {w: i for i, w in enumerate(words)}
    blank = sparse.lil_matrix((len(words), len(words)))
    for a, b in edges:
        for i in np.flatnonzero(bits[:, a] != bits[:, b]):
            j = ix[tuple(sorted(set(words[i]) ^ {a, b}))]
            blank[i, i] += 8.5
            blank[i, j] -= 0.5
    reduced = h - sparse.kron(blank, sparse.eye(4), format="csr")
    assert (reduced != reduced.T).nnz == 0
    return reduced, words, bits


def material():
    board, seam, dist, origins = m.fp.make_board()
    edges = clock_edges(board, seam, dist, origins)
    preps = m.preparations(seam, origins)
    radial = {}
    baselines = {}
    for n in (0, 1, 2):
        h, words, bits = material_h(board, seam, n, edges)
        for (source, extra), keys in preps.items():
            if source + int(extra is not None) != n:
                continue
            indices = [4 * words.index(w) + c for w, c in keys]
            ops = m.read_operators(h, bits, indices)
            radial[source, extra] = np.einsum("tpmk,p->tmk", ops["number"], dist)
            baselines[source, extra] = h[indices, :][:, indices].toarray()
    results = {}

    def read(arms, extra, ops):
        return sum(m.trace(r, ops[n, extra]) for n, r in enumerate(arms))

    arms = m.stocks()
    low, high = arms["balanced_mass"]
    # The 19/15 gain rescales only the probability of the complete stock;
    # failed launches stay in the high cohort.
    admitted = high[1] * 19
    raw_gain = float(
        read([np.zeros((2, 2)), admitted], None, baselines) - read(low, None, baselines)
    )
    probability = (17 / 8) / raw_gain
    arms["energy_rematched_mass"] = [
        low,
        [low[0] * (1 - probability), admitted * probability],
    ]
    for name in (
        "target",
        "phase_feedback",
        "exchange_1_1",
        "balanced_mass",
        "energy_rematched_mass",
    ):
        lo, hi = arms[name]
        none = read(hi, None, radial) - read(lo, None, radial)
        gain = float(read(hi, None, baselines) - read(lo, None, baselines))
        rows = []
        for p in origins:
            x = read(hi, p, radial) - read(lo, p, radial) - none
            rows.append(dict(port=p, radial=x.tolist()))
        results[name] = dict(initial_material_energy_increment=gain, launches=rows)
    return dict(
        seam=seam,
        ordinary_raw_increment=raw_gain,
        ordinary_rematched_probability=probability,
        clock_edges_near_far_by_owner=edges,
        times=[*m.fp.TIMES, "infinite_time_mean"],
        standards=results,
    )


def readings(clock, material):
    rng = np.random.default_rng(492105)
    shots = 2**20
    stocks = ("paired", "erased", "phase_feedback", "energy_rematched_mass")
    records = {}
    family = len(stocks) * 4 * 2 * len(clock["times"])
    for stock in stocks:
        rows = []
        for edge in material["clock_edges_near_far_by_owner"]:
            counts = {}
            for key, probs in clock["hand_calibration"]["8"].items():
                yes = rng.binomial(shots, np.clip(probs, 0, 1))
                counts[key] = np.column_stack([yes, shots - yes]).tolist()
            rows.append(dict(edge=edge, counts=counts))
        records[stock] = rows
    return dict(
        shots=shots,
        seed=492105,
        family=family,
        individual_bin_simultaneous_99_half_width=float(
            np.sqrt(np.log(2 * family / 0.01) / (2 * shots))
        ),
        records=records,
    )


def run():
    clock = clock_register()
    board = material()
    return dict(
        protocol_sha256=sha256((m.HERE / "NEUTRAL_CLOCK.md").read_bytes()).hexdigest(),
        clock=clock,
        material=board,
        pointer_counts=readings(clock, board),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    path = m.HERE / "neutral_clock.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        m.fp.check(json.loads(path.read_text()), json.loads(json.dumps(result)))
    for name, item in result["material"]["standards"].items():
        print(name, item)
    print(
        result["material"]["clock_edges_near_far_by_owner"],
        result["clock"]["hand_calibration"]["8"],
    )
