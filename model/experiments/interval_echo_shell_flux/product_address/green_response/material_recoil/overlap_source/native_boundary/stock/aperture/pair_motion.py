"""Free passive response of two identical ordinary units to stored pair records."""

import argparse
from hashlib import sha256
import json
import math

import numpy as np
from scipy import sparse
from scipy.special import jv

from . import cold_release as cold
from . import contact_timing as timing

HERE = cold.HERE
K = 8
TIMES = (1, 4, 16, 64)
RADII = (2, 3, 4, 6)
SHOTS = 2**50
VALUES = (0, -32, -33)


def hamiltonian(board, edges):
    """Two indistinguishable singles with the same permanent kinetic cuts on both arms."""
    i, j = np.triu_indices(board.n, 1)
    index = np.full((board.n, board.n), -1, dtype=int)
    index[i, j] = index[j, i] = np.arange(len(i))
    cut = {tuple(sorted(e)) for e in edges}
    adjacency = board.adjacency().tolil()
    for u, v in cut:
        adjacency[u, v] = adjacency[v, u] = 0
    adj = adjacency.tocsr()
    marked = np.bincount(np.asarray(edges, dtype=int).ravel(), minlength=board.n)
    u, v = adj.nonzero()
    rows, cols = [], []
    all_sites = np.arange(board.n)
    for a, b in zip(u, v, strict=True):
        spectator = all_sites[(all_sites != a) & (all_sites != b)]
        rows.append(index[a, spectator])
        cols.append(index[b, spectator])
    row, col = np.concatenate(rows), np.concatenate(cols)
    exchange = sparse.csr_matrix(
        (-np.ones(len(row)) / 2, (row, col)), shape=(len(i), len(i))
    )
    degree = np.bincount(row, minlength=len(i))
    near = np.asarray(adj[i, j]).ravel()
    h = exchange + sparse.diags(
        degree / 2 - K * (marked[i] + marked[j] + 2 * near), format="csr"
    )
    assert (h - h.T).nnz == 0
    return h, [index[a, b] for a, b in edges], i, j, index


def panel(board, edges):
    h, targets, i, j, _ = hamiltonian(board, edges)
    rad = cold.distances(board, sorted({x for e in edges for x in e}))
    meters = np.array(
        [((rad[i] <= r).astype(float) + (rad[j] <= r)) / 2 for r in RADII]
    )
    return h, targets, meters


def chebyshev(h, targets, coeffs, psi, times):
    """Simultaneously propagate all sectors using a common real Chebyshev recurrence."""
    coeffs = np.asarray(coeffs, dtype=float)
    row_radius = np.asarray(abs(h).sum(axis=1)).ravel() - abs(h.diagonal())
    low, high = h.diagonal() - row_radius, h.diagonal() + row_radius
    for k, z in enumerate(targets):
        low[z] += coeffs[:, k].min()
        high[z] += coeffs[:, k].max()
    left, right = math.floor(min(low)), math.ceil(max(high))
    mid, scale = (left + right) / 2, (right - left) / 2
    A = (h - sparse.eye(h.shape[0], format="csr") * mid) / scale

    def apply(w):
        out = A @ w
        out[targets] += coeffs.T * w[targets] / scale
        return out

    times = np.asarray(times)
    orders = np.ceil(math.e * scale * times / 2 + 80).astype(int)
    # For m first omitted, |J_m(x)| <= (x/2)^m/m! for real x.
    # Continue the factorial upper series geometrically and use m! >= (m/e)^m.
    m = orders + 1
    remainder = (
        2
        * (math.e * scale * times / (2 * m)) ** m
        / (1 - scale * times / (2 * (m + 1)))
    )
    assert max(remainder) < 1e-28
    maxorder = int(orders.max())
    coefficients = 2 * jv(np.arange(maxorder + 1)[:, None], scale * times[None, :])
    previous = np.repeat(psi[:, None], len(coeffs), axis=1)
    current = apply(previous)
    result = [
        (coefficients[0, k] / 2 * previous - 1j * coefficients[1, k] * current)
        for k in range(len(times))
    ]
    for n in range(2, maxorder + 1):
        following = 2 * apply(current) - previous
        sign = (1, -1j, -1, 1j)[n % 4]
        for k in np.flatnonzero(n <= orders):
            # The absolute propagator phase is irrelevant to the stored positions.
            if abs(coefficients[n, k]) > 0:
                result[k] += (sign * coefficients[n, k]) * following
        previous, current = current, following
    return result, dict(
        spectral_bounds=[left, right],
        orders=orders.tolist(),
        truncation_bound=max(remainder).item(),
    )


def responses(board, edges, times=TIMES):
    h, z, meters = panel(board, edges)
    cases = (
        [(a,) for a in VALUES]
        if len(z) == 1
        else [(a, b) for a in VALUES for b in VALUES]
    )
    n = h.shape[0]
    psi = np.ones(n) / np.sqrt(n)
    waves, proof = chebyshev(h, z, cases, psi, times)
    values = {}
    norm_error = 0.0
    for t, w in zip(times, waves, strict=True):
        prob = abs(w) ** 2
        norm_error = max(norm_error, float(max(abs(prob.sum(axis=0) - 1))))
        data = np.r_[prob[z], meters @ prob].T
        values[str(t)] = {
            ",".join(map(str, case)): readings.tolist()
            for case, readings in zip(cases, data, strict=True)
        }
    assert norm_error < 5e-11
    # Identity on the fixed material-number sector commutes with the complete
    # generator, including each head-dependent diagonal projector.
    identity = sparse.eye(n, format="csr")
    assert (identity @ h - h @ identity).nnz == 0
    invariant = np.r_[np.ones(len(z)) / n, meters.mean(axis=1)]
    return dict(
        edges=edges,
        pair_configurations=n,
        initial_endpoint_probability=1 / n,
        blank_parked_energies=h.diagonal()[z].tolist(),
        clock=proof,
        norm_error=norm_error,
        sector_readings=values,
        invariant_readings=invariant.tolist(),
    )


def binary_responses(row):
    def conditional(a, b):
        A = [(0, 1.0)] if not a else [(0, 0.5), (-32, 0.25), (-33, 0.25)]
        B = [(0, 1.0)] if not b else [(0, 0.5), (-32, 0.25), (-33, 0.25)]
        return sum(wa * wb * np.array(row[f"{x},{y}"]) for x, wa in A for y, wb in B)

    return np.array([conditional(a, b) for a, b in ((0, 0), (0, 1), (1, 0), (1, 1))])


def count_panel(output, joint, product, rng):
    banks, decoded = {}, {}
    err = 2 * np.sqrt(timing.lab.lab.RISK / (2 * SHOTS))
    for label, experiment in output.items():
        banks[label], decoded[label] = {}, {}
        for t, row in experiment["sector_readings"].items():
            if label == "one":
                probs = {tag: np.array(row[tag]) for tag in ("0", "-32", "-33")}
                control = "0"
            else:
                cond = binary_responses(row)
                probs = {"native": joint @ cond, "independent": product @ cond}
                diff = cond[0] - cond[1] - cond[2] + cond[3]
                np.testing.assert_allclose(
                    probs["native"] - probs["independent"],
                    (joint[3] - product[3]) * diff,
                    rtol=0,
                    atol=3e-15,
                )
                control = "native"
            # An independent bank for each endpoint coincidence and enclosed single.
            counts = {tag: rng.binomial(SHOTS, p).tolist() for tag, p in probs.items()}
            banks[label][t] = counts
            bands = {}
            for tag, readings in counts.items():
                if tag == control:
                    continue
                delta = (np.array(readings) - counts[control]) / SHOTS
                target = probs[tag] - probs[control]
                assert np.max(abs(delta - target)) < err
                bands[tag] = np.column_stack((delta - err, delta + err)).tolist()
            decoded[label][t] = bands
    return banks, decoded


def run():
    from . import pair_contact

    board = cold.federation(0)
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    marks, _ = cold.read_map(board, tape["spectral_tape"])
    edges, counts, _, _ = pair_contact.source_counts()
    joint = counts.ravel() / counts.sum()
    product = np.outer(counts.sum(1), counts.sum(0)).ravel() / counts.sum() ** 2
    experiments = {"one": responses(board, [marks[0]]), "two": responses(board, edges)}
    banks, decoded = count_panel(
        experiments, joint, product, np.random.default_rng(548102)
    )
    invariant = {}
    rng = np.random.default_rng(548103)
    err = 2 * np.sqrt(timing.lab.lab.RISK / (2 * SHOTS))
    for arm, experiment in experiments.items():
        p = experiment["invariant_readings"]
        matched = rng.binomial(SHOTS, p, size=(2, len(TIMES), len(p)))
        change = (matched[1] - matched[0]) / SHOTS
        assert max(abs(change.ravel())) <= err
        invariant[arm] = dict(
            counts=matched.tolist(), known_change=0, change_radius=err
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "PAIR_CONTACT.md",
                HERE / "mediated_release.json",
                HERE / "pair_contact.json",
            )
        },
        quantities=["joint_on_each_written_seam", "random_unit_inside_2_3_4_6_hops"],
        shots=SHOTS,
        experiments=experiments,
        source_banks=banks,
        incoherent_control=invariant,
        decoded=decoded,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "pair_motion.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        timing.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["decoded"])
