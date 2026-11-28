"""Cold release and subthreshold binding from a source reader's contact map."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.candidate_laws.interval_echo.openwave.photon_identity.native_emission import (
    measure as emission,
)
from . import measure as apparatus
from . import source_reader

HERE = Path(__file__).resolve().parent
TIMES = (1 / 32, 1 / 8, 0.5, 1.0, 4.0)
STARTS = 2**30
SEED = 501927


def check_contact(k, load, magnitude):
    """Raw integer columns on vacuum, left/right launch, every enrolled head."""
    law = emission.LAW
    four = emission.four_h(k)
    words = [law.SRC.index(v) for v in ((0, 0), (load, 0), (0, load))]
    labels = [i + 1 for i, (_, d) in enumerate(law.LABELS) if abs(d) == magnitude]
    vacuum = 9 * np.array(labels)
    head = four[np.ix_(vacuum, vacuum)]
    for word in words:
        cols = vacuum + word
        target = np.zeros((99, len(cols)), dtype=np.int64)
        target[cols] = head
        np.testing.assert_array_equal(four[:, cols], target)
    # A real write on the specified incident loads reaches this word magnitude.
    before = law.SRC.index((magnitude, 0))
    support = np.flatnonzero(law._write_pattern()[:, before])
    assert all(abs(law.LABELS[i // 9][1]) == magnitude for i in support)
    # On a blank contact the single travelling load has its native hopping.
    cols = words[1:]
    target = np.zeros((99, 2), dtype=np.int64)
    target[cols] = np.array([[2, -2], [-2, 2]]) + 4 * k * load**2 * np.eye(2, dtype=int)
    np.testing.assert_array_equal(four[:, cols], target)


def source_tape(board, path, magnitude, k, rng):
    seam_index = {e: i for i, e in enumerate(board.seams)}
    edges = sorted(
        seam_index[tuple(sorted(e))] for e in zip(path[:-1], path[1:], strict=True)
    )
    # A sharp native write after evacuation has equal reversal parities.
    law = emission.LAW
    state = np.eye(11)[1 + law.LABELS.index((0, magnitude))]
    spectrum = np.einsum("i,kij,j->k", state, source_reader.spectrometer(), state)
    ref_prob, _ = source_reader.probabilities(k)
    ref = rng.multinomial(2**26, ref_prob)
    per_edge = 2**14
    counts = [
        [i, rng.multinomial(per_edge, spectrum / spectrum.sum()).tolist()]
        for i in edges
    ]
    return dict(
        seams=len(board.seams),
        blank_outcomes_per_other_seam=per_edge,
        reference=ref.tolist(),
        written=counts,
    )


def read_map(board, tape):
    ref = np.array(tape["reference"]) / sum(tape["reference"])
    ref_error = np.sqrt(np.log(2 * 10000 / 0.01) / (2 * sum(tape["reference"])))
    gains = np.log((ref[0] / source_reader.MULT[0]) / (ref / source_reader.MULT))
    low = np.log(
        ((ref[0] - ref_error) / source_reader.MULT[0])
        / ((ref + ref_error) / source_reader.MULT)
    )
    high = np.log(
        ((ref[0] + ref_error) / source_reader.MULT[0])
        / ((ref - ref_error) / source_reader.MULT)
    )
    assert tape["seams"] == len(board.seams)
    mask = []
    count = np.zeros(6, dtype=np.int64)
    for i, row in tape["written"]:
        assert (
            len(row) == 6
            and all(type(x) is int and x >= 0 for x in row)
            and 0 <= i < tape["seams"]
            and sum(row) == tape["blank_outcomes_per_other_seam"]
        )
        assert row[0] == 0 and row[1] == 0
        mask.append(board.seams[i])
        count += row
    assert mask and len(mask) == len(set(mask))
    pair = [2, 3] if count[2:4].sum() else [4, 5]
    assert count[pair].sum() == count.sum()
    odd = count[pair[1]] / count.sum()
    err = np.sqrt(np.log(2 * 10000 / 0.01) / (2 * count.sum()))
    interval = [
        min((1 - f) * low[pair[0]] + f * low[pair[1]] for f in (odd - err, odd + err)),
        max(
            (1 - f) * high[pair[0]] + f * high[pair[1]] for f in (odd - err, odd + err)
        ),
    ]
    return mask, dict(
        per_seam=count @ gains / count.sum(),
        per_seam_interval=interval,
        seams=len(mask),
        total=count @ gains / tape["blank_outcomes_per_other_seam"],
    )


def motion(board, mask, load, k):
    mask = set(mask)
    incidence = np.bincount(
        np.asarray(sorted(mask), dtype=int).ravel(), minlength=board.n
    )
    u, v = np.asarray([e for e in board.seams if e not in mask]).T
    adj = sparse.csr_matrix(
        (np.ones(2 * len(u)), (np.r_[u, v], np.r_[v, u])), shape=(board.n, board.n)
    )
    deg = np.asarray(adj.sum(axis=1)).ravel()
    lap = sparse.diags(deg, format="csr") - adj
    h = lap / 2 - sparse.diags(k * load**2 * incidence, dtype=float, format="csr")
    return h.tocsr(), lap, incidence


def distances(board, sites):
    d = nx.multi_source_dijkstra_path_length(nx.Graph(board.seams), list(sites))
    return np.array([d[i] for i in range(board.n)])


def radial_law(h, psi, dist):
    result = []
    for t in TIMES:
        wave = expm_multiply(-1j * t * h, psi)
        p = np.bincount(dist, weights=abs(wave) ** 2)
        assert abs(p.sum() - 1) < 2e-11
        result.append(p / p.sum())
    return result


def observe(h, psi, dist, rng):
    rows = radial_law(h, psi, dist)
    return [rng.multinomial(STARTS, p).tolist() for p in rows]


def release(board, mask, q, k, rng):
    h, lap, m = motion(board, mask, q, k)
    dist = distances(board, np.flatnonzero(m))
    n = board.n
    psi = np.ones(n) / np.sqrt(n)
    law = k * q**2 * (lap @ m) / n
    hp = h @ psi
    direct = 2 * (hp**2 - psi * (h @ hp))
    np.testing.assert_allclose(law, direct, rtol=1e-10, atol=1e-13)
    shell = np.bincount(dist, weights=law)
    # Positive current across each surface points toward the source.
    flux = np.cumsum(shell)
    assert np.max(abs(flux[1:])) < 1e-12
    cold = observe(h, psi, dist, rng)
    control = observe(motion(board, [], q, k)[0], psi, dist, rng)
    result = dict(
        mean_distance_acceleration=float(dist @ law),
        inward_current_derivative=flux.tolist(),
        initial_shell_counts=np.bincount(dist).tolist(),
        cold=cold,
        blank=control,
    )
    if k:
        # Read the exchange energy of a localized load on each live blank
        # seam incident on the chosen port. Private and contrast terms are sharp.
        site = int(np.argmax(m))
        free_degree = int(board.degree[site] - m[site])
        exchanges = rng.binomial(2**18, 0.5, free_degree)
        cost = exchanges.sum() / 2**18
        err = np.sqrt(np.log(2 * 10000 / 0.01) * free_degree / (2 * 2**18))
        depth = k * q**2 * m[site]
        pos = np.eye(1, n, site).ravel()
        assert abs(h[site, site] - (free_degree / 2 - depth)) < 1e-12
        bound = (depth - cost - err) / depth
        radial = radial_law(h, pos, dist)
        assert all(p[0] >= 1 - free_degree / (2 * depth) - 1e-13 for p in radial)
        result["localized"] = dict(
            port=site,
            incidence=int(m[site]),
            blank_degree=free_degree,
            exchange_counts=exchanges.tolist(),
            shots_per_edge=2**18,
            insertion_below_free_band=[depth - cost - err, depth - cost + err],
            all_time_source_probability_lower=bound,
            shots_per_time=STARTS,
            source_counts=[rng.binomial(STARTS, min(1, p[0])) for p in radial],
        )
    return result


def probability_band(row, n):
    p = np.asarray(row) / n
    a = np.log(2 * 10000 / 0.01) / n
    width = np.sqrt(2 * a * (2 * p + 4 * a)) + 2 * a / 3 + 1e-10
    return np.stack((np.maximum(0, p - width), np.minimum(1, p + width)), axis=-1)


def difference(a, b):
    return np.stack((a[..., 0] - b[..., 1], a[..., 1] - b[..., 0]), axis=-1)


def decode_release_counts(record):
    """Use the public contact chart, prepared load, stiffness and count banks."""
    board = federation(1)
    shots = record["starts_per_time"]
    result = {}
    for name, case in record["rows"].items():
        stiffness = int(name.rsplit("_", 1)[1])
        marked, reading = read_map(board, case["tape"])
        incidence = np.bincount(np.asarray(marked).ravel(), minlength=board.n)
        # The acceleration prediction uses only the source readings.
        boundary = sum(
            int(incidence[a] + incidence[b])
            for a, b in board.seams
            if bool(incidence[a]) != bool(incidence[b])
        )
        result[name] = {}
        for load, probe in case["probes"].items():
            gain = stiffness * int(load) ** 2
            accel = -gain * boundary / board.n
            # Read the independent position bank only after computing the prediction.
            for records in (probe["cold"], probe["blank"]):
                assert len(records) == len(record["times"])
                assert all(min(row) >= 0 and sum(row) == shots for row in records)
            support = difference(
                probability_band(probe["cold"], shots),
                probability_band(probe["blank"], shots),
            )[:, 0]
            out = dict(
                marked_contacts=reading["seams"],
                boundary_weight=boundary,
                predicted_acceleration=accel,
                support_excess=support.tolist(),
            )
            if stiffness:
                read = probe["localized"]
                inc = int(incidence[read["port"]])
                exchanges = read["exchange_counts"]
                n = read["shots_per_edge"]
                assert len(exchanges) == board.degree[read["port"]] - inc
                assert all(0 <= x <= n for x in exchanges)
                error = np.sqrt(np.log(2 * 10000 / 0.01) * len(exchanges) / (2 * n))
                deficit = gain * inc - sum(exchanges) / n - error
                out["all_time_support_lower"] = deficit / (gain * max(incidence))
            result[name][load] = out
    return result


def run():
    archive = HERE / "results.json"
    original = next(
        b for b in json.loads(archive.read_bytes())["boards"] if b["owners"] == 80
    )
    board = federation(1)
    rng = np.random.default_rng(SEED)
    rows = {}
    for k in (8, 0):
        for q, mag in ((2, 1), (2, 2), (1, 2)):
            check_contact(k, q, mag)
    for name, key in (("outward", "geodesic"), ("folded", "alternative")):
        path = original[key]["path"]
        for mag, k in ((1, 8), (2, 8), (1, 0)):
            tape = source_tape(board, path, mag, k, rng)
            mask, reading = read_map(board, tape)
            probes = {
                str(q): release(board, mask, q, k, rng)
                for q in ((1, 2) if mag == 2 else (2,))
            }
            rows[f"{name}_{mag}_{k}"] = dict(
                tape=tape, modular_read=reading, probes=probes
            )
    return dict(
        protocol_sha256=sha256((HERE / "COLD_RELEASE.md").read_bytes()).hexdigest(),
        path_sha256=sha256(archive.read_bytes()).hexdigest(),
        seed=SEED,
        starts_per_time=STARTS,
        times=TIMES,
        rows=rows,
    )


if __name__ == "__main__":
    cmd = argparse.ArgumentParser(description=__doc__)
    modes = cmd.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cmd.parse_args()
    path = HERE / "cold_release.json"
    if args.counts_only:
        output = json.loads(path.read_bytes())
        apparatus.check(output["decoded"], decode_release_counts(output))
    else:
        output = run()
        output["decoded"] = decode_release_counts(output)
        if args.write:
            path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        else:
            apparatus.check(
                json.loads(path.read_bytes()), json.loads(json.dumps(output))
            )
    for name, row in output["decoded"].items():
        print(name, row["2"]["predicted_acceleration"], row["2"]["support_excess"][0])
