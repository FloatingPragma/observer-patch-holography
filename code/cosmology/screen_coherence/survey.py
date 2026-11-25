"""Survey strict-write histories through finite reciprocal material echoes."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from measure import HERE, ROOT, digest, face_basis, federation
from model.candidate_laws.interval_echo.openwave.row_28_acceleration import (
    experiment as repair,
)


def counts(n, edges):
    # Sparse rows retain all reflector identities out to three links from each source.
    x, y = edges.T
    graph = csr_matrix(
        (np.ones(2 * len(x), bool), (np.r_[x, y], np.r_[y, x])), shape=(n, n)
    )
    graph.setdiag(True)
    one = graph[::12]
    two = one @ graph
    three = two @ graph
    return np.stack([two.getnnz(axis=1), three.getnnz(axis=1)])


def spatial_masks(level):
    board = federation(level)
    x, y = np.asarray(board.external).T // 12
    g = csr_matrix((np.ones(len(x)), (x, y)), shape=(board.cells, board.cells))
    hop = shortest_path(g, directed=False, unweighted=True)
    return [hop == 0, hop == 1, hop == 2, 2 * hop >= hop.max()]


def statistics(q, masks):
    # q has independent preparation, repeated schedule and site axes.
    # Randomizing the two schedule labels leaves their expectation invariant.
    mid = q.mean(axis=1)
    mid -= mid.mean(axis=0)
    n = len(q)
    diff = (q[:, 0] - q[:, 1]) / np.sqrt(2.0)
    between_means = mid.T @ mid / (n - 1)
    within = diff.T @ diff / n
    shared = between_means - within / 2
    total = shared + within
    return np.array([[a[m].mean() for m in masks] for a in (total, shared, within)])


def analyze(fields, blank, masks):
    values, vectors = face_basis()
    n = len(fields)
    out = []
    for aperture, ref in enumerate(blank):
        q = np.log(fields[:, :, aperture, :] / ref) / 3
        stat = statistics(q, masks)
        leave = np.stack(
            [statistics(np.delete(q, i, axis=0), masks)[0] for i in range(n)]
        )
        error = np.sqrt((n - 1) / n * ((leave - leave.mean(axis=0)) ** 2).sum(axis=0))
        # Normalize on the 20 equal-weight parent faces. Center each site's ensemble mean.
        face = q.reshape(n, 2, 20, -1).mean(axis=3)
        mode = face @ vectors
        # The same preparation can induce dependence between schedules.
        # Average their two unbiased across-preparation sample variances.
        power = mode.var(axis=0, ddof=1).mean(axis=0)
        out.append(
            dict(
                return_ticks=4 + 2 * aperture,
                mean_log_scale=float(q.mean()),
                total_cov=stat[0].tolist(),
                shared_prep_cov=stat[1].tolist(),
                schedule_cov=stat[2].tolist(),
                total_jackknife_se=error.tolist(),
                projected_face_power=[
                    float(power[values == b].mean()) for b in np.unique(values)[2:]
                ],
                fraction_changed=float((fields[:, :, aperture, :] != ref).mean()),
            )
        )
    return out


def run(dest):
    repair.calibrate()  # Check both transport laws on both written signs using native matrices.
    raw, results = {}, {}
    for level in (0, 1, 2):
        n, edges, ranks, *_ = repair.layout(level)
        blank = counts(n, edges)
        masks = spatial_masks(level)
        raw[f"L{level}_blank"] = blank
        fields = np.empty((32, 2, 4, 2, n // 12), np.uint8)
        last = np.empty((32, 2), np.int16)
        for i in range(32):
            prep = np.random.default_rng([20261125, level, i, 0])
            loads = 2 * prep.integers(0, 2, n, dtype=np.int8)
            for j in range(2):
                offers = np.random.default_rng([20261125, level, i, j + 1])
                commit, _, end = repair.history(n, edges, ranks, loads, 0.25, offers)
                last[i, j] = end
                for k, t in enumerate((0, 4, 16, end)):
                    active = (commit == 0) | (commit > t)
                    read = counts(n, edges[active]) if t else blank
                    assert np.all((read > 0) & (read <= blank)) and blank.max() <= 255
                    fields[i, j, k] = read
        raw[f"L{level}_counts"] = fields
        raw[f"L{level}_last"] = last
        results[level] = {
            name: analyze(fields[:, :, k], blank, masks)
            for k, name in enumerate(("initial", "round4", "round16", "settled"))
        }
    dest.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest / "survey_counts.npz", **raw)
    receipt = dict(
        preparations_per_size=32,
        schedules_per_preparation=2,
        distance_bins=[
            "same site",
            "one hop",
            "two hops",
            "at least half the diameter",
        ],
        mode_eigenvalues=np.unique(face_basis()[0])[2:].tolist(),
        measurements=results,
        counts_sha256=digest(dest / "survey_counts.npz"),
        protocol_sha256=digest(HERE / "PROTOCOL.md"),
        producer_sha256=digest(Path(__file__)),
        source_sha256={
            str(p.relative_to(ROOT)): digest(p)
            for p in (
                Path(repair.__file__),
                Path(repair.priority.__code__.co_filename),
                Path(repair.candidate.native_law.__file__),
            )
        },
    )
    (dest / "survey.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE)
    run(parser.parse_args().out)
