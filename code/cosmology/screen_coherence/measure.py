"""Measure the finite stationary screen and its response to a sharp regional preparation."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from model.candidate_laws.interval_echo import observations
from model.executable_model import federation

APPOINTMENTS = (0, 1, 4, 16, 64, 128)
DENSITY = 0.16


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def face_basis():
    edges = np.asarray(federation(0).external) // 12
    lap = 3 * np.eye(20)
    i, j = edges.T
    lap[i, j] = lap[j, i] = -1
    values, vectors = np.linalg.eigh(lap)
    return values.round(8), vectors


def face_power(matter, basis):
    values, vectors = basis
    cells_per_face = len(matter) // 20
    face = matter.reshape(20, cells_per_face).sum(axis=1) - cells_per_face * DENSITY
    coeff = (face / np.sqrt(cells_per_face)) @ vectors
    return [float(np.mean(coeff[values == value] ** 2)) for value in np.unique(values)]


def edge_cov(first, last, x, y, center):
    # The contact graph is regular, so the end means equal the global site means.
    products = np.count_nonzero(first[x] & last[y]) + np.count_nonzero(
        first[y] & last[x]
    )
    return float(
        products / (2 * len(x)) - center * (first.mean() + last.mean()) + center**2
    )


def archive(folder, level, basis):
    path = folder / f"L{level}"
    manifest = json.loads((path / "universe.json").read_text())
    cells, duration = manifest["cells"], manifest["services"]
    with np.load(path / "universe.npz") as z:
        external, matter, primary = z["external"], z["double"], z["primary"]
    x, y = (external // 12).T
    assert cells == 20 * 4**level and len(x) == 3 * cells // 2
    assert np.all(
        np.bincount(x, minlength=cells) + np.bincount(y, minlength=cells) == 3
    )
    # Verify the cell ordering really lifts the native 20-face contact graph.
    parents = (external // 12) // (cells // 20)
    inter = parents[parents[:, 0] != parents[:, 1]]
    base_edges = np.sort(np.asarray(federation(0).external) // 12, axis=1)
    assert set(map(tuple, np.sort(inter, axis=1))) == set(map(tuple, base_edges))
    assert primary.shape == matter.shape == (duration + 1, cells // 8)
    bits = np.unpackbits(np.arange(256, dtype=np.uint8)[:, None], axis=1).sum(axis=1)
    total = manifest["screen_summary"]["initial_doubles"]
    assert all(int(bits[row].sum()) == total for row in matter)
    assert (
        sum(int(bits[row].sum()) for row in primary[1:])
        == manifest["screen_summary"]["positive_marks"]
    )
    first = np.unpackbits(matter[0])
    rows = []
    for t in sorted(set((*APPOINTMENTS, duration))):
        last = np.unpackbits(matter[t])
        row = dict(
            service=t,
            matter_mean=float(last.mean()),
            initial_onsite_cov=float(
                np.mean(first & last)
                - DENSITY * (first.mean() + last.mean())
                + DENSITY**2
            ),
            initial_neighbor_cov=edge_cov(first, last, x, y, DENSITY),
            equal_neighbor_cov=edge_cov(last, last, x, y, DENSITY),
            matter_face_mode_power=face_power(last, basis),
        )
        if t:
            mark = np.unpackbits(primary[t])
            row["primary_neighbor_cov"] = edge_cov(mark, mark, x, y, 0.5)
        rows.append(row)
    return external, dict(
        cells=cells,
        seed=manifest["seed"],
        services=duration,
        initial_doubles=total,
        rer_commit=manifest["rer_commit"],
        acquisition_sources_modified=manifest["rer_sources_modified_in_checkout"],
        acquisition_source_sha256=manifest["rer_sources_sha256"],
        npz_sha256=digest(path / "universe.npz"),
        metadata_sha256=digest(path / "universe.json"),
        rows=rows,
    )


def near_half(cells, groups):
    edges = np.concatenate([np.stack((x, y)) for x, y in groups], axis=1)
    graph = coo_matrix(
        (np.ones(edges.shape[1], dtype=np.int8), edges), shape=(cells, cells)
    ).tocsr()
    order, parent = breadth_first_order(graph, 0, directed=False)
    assert len(order) == cells
    distance = np.zeros(cells, dtype=np.int32)
    for cell in order[1:]:
        distance[cell] = distance[parent[cell]] + 1
    return np.argsort(distance, kind="stable")[: cells // 2], int(distance.max())


def apply_mean(field, x, y, p):
    delta = p * (field[y] - field[x])
    field[x] += delta
    field[y] -= delta


def paired_moves(field, x, y, swap):
    a, b = field[..., x], field[..., y]
    field[..., x] = np.where(swap, b, a)
    field[..., y] = np.where(swap, a, b)


def mixing(level, groups, p, trajectories):
    cells = 20 * 4**level
    small = level <= 2
    near, radius = near_half(cells, groups)
    mean = np.zeros(cells)
    mean[near] = 1
    impulse = np.zeros(cells)
    impulse[0] = 1
    reach = impulse.astype(bool)
    half_time = reach_time = None
    if small:
        copies = np.tile(mean.astype(np.uint8), (128, 1))
        rng = np.random.default_rng(np.random.SeedSequence([20261125, level]))
    rows = []
    for t in range((1024 if small else 512) + 1):
        if t:
            x, y = groups[(t - 1) % len(groups)]
            apply_mean(mean, x, y, p)
            apply_mean(impulse, x, y, p)
            reach[x] = reach[y] = reach[x] | reach[y]
            if small:
                paired_moves(copies, x, y, rng.random((len(copies), len(x))) < p)
        contrast = 2 * mean[near].mean() - 1
        if contrast <= 0.5 and half_time is None:
            half_time = t
        if reach.all() and reach_time is None:
            reach_time = t
        if t and t & (t - 1):
            continue
        assert np.isclose(mean.sum(), cells / 2) and np.isclose(impulse.sum(), 1)
        assert not np.any(impulse[~reach])
        row = dict(
            service=t,
            contrast=float(contrast),
            frozen_exchange_contrast=1,
            reachable=int(reach.sum()),
            impulse_home=float(impulse[0]),
            impulse_outside_support=float(impulse[~reach].sum()),
        )
        if small:
            assert (copies.sum(axis=1) == cells // 2).all()
            samples = 2 * copies[:, near].sum(axis=1).astype(np.int32) - cells // 2
            trajectories[f"L{level}_t{t}"] = samples
            sample = samples / (cells // 2)
            row.update(
                sample_mean=float(sample.mean()),
                sample_se=float(sample.std(ddof=1) / np.sqrt(len(sample))),
            )
        rows.append(row)
    return dict(
        cells=cells,
        root_eccentricity=radius,
        first_half_contrast_service=half_time,
        first_full_support_service=reach_time,
        matchings=len(groups),
        rows=rows,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sim", type=Path, default=ROOT.parent / "oph-model-sim")
    parser.add_argument("--out", type=Path, default=HERE)
    parser.add_argument(
        "--small-only",
        action="store_true",
        help="Run the three native intervention boards without the external archive.",
    )
    args = parser.parse_args()
    if args.small_only and args.out.resolve() == HERE:
        parser.error("Choose a separate --out directory for the small-only receipt.")
    sim = args.sim.resolve()
    if not args.small_only:
        sys.path.insert(0, str(sim))
        from oph_model_sim import engine
    args.out.mkdir(parents=True, exist_ok=True)
    tables = observations.instruments()
    p = float(tables[1][0, 2])
    assert np.allclose(tables[1][[0, 2], [2, 0]], 0.5, atol=1e-13, rtol=0)
    basis = face_basis()
    bands, dims = np.unique(basis[0], return_counts=True)
    result = dict(
        protocol_sha256=digest(HERE / "PROTOCOL.md"),
        producer_sha256=digest(Path(__file__)),
        source_sha256={
            str(path.relative_to(ROOT)): digest(path)
            for path in [
                Path(observations.__file__),
                Path(observations.native_law.__file__),
                ROOT / "model/executable_model/carrier.py",
            ]
        },
        schedule_helper_sha256=None
        if args.small_only
        else digest(Path(engine.__file__)),
        exchange_probability=p,
        mode_eigenvalues=bands.tolist(),
        mode_dimensions=dims.tolist(),
        archive={},
        interventions={},
    )
    trajectories = {}
    for level in (0, 1, 2) if args.small_only else (0, 1, 2, 6, 8):
        if not args.small_only and level in (2, 6, 8):
            external, result["archive"][level] = archive(sim / "runs", level, basis)
            groups = engine.schedule(external, 20 * 4**level)
        if level <= 2:
            native = [
                (g[0], g[1]) for g in observations.schedules(federation(level))[0]
            ]
            if level == 2 and not args.small_only:
                assert len(native) == len(groups) and all(
                    np.array_equal(a, b) for a, b in zip(native, groups)
                )
            groups = native
        result["interventions"][level] = mixing(level, groups, p, trajectories)
    np.savez_compressed(args.out / "regional_contrasts.npz", **trajectories)
    result["trajectory_sha256"] = digest(args.out / "regional_contrasts.npz")
    (args.out / "measurements.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
