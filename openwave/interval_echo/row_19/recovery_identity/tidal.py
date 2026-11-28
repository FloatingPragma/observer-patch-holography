"""Free two-ray release and delayed radial response on the measured overlap source."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from matched_einstein_limit.source_implication.common_front import measure as front
from . import transfer

HERE = transfer.HERE
RECEIPT = HERE / "tidal.json"
STOCKS = ("ordinary_minus", "ordinary_plus", "held_minus", "held_plus")
CENTERS = (1, 4, 8, 16)
DT = 0.25
SHOTS = 2**60
SEED = 191129


def chart(board, route):
    radii = (
        shortest_path(
            board.adjacency(), unweighted=True, directed=False, indices=list(route[1:3])
        )
        .min(axis=0)
        .astype(int)
    )

    def extensions(path):
        if len(path) == 4:
            return [path]
        return [
            chain
            for j in sorted(board.neighbours[path[-1]])
            if radii[j] == len(path)
            for chain in extensions([*path, j])
        ]

    paths = [
        min(extensions([a, b])) for a, b in ((route[1], route[0]), (route[2], route[3]))
    ]
    rails = [pair for path in paths for pair in zip(path, path[1:])]
    positions = np.zeros((board.n, len(rails)))
    for i, rail in enumerate(rails):
        positions[list(rail), i] = 1 / np.sqrt(2)
    return radii, paths, positions


def check_contacts(integer_h, route):
    """Check disjoint support and one-head factorization on every integer port block."""
    blocks = integer_h.tobsr(blocksize=(4, 4))
    row = np.repeat(np.arange(len(blocks.indptr) - 1), np.diff(blocks.indptr))
    col = blocks.indices
    left = np.isin(row, route[:2]) & np.isin(col, route[:2])
    right = np.isin(row, route[2:]) & np.isin(col, route[2:])
    if np.any(left & right):
        raise ValueError("Marked contacts overlap")
    eye = np.eye(2, dtype=int)
    for mask, which in ((left, "left"), (right, "right"), (~(left | right), "scalar")):
        values = blocks.data[mask]
        if which == "left":
            rebuilt = np.einsum("nab,ij->naibj", values[:, ::2, ::2], eye)
        elif which == "right":
            rebuilt = np.einsum("ab,nij->naibj", eye, values[:, :2, :2])
        else:
            rebuilt = values[:, 0, 0, None, None] * np.eye(4, dtype=int)
        np.testing.assert_array_equal(values, rebuilt.reshape(-1, 4, 4))
    # Head-trace control uses the same spatial blocks and parity weights.
    return sparse.csr_matrix(
        (np.trace(blocks.data, axis1=1, axis2=2) / 8, (row, col)),
        shape=(len(blocks.indptr) - 1,) * 2,
    )


def spatial_reads(h, initial, densities, radii):
    n, rays = initial.shape
    dim = densities.shape[-1]
    launch = np.kron(initial, np.eye(dim))

    def position_mean(left, right):
        a = left.reshape(n, dim, rays, dim)
        b = right.reshape(n, dim, rays, dim)
        return np.einsum("pkra,sab,pkrb->srp", a, densities, b.conj(), optimize=True)

    v = h @ launch
    acceleration = 2 * (position_mean(v, v) - position_mean(h @ v, launch)).real @ radii
    velocity = 2 * position_mean(-1j * v, launch).real @ radii
    waves = expm_multiply(
        -1j * h, launch, start=0, stop=17 * DT, num=18, traceA=-1j * h.diagonal().sum()
    )
    shell = np.eye(radii.max() + 1)[radii]
    distributions = np.stack([position_mean(w, w).real @ shell for w in waves], axis=2)
    return acceleration, velocity, distributions


def two_rays(means):
    """Outer minus inner on each successive pair within a side."""
    return means[:, [1, 2, 4, 5]] - means[:, [0, 1, 3, 4]]


def compare(relative, steps, source_ratio, uncertainty=0.0):
    slopes = [
        (relative[i + 1] - relative[i]) / (2 * steps[key])
        for i, key in ((0, "ordinary"), (2, "held"))
    ]
    per_slope = [uncertainty / (2 * steps[k]) for k in ("ordinary", "held")]
    out = []
    for ordinary, held in zip(*slopes, strict=True):
        prediction = sorted(
            x * y
            for x in source_ratio
            for y in (ordinary - per_slope[0], ordinary + per_slope[0])
        )
        resolved = abs(ordinary) > max(per_slope[0], 1e-10)
        out.append(
            dict(
                ordinary=float(ordinary),
                held=float(held),
                ratio=float(held / ordinary) if resolved else None,
                ordinary_radius=per_slope[0],
                held_radius=per_slope[1],
                predicted=[float(prediction[0]), float(prediction[-1])],
                compatible=(
                    bool(
                        held + per_slope[1] >= prediction[0]
                        and held - per_slope[1] <= prediction[-1]
                    )
                    if resolved
                    else None
                ),
            )
        )
    return out


def timed_comparisons(means, steps, source_ratio, bound):
    relative = two_rays(means)
    result = {}
    for j in CENTERS:
        accel = (
            relative[:, :, j + 1] - 2 * relative[:, :, j] + relative[:, :, j - 1]
        ) / DT**2
        # Two rays, two endpoints, and sum of absolute stencil weights four.
        result[str(j * DT)] = compare(accel, steps, source_ratio, 16 * bound / DT**2)
    return result


def counted_comparisons(histograms, steps, source_ratio):
    if (
        histograms.shape[:3] != (4, 6, 18)
        or histograms.ndim != 4
        or histograms.dtype != np.int64
        or histograms.min() < 0
        or histograms.max() > SHOTS
        or not np.all(histograms.sum(axis=-1, dtype=object) == SHOTS)
    ):
        raise ValueError("Incomplete shell counts")
    shell = np.arange(histograms.shape[-1])
    total_means = 2 * np.prod(histograms.shape[:-1])
    radius = float(shell[-1] * np.sqrt(np.log(2 * total_means / 0.01) / (2 * SHOTS)))
    return radius, timed_comparisons(
        (histograms / SHOTS) @ shell, steps, source_ratio, radius
    )


def inspect_bank(laws, steps, source_ratio, rng):
    shell = np.arange(laws.shape[-1])
    np.testing.assert_allclose(laws.sum(axis=-1), 1, atol=5e-12, rtol=0)
    if laws.min() < 0:
        raise ValueError("Negative shell probability")
    histograms = np.array(
        [
            rng.multinomial(SHOTS, row / row.sum())
            for row in laws.reshape(-1, len(shell))
        ]
    ).reshape(laws.shape)
    radius, sampled = counted_comparisons(histograms, steps, source_ratio)
    return dict(
        counts=histograms.tolist(),
        single_mean_radius=radius,
        sampled=sampled,
        population=timed_comparisons(laws @ shell, steps, source_ratio, 0.0),
    )


def run():
    stocks, steps, predicted, source = transfer.inputs()
    source_ratio = float(predicted.a), float(predicted.b)
    formation = front.optics.ranging.formation
    board = formation.federation(1)
    route = formation.region(board)
    radii, paths, launch = chart(board, route)
    densities = np.array([stocks[key] for key in STOCKS])
    rng = np.random.default_rng(SEED)
    banks = {}
    for q in (-1, 1):
        active, control = [], []
        for parity in (0, 1):
            integer_h = front.native.generator(board, route, q, bool(parity))
            blind = check_contacts(integer_h, route)
            states = densities[:, parity]
            active.append(spatial_reads(integer_h / 2, launch, states, radii))
            weights = states.trace(axis1=1, axis2=2)[:, None, None]
            control.append(spatial_reads(blind, launch, weights, radii))
        acc, velocity, laws = map(sum, zip(*active, strict=True))
        _, _, blind_laws = map(sum, zip(*control, strict=True))
        control_difference = float(np.max(abs(blind_laws[3] - blind_laws[2])))
        recorded = inspect_bank(laws, steps, source_ratio, rng)
        recorded.update(
            release=compare(two_rays(acc), steps, source_ratio),
            release_velocities=velocity.tolist(),
            parity_trace_control_max_held_probability_change=control_difference,
        )
        banks[str(q)] = recorded
    return dict(
        protocol_sha256=sha256((HERE / "TIDAL.md").read_bytes()).hexdigest(),
        source_receipt_sha256=sha256(transfer.RECEIPT.read_bytes()).hexdigest(),
        shot_count=SHOTS,
        seed=SEED,
        sources=STOCKS,
        source=source,
        paths=paths,
        largest_shell=int(radii.max()),
        half_steps=steps,
        banks=banks,
    )


def verify(actual, expected):
    if isinstance(actual, dict):
        if set(actual) != set(expected):
            raise ValueError("Receipt fields differ")
        for key in actual:
            verify(actual[key], expected[key])
    elif isinstance(actual, (list, tuple)):
        for a, b in zip(actual, expected, strict=True):
            verify(a, b)
    elif isinstance(actual, float):
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-10)
    elif actual != expected:
        raise ValueError("Receipt value differs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    if args.counts_only:
        result = json.loads(RECEIPT.read_text())
        src = json.loads(transfer.RECEIPT.read_text())
        for key, val in (
            ("protocol_sha256", sha256((HERE / "TIDAL.md").read_bytes()).hexdigest()),
            (
                "source_receipt_sha256",
                sha256(transfer.RECEIPT.read_bytes()).hexdigest(),
            ),
            ("shot_count", SHOTS),
            ("half_steps", src["half_steps"]),
            ("source", src["source"]),
            ("sources", STOCKS),
        ):
            verify(val, result[key])
        if set(result["banks"]) != {"-1", "1"}:
            raise ValueError("Missing charge")
        interval = src["source"]["predicted_ratio_interval"]
        lower, upper = map(float, interval.strip("[]").split(","))
        for bank in result["banks"].values():
            raw = np.asarray(bank["counts"])
            if raw.shape[-1] != result["largest_shell"] + 1:
                raise ValueError("Missing terminal shell")
            radius, reports = counted_comparisons(
                raw, result["half_steps"], (lower, upper)
            )
            verify(radius, bank["single_mean_radius"])
            verify(reports, bank["sampled"])
    else:
        result = run()
    if args.write:
        with RECEIPT.open("x") as output:
            json.dump(result, output, allow_nan=False, separators=(",", ":"))
            output.write("\n")
    elif not args.counts_only:
        verify(result, json.loads(RECEIPT.read_text()))
    for q, bank in result["banks"].items():
        if not args.counts_only:
            print(q, "release", [(r["ordinary"], r["held"]) for r in bank["release"]])
        for t, rows in bank["sampled"].items():
            print(q, t, [(r["ratio"], r["compatible"]) for r in rows])
