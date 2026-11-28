"""Compare the original 90% selector across the complete transport menu."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, norm

from model.executable_model import federation
from ...row_30_arrows.experiment import ARMS, local_check, sector
from ..graph import DONOR


HERE = Path(__file__).resolve().parent
THRESHOLD = 0.9
# Cluster more broadly than the arithmetic allowance. A cluster may contain
# more than one exact level. Ruling out the entire cluster rules out each level.
GROUP_WIDTH = 1e-6
RESIDUAL = 2e-10
ORTHOGONALITY = 1e-8
SUBSPACE_ERROR = 0.000202
BLOCK_PRODUCT_ALLOWANCE = 1e-7


def verify_frame(operator, rays, spectrum):
    # Gershgorin and the explicit integer sparse coefficients control each
    # residual evaluation. The array size bounds dot/sum rounding separately.
    size = len(spectrum)
    assert size <= 2730 and np.max(np.abs(rays)) <= 1
    assert np.max(np.diff(operator.indptr)) <= 7
    assert np.max(np.asarray(abs(operator).sum(axis=1))) <= 64
    assert np.max(abs(spectrum)) <= 64
    assert norm(operator @ rays - rays * spectrum, "fro") < 1e-10
    assert norm(rays.T @ rays - np.eye(size), "fro") < 5e-9


def exact_weight_upper(approximate_amplitude):
    return (
        (1 + ORTHOGONALITY) * (approximate_amplitude + BLOCK_PRODUCT_ALLOWANCE)
        + 2 * SUBSPACE_ERROR
        + SUBSPACE_ERROR**2
    ) ** 2


def labels(values, resolution):
    return np.split(
        np.arange(len(values)), np.where(np.diff(values) > resolution)[0] + 1
    )


def measure_overlap(
    temporal, spatial, energies, momenta, resolution=GROUP_WIDTH, enclose=True
):
    e_groups = labels(energies, resolution)
    p_groups = labels(momenta, resolution)
    overlap = temporal.T @ spatial
    squared = overlap**2
    row_starts = [group[0] for group in e_groups]
    col_starts = [group[0] for group in p_groups]
    block_bound = np.add.reduceat(
        np.add.reduceat(squared, row_starts, axis=0), col_starts, axis=1
    )
    rows, cols = np.where(block_bound >= 0.88)
    matches = []
    reject_block_amplitude = float(
        np.sqrt(np.max(block_bound[block_bound < 0.88], initial=0))
    )
    rejected = reject_block_amplitude**2
    for r, c in zip(rows, cols, strict=True):
        i, j = e_groups[r], p_groups[c]
        block = overlap[np.ix_(i, j)]
        left, singular, right = np.linalg.svd(block, full_matrices=False)
        assert norm((left * singular) @ right - block, "fro") < 1e-10
        assert norm(left.T @ left - np.eye(len(singular)), "fro") < 1e-10
        assert norm(right @ right.T - np.eye(len(singular)), "fro") < 1e-10
        gain = float(singular[0] ** 2)
        if gain >= THRESHOLD:
            assert gain > THRESHOLD + 1e-5
            matches.append(
                {
                    "energy": round(float(energies[i[0]]), 8),
                    "energy_index_range": [int(i[0]), int(i[-1])],
                    "energy_dimension": len(i),
                    "wave_number_squared": round(float(momenta[j[0]]), 8),
                    "wave_dimension": len(j),
                    "maximum_joint_weight": round(gain, 9),
                }
            )
        else:
            assert gain < THRESHOLD - 1e-5
            rejected = max(rejected, gain)
    result = {
        "distinct_temporal_groups": len(e_groups),
        "distinct_spatial_groups": len(p_groups),
        "passing_pairs": matches,
    }
    if enclose:
        reject_upper = exact_weight_upper(np.sqrt(rejected))
        assert reject_upper < THRESHOLD
        result["upper_on_exact_joint_weight_outside_passing_clusters"] = round(
            float(reject_upper) + 1e-8, 8
        )
    return result


def native_leaf_modes(h, lap, sites, board):
    # The integer difference of two leaf children of one parent is a joint
    # eigenvector whenever their degree-reference terms agree.
    result = []
    graph = lap.copy()
    graph.setdiag(0)
    graph.eliminate_zeros()
    for parent in range(len(sites)):
        row = graph.indices[graph.indptr[parent] : graph.indptr[parent + 1]]
        for degree in (5, 6):
            leaves = [
                j
                for j in row
                if graph.indptr[j + 1] - graph.indptr[j] == 1
                and j != 0
                and int(board.degree[sites[j]]) == degree
            ]
            if len(leaves) < 2:
                continue
            energy = 8 * (6 - degree) - 4
            for other in leaves[1:]:
                d = np.zeros(len(sites), int)
                d[leaves[0]], d[other] = (1, -1)
                np.testing.assert_array_equal(h @ d, energy * d)
                np.testing.assert_array_equal(lap @ d, d)
            result.append(
                {
                    "parent_port": int(sites[parent]),
                    "leaf_ports": [int(sites[i]) for i in leaves],
                    "joint_dimension_contributed": len(leaves) - 1,
                    "energy": energy,
                    "wave_number_squared": 1,
                }
            )
    return result


def run():
    assert (
        ORTHOGONALITY
        + np.sqrt(1 + ORTHOGONALITY)
        * RESIDUAL
        / (GROUP_WIDTH - RESIDUAL / np.sqrt(1 - ORTHOGONALITY))
        < SUBSPACE_ERROR
    )
    local_check()
    board = federation(0)
    addresses = list(range(board.n))
    orders = [
        addresses,
        sorted(addresses, key=lambda q: (q // 12, (-q) % 12)),
        addresses[::-1],
    ]
    comparison = []
    for orientation, order in zip(ARMS, orders, strict=True):
        h, sites, census = sector(board, order)
        assert census["omitted_states"] == 0 and h.shape[0] == len(sites)
        # Every off-diagonal entry of the temporal generator is -4 on
        # an allowed joint traversal. Drop only the native degree diagonal.
        graph = h.copy() / 4
        graph.setdiag(0)
        graph.eliminate_zeros()
        lap = graph.copy()
        lap.setdiag(-np.asarray(graph.sum(axis=1)).ravel())
        x, v = eigh(h.toarray())
        y, z = eigh(lap.toarray())
        verify_frame(h, v, x)
        verify_frame(lap, z, y)
        depths = np.zeros(len(sites), dtype=int)
        for i in range(1, len(sites)):
            row = graph.indices[graph.indptr[i] : graph.indptr[i + 1]]
            (parent,) = row[row < i]
            depths[i] = depths[parent] + 1
        marked = (-1.0) ** depths
        plain = measure_overlap(v, z, x, y)
        staggered = measure_overlap(v, marked[:, None] * z, x, y)
        fine_staggered = measure_overlap(
            v, marked[:, None] * z, x, y, resolution=1e-8, enclose=False
        )
        staggered["historical_resolution"] = 1e-8
        staggered["historical_resolution_passing_group_count"] = len(
            fine_staggered["passing_pairs"]
        )
        # Each observed pass sits well away from both the resolution and the
        # eigenspace overlap threshold. Keep the exact algebra only for leaves.
        if orientation == "ascending":
            archived = json.loads((DONOR / "results.json").read_text())[
                "spectral_audit"
            ]["path"]
            observed = [
                (round(row["energy"], 6), round(row["wave_number_squared"], 6))
                for row in plain["passing_pairs"]
            ]
            expect = [
                (round(row["energy"], 6), round(k, 6))
                for row in archived
                for k in row["compatible_k2"]
            ]
            assert observed == expect
            stored_staggered = json.loads((DONOR / "staggered.json").read_text())[
                "eigenspaces"
            ]
            archived_staggered = [
                (round(row["energy"], 6), round(row["best_k2"], 6))
                for row in stored_staggered
                if row["best_weight"] >= THRESHOLD
            ]
            observed_staggered = [
                (round(row["energy"], 6), round(row["wave_number_squared"], 6))
                for row in staggered["passing_pairs"]
            ]
            assert observed_staggered == archived_staggered
        comparison.append(
            {
                "arrow_order": orientation,
                "histories": len(sites),
                "reachable_physical_ports": census["reachable_ports"],
                "plain": plain,
                "staggered": staggered,
                "exact_equal_degree_leaf_spaces": native_leaf_modes(
                    h, lap, sites, board
                ),
            }
        )
    return {
        "scope": "20-cell spectrum with excluded-group bounds; secondary passing clusters are diagnostic",
        "joint_selector_threshold": THRESHOLD,
        "spectral_cluster_width": GROUP_WIDTH,
        "residual_norm_bound": RESIDUAL,
        "orthogonality_norm_bound": ORTHOGONALITY,
        "exact_eigenvector_distance_allowance": SUBSPACE_ERROR,
        "families": comparison,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.record:
        with args.record.open("x") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads((HERE / "menu.json").read_text())
    print(
        [
            (
                r["arrow_order"],
                r["histories"],
                len(r["plain"]["passing_pairs"]),
                len(r["staggered"]["passing_pairs"]),
            )
            for r in result["families"]
        ]
    )


if __name__ == "__main__":
    main()
