"""Recover a patch's angular frame and test which receiver subsets span its null design."""

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo import candidate
from model.evaluation.patch_spatial_readout import transitions
from model.executable_model import cell


HERE = Path(__file__).resolve().parent


def reconstruct(counts):
    symmetric = counts + counts.T
    upper = symmetric[np.triu_indices(len(counts), k=1)]
    ordered = np.sort(upper)
    gaps = np.diff(ordered)
    largest = np.flatnonzero(gaps == gaps.max())
    if len(largest) != 1 or gaps.max() == 0:
        return None
    neighbours = symmetric > ordered[largest[0]]
    np.fill_diagonal(neighbours, False)
    degree = neighbours.sum(axis=1)
    if len(set(degree)) != 1 or connected_components(neighbours, directed=False)[0] != 1:
        return None
    eigenvalues, eigenvectors = np.linalg.eigh(np.diag(degree) - neighbours)
    modes = eigenvectors[:, np.isclose(eigenvalues, eigenvalues[1], atol=1e-9, rtol=0)]
    lengths = np.linalg.norm(modes, axis=1)
    if lengths.min() < 1e-9:
        return None
    frame = modes / lengths[:, None]
    return neighbours, frame


def angular_coverage(frame, labels):
    gram = frame @ frame.T
    kernel = (1 + gram)**2
    # Orthonormal coordinates of k k^T in the Euclidean Frobenius inner product.
    k = np.column_stack((np.ones(len(gram)), frame))
    features = np.column_stack([k[:, i]**2 for i in range(k.shape[1])] + [
        np.sqrt(2) * k[:, i] * k[:, j] for i, j in combinations(range(k.shape[1]), 2)
    ])
    np.testing.assert_allclose(features @ features.T, kernel, atol=1e-10)
    coverage = {}
    least = []
    best_score = -1.
    best_addresses = None
    for size in (9, 10, 11, 12):
        ranks = Counter()
        for group in combinations(range(len(gram)), size):
            spectrum = np.linalg.eigvalsh(kernel[np.ix_(group, group)])
            direct = np.linalg.svd(features[list(group)], compute_uv=False)
            rank = int(np.count_nonzero(spectrum > 1e-9))
            if rank != np.count_nonzero(direct > 1e-9):
                raise ArithmeticError("independent design factorizations disagree")
            ranks[str(rank)] += 1
            if size == 9 and rank == 9:
                least.append(direct[8])
                score = round(float(direct[8]), 9)
                addresses = sorted(int(labels[i]) for i in group)
                if score > best_score or (score == best_score and addresses < best_addresses):
                    best_score, best_addresses = score, addresses
        coverage[str(size)] = dict(sorted(ranks.items()))
    antipodes = np.argmin(gram, axis=1)
    if not np.allclose(gram[np.arange(len(gram)), antipodes], -1):
        raise ArithmeticError("the recovered frame has no matched antipodes")
    even = (features + features[antipodes]) / 2
    return dict(frame_rank=frame.shape[1], subset_size_to_rank_counts=coverage,
                nine_direction_minimum_singular_value_range=[round(float(f(least)), 9) for f in (min, max)],
                best_nine_example_native_ports=best_addresses,
                antipodal_average_rank=int(np.linalg.matrix_rank(even, tol=1e-9)))


def run():
    car = cell()
    generator = candidate.generator(car).toarray()
    probabilities = transitions(generator, 1/4)
    truth = car.adjacency().toarray() > 0
    receipts = []
    design = None
    for budget in (128, 65536):
        counts_by_outcome = Counter()
        ranks = Counter()
        for repeat in range(32):
            rng = np.random.default_rng([74291, budget, repeat])
            labels = rng.permutation(car.n)
            renamed = probabilities[np.ix_(labels, labels)]
            counts = np.column_stack([rng.multinomial(budget, column) for column in renamed.T])
            rebuilt = reconstruct(counts)
            if rebuilt is None:
                counts_by_outcome["abstained"] += 1
                continue
            neighbours, frame = rebuilt
            ranks[str(frame.shape[1])] += 1
            correct = np.array_equal(neighbours, truth[np.ix_(labels, labels)])
            counts_by_outcome["correct" if correct else "incorrect"] += 1
            if budget == 65536:
                measured = angular_coverage(frame, labels)
                if design is not None and design != measured:
                    raise ArithmeticError("direction designs differ between held tables")
                design = measured
        receipts.append(dict(preparations_per_launch=budget, tables=32, outcomes=dict(counts_by_outcome), ranks=dict(ranks)))
    return dict(exposure=0.25, count_tables=receipts, angular_design=design)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a fresh receipt")
    output = cli.parse_args().output
    result = json.dumps(run(), indent=2) + "\n"
    if output:
        with output.open("x") as handle:
            handle.write(result)
    elif (HERE / "results.json").read_text() != result:
        raise SystemExit("angular measurement differs from the deposited receipt")
    print(result, end="")
