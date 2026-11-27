"""Check the incidence Gram projector and the finite operators transporting it for each candidate."""
from itertools import product
import json

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast
from model.executable_model import cell, laws
from .evaluate import dense, rank_check


def measure(candidate):
    edges = list(candidate.LABELLING)
    incidence = np.zeros((len(edges), 12))
    for i, (a, b) in enumerate(edges):
        incidence[i, [a, b]] = 1
    adjacency = incidence.T @ incidence - 5 * np.eye(12)
    distances = shortest_path(adjacency, directed=False, unweighted=True).astype(int)
    gram = np.array([1, 1 / np.sqrt(5), -1 / np.sqrt(5), -1])[distances]
    projector = gram / 4
    kinetic = (5 * np.eye(12) - adjacency) / 2
    slow = (5 - np.sqrt(5)) / 2

    # With each internal seam attempted at unit rate, the fair binary write has
    # average outgoing load (a+b)/2 on each endpoint of the addressed seam.
    conditional = {}
    for a, b in product((0, 1), repeat=2):
        outgoing = np.mean([contrast.law.native_event(a, b, coin) for coin in (0, 1)], axis=0)
        if not np.array_equal(outgoing, [(a + b) / 2] * 2):
            raise ValueError("the native binary write has a different conditional mean")
        conditional[a, b] = outgoing - np.array([a, b])
    mean_generator = np.zeros((12, 12))
    for a, b in edges:
        mean_generator[[a, b], a] += conditional[1, 0]
        mean_generator[[a, b], b] += conditional[0, 1]
    board = cell()
    generator = dense(candidate.generator(board))
    response = dense(candidate.response(board))
    labels = candidate.basis(board)
    if not np.isfinite(generator).all() or not np.isfinite(response).all():
        raise ArithmeticError(f"{candidate.NAME}: nonfinite candidate operator")
    check = rank_check(response, generator, labels)
    if check["status"] != "pass":
        raise ArithmeticError(f"{candidate.NAME}: {check['detail']}")
    is_seam_response = response.shape != generator.shape
    if is_seam_response:
        # The scalar lives on seams. The endpoint sum maps port amplitudes to seam amplitudes.
        seams = [where for kind, where in labels if kind == "seam"]
        embedding = np.zeros((len(seams), 12))
        for i, (a, b) in enumerate(seams):
            embedding[i, [a, b]] = 1
        offset = 0
    else:
        if set(labels) != {("port", i) for i in range(12)} or len(labels) != 12:
            raise ValueError(f"{candidate.NAME}: no port-to-response map for this basis")
        embedding = np.eye(12)[[v for _, v in labels]]
        offset = float(response[0].sum().real)
    centered = response - offset * np.eye(len(response))
    lifted = embedding @ projector
    errors = {
        "binary_unit_rate_mean": np.linalg.norm(mean_generator + kinetic),
        "projector_idempotence": np.linalg.norm(projector @ projector - projector),
        "gram_eigenrelation": np.linalg.norm(kinetic @ projector - slow * projector),
        "full_operator_intertwining": np.linalg.norm(centered @ embedding - embedding @ kinetic),
        "candidate_band_eigenrelation": np.linalg.norm(centered @ lifted - slow * lifted),
    }
    rank = np.linalg.matrix_rank(lifted)
    if rank != 3 or max(errors.values()) > 1e-10:
        raise ArithmeticError(f"{candidate.NAME}: incompatible generator or Gram frame")
    return {
        "port_to_candidate_map": ("unsigned seam-port incidence" if is_seam_response else
                                  "identity" if np.array_equal(embedding, np.eye(12)) else "registered port order"),
        "response": "squared frequency minus squared rest level" if is_seam_response else "one-load Hamiltonian",
        "one_load_scalar_offset": offset,
        "lifted_rank": int(rank),
        "frobenius_errors": {name: float(value) for name, value in errors.items()},
    }


def run():
    return {name: measure(laws.load(name)) for name in laws.names()}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
