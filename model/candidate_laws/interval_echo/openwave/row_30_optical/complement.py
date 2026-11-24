"""Remove the stationary wave by spatial preparation and read every frequency bin."""

import argparse
import json

import numpy as np
from scipy import linalg

from model.executable_model import federation
from .experiment import HERE, counts, exact_stationary, preparations, spatial_modes
from ..row_30_formed.experiment import synthesize


def panel(level):
    board = federation(level)
    axes, edges, profile, norm, lcell = preparations(board)
    _, project = exact_stationary(board, axes, edges)
    h = board.laplacian() / 2
    energies, vectors = linalg.eigh(h.toarray(), driver="evd")
    # The uniform load provides the zero. The relative spectrum lies in [0,6].
    # Use 128 release-time arms, with delay step 2 pi / 8, as a Fourier meter.
    q, period = 128, 8
    bins = period * np.arange(q) / q
    diff = (energies[:, None] - bins) / period
    response = (np.sinc(q * diff) / np.sinc(diff)) ** 2
    assert np.allclose(response.sum(axis=1), 1, rtol=0, atol=1e-12)
    rng = np.random.default_rng(303442 + level)
    root = int(np.flatnonzero(board.degree == 5)[0])
    rows = []
    for i, (_, envelope) in enumerate(spatial_modes(lcell, norm)):
        if i == 0:
            continue
        original = (envelope[:, None] * profile).ravel()
        complement = original - project(original)
        norm2 = np.dot(complement, complement)
        source, _ = synthesize(
            complement / np.sqrt(norm2), board.adjacency(), board.degree, root
        )
        weights = abs(vectors.T @ source) ** 2
        odd = abs(source[edges[:, 0]] - source[edges[:, 1]]) ** 2 / 2
        observed = counts(rng, weights @ response)
        rows.append(
            dict(
                source_mode=i,
                source_overlap_squared=round(float(norm2), 10),
                odd=int(rng.binomial(sum(observed), odd.mean())),
                clock=observed,
            )
        )
    return dict(cells=board.cells, bins=bins.tolist(), rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    result = [panel(i) for i in range(3)]
    path = HERE / "complement.json"
    if args.record:
        with path.open("x") as f:
            json.dump(result, f, indent=1)
            f.write("\n")
    elif json.loads(path.read_text()) != result:
        raise SystemExit(
            "Complement source observations differ from the stored receipt"
        )
    else:
        print("Verified all 15 complement-source clock histograms and bond counts.")


if __name__ == "__main__":
    main()
