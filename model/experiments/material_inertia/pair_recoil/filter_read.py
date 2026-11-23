"""Independent calibration of the two-address bound-band preparation."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from . import measure


def run(primary):
    result = []
    for index, source in zip(measure.PLACEMENTS, primary["placements"], strict=True):
        gram = np.array(source["filter_gram"])
        vectors = np.array(
            [[1, 0], [0, 1], [1, 1], [1, -1], [1, 1j], [1, -1j]], complex
        )
        vectors /= np.linalg.norm(vectors, axis=1)[:, None]
        probabilities = np.einsum("bi,ij,bj->b", vectors.conj(), gram, vectors).real
        data = [
            measure.read_ramsey(p, 2**52, [measure.SEED, 8, index, i])
            for i, p in enumerate(probabilities)
        ]
        p = np.array([r["click"] / 2**52 for r in data])
        sigma = probabilities * (1 - probabilities) / 2**52
        result.append(
            dict(
                bridge=source["bridge"],
                counts=data,
                measured_diagonal=p[:2].tolist(),
                measured_cross_real=float((p[2] - p[3]) / 2),
                measured_cross_imag=float((p[5] - p[4]) / 2),
                cross_standard_errors=(
                    np.sqrt([sigma[2] + sigma[3], sigma[4] + sigma[5]]) / 2
                ).tolist(),
            )
        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (measure.HERE / "FILTER_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        placements=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(json.loads((measure.HERE / "results.json").read_text()))
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
