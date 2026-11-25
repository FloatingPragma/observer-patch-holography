"""Check a Born-free sharp-source forecast against a second native reader."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import circuit as c
from . import measure

HERE = Path(__file__).parent


def run():
    rng = np.random.default_rng(113031208)
    # Words 0 and 1 differ at j=1. Word 1 holds load 2 there.
    fix = c.label_phase(1, [0, 0, np.pi / 2])
    balanced = c.word_rotation(0, np.pi / 4, -np.pi / 2) + fix
    swap = [c.contact(1, 2, np.pi)]
    coin, _ = measure.committed_coin()
    bits = np.einsum("i,j,k->ijk", coin, coin, coin).ravel()
    data = {}
    for name, pulses in (
        ("balanced", balanced),
        ("swapped", balanced + swap),
        ("unbalanced", fix),
    ):
        banks = {}
        for src in ["zero", "one"]:
            v = c.play(pulses, c.prepare(src, c.K), c.K)
            if name == "balanced":
                w = c.play(swap, v, c.K)
                sign = 1 if src == "zero" else -1
                np.testing.assert_allclose(w, sign * v, rtol=0, atol=3e-12)
            p = np.kron(c.read(v), bits)
            banks[src] = measure.counts(rng, measure.N[src], p).tolist()
        data[name] = banks
    return {
        "seed": 113031208,
        "sampling_grid_bits": measure.GRID_BITS,
        "banks": data,
        "exact_source_forecast": {"fine": [4, 4, 8], "denominator": 16},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "symmetry.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(path.read_bytes()):
        raise ValueError("Symmetry receipt changed")
    print(
        "Verified both sharp-source swap eigenrays and the full stochastic read banks."
    )
