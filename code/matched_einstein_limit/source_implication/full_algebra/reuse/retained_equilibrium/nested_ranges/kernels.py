"""Native cap-2, 4 and 6 effects from the same four integer propagators."""

import argparse
from pathlib import Path

import numpy as np

from .... import native_precision as native

HERE = Path(__file__).parent
GRID = 2**56
CAPS = (2, 4, 6)


def build():
    board = native.iso.pairing.formation.federation(1)
    ports = native.iso.pairing.formation.region(board)
    ix = np.asarray([4 * x + j for x in ports for j in range(4)])
    coeff = native.prior.weights()
    bank = np.zeros((3, 12, 2, 2, 4, 4), dtype=object)
    for parity in (0, 1):
        for q in (-1, 1):
            h = native.generator(board, ports, q, bool(parity))
            u = native.kernel_panels(h, ix, coeff)
            for i, cap in enumerate(CAPS):
                bank[i, :, parity] += np.array(
                    native.prior.direction_effects(u, 2 * cap), dtype=object
                )
    den = 8 * native.prior.FRAC
    # The retained constant cI/2 is added after reading the quantized deficit.
    bank = (bank * GRID + den // 2) // den
    assert native.prior.error_effect() < 7.82e-16
    return bank.astype(np.int64)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    archive = HERE / "deficits.npz"
    fresh = build()
    if args.write:
        with archive.open("xb") as stream:
            np.savez_compressed(stream, numerators=fresh)
    else:
        with np.load(archive) as saved:
            np.testing.assert_array_equal(saved["numerators"], fresh)
        print("Verified all 72 pooled, parity-resolved clock effects.")


if __name__ == "__main__":
    main()
