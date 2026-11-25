"""Region of one station in the native marked one-mobile register, after two periods."""

import numpy as np

from model.experiments.seam_program.clocked_apparatus import patches
from model.experiments.seam_program.compiler import Compiler


def read():
    _, data, _, _, _ = patches(0)
    comp = Compiler(data, 0, (2, 2, 2, 1))
    coords = np.indices((2, 2, 2)).reshape(3, -1).T
    # A pi wave along x or along y. Both have a nonzero and identical scalar frequency.
    launches = np.repeat(1 - 2*coords[:, :2], 8, axis=0) / 8
    records = {}
    for member in ("nominal", "shift_only"):
        step = comp.columns(comp.build(member))[comp.slots]
        final = np.linalg.matrix_power(step, 2) @ launches
        inside = final[:8]  # station at the origin, all eight slots retained
        gram = inside.conj().T @ inside
        p = gram.diagonal().real
        overlap = abs(gram[0, 1]) / np.sqrt(p.prod())
        # The restricted state has one vacuum component and the one-mobile Gram matrix.
        target = np.array([0.125, 0.125])
        np.testing.assert_allclose(p, target, rtol=0, atol=2e-13)
        records[member] = (float(p[0]), float(overlap))
    return records


if __name__ == "__main__":
    for member, pair in read().items():
        print(f"{member}: per-station probability {pair[0]:.8g}, local ray overlap {pair[1]:.8g}")
