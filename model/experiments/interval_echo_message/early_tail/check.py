"""Check the two addressed matrix elements behind the short-exposure signal."""
from fractions import Fraction

import numpy as np
from scipy.sparse import eye

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import channel, headline


def check(station):
    board = federation(station["level"])
    h = candidate.generator(board).tocsr()
    port, sender = station["contact"]
    near, far = station["ports"]
    if near != sender or port >= 12 or min(near, far) < 12:
        raise ValueError("the receiver must begin empty on both launches")
    launch = channel.encode(board, station)
    np.testing.assert_allclose(np.abs(launch[[near, far], [0, 1]]), 1, atol=1e-12, rtol=0)

    # Reference seams give a uniform 6*kappa shift. All row entries below
    # are integers or half-integers, exactly representable in this sector.
    k = h - (6 * candidate.KAPPA + 3) * eye(board.n, format="csr")
    skew = k - k.T.conj()
    if skew.count_nonzero() or np.max(np.asarray(abs(k).sum(axis=1))) > 3:
        raise ValueError("centered Hermitian generator exceeds the norm bound")
    if k[port, near] != -0.5 or k[port, far] != 0:
        raise ValueError("the leading receiver amplitudes changed")

    for exposure in (Fraction(1, 1000), Fraction(1, 100), Fraction(1, 36)):
        t = float(exposure)
        floor = exposure**2 * (Fraction(1, 4) - Fraction(9, 2) * exposure)
        p, _ = channel.readout(h, launch, t)
        np.testing.assert_allclose(p, headline.rebuild(station, t), rtol=0, atol=1e-12)
        delta = p[0, port] - p[1, port]
        if not delta >= float(floor) > 0:
            raise ValueError("numerical receiver probabilities missed the analytic floor")
        print(f"{board.cells} cells, t={exposure}: click difference {delta:.9g}, floor {float(floor):.9g}")


if __name__ == "__main__":
    for station in channel.SPEC["stations"]:
        check(station)
