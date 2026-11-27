"""Ordered modular-flip experiment on the same two received quantum record slots."""
from math import log, pi

import numpy as np


def phase_read(vacuum, records, time):
    """Ramsey imaginary quadrature of X U X U†; U is calibrated on the blank."""
    vacuum = np.asarray(vacuum, dtype=float)
    records = np.asarray(records, dtype=float)
    if len(vacuum) != 4 or len(records) != 4 or np.min(vacuum)<=0 or np.min(records)<0:
        raise ValueError("two received signs and a faithful four-word blank are required")
    if not np.isclose(vacuum.sum(), 1) or not np.isclose(records.sum(), 1):
        raise ValueError("the two word distributions must be normalized")
    u = np.diag(np.exp(-1j*time*np.log(vacuum)))
    flip = np.kron([[0, 1], [1, 0]], np.eye(2))
    loop = flip @ u @ flip @ u.conj().T
    # A controlled loop on the retained slots multiplies one coherence of
    # the path qubit by tr(rho V). Measuring its Y quadrature reads Im tr(rho V).
    return float(np.trace(np.diag(records) @ loop).imag)


def verify():
    from .census import census

    blank = (3/8, 1/8, 1/8, 3/8)
    gap = log(blank[0]/blank[1])
    t = pi/(2*gap)
    _, rows = census()
    for c in (0.48, -0.125, *(float(row["correlation"]) for row in rows.values())):
        held = ((1+c)/4, (1-c)/4, (1-c)/4, (1+c)/4)
        contrast = phase_read(blank, held, t)-phase_read(blank, blank, t)
        energy = -sum((p-q)*log(q) for p, q in zip(held, blank, strict=True))
        if abs(energy+gap*contrast/2) > 3e-15:
            raise ArithmeticError("the ordered modular clock does not give the regional modular-energy change")
        if abs(phase_read((1/4,)*4, held, t)) > 3e-15:
            raise ArithmeticError("a flat reference produced a modular frequency")
    return gap
