"""Check conservation of charge dressed by the recorded signed delivery on all echo rows."""

import json

import numpy as np

from .apparatus import K, LAW
from .code_frame import report
from .exact_bands import certify


def delivery(c, d):
    return (d + (1 - 2 * c) * (d % 2)) // 2


def local_certificate():
    w = LAW._write_pattern()
    # Nonzero supports of W and (I-WW*) are exact integer matrices after
    # clearing sqrt(2) and 2 separately, which suffices for diagonal charges.
    tail = 2 * np.eye(90, dtype=int) - w @ w.T
    support = np.block([[np.zeros((9, 9), dtype=int), w.T], [w, tail]])
    a, b = map(np.array, zip(*LAW.SRC))
    currents = np.repeat([0] + [delivery(c, d) for c, d in LAW.LABELS], 9)
    qa, qb = np.tile(a, 11) + currents, np.tile(b, 11) - currents
    assert np.all((qa[:, None] - qa[None, :]) * support == 0)
    assert np.all((qb[:, None] - qb[None, :]) * support == 0)
    for c, d in LAW.LABELS:
        assert delivery(c, d) - delivery(c, -d) == d
    for a, b in LAW.SRC:
        for c in (0, 1):
            x, y = LAW.native_event(a, b, c)
            assert delivery(c, a - b) == a - x == y - b
    Ua = LAW.echo()
    ga, gb = LAW.generator_parts()
    hx, hc = LAW.local_term_parts()
    u = np.exp(1j * (0.37 * qa - 0.83 * qb))

    def rephase(matrix):
        return u[:, None] * matrix * u.conj()[None, :]

    err = np.linalg.norm(rephase(hx + K * hc) - Ua @ rephase(ga + K * gb) @ Ua.T, 2)
    assert err < 2e-13
    return dict(
        echo_rows=99,
        exact_dressed_charge_commutators_zero=True,
        signed_contrast_witnesses=18,
        compensated_hamiltonian_check_error=float(err),
    )


if __name__ == "__main__":
    print(
        json.dumps(
            dict(echo=local_certificate(), code=report(), exact_currents=certify()),
            indent=2,
        )
    )
