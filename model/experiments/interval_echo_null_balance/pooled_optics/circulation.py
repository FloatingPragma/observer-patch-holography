"""Exterior phase reversal measures steady internal current in native exposure units."""

import json

import numpy as np
import sympy as sp

from model.executable_model import cell, federation
from ..active_current import EDGE, TIMES, pulsed_records
from ..hidden_current import compile_wave, exact_modes, read_current


def finite_replay(board, bridges, u, v):
    # A mixed source which includes a real coherence and both unequal populations.
    pure = [(2 * u + (1 + sp.I) * v) / sp.sqrt(6), (u - sp.I * v) / sp.sqrt(2)]
    weights = np.array([0.8, 0.2])
    waves = [compile_wave(cell(), np.asarray(w, dtype=complex).ravel()) for w in pure]
    initial = np.zeros((board.n, len(pure)), dtype=complex)
    initial[:12] = np.column_stack(waves)
    current = sum(
        p * read_current(cell(), w, EDGE)["kinetic_current"]
        for p, w in zip(weights, waves, strict=True)
    )
    pairs = {}
    for strength, name in ((8, "native"), (0, "exchange_only")):
        forward, reverse = [
            pulsed_records(board, initial, bridges, strength, angle) @ weights
            for angle in (np.pi / 3, 5 * np.pi / 3)
        ]
        contrast = (forward[0] - reverse[0]) / TIMES[0] ** 4
        pairs[name] = contrast.tolist()
    return dict(
        dwell=TIMES[0],
        independent_current=float(current),
        probabilities_difference_over_t4=pairs,
    )


def run():
    # Both settings use forward exposures. The native binary contact spectrum is integral.
    theta = sp.pi / 3
    even, odd = (
        sp.exp(-8 * sp.I * theta).expand(complex=True),
        sp.exp(-9 * sp.I * theta),
    )
    diag = (even + odd) / 2 - 1
    off = (even - odd) / 2
    kick = sp.Matrix([[diag, off], [off, diag]])
    # Read kinetic energy above the homogeneous offset, with E=3 on this stationary sector.
    current = 3 * sp.Matrix([[0, -sp.I / 2], [sp.I / 2, 0]])
    board_rows = []
    for level in range(3):
        board = federation(level)
        bridges = sorted((a, b) for a, b in board.external if a < 12)
        gateways = [a for a, _ in bridges]
        u, v = exact_modes(cell(), gateways)
        assert (
            5 * sp.eye(12) - sp.Matrix(cell().adjacency().toarray().astype(int))
        ) * u == 6 * u
        ratios = []
        for g, r in bridges:
            assert set(board.neighbours[r]).intersection(range(12)) == {g}
            assert not set(EDGE).intersection([g])
            path = sp.Matrix([[int(p in board.neighbours[g]) for p in EDGE]])
            leading = path * kick / 8
            # The reverse pulse conjugates the kick. Global second-order path sign cancels.
            effect = leading.H * leading
            difference = sp.simplify(effect - effect.conjugate())
            ratio = sp.simplify(difference[0, 1] / current[0, 1])
            assert sp.simplify(difference - ratio * current) == sp.zeros(2)
            ratios.append(str(ratio))
        # Every stationary internal number current has zero divergence at every port.
        psi = (u + sp.I * v) / sp.sqrt(2)
        div = [sp.S.Zero] * 12
        for a, b in cell().seams:
            j = 3 * sp.im(sp.conjugate(psi[a]) * psi[b]).simplify()
            div[a] += j
            div[b] -= j
        assert all(sp.simplify(d) == 0 for d in div)
        board_rows.append(
            dict(
                cells=board.cells,
                gateways=gateways,
                receivers=[b for _, b in bridges],
                currents_to_optical_gains=ratios,
            )
        )
    return dict(
        theta="pi/3",
        alternate="5*pi/3",
        seam=EDGE,
        boards=board_rows,
        finite_mixed_320=finite_replay(board, bridges, u, v),
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
