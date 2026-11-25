"""Native four-tick force coefficients in the 18 fixed-marginal source coordinates."""

import numpy as np
import sympy as sp

from ...full_algebra import certificate as cert
from .. import bridge
from ..traction.certificate import four_tick_columns
from ..translation import certificate as regional
from ..translation import measure


def force_effects():
    """Return ((C_v,even,C_v,odd),(C_e,even,C_e,odd)); each block has operator error <2e-15."""
    atlas = measure.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    launch = bridge.beam.indices(sum(bridge.beam.screens(board, route), []))
    patch = measure.chosen(atlas)
    weights = bridge.native.prior.weights()[7]
    pairs = []
    for q in (-1, 1):
        pair = []
        for parity in (False, True):
            h, entries = measure.op_table(board, route, patch, q, parity)
            numerator, unit = entries["marked_force"]
            assert unit == 4 and abs(numerator).sum(axis=1).max() < 200 * unit
            idx = np.unique(numerator.nonzero()[0])
            w = four_tick_columns(h, idx, launch, weights)
            block = regional.compress(w, numerator[idx, :][:, idx], unit)
            pair.append(block)
        pairs.append(tuple(pair))
    return tuple(pairs)


def force_gradients():
    """Return rational central 2x18 gradients (vacancy, excess) and an L2 error per row."""
    basis = [cert.prior.sym_matrix(a) for a in cert.iso.basis(True)]
    rows = [
        [cert.inner(block, a) for block in pair for a in basis]
        for pair in force_effects()
    ]
    # The independent regional generator certificate gives <=2e-15 operator
    # error per 4x4 head block. Projection on orthonormal B_i is contractive.
    return sp.Matrix(rows), sp.Rational("6e-15")
