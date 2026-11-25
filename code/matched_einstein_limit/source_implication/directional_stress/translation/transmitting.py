"""Smallest cell-local reducing sector containing all gateway channels."""

import argparse
import json
from collections import Counter
from functools import lru_cache

import numpy as np
import sympy as sp

from .. import bridge
from . import measure


@lru_cache(maxsize=None)
def dark_projector(action_bytes, gates):
    h = sp.Matrix(np.frombuffer(action_bytes, dtype=np.int64).reshape(48, 48))
    pick = [4 * p + i for p in gates for i in range(4)]
    old = sp.eye(48)[:, pick]
    added = old
    while True:
        offered = old.row_join(h * added)
        piv = offered.to_DM().rref()[1]
        extra = [j for j in piv if j >= old.cols]
        if not extra:
            break
        added = offered[:, extra]
        old = old.row_join(added)
    dark = old.T.nullspace()
    proj = sp.zeros(48)
    if dark:
        box = sp.Matrix.hstack(*dark)
        proj = box * (box.T * box).inv() * box.T
    assert proj * proj == proj and proj.T == proj and h * proj == proj * h
    assert proj[:, pick] == sp.zeros(48, len(pick))
    return proj


def run():
    atlas = measure.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    patch = measure.chosen(atlas)
    launches = sum((list(v) for v in bridge.beam.screens(board, route)), [])
    states = bridge.action.preparations()
    rows = {}
    s = measure.map_matrix(board, patch)
    p = route[1]
    dest = int(s[:, 4 * p].nonzero()[0][0]) // 4
    assert board.degree[p] == board.degree[dest] == 6
    pair = [4 * p, 4 * dest]
    shift = (s - s.T)[pair, :][:, pair].toarray() / (2j)
    for sign in (-1, 1):
        amp = np.array([1, sign * 1j]) / np.sqrt(2)
        np.testing.assert_allclose(np.vdot(amp, shift @ amp), -sign / 2, atol=2e-15)
    sector = []
    for q in (-1, 1):
        effects = []
        for odd in (False, True):
            h, ops = measure.op_table(board, route, patch, q, odd)
            force = ops["marked_force"][0].tocoo()
            assert np.all(force.row // 48 != force.col // 48)
            sizes = Counter()
            e = sp.zeros(4)
            projectors = []
            for c in range(board.cells):
                loc = h[48 * c : 48 * c + 48, 48 * c : 48 * c + 48].toarray()
                gates = tuple(
                    int(v)
                    for v in np.flatnonzero(board.degree[12 * c : 12 * c + 12] == 6)
                )
                dark = dark_projector(loc.tobytes(), gates)
                projectors.append(dark)
                sizes[48 - int(dark.trace())] += 1
                for point in launches:
                    if point // 12 == c:
                        i = 4 * (point % 12)
                        e += sp.eye(4) - dark[i : i + 4, i : i + 4]
            for c, target, perm in zip(patch["cells"], patch["images"], patch["lifts"]):
                local_t = sp.kronecker_product(sp.eye(12)[:, list(perm)], sp.eye(4))
                assert projectors[target] * local_t == local_t * projectors[c]
            for point in route:
                d = projectors[point // 12]
                assert d[:, 4 * (point % 12) : 4 * (point % 12) + 4] == sp.zeros(48, 4)
            assert e / 10 == sp.Rational(4, 5) * sp.eye(4)
            effects.append(np.asarray(e / 10, dtype=float))
            sector.append(
                dict(
                    charge=q,
                    odd=odd,
                    cell_transmitting_dimensions=dict(sorted(sizes.items())),
                    transmitting_dimensions=sum(k * v for k, v in sizes.items()),
                    full_dimensions=h.shape[0],
                )
            )
        rows[str(q)] = bridge.evaluate(effects, states)
    return dict(
        shift_witness_gateway_ports=[p, dest],
        shift_contrast=1,
        energy_current_contrast=0,
        native_trace_H7_force_restricted_to_transmitting_sector="1/8",
        registered_transmitting_weight_exact="4/5",
        boundary_and_marked_force_supported_on_transmitting_sector=True,
        sectors=sector,
        conserved_transmitting_weight=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "transmitting.json").write_text(data)
    else:
        print(data, end="")
