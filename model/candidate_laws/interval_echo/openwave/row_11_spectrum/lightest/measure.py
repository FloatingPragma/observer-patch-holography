"""Certify the lowest net-minus-one formation line on each deposited blank board."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components

from model.executable_model import federation

HERE = Path(__file__).resolve().parent
TICKS = 10**8


def double_h_times(board, v):
    """Exact integer action of twice the host-relative vacancy Hamiltonian."""
    return [
        (16 * (2 * len(ns) - 6) + len(ns)) * v[i] - sum(v[j] for j in ns)
        for i, ns in enumerate(board.neighbours)
    ]


def make_witness(board):
    """Floating-point discovery only. The independent check uses integral arithmetic."""
    from scipy.sparse import diags, eye
    from scipy.sparse.linalg import eigsh, spsolve

    degree = board.degree
    h = diags(8 * (2 * degree - 6) + degree / 2) - board.adjacency() / 2
    guess = eigsh(
        h, k=1, which="SA", tol=1e-11, v0=np.ones(board.n), return_eigenvectors=False
    )[0]
    shift = (np.floor(guess * TICKS) - 1) / TICKS
    # The positive resolvent of 1 gives each component an absolute slack of 1,
    # even where ground amplitudes are extremely small. Then round once.
    trial = spsolve((h - shift * eye(board.n)).tocsc(), np.ones(board.n))
    return np.rint(128 * trial).astype(np.int64)


def certify(board, witness):
    assert witness.shape == (board.n,) and witness.dtype == np.int64
    assert connected_components(board.adjacency(), directed=False)[0] == 1
    v = witness.tolist()  # Python integers for all products and dot products.
    assert min(v) > 0
    twice_hv = double_h_times(board, v)
    lower = min(Fraction(y, 2 * x) for x, y in zip(v, twice_hv, strict=True))
    upper = Fraction(
        sum(x * y for x, y in zip(v, twice_hv, strict=True)), 2 * sum(x * x for x in v)
    )
    assert lower <= upper < Fraction(69, 2)
    low_tick = (TICKS * lower).__floor__()
    high_tick = (TICKS * upper).__ceil__()
    assert high_tick - low_tick <= 3
    return dict(
        cells=board.cells,
        ports=board.n,
        exact_componentwise_lower=str(lower),
        exact_rayleigh_upper=str(upper),
        ground_interval_ticks=[low_tick, high_tick],
        compact_minus_ground_interval_ticks=[
            35 * TICKS - high_tick,
            35 * TICKS - low_tick,
        ],
        port_probability_support=board.n,
        multiplicity=1,
    )


def run():
    # Every other composition is excluded by the parent's exact 3^12
    # local-word certificate and the PSD exchange inequality for D >= 9.
    prior = json.loads((HERE.parent / "bounds.json").read_text())
    floors = prior["few_defect_potential_floors"]
    few = {d: 8 * floors[str(d)] for d in (3, 5, 7)}
    many = 6 * 9 - 12
    assert min(*few.values(), many) > Fraction(69, 2)
    with np.load(HERE / "positive_witnesses.npz", allow_pickle=False) as witnesses:
        rows = [certify(federation(i), witnesses[str(i)]) for i in range(3)]
    return dict(
        source="blank memories; unit-filled host; total signed load -1; all contacts free",
        native_units_per_tick=f"1/{TICKS}",
        other_composition_floors={
            "3": few[3],
            "5": few[5],
            "7": few[7],
            "9_or_more": many,
        },
        boards=rows,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="rediscover and replace the integer witnesses",
    )
    args = parser.parse_args()
    if args.write:
        np.savez_compressed(
            HERE / "positive_witnesses.npz",
            **{str(i): make_witness(federation(i)) for i in range(3)},
        )
    results = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(results, indent=2) + "\n")
    else:
        assert results == json.loads(path.read_text())
    for row in results["boards"]:
        print(
            f"{row['cells']} cells: {row['ground_interval_ticks']} in units of 1/{TICKS}"
        )


if __name__ == "__main__":
    main()
