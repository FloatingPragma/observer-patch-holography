"""Compare the transported modular current to the separate pooled area differential."""

import json
from pathlib import Path

from .. import signed_transfers
from ..dressed_charge import measure as charge
from ..full_algebra import station
from .measure import preparations


def run():
    states = preparations()
    chi = states["equilibrium"]
    product = states["pairing"]
    total_k = -station.log_blocks(chi)
    marginal_k = -station.log_blocks(product)
    collar_k = total_k - marginal_k
    eta = json.loads(Path(station.__file__).with_name("results.json").read_text())[
        "ordinary_eta"
    ]
    arms = json.loads(
        Path(signed_transfers.__file__).with_name("signed_transfers.json").read_text()
    )["tests"]
    board = station.iso.pairing.formation.federation(1)
    route = station.iso.pairing.formation.region(board)
    _, _, twice_m = charge.dark_cell(board, route)
    right_turn = (states["plus"] - states["minus"]) / 2
    times = (1 / 128, 1 / 256, 1 / 512)
    by_sign = {}
    for q in (-1, 1):
        hamiltonians = [
            charge.native_precision.generator(board, route, q, b) for b in (False, True)
        ]
        flights = [[charge.wave_at(h, 5, t) for h in hamiltonians] for t in times]
        out = {}
        for arm in arms:
            key = "ordinary" if arm["source"] == "ordinary" else "pairing"
            tangent = states[key] - chi + arm["angle_per_unit"] * right_turn
            curr = []
            for t, waves in zip(times, flights):
                edge = 0
                for b, (h, wave) in enumerate(zip(hamiltonians, waves)):
                    _, _, j, _ = charge.read_flux(
                        h, tangent[b], total_k[b], twice_m / 2, wave, board
                    )
                    edge += j[4, 5]
                curr.append(float(edge / t))
            predicted = float(3 * eta * arm["matched"][0] / 8)
            matter = -3 * float(station.pair(tangent, total_k)) / 8
            out[key] = dict(
                area_prediction=predicted,
                source_prediction=matter,
                marginal_contribution=-3 * float(station.pair(tangent, marginal_k)) / 8,
                collar_contribution=-3 * float(station.pair(tangent, collar_k)) / 8,
                population_current_per_time=curr,
            )
        by_sign[str(q)] = out
    return dict(
        launch_port=5, current_edge=[4, 5], times=list(times), eta=eta, charges=by_sign
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
