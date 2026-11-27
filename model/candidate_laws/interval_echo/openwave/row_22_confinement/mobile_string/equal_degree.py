"""Static Ramsey work on shifted dimer pairings of the deposited material routes."""
import argparse
from pathlib import Path

import numpy as np

from . import commission, material as mat, measure

RECEIVERS = Path(__file__).with_name("equal_degree.npz")
SEED = 220032


def ports(level):
    source = mat.route(level)
    board = mat.federation(level)
    external = dict((a, b) for edge in board.external for a, b in (edge, edge[::-1]))
    sites = source[1:] + [external[source[-1]]]
    if (len(set(sites)) != 32 or
            any(b not in board.neighbours[a] for a, b in zip(sites, sites[1:])) or
            any(board.degree[a] != board.degree[b] for a, b in zip(sites[::2], sites[1::2]))):
        raise ValueError("shifted route is not a simple native chain of 16 equal-degree dimers")
    return sites


def pinned_matrix(board, sites, length, pin, calibration):
    """Sum connector and private-port energies on physical occupations, then add dimer swaps."""
    contact, gap, private = calibration
    dimers = np.asarray(sites).reshape(-1, 2)
    first = (len(dimers) - length) // 2
    if length < 2 or first < 1 or first + length >= len(dimers):
        raise ValueError("both outer pins must lie on the material")
    block = dimers[first - 1:first + length + 1]
    inside = mat.bits(length - 2)
    zeros = np.zeros(len(inside), dtype=int)
    pins = np.full(len(inside), pin)
    upper = np.column_stack((zeros, pins, inside, pins, zeros))
    # The upper member is on the right of even dimers and the left of odd ones.
    on_left = np.where(np.arange(first - 1, first + length + 1) % 2, upper, 1 - upper)
    on_right = 1 - on_left
    energy = contact * (1 - on_right[:, :-1] * on_left[:, 1:]).sum(axis=1)
    energy += private * (on_left @ (6 - board.degree[block[:, 0]]) +
                         on_right @ (6 - board.degree[block[:, 1]]))
    return mat.generate(energy, length - 2, gap / 2, False)


def acquire(path):
    rng = np.random.default_rng(SEED)
    archive = {}
    for level in range(3):
        board = mat.federation(level)
        geometry = ports(level)
        stock = commission.material_stock(level, 16, geometry)
        archive[f"ports_{level}"] = geometry
        archive[f"stock_sha_{level}"] = stock["ordered_address_sha256"]
        for kappa in (8, 0):
            cal = mat.calibration_energies(kappa)
            archive[f"cal_{level}_{kappa}"] = measure.sample_ports(rng, measure.quadratures(cal * measure.CAL_HOLDS))
            archive[f"clock_{level}_{kappa}"] = measure.sample_ports(rng, measure.quadratures(cal[1] * measure.TICK))
            for length in mat.LENGTHS:
                h0, h1 = (pinned_matrix(board, geometry, length, pin, cal) for pin in (0, 1))
                _, panels = mat.pinned_laws(h0, h1)
                # Intact and screened preparations, each including both quadratures and all holds.
                archive[f"static_{level}_{kappa}_{length}"] = np.stack(
                    [measure.sample_ports(rng, signal) for signal, _ in panels]
                )
    np.savez_compressed(path, **archive)


def decode(path=RECEIVERS):
    out = {}
    with np.load(path, allow_pickle=False) as archive:
        for level in range(3):
            for kappa in (8, 0):
                clock = np.polyfit(measure.TICK, np.unwrap(measure.port_phase(archive[f"clock_{level}_{kappa}"])), 1)[0]
                cal = measure.port_phase(archive[f"cal_{level}_{kappa}"]) / (clock * measure.CAL_HOLDS)
                readings = {}
                for length in mat.LENGTHS:
                    counts = archive[f"static_{level}_{kappa}_{length}"]
                    readings[length] = {name: measure.ramsey(panel, clock) for name, panel in
                                        zip(("sharp_intact", "sharp_screened"), counts, strict=True)}
                out[level, kappa] = cal, readings
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", metavar="NPZ", type=Path, help="write a fresh receiver bank")
    args = parser.parse_args()
    if args.acquire:
        acquire(args.acquire)
    for (level, kappa), (cal, readings) in decode(args.acquire or RECEIVERS).items():
        print(f"level {level}, κ={kappa}, contact/gap/private {cal.round(4).tolist()}")
        print("length, preparation, status, accepted line, attempted line, minimum visibility")
        for length, pair in readings.items():
            for name, word in pair.items():
                print(length, name, word["status"], word["line"], word.get("attempted_line"), word["min_visibility"])
