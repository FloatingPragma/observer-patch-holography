"""Free release generated from local record changes on the registered M10.3 board."""
import argparse
import base64
import json
import zlib
from math import acos, cos, pi, sqrt
from pathlib import Path

import numpy as np
from scipy import sparse

from model.candidate_laws.interval_echo.openwave.submission import row_04_check as pair
from model.executable_model import carrier, cell, federation
from model.experiments.interval_echo_encounter.check_transport import checked_steps
from .record_action import addresses, changes, signed_table


def source(arm):
    board = federation(0)
    edges, index, exchange, _ = pair.screen(carrier.HERE / "carriers/FEDERATION_LEVEL0.npz")
    if {tuple(e) for e in edges} != set(board.seams):
        raise ValueError("The M10.3 fixture differs from the native carrier")
    words = signed_table()
    primitive = cell()
    tile_edge = primitive.seams[0]

    def changed_count(layout, edge, transport):
        state = {v: int(v == edge[0]) for f in addresses(layout, edge) for v in f}
        return np.count_nonzero(changes(addresses(layout, edge), edge, state, words, transport))

    if arm == "assembled":
        weights = [changed_count(board, tuple(e), True) for e in edges]
    elif arm in ("primitive", "off", "double", "untransported"):
        weight = changed_count(primitive, tile_edge, arm != "untransported")
        weights = [dict(off=0, double=2*weight).get(arm, weight)] * len(edges)
    else:
        raise ValueError("Unknown record-change circuit")

    positions = np.array(np.triu_indices(board.n, 1)).T
    i, j = positions.T
    if not np.array_equal(index[i, j], np.arange(len(positions))):
        raise ArithmeticError("Pair ordering differs from the registered apparatus")
    onsite = np.zeros(board.n, dtype=int)
    diagonal = np.zeros(len(i), dtype=int)
    contact = np.zeros(len(i), dtype=int)
    for (a, b), w in zip(edges, weights, strict=True):
        onsite[[a, b]] += w
        occupied_a, occupied_b = (i == a) | (j == a), (i == b) | (j == b)
        diagonal += w * (occupied_a ^ occupied_b)
        contact[index[a, b]] = w
    reference = onsite.max() - onsite
    diagonal += reference[i] + reference[j] - 2*onsite.max()
    if not np.array_equal(diagonal, -2*contact):
        raise ArithmeticError("The one-load reference fails on binary pairs")
    h = exchange + sparse.diags(diagonal, format="csr")
    home = edges[np.all(edges//12 == 13, axis=1)]
    gateways = np.unique(edges[edges[:, 0]//12 != edges[:, 1]//12])
    launch = tuple(home[np.isin(home, gateways).any(axis=1)][0])
    if launch != (156, 157) or len(home) != 30:
        raise ValueError("The enrolled source or home observer changed")
    return h, index[launch], index[home[:, 0], home[:, 1]], reference


def release(arm):
    h, launch, home, reference = source(arm)
    start = np.zeros(h.shape[0], complex)
    start[launch] = 1
    speed = np.linalg.norm(h @ start - h[launch, launch]*start)
    grid = range(1, 297)
    samples = checked_steps(h/8, start, grid, degree=16)
    floor, appointments = 1., {}
    for wanted, (tick, psi, error) in zip(grid, samples, strict=True):
        if tick != wanted or psi.shape != start.shape or not 0 <= error < 1e-5:
            raise ArithmeticError("Incomplete grid or excessive propagation budget")
        norm = float(np.vdot(psi, psi).real)
        if not abs(norm-1) <= error:
            raise ArithmeticError("Norm outside the numerical enclosure")
        local = abs(psi[home])**2
        floor = min(floor, float(local.sum())-error)
        if tick in (8, 56, 296):
            appointments[tick//8] = [*local.tolist(), float(norm-local.sum())]
    angle = acos(sqrt(max(0., floor))) + speed/16
    continuous = float(cos(angle)**2) if angle < pi/2 else 0.
    result = dict(arm=arm, continuous_floor=continuous, threshold_met=continuous > .94,
                  variance=float(speed**2), reference_span=int(reference.max()), appointments=appointments)
    if arm in ("primitive", "off"):
        root = Path(pair.__file__).parent
        data = json.loads((root/"row_04_receipt.json").read_text())
        arm_index = int(arm == "off")
        archived = data["arms"][arm_index]["readings"]
        discrepancies = []
        for j, (time, row) in enumerate(zip(pair.EXPOSURES, archived, strict=True)):
            p = np.array(appointments[time])
            if row["exposure"] != time or np.min(p) < 0:
                raise ValueError("Invalid terminal law")
            discrepancies.append(float(np.max(abs(p-row["probabilities"]))))
            rng = np.random.default_rng([917421, arm_index, j])
            draws = rng.choice(len(p), size=pair.SHOTS, p=p/p.sum()).astype(np.uint8)
            if draws.tobytes() != zlib.decompress(base64.b64decode(row["trials_zlib_base64"], validate=True)):
                raise ValueError("Changed registered elementary outcomes")
        if max(discrepancies) >= 5.1e-10:
            raise ArithmeticError("Complete distributions differ beyond archive rounding")
        result["max_archived_difference"] = max(discrepancies)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arm", choices=("primitive", "assembled", "off", "double", "untransported"))
    print(json.dumps(release(parser.parse_args().arm), allow_nan=False))
