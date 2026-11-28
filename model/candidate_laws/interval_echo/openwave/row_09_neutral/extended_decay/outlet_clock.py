"""Read every positive low line after a prompt or pre-exposed inlet."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from . import outlet as base
from .measure import HERE
from ..passive_decay.measure import compare_records

LABELS = 32768


def ready_inlet(board, path, cell, edges, sectors):
    ports = tuple(sorted((path[1], *cell)))
    _, h = base.binary_generator(board, ports, [*edges, (path[1], path[2])], 1)
    values, v = eigh(h)
    initial = np.zeros(len(ports))
    for i, p in enumerate(cell):
        initial[ports.index(p)] = sectors[1][2][i, 0]
    phases = np.exp(-2j * np.pi * values[:, None] * np.arange(LABELS) / base.PERIOD)
    waves = (v.T @ initial)[:, None] * phases
    rho = v @ (waves @ waves.conj().T / LABELS) @ v.T
    weights, states = eigh(rho)
    assert min(weights) > -1e-13
    mixture = states * np.sqrt(np.maximum(0, weights))
    np.testing.assert_allclose(mixture @ mixture.conj().T, rho, atol=1e-13)
    return ports, {"prompt": initial[:, None], "pre_exposed": mixture}


def acquire():
    board, dimers, path, cell, edges = base.outlet_geometry()
    sectors = base.receiving_cell(board, cell, edges)
    ports, preps = ready_inlet(board, path, cell, edges, sectors)
    stocks, gap = base.source_stocks()
    _, one, v = sectors[1]
    response = base.clock(
        one, one[0], gap / 10, labels=LABELS, low_edge=base.PERIOD / LABELS
    )
    allowed = (one - one[0] > 1e-8) & (one - one[0] < gap / 10)
    groups = np.r_[0, np.flatnonzero(np.diff(one) > 1e-8) + 1]
    result = dict(
        protocol_sha256=sha256((HERE / "OUTLET_CLOCK.md").read_bytes()).hexdigest(),
        parent_gap=gap,
        labels=LABELS,
        period=base.PERIOD,
        waits=list(base.WAITS),
        false_light_ceiling=float(response[~allowed, 1].max()),
        low_line_acceptance_floor=float(response[allowed, 1].min()),
        sharp_address_clock_other_light_heavy=((v**2) @ response).round(13).tolist(),
        truly_light_line_indices=np.flatnonzero(allowed[groups]).tolist(),
        panels={},
    )
    for name, ready in preps.items():
        bridge = (abs(ready[ports.index(path[1])]) ** 2).sum()
        ebit = int(path[0] == dimers[-1][0])
        change = sum(
            sign * (abs(stocks["native", arm][:, :, :, ebit]) ** 2).sum()
            for sign, arm in ((1, "pulse"), (-1, "idle"))
        )
        panels = {}
        result["panels"][name] = dict(
            ready_load_in_cell=float(1 - bridge),
            first_shutter_work_contrast=float(8.5 * (1 - 2 * bridge) * change),
            settings=panels,
        )
        for closed in (False, True):
            words, h = base.material(board, dimers, path, cell, edges, closed)
            values, vectors = eigh(h)
            readers = base.read_groups(words, sectors, cell)
            for (setting, arm), guide in stocks.items():
                if closed and setting != "native":
                    continue
                detector = "outlet_cut" if closed else setting
                spectral = vectors.T @ base.inject(guide, words, dimers, ports, ready)
                rows = []
                for t in base.WAITS:
                    state = vectors @ (np.exp(-1j * t * values)[:, None] * spectral)
                    rows.append(base.flight_read(state, readers, sectors, response))
                panels.setdefault(detector, {})[arm] = rows
    return result


def decode(data):
    for name, bank in data["panels"].items():
        base.decode_panels(bank["settings"], data["waits"])
        native = bank["settings"]["native"]
        i = data["waits"].index(32)
        good = data["truly_light_line_indices"]

        def line(arm):
            return np.asarray(native[arm][i]["core_by_sector_energy"]["1"])[
                :, good
            ].sum()

        print(
            name,
            "ready in cell",
            bank["ready_load_in_cell"],
            "true-low excess",
            line("pulse") - line("idle"),
            "differential entrance work",
            bank["first_shutter_work_contrast"],
        )
    print(
        "heavy/ground false-light ceiling",
        data["false_light_ceiling"],
        "true low response floor",
        data["low_line_acceptance_floor"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    result = acquire()
    decode(result)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "outlet_clock.json").read_text()), result)
