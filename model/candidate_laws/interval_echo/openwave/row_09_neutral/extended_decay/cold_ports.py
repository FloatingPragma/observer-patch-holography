"""Qualify both open interfaces on an externally measured cold receiving cell."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

from . import open_ports, outlet, transmission
from .measure import HERE
from .outlet_clock import LABELS
from ..passive_decay.measure import compare_records, hamiltonians


def zero_bin(groups, sectors):
    _, values, v = sectors[1]
    phases = (values - values[0]) / outlet.PERIOD
    amplitude = (
        np.exp(-1j * np.pi * (LABELS - 1) * phases)
        * np.sinc(LABELS * phases)
        / np.sinc(phases)
    )
    cell_kraus = (v * amplitude) @ v.T
    size = sum(len(rows) for _, rows, _ in groups)
    out = np.zeros((size, size), complex)
    for n, rows, _ in groups:
        if n == 1:
            out[np.ix_(rows, rows)] = cell_kraus
    return out


def acquire(
    protocol="COLD_PORTS.md", waits=outlet.WAITS, geometries=("single", "complete")
):
    board, dimers, cell, available, sites, contacts, cut = open_ports.aperture()
    _, _, path, _, internal = outlet.outlet_geometry()
    sectors = outlet.receiving_cell(board, cell, internal)
    _, one, cv = sectors[1]
    core, wing, *_ = hamiltonians(4)
    limit = np.ptp(eigh(core, eigvals_only=True)) / 10
    response = outlet.clock(
        one, one[0], limit, labels=LABELS, low_edge=outlet.PERIOD / LABELS
    )
    soft = (one - one[0] > 1e-8) & (one - one[0] < limit)
    ceiling = float(response[~soft, 1].max())
    _, v = eigh(wing)
    guide = np.zeros((4, 2, 2, 2))
    guide[0, 0] = v[:, 0].reshape(2, 2)
    result = dict(
        protocol_sha256=sha256((HERE / protocol).read_bytes()).hexdigest(),
        waits=list(waits),
        clock_false_light=ceiling,
        mounts={},
    )
    for geometry in geometries:
        settings = {}
        for name in ("native", "exit_cut", "no_reversal"):
            if geometry == "single":
                words, h = outlet.material(
                    board, dimers, path, cell, internal, name == "exit_cut"
                )
            else:
                kept = [e for e in contacts if name != "exit_cut" or e not in cut]
                words, h = outlet.connected_wing(board, dimers, sites, kept)
            if name == "no_reversal":
                h += np.diag(
                    [int(d == 1 and dimers[-1][0] in w) for d in (0, 1) for w in words]
                )
            groups = outlet.read_groups(words, sectors, cell)
            effects = transmission.outlet_effects(groups, sectors, response)
            initial = outlet.inject(guide, words, dimers, cell, cv[:, :1])[:, 0]
            settings[name] = open_ports.analyze_clock(
                h, initial, effects, limit, zero_bin(groups, sectors), waits=waits
            )
        result["mounts"][geometry] = settings
    return result


def verdict(data):
    passing = []
    idx = data["waits"].index(32)
    for geometry, settings in data["mounts"].items():
        native = settings["native"]
        admission = native["cold_preparation"]["admission"]
        cut = settings["exit_cut"]["fixed"]
        for name in open_ports.NAMES[1:]:
            base = native["fixed"]["idle"]
            pulse = native["fixed"][name]
            excess = (
                np.asarray(pulse["panels"])[:, 2] - np.asarray(base["panels"])[:, 2]
            )
            cut_excess = (
                np.array(cut[name]["panels"])[:, 2]
                - np.array(cut["idle"]["panels"])[:, 2]
            )
            np.testing.assert_allclose(cut_excess, 0, atol=3e-12)
            _, _, light, heavy = pulse["unresolved_loss_light_gain_heavy_gain"]
            if (
                admission * (excess[idx] - data["clock_false_light"]) > 0.01
                and pulse["work"] > 0
                and light > heavy
            ):
                passing.append((geometry, name))
            print(
                geometry,
                name,
                "admission",
                admission,
                "original launch excess",
                np.round(admission * excess, 7),
                "input work",
                pulse["work"],
            )
        print(
            geometry,
            "all-unitary maxima per launch",
            np.round(
                admission * np.array(native["all_D_unitaries_max_light_excess"]), 7
            ),
        )
        print(
            "joined ground one-load inventory",
            native["cold_preparation"]["static_ground_inventory_one"],
        )
    print("fixed qualifying controls", passing)
    return passing


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    result = acquire()
    verdict(result)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "cold_ports.json").read_text()), result)
