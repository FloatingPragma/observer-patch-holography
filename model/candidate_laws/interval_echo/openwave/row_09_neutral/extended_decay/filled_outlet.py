"""Test the heavy source on the unit-filled blank host and bound bare neutral defects."""

import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from . import joined_outlet as joined, open_ports, outlet, wing
from .measure import HERE
from ..passive_decay.material import columns
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.readers import Z, circuit
from ..joined_decay.measure import core_pulse

KAPPA = candidate.KAPPA


def pair_standard(board, cell, internal, light_limit):
    """All 132 ternary vacancy-excess basis words in the cold receiver."""
    pairs = list(permutations(cell, 2))
    lookup = {w: i for i, w in enumerate(pairs)}
    h = np.zeros((len(pairs), len(pairs)))
    potential = []
    for j, (vacancy, excess) in enumerate(pairs):
        loads = {vacancy: 0, excess: 2}
        contrast = 0
        for a, b in internal:
            x, y = loads.get(a, 1), loads.get(b, 1)
            contrast += KAPPA * (x - y) ** 2
            if x != y:
                swapped = tuple(
                    b if p == a else a if p == b else p for p in (vacancy, excess)
                )
                h[lookup[swapped], j] -= 0.5
                h[j, j] += 0.5
        cost = contrast + KAPPA * (
            -(6 - board.degree[vacancy]) + 3 * (6 - board.degree[excess])
        )
        expression = KAPPA * (
            10
            + 2 * int(excess in board.neighbours[vacancy])
            - (6 - board.degree[vacancy])
            + 3 * (6 - board.degree[excess])
        )
        assert cost == expression
        potential.append(cost)
    exchange_min = eigh(h, eigvals_only=True)[0]
    assert exchange_min >= -1e-13
    values, v = eigh(h + np.diag(potential))
    period, labels = 512, 32768
    bins = np.fft.fftfreq(labels) * period
    clock = (
        abs(
            np.fft.ifft(
                np.exp(-2j * np.pi * values[:, None] * np.arange(labels) / period),
                axis=1,
            )
        )
        ** 2
    )
    low = (bins >= period / labels) & (bins < light_limit)
    heavy = bins >= light_limit
    np.testing.assert_allclose(clock.sum(axis=1), 1, rtol=0, atol=2e-12)
    assert values.max() < period / 2
    sharp = (v**2) @ np.column_stack([clock[:, low].sum(1), clock[:, heavy].sum(1)])
    boards = []
    for level in (0, 1, 2):
        b = federation(level)
        assert min(b.degree) == 5 and max(b.degree) == 6
        hole = next(p for p in range(b.n) if b.degree[p] == 5)
        surplus = next(
            p for p in range(b.n) if b.degree[p] == 6 and p not in b.neighbours[hole]
        )
        # On the full native graph the potential difference from all unit loads
        # in this sector is k[12+2(d_v-d_e)+2*adj(v,e)].
        boards.append(
            dict(
                cells=b.cells,
                ports=b.n,
                full_contrast_min=10 * KAPPA,
                vacancy_excess_witness=[hole, surplus],
            )
        )
    return dict(
        pairs=len(pairs),
        isolated_potential_floor=min(potential),
        isolated_ground=float(values[0]),
        isolated_ceiling=float(values[-1]),
        all_line_false_light_upper=float(clock[:, low].sum(1).max()),
        smallest_sharp_heavy_probability=float(sharp[:, 1].min()),
        global_boards=boards,
    )


def system(name):
    b, d, cell, _, sites, seams, entrance = open_ports.aperture()
    if name == "bridge_cut":
        seams = [e for e in seams if e != (80, 192)]
    elif name == "exit_cut":
        seams = [e for e in seams if e not in entrance]
    words, tail = outlet.connected_wing(b, d, sites, seams, number=len(sites) - 1)
    holds = wing.installation()[2]
    h = joined.generator(
        holds[
            {"no_reversal": "no_reversal", "source_cut": "boundary_cut"}.get(
                name, "native"
            )
        ],
        tail,
        words,
        d,
        name if name != "bridge_cut" else "native",
    )
    hole = [next(p for p in sites if p not in w) for w in words]
    inside = np.isin(hole, cell)
    return h, d, np.asarray(hole), inside, words


def acquire():
    b, _, cell, _, _, seams, _ = open_ports.aperture()
    # Every blank-bond column is checked in binary_generator. These additional
    # written heads must also reduce the same binary subspace in each law.
    for comparison in (False, True):
        columns(comparison)
    occupied = tuple(cell)
    sectors = {
        len(cell) - 1: ([set(cell) - {vacancy} for vacancy in cell], None, None),
        len(cell): ([occupied], None, None),
    }
    _, pointers = outlet.relative_inventory_reader(
        b, cell, sectors, labels=16, reference=len(cell)
    )
    assert pointers == [15, 0]
    core = hamiltonians(4)[0]
    energies, frame = eigh(core)
    parent_gap = float(np.ptp(energies))
    read = circuit(frame, 32)
    pulse = circuit(core_pulse(frame), 32)
    heads = wing.installation()[2]
    _, cold = eigh(heads["native"])
    settings = {}
    for name in joined.SETTINGS:
        h, d, hole, inside, words = system(name)
        n = len(words)
        k = np.diag(np.tile(~inside, 2).astype(float))
        stock = np.zeros((4, 2, 2, n))
        first = cold[:, 0].reshape(4, 2, 2, 2)
        for e in (0, 1):
            slot = np.flatnonzero(hole == d[-1][e])[0]
            stock[:, :, :, slot] = first[:, :, :, e]
        ev, v, prepared, receipt = joined.prepare(h, stock.ravel(), k)
        rank = prepared.shape[1]
        arms = {}
        for arm in ("idle", "D_phase", "core_pulse"):
            launch = prepared.reshape(4, 2, 2, n, rank)
            if arm == "core_pulse":
                launch = np.einsum("ba,acdmr->bcdmr", pulse, launch)
            elif arm == "D_phase":
                launch = launch * Z.diagonal()[None, None, :, None, None]
            amplitudes = v.T @ launch.reshape(len(h), rank)
            panels = []
            for time in joined.WAITS:
                psi = (v @ (np.exp(-1j * time * ev)[:, None] * amplitudes)).reshape(
                    4, 2, 2, n, rank
                )
                p = abs(np.einsum("ab,acdmr->bcdmr", read, psi)) ** 2
                p = p.sum(axis=(1, 2, 4))
                # Each unchanged-inventory word is the unique all-one word on
                # the receiver. Its isolated clock has only the zero bin.
                panels.append(
                    np.column_stack([p[:, inside].sum(1), p[:, ~inside].sum(1)])
                    .round(13)
                    .tolist()
                )
            arms[arm] = panels
        settings[name] = dict(preparation=receipt, arms=arms)
    controls = {}
    waits = np.r_[0, 2 * np.pi / outlet.PERIOD, joined.WAITS[1:]]
    for name in ("native", "bridge_cut"):
        h, _, hole, inside, words = system(name)
        ev, v = eigh(h)
        sharp = np.zeros(len(h))
        sharp[15 * len(words) + np.flatnonzero(hole == 80)[0]] = 1
        stocks = v.T @ sharp
        trials = []
        for time in waits:
            psi = (v @ (np.exp(-1j * time * ev) * stocks)).reshape(16, len(words))
            trials.append(round(float((abs(psi[:, inside]) ** 2).sum()), 13))
        controls[name] = trials
    return dict(
        protocol_sha256=sha256((HERE / "FILLED_OUTLET.md").read_bytes()).hexdigest(),
        parent_gap=parent_gap,
        waits=joined.WAITS.tolist(),
        inventory_outcomes=["one_vacancy", "unchanged_filled_ground"],
        vacancy_calibration_waits=waits.tolist(),
        vacancy_calibration=controls,
        settings=settings,
        balanced_pair=pair_standard(
            b, cell, [e for e in seams if set(e) <= set(cell)], parent_gap / 10
        ),
    )


def decode(data):
    print("balanced comparison", data["balanced_pair"])
    print("vacancy arrival", data["vacancy_calibration"])
    for name, setting in data["settings"].items():
        p, i = (np.array(setting["arms"][a]) for a in ("core_pulse", "idle"))
        np.testing.assert_allclose(
            p.sum(axis=(1, 2)), setting["preparation"]["admission"], atol=2e-10
        )
        print(
            name,
            "vacancy excess",
            (p - i)[:, :, 0].sum(1),
            "lower core and vacancy",
            p[:, 1:3, 0].sum(1),
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
        compare_records(json.loads((HERE / "filled_outlet.json").read_text()), result)
