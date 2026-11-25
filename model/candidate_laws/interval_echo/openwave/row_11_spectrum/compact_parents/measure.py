"""Native spectroscopy and spatial phase impulses on charged compact product sources."""

from argparse import ArgumentParser
from functools import cache
from math import pi
import json

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.experiments.interval_echo_lepton_register import measure as rails
from ...row_04_compact import material_reference as store
from . import spectrum as spec
from . import metrology

SHOTS = 2**18
WAITS = (0, 0.5, 1, 2, 4, 8)


@cache
def wave(gates, sign, tag, sharp):
    """Native routed splits, followed by the signed private phase on each occupied interior port."""
    target = spec.orb(gates, tag)
    ports = np.flatnonzero(abs(target) > 1e-12)
    source = int(ports[0])
    u = np.eye(12, dtype=complex)
    law = store.LAW
    idx = [law.SRC.index(q) for q in ((1 + sign, 1), (1, 1 + sign))]
    block = store.anchor.free_seam(store.K)[np.ix_(idx, idx)]
    np.testing.assert_allclose(block, [[8.5, -0.5], [-0.5, 8.5]])

    def contact(a, b, t):
        u[[a, b]] = expm(-1j * t * block) @ u[[a, b]]

    trail = rails.route(sharp, source)
    for a, b in zip(trail, trail[1:]):
        contact(a, b, pi)
    for dest in ports[1:]:
        trail = rails.route(int(dest), source)
        swaps = list(zip(trail[:-2], trail[1:-1]))
        for a, b in swaps:
            contact(a, b, pi)
        split = abs(target[dest]) / abs(u[source, sharp])
        # The final beam can use up the entire amplitude to roundoff.
        if split > 1 + 1e-11:
            raise ArithmeticError("Unavailable beam weight")
        contact(source, trail[-2], 2 * np.arcsin(min(1.0, split)))
        for a, b in reversed(swaps):
            contact(a, b, pi)
    rate = store.K * ((1 + sign) ** 2 - 1)
    for j in ports:
        assert j not in gates
        duration = ((np.angle(u[j, sharp]) - np.angle(target[j])) / rate) % (
            2 * pi / abs(rate)
        )
        u[j] *= np.exp(-1j * rate * duration)
    np.testing.assert_allclose(u[:, sharp], target, rtol=0, atol=2e-12)
    return u


def layout(board, charts, candidate):
    special = {i for i, g in enumerate(charts) if 1 in spec.modes(g)}
    # Choose an existing seam with two five-mode owners. Their two gateway
    # endpoints are the splitter and recombiner, so no spectator lies on the route.
    endpoints = next(
        (a, b) for a, b in board.external if a // 12 in special and b // 12 in special
    )
    home, remote = (x // 12 for x in endpoints)
    others = sorted(special - {home, remote})
    spectators = []
    for tag, key in [(1, "plus_branch"), (-1, "minus_branch")]:
        k = candidate[key]
        spectators.extend([(i, tag) for i in others[:k]])
        others = others[k:]
    available = sorted(
        set(range(board.cells)) - {home, remote} - {i for i, _ in spectators}
    )
    spectators.extend((i, 0) for i in available[: candidate["ordinary_branch"] - 1])
    assert len(spectators) == 2 * candidate["pairs"]
    records = []
    for j, (i, tag) in enumerate(spectators):
        sign = (-1) ** j
        gate = charts[i][0]
        wave(charts[i], sign, tag, gate)
        records.append(
            dict(
                owner=i, signed_load=sign, compact_branch=tag, sharp_port=12 * i + gate
            )
        )
    return endpoints, records


def native_home(gates, sign):
    _, lap, _ = spec.compact.blocks(gates, "single", 0)
    degree = np.diag(lap)
    potential = store.K * (degree + ((1 + sign) ** 2 - 1) * (6 - degree))
    return lap / 2 + np.diag(potential)


def observe(board, charts, ends):
    rng = np.random.default_rng(171129)
    a, b = ends
    home, remote = a // 12, b // 12
    split, combiner = store.optics(-1)
    source = [a % 12, b % 12]
    records = {}
    for tag in (0, -1, 1):
        families = [0, tag]
        us = [
            wave(charts[i], -1, k, j)
            for i, k, j in zip((home, remote), families, source, strict=True)
        ]
        read = {}
        phased = {}
        for t in WAITS:
            arms = []
            for start, i, j, u in zip(split, (home, remote), source, us, strict=True):
                evolved = expm(-1j * t * native_home(charts[i], -1)) @ u[:, j]
                assert np.linalg.norm(evolved[list(charts[i])]) < 3e-12
                arms.append(u.conj().T @ (start * evolved))
            q = []
            for angle in (0, pi / 2, pi):
                # A uniform phase on the second owner's occupied interior
                # rails is a native signed-load impulse. On these words it
                # commutes with that owner's free restricted Hamiltonian.
                phase = np.exp(1j * store.K * (angle / store.K))
                pair = combiner @ np.array(
                    [arms[0][source[0]], phase * arms[1][source[1]]]
                )
                escaped = sum(
                    np.vdot(z, z).real - abs(z[j]) ** 2
                    for z, j in zip(arms, source, strict=True)
                )
                prob = np.maximum(np.r_[abs(pair) ** 2, escaped], 0)
                np.testing.assert_allclose(sum(prob), 1, atol=5e-12)
                q.append(
                    rng.multinomial(SHOTS, np.round(prob / prob.sum(), 12)).tolist()
                )
            read[str(t)] = q
            if t == 1:
                for k, arr in enumerate(q[:2]):
                    phased[k] = (arr[0] - arr[1]) / SHOTS
        # X and pi/2 reads give cos(delta t), sin(delta t) in this port convention.
        signed_gap = float(np.arctan2(phased[1], phased[0]))
        predicted = 0 if tag == 0 else (-1 + tag * np.sqrt(5)) / 2
        assert abs(signed_gap - predicted) < 0.012
        records[str(tag)] = dict(
            predicted_difference=predicted,
            signed_gap_from_time_one=signed_gap,
            bright_dark_else_by_wait=read,
        )
    return records


def run():
    ground, targets, uncertainty = spec.datasets()
    floors = []
    bs = None
    for level, g in enumerate(ground):
        board = federation(level)
        charts = spec.compact.gate_profiles(board)
        count = sum(1 in spec.modes(c) for c in charts)
        floors.append(
            spec.scan(board.cells, count, g["ground_interval_ticks"], targets)
        )
        if board.cells == 320:
            bs = board, charts
    board, charts = bs
    selected = floors[-1]["comparisons"]["ordinary_compact"]["measured_muon"]
    ends, spectators = layout(board, charts, selected)
    ins = observe(board, charts, ends)
    parent = spectators + [
        dict(owner=ends[0] // 12, signed_load=-1, compact_branch=0, sharp_port=ends[0])
    ]
    direct = metrology.measure(board, charts, parent, ends[1], native_home, wave)
    ci = direct["simultaneous_interval_99"]
    assert ci[0] <= selected["energy"] <= ci[1]

    n = selected["pairs"]
    plus = selected["plus_branch"]
    minus = selected["minus_branch"]
    baseline = json.loads(
        (store.anchor.experiment.HERE / "material_reference.json").read_text()
    )["native"][0]["panels"]["6"]
    inferred = (
        baseline["vacancy"]["frequency_per_clock_gap"]
        + n * baseline["opposite_signs"]["frequency_per_clock_gap"]
    )
    for k, count in [(1, plus), (-1, minus)]:
        inferred += count * ins[str(k)]["signed_gap_from_time_one"]

    def count_foreign():
        used = {r["owner"] for r in spectators} | {ends[0] // 12}
        return sum(a // 12 in used and b // 12 in used for a, b in board.external)

    return dict(
        material="filled blank host, kappa=8; signed load minus one",
        targets=targets,
        measured_muon_ratio_standard_uncertainty=uncertainty,
        levels=floors,
        held_source=dict(
            selection="nearest measured muon ratio / ordinary compact vacancy on 320 owners",
            declared_arms=list(map(int, ends)),
            spectator_recipe=spectators,
            occupied_neighbor_seams=count_foreign(),
            independent_inferred_formation=inferred,
            storage_baselines={
                k: baseline[k]["frequency_per_clock_gap"]
                for k in ("vacancy", "opposite_signs")
            },
        ),
        complete_storage_scan=direct,
        shots_per_apparatus=SHOTS,
        owner_phase_impulses=[0, "pi/2", "pi"],
        spectroscopy=ins,
    )


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = run()
    out = json.dumps(data, indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (spec.HERE / "results.json").write_text(out)
    else:
        print(out, end="")


if __name__ == "__main__":
    main()
