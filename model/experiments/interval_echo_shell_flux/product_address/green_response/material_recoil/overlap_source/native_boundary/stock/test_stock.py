"""Independent checks of the conservation law, stationary classes and priced trail."""

from fractions import Fraction as F
from pathlib import Path
import json

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation

HERE = Path(__file__).resolve().parent


def test_fuel_and_all_passive_magnitude_projectors():
    mat = law.echo()
    h0, h1 = law.local_term_parts()
    mem = np.repeat([0, *[abs(d) for _, d in law.LABELS]], 9)
    for value in [0, 2]:
        fuel = np.tile([p.count(value) for p in law.SRC], 11) + (mem == 2)
        off = fuel[:, None] != fuel
        for op in [mat, h0, h1]:
            assert np.max(abs(op[off])) < 1e-12
    for m in [0, 1, 2]:
        off = (mem[:, None] == m) != (mem == m)
        for op in [h0, h1]:
            assert np.max(abs(op[off])) < 1e-12


def test_known_stationary_resources_and_parked_full_board():
    summary = json.loads((HERE / "results.json").read_text())["panels"]
    p = json.loads((HERE / "parked.json").read_text())
    assert len(p["all81_parked_inputs"]) == 81
    for name, n in [("echo", 81), ("echo_plus_exchange", 15), ("exchange_only", 15)]:
        assert len(summary[name]["classes"]) == n
    mixed = summary["echo_plus_exchange"]
    for c in mixed["classes"]:
        assert c["stationary_for_direct_path_free_hold"]
        for flow in c["repair_flows"].values():
            assert F(flow["birth_per_cycle"]) == F(flow["return_per_cycle"])
        assert c["extrema"]["2"][1] <= min(c["fuel"])
    c = mixed["classes"][mixed["initial_to_class"]["1100"]]
    assert c["size"] == 666 and F(c["resident"]["1"]) == F(54, 37)
    r = p["all81_parked_inputs"]["1100"]["echo_plus_exchange"]
    assert F(r["work_from_original"]) == F(476, 111)
    r = p["all81_parked_inputs"]["2200"]["echo_plus_exchange"]
    assert (
        F(r["actual_minus_matched_blank"]) > 21
        and F(r["relative_to_same_fuel_exchange_only"]) < -74
    )


def test_actual_native_path_energy_and_cost_of_read():
    data = json.loads((HERE / "trail.json").read_text())
    h0, h1 = law.local_term_parts()
    H = h0 + 8 * h1
    u = law.echo()
    for name in ["forward", "reverse"]:
        cfg = data[name]
        ids = cfg["basis"]
        visits = [F(x) for x in cfg["visits"]]
        out = []
        for x in ids[:2]:
            psi = u[:, x]
            out.append((abs(psi) ** 2) @ H.diagonal() - psi @ H @ psi)
        assert abs(np.array(list(map(float, visits))) @ out - 4.25) < 1e-12
    for n, row in enumerate(data["boards"]):
        route = row["path"]
        car = federation(n)
        keys = {tuple(sorted((a, b))): (a, b) for a, b in zip(route, route[1:])}
        assert len(keys) == row["live_seams"] and set(keys) <= set(car.seams)
        ev = 0.0
        for edge in car.seams:
            a, b = edge
            load = (int(a == route[-1]), int(b == route[-1]))
            label = 0
            if edge in keys:
                left, right = keys[edge]
                coin, d = (0, 1) if left < right else (1, -1)
                label = 1 + law.LABELS.index((coin, d))
            x = label * 9 + law.SRC.index(load)
            ev += H[x, x]
        ev += 8 * (6 - car.degree[route[-1]])
        assert abs(ev - float(F(row["simultaneous_terminal_energy"]))) < 1e-9
        assert F(row["net_work"]) == F(row["work_from_sharp_reads"]) + F(
            row["work_from_addressed_echoes"]
        )
        assert row["passive_excess_vs_any_one_load_blank"] == 8 * (len(keys) - 2) - 18
