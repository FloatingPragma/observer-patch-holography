"""Test the measured rate of whole-cell energy imbalance against the held camera."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from mpmath import iv
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from . import matched_readings, measure, observer_records as camera

HERE = Path(__file__).parent
RECEIPT = HERE / "regional_rate.json"
STOCKS = ("ordinary", "matched", "exchange")
SIGNS = ("minus", "plus")
CENTER = 4.25
SHOTS, CLOCK_SHOTS = 2**78, 2**68
SCALE = 96 * 32


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def preparations():
    ref, prior, *_ = camera.BEAM.aperture.source_and_changes()
    matched = json.loads(matched_readings.RECEIPT.read_text())["source"]
    basis = measure.joint_basis()
    # Odd-sector X.X + Y.Y exchanges |01> and |10>.
    exchange = (basis[9] + basis[13]) / np.sqrt(2)
    np.testing.assert_allclose(measure.marginals(exchange), 0, atol=1e-14)
    directions = {
        "ordinary": prior["ordinary"],
        "matched": np.einsum("i,ibjk->bjk", matched["selected_coefficients"], basis),
        "exchange": exchange,
    }
    steps = {"ordinary": 1 / 2048, "matched": matched["half_step"], "exchange": 1 / 128}
    recipes, states = {}, {}
    for name, direction in directions.items():
        for _ in range(21):
            endpoints = [ref + sign * steps[name] * direction for sign in (-1, 1)]
            if min(np.linalg.eigvalsh(x).min() for x in endpoints) > 0:
                prepared = [measure.station.preparation(x) for x in endpoints]
                if min(min(p["min_kept"], p["min_dephased"]) for p in prepared) > 0:
                    break
            if name != "exchange":
                raise ValueError("Previously frozen endpoint failed")
            steps[name] /= 2
        else:
            raise ValueError("No physical exchange endpoints")
        for sign, x, recipe in zip(SIGNS, endpoints, prepared, strict=True):
            if recipe["max_error"] > 1e-11:
                raise ValueError("Endpoint has no accurate native realization")
            key = f"{name}_{sign}"
            recipes[key], states[key] = recipe, x
    return states, steps, recipes


def regional_difference(board, route, odd):
    """Return (E_right - E_left)/2 and audit each full-cell balance in integer units."""
    bridge = measure.bridge
    cells = bridge.operators(board, route, camera.Q, odd)
    generator = bridge.native.generator(board, route, camera.Q, odd)
    regions4 = []
    slots = 0
    for port in route[1:3]:
        owner = port // 12
        edges = [
            (a, b) if b // 12 == owner else (b, a)
            for a, b in board.external
            if owner in (a // 12, b // 12)
        ]
        slots += (
            12 + sum(a // 12 == b // 12 == owner for a, b in board.seams) + len(edges)
        )
        region4 = 2 * bridge.masked(cells, [owner]) + bridge.bridges(board, edges)
        incoming = sum(bridge.current(board, cells, edge) for edge in edges)
        if (generator @ region4 - region4 @ generator - incoming).count_nonzero():
            raise ValueError("Whole-cell continuity equation failed")
        regions4.append(region4)
    if slots > 96:
        raise ValueError("Regional constituents exceed registered meter slots")
    return (regions4[1] - regions4[0]) / 8


def energy_effect(board, route, rings, time, strength, operators):
    src = camera.BEAM.indices(rings[0] + rings[1])
    effects = []
    for odd, operator in enumerate(operators):
        if strength == 8:
            h = measure.bridge.native.generator(board, route, camera.Q, odd) / 2
        else:
            h = camera.BEAM.front.optics.ranging.generator(
                board, route, camera.Q, strength, odd
            )
        launch = np.zeros((h.shape[0], len(src)))
        launch[src] = np.eye(len(src))
        ht = -1j * time * h
        flight = expm_multiply(ht, launch, traceA=ht.diagonal().sum())
        effects.append(
            np.mean(
                [v.conj().T @ (operator @ v) for v in np.hsplit(flight, 10)], axis=0
            )
        )
    return np.asarray(effects)


def acquire(states, steps, board, route, rings, operators, strength, seed):
    rng = np.random.default_rng(seed)
    panels = {key: {"forward": [], "reverse": [], "energy_plus": []} for key in states}
    clocks = []
    for offset in (-camera.OFFSET, 0, camera.OFFSET):
        t = CENTER + offset
        clock = expm(-1j * t * camera.blank_clock(1, 8))
        clocks.append(camera.precision_draw(rng, CLOCK_SHOTS, abs(clock[1, 0]) ** 2))
        for tag, screen in (("forward", rings), ("reverse", rings[::-1])):
            table = camera.BEAM.effects(board, route, screen, camera.Q, t, strength)
            for key, state in states.items():
                hits = camera.BEAM.probabilities(table, state)
                panels[key][tag].append(camera.precision_draw(rng, SHOTS, hits))
        if offset:
            effect = energy_effect(board, route, rings, t, strength, operators)
            for key, state in states.items():
                reading = measure.station.pair(effect, state)
                panels[key]["energy_plus"].append(
                    camera.precision_draw(rng, SHOTS, (1 + reading / SCALE) / 2)
                )
    return {
        "strength": strength,
        "seed": seed,
        "half_steps": steps,
        "clock_plus_counts": clocks,
        "shots_per_setting": SHOTS,
        "shots_per_clock": CLOCK_SHOTS,
        "stocks": panels,
    }


def decode(record, *, energy_shots=SHOTS):
    steps = record["half_steps"]
    panels = record["stocks"]
    if set(panels) != {f"{name}_{sign}" for name in STOCKS for sign in SIGNS}:
        raise ValueError("Missing or unexpected prepared stock")
    if record["shots_per_setting"] != SHOTS or record["shots_per_clock"] != CLOCK_SHOTS:
        raise ValueError("Unregistered measurement budget")
    counts = record["clock_plus_counts"]
    if len(counts) != 3 or not all(0 <= n <= CLOCK_SHOTS for n in counts):
        raise ValueError("Invalid clock record")
    clocks = camera.recorded_clocks(counts, CLOCK_SHOTS)
    # Two energy settings per prepared stock, each sign lies in {-1,1}.
    if energy_shots not in (SHOTS, 2**94):
        raise ValueError("Unregistered energy budget")
    error = SCALE * iv.sqrt(
        2 * iv.ln(iv.mpf(4 * len(panels)) / iv.mpf(".003")) / energy_shots
    )
    rates, focus = {}, {}
    for key, panel in panels.items():
        counts = panel["energy_plus"]
        if len(counts) != 2 or not all(0 <= n <= energy_shots for n in counts):
            raise ValueError("Incomplete regional energy record")
        delta = 2 * SCALE * iv.mpf(counts[1] - counts[0]) / energy_shots
        rates[key] = (delta + iv.mpf([-2, 2]) * error) / (clocks[2] - clocks[0])
        readings = []
        for tag in ("forward", "reverse"):
            if len(panel[tag]) != 3:
                raise ValueError("Missing camera appointment")
            images = [
                camera.pixel_image(rows, SHOTS, stock_count=len(panels))
                for rows in panel[tag]
            ]
            try:
                readings.append(camera.focusing(images, clocks))
            except AssertionError:
                return {
                    "decision": "unresolved clock or camera inversion",
                    "stock": key,
                    "tag": tag,
                }
        focus[key] = readings[0] - readings[1]

    def contrast(values, name):
        return (values[f"{name}_plus"] - values[f"{name}_minus"]) / (2 * steps[name])

    ordinary_f = contrast(focus, "ordinary")
    ordinary_rate = contrast(rates, "ordinary")
    out = {"ordinary_rate": str(ordinary_rate), "ordinary_focusing": str(ordinary_f)}
    if 0 in ordinary_rate:
        return {**out, "decision": "unidentified ordinary calibration"}
    coupling = ordinary_f / ordinary_rate
    out["ordinary_multiplier"] = str(coupling)
    out["comparisons"] = {}
    for name in STOCKS[1:]:
        rate = contrast(rates, name)
        forecast, observed = coupling * rate, contrast(focus, name)
        residual = observed - forecast
        tol = iv.mpf(".01") * abs(forecast)
        out["comparisons"][name] = {
            "regional_rate": str(rate),
            "predicted_focusing": str(forecast),
            "measured_focusing": str(observed),
            "measured_minus_predicted": str(residual),
            "equality_rejected": 0 not in residual,
            "one_percent_rejected": bool(residual.a > tol.b or residual.b < -tol.b),
            "one_percent_established": bool(abs(residual).b <= tol.a),
        }
    return out


def run():
    states, steps, recipes = preparations()
    formation = camera.BEAM.front.optics.pairing.formation
    board = formation.federation(1)
    route = formation.region(board)
    rings = camera.BEAM.screens(board, route)
    operators = [regional_difference(board, route, odd) for odd in (False, True)]
    banks = {}
    for strength, seed in ((8, 190821), (0, 190822)):
        bank = acquire(states, steps, board, route, rings, operators, strength, seed)
        banks[str(strength)] = {"counts": bank, "decoded": decode(bank)}
    return {
        "protocol_sha256": digest(HERE / "REGIONAL_RATE.md"),
        "matched_source_sha256": digest(matched_readings.RECEIPT),
        "endpoint_recipes": recipes,
        "banks": banks,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    if args.counts_only:
        data = json.loads(RECEIPT.read_text())
        if data["protocol_sha256"] != digest(HERE / "REGIONAL_RATE.md") or data[
            "matched_source_sha256"
        ] != digest(matched_readings.RECEIPT):
            raise ValueError("Protocol or source differs from the acquisition")
        if any(
            decode(bank["counts"]) != bank["decoded"] for bank in data["banks"].values()
        ):
            raise ValueError("Counted intervals differ from acquisition")
    else:
        data = run()
        if args.write:
            with RECEIPT.open("x") as handle:
                json.dump(data, handle, indent=2, allow_nan=False)
                handle.write("\n")
        elif data != json.loads(RECEIPT.read_text()):
            raise ValueError("Regional energy acquisition differs from archive")
    print(json.dumps({k: b["decoded"] for k, b in data["banks"].items()}, indent=2))


if __name__ == "__main__":
    main()
