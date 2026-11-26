"""Calibrated exterior charge and signed spatial moments after opening the board."""

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np

from model.executable_model import federation
from model.experiments.count_intervals import clopper_pearson

from . import reader, words

HERE = Path(__file__).resolve().parent
PORTS = HERE.parent / "records.json"
SCATTER = HERE.parent / "scattering" / "results.json"
SOURCE_GATEWAYS = (1, 3, 10)
STARTS = 2**36
GRID = 2**40
SEED = 140526
WAITS = (F(1, 8192), F(1, 4096), F(1, 2048), F(1, 1024))
STRENGTHS = (F(1), F(1, 2), F(0))


def audit_board(level):
    board = federation(level)
    edge = [(a, b) for a, b in board.external if a < 12]
    if tuple(a for a, _ in edge) != SOURCE_GATEWAYS or len({b for _, b in edge}) != 3:
        raise ValueError("Expected three distinct exits from the sender")
    if set(board.degree.tolist()) != {5, 6}:
        raise ValueError("Width bound requires port degrees five or six")
    endpoints = [p for seam in board.external for p in seam]
    if len(endpoints) != len(set(endpoints)):
        raise ValueError(
            "The locality bound requires at most one external seam per port"
        )
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    distance = nx.single_source_shortest_path_length(owners, 0)
    if len(distance) != board.cells:
        raise ValueError("Owner chart incomplete")
    return dict(owners=board.cells, exits=edge, maximum_distance=max(distance.values()))


def sharp_source(record):
    panel = record["panels"]["neutral"]["native"]["profile"]
    out = []
    for gate in SOURCE_GATEWAYS:
        p = panel[f"gateway_{gate}"]
        counts = [
            sum(c for key, c in zip(p["keys"], p["counts"], strict=True) if key[0] == q)
            for q in (-1, 0, 1)
        ]
        n = sum(counts)
        if n != p["accepted"] or n + p["rejected"] != record["starts"] or n < 4_000_000:
            raise ValueError("Incomplete or insufficient source calibration")
        mean = F(counts[2] - counts[0], n)
        out.append((counts, (mean - F(1, 500), mean + F(1, 500))))
    return out


def exact_bounds(t, exchange, diameter, intervals):
    """Rational enclosures on the finite signed means divided by t squared."""
    # Relative contrast lies in [16,36]*8, and exchange is in [0,24*exchange].
    width = 80 + 12 * exchange
    cross = F(866026, 10**6) * exchange  # greater than sqrt(3)*exchange/2
    remainder = cross * width * t  # norm error after division by t
    square = exchange**2 / 4

    def enclose(source, leading_norm, meter_norm):
        err = 2 * leading_norm * remainder + meter_norm * remainder**2
        return [square * source[0] - err, square * source[1] + err]

    source = [sum(row[j] for row in intervals) for j in range(2)]
    coarse = enclose(source, cross, 2 * diameter**2)
    exterior = enclose(source, cross, 2)
    # Once a load has crossed, reaching owner shell j >= 2 costs at least
    # 2j-2 more contacts. The farther shells have summable reciprocal-factorial tails.
    z = width * t
    if z >= 1:
        raise ValueError("This positive-time bound uses z < 1")
    far = 6 * (cross * z**2 / ((1 - z) * 6)) ** 2 / (1 - z**4)
    radius = [max(coarse[0], exterior[0] - far), min(coarse[1], exterior[1] + far)]
    return {
        "three_first_ports": [enclose(row, exchange / 2, 1) for row in intervals],
        "net_exterior_charge": exterior,
        "diameter_only_squared_radius": coarse,
        "signed_squared_radius": radius,
        "slope_in_squared_wave_number": [-radius[1] / 2, -radius[0] / 2],
    }


def display(x):
    """Round lower/upper rational endpoints outwards, to six decimal places."""
    if isinstance(x, dict):
        return {k: display(v) for k, v in x.items()}
    if isinstance(x[0], list):
        return [display(row) for row in x]
    return [math.floor(x[0] * 10**6) / 10**6, math.ceil(x[1] * 10**6) / 10**6]


def scattering_source(record):
    # First and second moments from two separate detector settings at each gate.
    std = record["standards"]["native"]["2"]["probability"]
    real_gain = std[0] - std[1]
    gain_squared = 1 - 2 * std[2]
    p = record["panels"]["neutral"]
    accepted = p["acceptance"]
    rows = []
    for gate in SOURCE_GATEWAYS:
        panel = p["settings"][f"gateway_{gate}/center"]
        plus, minus, _ = np.reshape(
            panel["native_fringe"]["probability"][:-1], (2, 3)
        ).sum(axis=0)
        signed = (plus - minus) / (accepted * real_gain)
        occupied = sum(panel["native_elastic"]["probability"][:-1:2]) / (
            accepted * gain_squared
        )
        rows.append([(occupied - signed) / 2, 1 - occupied, (occupied + signed) / 2])
    return accepted, np.array(rows)


def transfer_probability(ratio, exchange):
    eighth_turn = 4 * ratio * exchange
    if eighth_turn.denominator != 1:
        raise ValueError("Unregistered clock angle")
    low = (2 - math.sqrt(2)) / 4
    return [0, low, 0.5, 1 - low, 1, 1 - low, 0.5, low][int(eighth_turn) % 8]


def draw(rng, source, accepted, transmission):
    vacancy, _, excess = source
    p = np.array(
        [
            accepted * vacancy * transmission,
            accepted * (1 - (vacancy + excess) * transmission),
            accepted * excess * transmission,
            1 - accepted,
        ]
    )
    if np.any(p < 0) or abs(p.sum() - 1) > 2e-13:
        raise ValueError("Invalid terminal probabilities")
    weights = np.rint(GRID * p).astype(np.int64)
    weights[np.argmax(p)] += GRID - weights.sum()
    counts = rng.multinomial(STARTS, weights / GRID)
    return dict(probabilities=p.tolist(), counts=counts.tolist(), grid=weights.tolist())


def run():
    places = [audit_board(i) for i in range(3)]
    sharp = sharp_source(json.loads(PORTS.read_text()))
    accepted, incident = scattering_source(json.loads(SCATTER.read_text()))
    for (counts, _), row in zip(sharp, incident, strict=True):
        if np.max(np.abs(np.array(counts) / sum(counts) - row)) > 0.002:
            raise ArithmeticError("Independent sharp and scattering sources disagree")
    bounds = {}
    for chart in places:
        arcs = {}
        for x in STRENGTHS:
            arcs[str(x)] = {
                str(t): display(
                    exact_bounds(
                        t, x, chart["maximum_distance"], [row[1] for row in sharp]
                    )
                )
                for t in WAITS
            }
        bounds[str(chart["owners"])] = arcs
    # A calibrated binary dimer angle and a first receiving port on separate starts.
    rng = np.random.default_rng(SEED)
    pulse = {}
    for x in STRENGTHS:
        for ratio in (F(0), F(1, 2), F(1), F(2)):
            transmission = transfer_probability(ratio, x)
            pulse[f"exchange={x}, duration_over_pi={ratio}"] = {
                key: draw(rng, values, admission, transmission)
                for key, values, admission in [
                    ("vacancy_standard", [1, 0, 0], 1),
                    ("excess_standard", [0, 0, 1], 1),
                    (
                        "uniformly_chosen_first_receiver",
                        incident.mean(axis=0),
                        accepted,
                    ),
                ]
            }
    source_estimate = [sum(float(v[1][j]) for v in sharp) / 3 for j in range(2)]
    forecasts = {}
    for setting, bank in pulse.items():
        negative = bank["vacancy_standard"]["counts"][0] / STARTS
        positive = bank["excess_standard"]["counts"][2] / STARTS
        gain = (positive + negative) / 2
        g_lo, g_hi = max(0, gain - 0.00002), min(1, gain + 0.00002)
        expected = [source_estimate[0] * g_lo, source_estimate[1] * g_hi]
        vacancy, unit, excess, _ = bank["uniformly_chosen_first_receiver"]["counts"]
        admitted = vacancy + unit + excess
        if admitted < 39_000_000:
            raise ValueError("Too few primary receipts for the fixed error budget")
        measured = (excess - vacancy) / admitted
        interval = [measured - 0.001, measured + 0.001]
        if expected[0] > interval[1] or expected[1] < interval[0]:
            raise ArithmeticError("Calibrated signed transfer missed the receiver")
        forecasts[setting] = {
            "from_clock_and_sharp_source": expected,
            "receiving_score_interval": interval,
        }
    clocks = {}
    for t in WAITS:
        # The local controller restricts this standard to its first half-period.
        # Inferring an unwrapped duration from one terminal bit would be ambiguous.
        probability = math.sin(float(t) / 2) ** 2
        bank = draw(rng, [0, 0, 1], 1, probability)
        n = bank["counts"][2]
        lo, hi = (float(v) for v in clopper_pearson(n, STARTS, 0.01 / 32))
        a = 2 * math.asin(math.sqrt(max(0.0, lo - 3 / GRID)))
        b = 2 * math.asin(math.sqrt(min(1.0, hi + 3 / GRID)))
        lower = F(math.floor(10**10 * a), 10**10)
        upper = F(math.ceil(10**10 * b), 10**10)
        clocks[str(t)] = dict(
            standard=bank,
            duration_interval=[float(lower), float(upper)],
            full_native_radius_per_squared_time=display(
                exact_bounds(
                    upper,
                    F(1),
                    places[0]["maximum_distance"],
                    [row[1] for row in sharp],
                )["signed_squared_radius"]
            ),
        )
    return dict(
        source_sha256={
            p.name if p == PORTS else "scattering/results.json": hashlib.sha256(
                p.read_bytes()
            ).hexdigest()
            for p in (PORTS, SCATTER)
        },
        terminal_reader=reader.census(),
        planted_full_flight=words.standards(),
        geometry=places,
        all_sharp_source_counts=[row[0] for row in sharp],
        source_charge_interval=display([list(row[1]) for row in sharp]),
        full_flight_bounds_per_squared_exchange_time=bounds,
        local_clock_checks=clocks,
        bridge_only=dict(
            starts=STARTS,
            grid=GRID,
            seed=SEED,
            source_acceptance=accepted,
            incoming_vacancy_unit_excess=incident.tolist(),
            outcomes=["vacancy", "unit", "excess", "source refused"],
            forecasts=forecasts,
            banks=pulse,
        ),
    )


def verify(got, archived):
    if isinstance(archived, dict):
        if got.keys() != archived.keys():
            raise ValueError("Receipt fields differ")
        for key in archived:
            verify(got[key], archived[key])
    elif isinstance(archived, list):
        for actual, expected in zip(got, archived, strict=True):
            verify(actual, expected)
    elif isinstance(archived, float):
        if not math.isclose(got, archived, rel_tol=5e-12, abs_tol=1e-15):
            raise ValueError("Numerical result differs")
    elif got != archived:
        raise ValueError("Detector count or setting differs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    target = parser.parse_args().output
    result = run()
    if target:
        with target.open("x") as stream:
            json.dump(result, stream, indent=1, allow_nan=False)
            stream.write("\n")
    else:
        verify(result, json.loads((HERE / "results.json").read_text()))
        print("All free-flight bounds and pulse/clock banks verified")
