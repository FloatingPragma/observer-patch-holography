"""Independent sharp local returns, with no intermediate station measurement."""

import argparse
import json
import numpy as np
from mpmath import iv

from .. import measure as formed
from . import chart
from model.experiments.interval_echo_einstein_station.collar_pairing.fringe_ruler.certificate import (
    contraction,
)
from matched_einstein_limit.source_implication.directional_stress.traction.certificate import (
    four_tick_columns,
)
from matched_einstein_limit.source_implication.full_geometry.station_balls.formed_optics import (
    certificate as optics,
)

SHOTS = 2**42
GRID = 2**48
SEED = 851203


def probabilities():
    iv.dps = 42
    states = formed.preparations()
    coins = {s: optics.quantize(x) for s, x in states.items()}
    n = formed.native
    board = n.iso.pairing.formation.federation(1)
    marks = n.iso.pairing.formation.region(board)
    ports = [*marks, 5]
    indices = np.array([4 * p + j for p in ports for j in range(4)])
    den = optics.DEN * n.prior.FRAC**2
    co = n.prior.weights()[7]
    tables = {}
    for charge, name in ((-1, "vacancy"), (1, "excess")):
        panels = {}
        for odd in (False, True):
            h = n.generator(board, marks, charge, odd)
            waves = four_tick_columns(h, indices, indices, co).reshape(2, 5, 4, 5, 4)
            for stock, c in coins.items():
                if odd and stock != "independent":
                    continue
                ps = []
                for i in range(5):
                    u = waves[:, i, :, i, :]
                    numer, imaginary = contraction(u, u, c)
                    assert imaginary == 0
                    prob = iv.mpf(numer) / den + iv.mpf(["-3e-18", "3e-18"])
                    if odd:
                        prob = (prob + iv.mpf(panels["independent"][i])) / 2
                    ps.append(chart.bounds(prob))
                panels["ordinary" if odd else stock] = ps
        tables[name] = panels
    return tables


def acquire(probs):
    rng = np.random.default_rng(SEED)
    ordinary_rng = np.random.default_rng(SEED + 1)
    out = {}
    for q, stocks in probs.items():
        out[q] = {}
        for name, rows in stocks.items():
            hits = []
            for lower, upper in rows:
                p = round(GRID * (lower + upper) / 2) / GRID
                assert max(abs(p - lower), abs(p - upper)) < 2 / GRID
                sampler = ordinary_rng if name == "ordinary" else rng
                hits.append(int(sampler.binomial(SHOTS, p)))
            out[q][name] = hits
    return dict(
        starts_each=SHOTS,
        seed=SEED,
        ordinary_seed=SEED + 1,
        ports=[0, 1, 38, 39, 5],
        hits=out,
    )


def report(geo, probs, raw):
    assert raw["starts_each"] == SHOTS
    iv.dps = 42
    # Thirty two-sided binomial frequency checks, including ticket rounding.
    width = iv.sqrt(iv.ln(6000) / (2 * SHOTS)) + iv.mpf(2) / GRID
    weights = [iv.mpf(w) for w in geo["pooled"]["weights"]]
    moment = iv.mpf(geo["pooled"]["second_moment"])

    def delta(data):
        def change(stock):
            step = [
                f - b for f, b in zip(data[stock], data["independent"], strict=True)
            ]
            return sum(w * x for w, x in zip(weights, step[1:], strict=True)) - step[0]

        val = change("formed")
        ordinary = 2 * change("ordinary")
        assert ordinary.a > 0
        return dict(
            second_difference=chart.bounds(val),
            per_squared_radar_unit=chart.bounds(2 * val / moment),
            ordinary_per_energy_unit=chart.bounds(ordinary),
            equivalent_exchange_energy=chart.bounds(val / ordinary),
        )

    results = {}
    for q, panels in probs.items():
        measured = {}
        native = {}
        for s, ps in panels.items():
            native[s] = [iv.mpf(r) for r in ps]
            measured[s] = []
            for p, h in zip(native[s], raw["hits"][q][s], strict=True):
                center = iv.mpf(h) / SHOTS
                interval = center + iv.mpf([-1, 1]) * width
                assert p.a >= interval.a and p.b <= interval.b
                measured[s].append(interval)
        results[q] = {"native": delta(native), "counts": delta(measured)}
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    geo = chart.analyze()
    probs = probabilities()
    path = chart.HERE / "results.json"
    if args.write:
        raw = acquire(probs)
    else:
        saved = json.loads(path.read_text())
        raw = saved["raw"]
    result = dict(
        geometry=geo, probabilities=probs, raw=raw, readout=report(geo, probs, raw)
    )
    if args.write:
        with path.open("x") as out:
            out.write(json.dumps(result, indent=2) + "\n")
    else:
        for q, model in result["readout"].items():
            for record, values in model.items():
                for k, endpoints in values.items():
                    np.testing.assert_allclose(
                        endpoints, saved["readout"][q][record][k], atol=1e-13, rtol=0
                    )
        print(
            "Verified five-station geometry, both formed-source returns and all thirty count banks."
        )


if __name__ == "__main__":
    main()
