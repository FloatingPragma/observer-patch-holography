"""Reproduce all ordered two-offer ports and the correlated four-offer bank."""

import argparse
import gzip
import hashlib
import json
from fractions import Fraction
from itertools import product
from math import pi, sqrt
from pathlib import Path

import numpy as np

from . import apparatus as app
from .chart import scalar_chart_report

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(x) for x in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, float):
        val = round(value, 11)
        return val if val else 0.0
    return value


def forecast(word, initial):
    out = initial.copy()
    for e in word:
        a, b = app.SEAMS[e]
        out[[a, b]] = (out[a] + out[b]) / 2
    return out


def series():
    a = 0
    b = next(
        k for k, e in enumerate(app.SEAMS) if k != a and set(e) & set(app.SEAMS[a])
    )
    c = next(
        k
        for k, e in enumerate(app.SEAMS)
        if set(e).isdisjoint(app.SEAMS[a] + app.SEAMS[b])
    )
    return {
        (a, a, b, b): Fraction(1, 6),
        (b, b, a, a): Fraction(1, 3),
        (a, b, a, b): Fraction(1, 10),
        (a, c, b, a): Fraction(2, 5),
    }


def two_steps():
    starts = np.eye(app.N)
    rows = []
    for w in product(range(len(app.SEAMS)), repeat=2):
        bins = app.observed_word(w, starts)
        raw = np.stack([bins[x] for x in product((0, 1), repeat=2)])
        records = np.rint(4 * raw).astype(int)
        np.testing.assert_allclose(records / 4, raw, atol=app.TOL)
        np.testing.assert_allclose(
            records.sum(axis=0) / 4, forecast(w, starts), atol=app.TOL
        )
        rows.append(records.tolist())
    return rows


def trials():
    weights = series()
    exact, calendar, nodes = app.schedule(weights)
    if any(p is None for p in calendar.values()):
        raise ValueError(
            "The registered observation needs a resolved dimer on every positive branch"
        )
    source = np.eye(app.N)[:, app.SEAMS[0][0]]
    coherent = (1j ** np.arange(app.N)) / sqrt(app.N)
    banks = {}
    for name, state in (("localized", source), ("quarter_phases", coherent)):
        initial = abs(state) ** 2
        position_law = []
        for i in range(4):
            prob = {}
            for w, m in weights.items():
                prob[w[i]] = prob.get(w[i], 0) + m
            position_law.append(prob)
        independent = initial.copy()
        for p in position_law:
            independent = sum(
                float(m) * forecast([edge], independent) for edge, m in p.items()
            )
        modes = {}
        for mode, t, read in (
            ("native", pi / 2, True),
            ("third_turn", pi / 3, True),
            ("unread", pi / 2, False),
        ):
            mixture = np.zeros(app.N)
            details = []
            for word, branch in calendar.items():
                rec = app.observed_word(word, state[:, None], t, read)
                terminal = sum(rec.values())[:, 0]
                if mode == "native":
                    np.testing.assert_allclose(
                        terminal, forecast(word, initial), atol=app.TOL
                    )
                np.testing.assert_allclose(terminal.sum(), 1, atol=app.TOL)
                mixture += branch * terminal
                details.append(
                    {
                        "word": word,
                        "paired_read_transcript_and_port": [
                            [list(key), p[:, 0]] for key, p in sorted(rec.items())
                        ],
                    }
                )
            modes[mode] = {
                "conditional_records": details,
                "all_launch_terminal": mixture,
            }
        banks[name] = {
            "modes": modes,
            "independent_positions": independent,
            "joint_minus_independent_TV": abs(
                modes["native"]["all_launch_terminal"] - independent
            ).sum()
            / 2,
        }
    return {
        "weights": [[w, str(p)] for w, p in weights.items()],
        "calendar_nodes": nodes,
        "exact_emitted_word_law": [[w, str(p)] for w, p in exact.items()],
        "returned_word_law": [[w, p] for w, p in calendar.items()],
        "calendar_controls": [
            [p, app.calendar_bit(p)] for p in (0.0, 1.0, sqrt(2) / 3)
        ],
        "banks": banks,
    }


def repeated_seam_response(records):
    """Infer the uniform repair mean from the observed repeated-offer terminal locations."""
    maps = [
        np.sum(records[e * len(app.SEAMS) + e], axis=0) / 4
        for e in range(len(app.SEAMS))
    ]
    integral = np.rint(2 * np.sum(maps, axis=0)).astype(int)
    graph = np.zeros((app.N, app.N), dtype=int)
    for a, b in app.SEAMS:
        graph[a, b] = graph[b, a] = 1
    laplacian = np.diag(graph.sum(axis=1)) - graph
    np.testing.assert_array_equal(60 * np.eye(app.N, dtype=int) - integral, laplacian)
    distinct = len({m.tobytes() for m in maps})
    return {
        "distinct_maps": distinct,
        "uniform_mean_times_60": integral,
        "eigenvalues": np.linalg.eigvalsh(integral / 60),
    }


def report():
    files = [
        HERE / "PROTOCOL.md",
        HERE / "apparatus.py",
        HERE / "chart.py",
        Path(__file__),
        Path(app.candidate.__file__),
        Path(app.LAW.__file__),
        ROOT / "model/executable_model/carrier.py",
        ROOT / "model/executable_model/carriers/FEDERATION_LEVEL0.npz",
    ]
    records = two_steps()
    return clean(
        {
            "schema": "oph.seam-word-meter.v1",
            "sources": {
                str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in files
            },
            "provenance": "Probabilities are deterministic evaluations of the native operators. All twelve ports are retained. No Monte Carlo counts.",
            "layout": app.layout(),
            "subsequent_scalar_chart_check": scalar_chart_report(),
            "local": app.local_effects(),
            "two_offers_axes": [
                "lexicographic seam pair",
                "first and second reads 00,01,10,11",
                "terminal port",
                "source port",
            ],
            "two_offers_unit": "Units of 1/4 per launch. The input total is one, so a load-two pointer outcome always has probability zero.",
            "two_offers": records,
            "repeated_seam_response": repeated_seam_response(records),
            "held": trials(),
        }
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    outfile = HERE / "records.json.gz"
    args = parser.parse_args()
    result = report()
    if args.write:
        raw = json.dumps(result, separators=(",", ":"), ensure_ascii=False).encode()
        outfile.write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))
    if result != json.loads(gzip.decompress(outfile.read_bytes())):
        raise ValueError("Stored measurements differ")
    print("900 seam pairs, 12 source ports per pair")
    print("phase holds", result["local"]["phase_inputs_initial_unread_read"])
    print(
        "joint vs independent TV",
        [
            (k, v["joint_minus_independent_TV"])
            for k, v in result["held"]["banks"].items()
        ],
    )
