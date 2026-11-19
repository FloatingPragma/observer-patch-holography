"""Command line of the simulator.

    python -m simulator laws
    python -m simulator check
    python -m simulator run continuum [--size 32] [--chart-time 3] [--law reflected-routing] [--out FILE]
    python -m simulator run journeys [--out FILE]
    python -m simulator app [RUN ...] [--out FILE]
"""
from __future__ import annotations

import argparse
import base64
import json
import sys
import time
from pathlib import Path

import numpy as np

from . import runs
from .experiments import clock_journey, continuum_run
from .law import LAWS

HERE = Path(__file__).resolve().parent
OUT = HERE / "runs"
STEP3_ROUTES = ((1, 0), (1, 1))


def _laws(_args) -> int:
    for law in LAWS.values():
        print(f"{law.key:28s} {law.status:12s} {law.title}")
    return 0


def _check(_args) -> int:
    from .test_simulator import CHECKS
    started = time.time()
    for name, check in CHECKS:
        print(f"{name}: {check()}")
    print(f"all checks passed in {time.time() - started:.1f} s")
    return 0


def _run(args) -> int:
    started = time.time()
    if args.experiment == "continuum":
        convergence = []
        for size in sorted({16, 24, args.size} - {0}):
            if size == args.size:
                continue
            r = continuum_run(size, args.chart_time, args.law, keep_frames=False)
            convergence.append({"size": size, "max_error": r["max_error"]})
        result = continuum_run(args.size, args.chart_time, args.law)
        convergence.append({"size": args.size, "max_error": result["max_error"]})
        convergence.sort(key=lambda row: row["size"])
        document = runs.continuum_document(result, convergence)
        default = OUT / f"continuum-{args.law}-N{args.size}.run.json.gz"
        for row in convergence:
            print(f"N = {row['size']:3d}: largest deviation from the continuum {row['max_error']:.6f}")
    else:
        journeys = [clock_journey(route) for route in STEP3_ROUTES]
        for j in journeys:
            print(f"route {j['route']}: reunion {j['reunion']} ticks, forecast {j['forecast']:.5f}, "
                  f"beat {j['beat_ratio']:.5f}, optical {j['optical_ratio']:.5f}")
        document = runs.journeys_document(journeys, time.time() - started)
        default = OUT / "journeys-reflected-routing.run.json.gz"
    path = runs.write(document, Path(args.out) if args.out else default)
    print(f"wrote {path} in {time.time() - started:.1f} s")
    return 0


def _app(args) -> int:
    template = (HERE / "app" / "index.html").read_text(encoding="utf-8")
    embedded = [base64.b64encode(Path(p).read_bytes()).decode() for p in args.runs]
    payload = json.dumps(embedded)
    marker = "/*EMBEDDED_RUNS*/[]"
    if marker not in template:
        raise SystemExit("the app template lacks its embedding marker")
    out = Path(args.out) if args.out else OUT / "replay.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(template.replace(marker, payload), encoding="utf-8")
    print(f"wrote {out} with {len(embedded)} embedded run(s)")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m simulator", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("laws", help="list the laws of the family and their standing").set_defaults(func=_laws)
    sub.add_parser("check", help="exactness and regression checks (under a minute)").set_defaults(func=_check)
    run = sub.add_parser("run", help="run an experiment and write a run file")
    run.add_argument("experiment", choices=("continuum", "journeys"))
    run.add_argument("--size", type=int, default=32)
    run.add_argument("--chart-time", type=float, default=3.0)
    run.add_argument("--law", choices=sorted(LAWS), default="reflected-routing")
    run.add_argument("--out")
    run.set_defaults(func=_run)
    app = sub.add_parser("app", help="write the replay app, optionally with embedded runs")
    app.add_argument("runs", nargs="*")
    app.add_argument("--out")
    app.set_defaults(func=_app)
    args = parser.parse_args(argv)
    np.seterr(over="raise", invalid="raise", divide="raise")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
