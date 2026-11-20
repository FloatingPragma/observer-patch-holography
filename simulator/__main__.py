"""Command line of the simulator.

    python -m simulator laws
    python -m simulator check
    python -m simulator run continuum [--size 32] [--chart-time 3] [--law oph] [--out FILE]
    python -m simulator run journeys [--out FILE]
    python -m simulator app [RUN ...] [--out DIR]
"""
from __future__ import annotations

import argparse
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


def _run_all(args) -> int:
    """Every successful experiment of steps 1 to 5 as a run file, and the app with all of them."""
    from .experiments import (dispersion_scan, journey_run, pair_correction_scan, race_run)
    from .law import LAWS
    paths = []

    def save(document, name, zero_line=False):
        if zero_line:
            document["zero_line"] = True
        path = runs.write(document, OUT / name)
        paths.append(path)
        print(f"wrote {path.name}")

    started = time.time()
    convergence = [{"size": n, "max_error": continuum_run(n, 3.0, keep_frames=False)["max_error"]} for n in (16, 24)]
    field = continuum_run(32, 3.0)
    convergence.append({"size": 32, "max_error": field["max_error"]})
    document = runs.continuum_document(field, convergence)
    document["step"] = "Spacetime"
    document["title"] = "The 3D scalar continuum, seen from an observer"
    save(document, "01-continuum.run.json.gz")

    scan = dispersion_scan()
    pts = scan["points"]
    curve = [[x, float(np.sqrt(1 - x * x))] for x in np.linspace(0, 0.995, 120)]
    save(runs.scan_document(
        "A clock carried by the travelling wave", "Clocks",
        "The beat of two private modes on a freely moving packet, divided by the same beat at rest, "
        "for launch momenta from long to short wavelengths. Speeds come from the packet's own envelope "
        "and the archived free-wave speed. At long wavelengths the clock follows special relativity; at "
        "short wavelengths the finite law makes it run slower.",
        "speed / free-wave speed", "moving / rest clock rate",
        [{"name": "clock on the travelling packet", "points": [[p["v_over_c"], p["ratio"]] for p in pts]}],
        [{"name": "special relativity, sqrt(1 - v^2/c^2)", "points": curve}],
        LAWS["oph"], scan["seconds"],
        [["launch order", "v / c", "forecast", "clock"]] +
        [[p["order"], round(p["v_over_c"], 5), round(p["forecast"], 5), round(p["ratio"], 5)] for p in pts]),
        "02-moving-clock-scan.run.json.gz")

    race = race_run()
    save(runs.race_document(race), "03-signal-race.run.json.gz")

    two = journey_run(("beat", "optical"), "oph", (1, 1))
    save(runs.journey_document(
        two, "Two clocks on one journey", "Clocks",
        "A private-mode beat clock and an optical clock, a pulse circling the third spatial axis, travel "
        "the forward diagonal and return. Their stationary twins stay at home. One calibration predicts "
        "both readings at reunion."), "04-two-clocks.run.json.gz")

    pair = journey_run(("rails",), "oph", (1, 1))
    save(runs.journey_document(
        pair, "A clock held together by an interaction", "Clocks",
        "Two constituents bound by the contact law travel as one pair. The relative phase of one "
        "constituent's two rails is the clock hand, and the rails differ only in their contact coupling."),
        "05-bound-pair-clock.run.json.gz")

    correction = pair_correction_scan()
    pts = correction["points"]
    save(runs.scan_document(
        "The finite-wavelength correction of the bound clock", "Clocks",
        "How far three forecasts of the bound pair's moving/rest rate lie from the rate the observer "
        "recorded, on the diagonal at three central wavelengths: the smooth metric with the "
        "measured speed, the finite law at the central wavelength, and the finite law over the whole "
        "packet read at the home address. The recorded rates come from a signed-return detector; "
        "the remaining 0.1 to 0.4 thousandths are the difference between that detector and the "
        "untagged home reading.",
        "central launch order", "forecast − recorded rate (×1000)",
        [{"name": "smooth metric", "points": [[p["order"], 1000 * (p["smooth"] - p["recorded"])] for p in pts]},
         {"name": "finite law, central wavelength", "points": [[p["order"], 1000 * (p["central"] - p["recorded"])] for p in pts]},
         {"name": "finite law, whole packet at home", "points": [[p["order"], 1000 * (p["packet"] - p["recorded"])] for p in pts]}],
        [], LAWS["oph"], correction["seconds"],
        [["order", "smooth", "central mode", "whole packet", "recorded"]] +
        [[p["order"], round(p["smooth"], 6), round(p["central"], 6), round(p["packet"], 6), p["recorded"]] for p in pts]),
        "06-bound-clock-correction.run.json.gz", zero_line=True)

    from .formation import formation_run
    formed = formation_run(size=96, updates=200, orders=(14, -10), separation=20, law_key="oph")
    control = formation_run(size=96, updates=200, orders=(14, -10), separation=20, law_key="oph", angle=0.0)
    save(runs.formation_document(formed, control), "07-formation.run.json.gz")

    from .row30 import matter_frames, row30_run
    frames, stride = matter_frames()
    save(runs.row30_document(row30_run(), frames, stride), "08-openwave-row-30.run.json.gz")

    from .row04 import row04_run
    save(runs.row04_document(row04_run()), "09-openwave-row-04.run.json.gz")

    args.runs = [str(p) for p in paths]
    args.out = None
    _app(args)
    print(f"all runs in {time.time() - started:.1f} s")
    return 0


def _run(args) -> int:
    if args.experiment == "all":
        return _run_all(args)
    if args.experiment == "row04":
        from .row04 import row04_run
        path = runs.write(runs.row04_document(row04_run()),
                          Path(args.out) if args.out else OUT / "09-openwave-row-04.run.json.gz")
        print(f"wrote {path}")
        return 0
    if args.experiment == "row30":
        from .row30 import matter_frames, row30_run
        frames, stride = matter_frames()
        path = runs.write(runs.row30_document(row30_run(), frames, stride),
                          Path(args.out) if args.out else OUT / "08-openwave-row-30.run.json.gz")
        print(f"wrote {path}")
        return 0
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
        default = OUT / "journeys.run.json.gz"
    path = runs.write(document, Path(args.out) if args.out else default)
    print(f"wrote {path} in {time.time() - started:.1f} s")
    return 0


def _app(args) -> int:
    from . import openwave_sync, webapp
    if openwave_sync.DEFAULT_REPORTS.is_dir():
        reported = openwave_sync.sync()
        waiting = [n for n, s in reported.items() if s != "ported"]
        print(f"OpenWave rows: ported {sorted(openwave_sync.PORTED)}; reported, not ported {waiting}")
    out = Path(args.out) if args.out else webapp.DEFAULT_OUT
    pages = webapp.build([Path(p) for p in args.runs], out)
    print(f"wrote the app to {out}: {len(pages)} pages, {len(args.runs)} embedded run(s); open {out / 'index.html'}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m simulator", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("laws", help="list the laws of the family and their standing").set_defaults(func=_laws)
    sub.add_parser("check", help="exactness and regression checks (under a minute)").set_defaults(func=_check)
    run = sub.add_parser("run", help="run an experiment and write a run file")
    run.add_argument("experiment", choices=("continuum", "journeys", "row30", "row04", "all"))
    run.add_argument("--size", type=int, default=32)
    run.add_argument("--chart-time", type=float, default=3.0)
    run.add_argument("--law", choices=sorted(LAWS), default="oph")
    run.add_argument("--out")
    run.set_defaults(func=_run)
    app = sub.add_parser("app", help="write the replay app, one page per run and per OpenWave row")
    app.add_argument("runs", nargs="*")
    app.add_argument("--out", help="folder of the site (default simulator/runs/site)")
    app.set_defaults(func=_app)
    args = parser.parse_args(argv)
    np.seterr(over="raise", invalid="raise", divide="raise")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
