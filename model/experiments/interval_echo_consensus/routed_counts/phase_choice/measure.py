"""Let a repair-count comparison select a physical analyzer, then read four pointers."""

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ...signed_effects.pair_return import measure as pointer
from ...signed_effects.sequential import grid_counts
from .. import material

HERE = Path(__file__).resolve().parent
ROOT = pointer.address.ROOT
ARCHIVE = ROOT / "code/thermodynamics/repair_current_orientation/repair_current_payload.v3.json"
SHOTS = 8192
ANGLES = tuple(map(Fraction, ("0", "1/2", "-1/2", "1", "1/6", "-1/6", "1/3", "-1/3")))
WINDOW = 20
EPISODES = 96


def orientation(table):
    """Lexicographically least absolute maximizer, with an explicit zero-gap refusal."""
    if len(table) != 8 or any(len(row) != 8 for row in table):
        raise ValueError("an eight-bucket table is required")
    if any(type(x) is not int or x < 0 for row in table for x in row):
        raise ValueError("counts must be nonnegative integers")
    # Repeated indices have gap zero. Each distinct triple has six ordered
    # forms of the same absolute gap; the increasing order is the least.
    best = 0
    triple = None
    for a, b, c in combinations(range(8), 3):
        gap = table[a][b] * table[b][c] * table[c][a] - table[a][c] * table[c][b] * table[b][a]
        if abs(gap) > abs(best):
            best, triple = gap, [a, b, c]
    return {"triple": triple, "gap": best, "sign": (best > 0) - (best < 0)}


def counts(words, *, reverse=False, across_resets=True):
    windows = []
    for word in words:
        if len(word) != 3 * WINDOW or set(word) - {"0", "1"}:
            raise ValueError("incomplete repair log")
        windows.append([min(7, word[j:j + WINDOW].count("1")) for j in range(0, len(word), WINDOW)])
    streams = [[v for row in windows for v in row]] if across_resets else windows
    table = [[0] * 8 for _ in range(8)]
    for stream in streams:
        if reverse:
            stream = stream[::-1]
        for a, b in zip(stream, stream[1:]):
            table[a][b] += 1
    return table


def episode(sites, pairs, board, routing, chords, rng):
    stocks = rng.permutation([0] * 9 + [2] * 9).tolist()
    menu = rng.integers(len(pairs), size=3 * WINDOW).tolist()
    result = {"preparation": "".join(map(str, stocks)), "offer_indices_hex": bytes(menu).hex()}
    for name, disabled in (("native", False), ("disabled", True)):
        loads = [1] * board.n
        for site, stock in zip(sites, stocks, strict=True):
            loads[site] = stock
        used = 0
        reports = []
        for k in menu:
            marked, _ = material.attempt(loads, routing, pairs[k], chords[used], disabled)
            used += int(marked)
            reports.append("1" if marked else "0")
        result[name] = "".join(reports)
    return result


def histories():
    board = federation(0)
    routing = material.Routing(board)
    sites = np.random.default_rng(129001).choice(np.arange(24, 216), size=18, replace=False).tolist()
    pairs = list(combinations(sites, 2))
    chords = [c for c in routing.chords if not set(c).intersection(sites)]
    reports = []
    for seed in (129002, 129003):
        rng = np.random.default_rng(seed)
        rows = [episode(sites, pairs, board, routing, chords, rng) for _ in range(EPISODES)]
        raw = [row["native"] for row in rows]
        fwd, rev = counts(raw), counts(raw, reverse=True)
        tables = {"forward": fwd, "backward": rev,
                  "pooled": [[x + y for x, y in zip(a, b, strict=True)] for a, b in zip(fwd, rev, strict=True)],
                  "disconnected": counts([row["disabled"] for row in rows]),
                  "within_episode": counts(raw, across_resets=False)}
        reports.append({"seed": seed, "episodes": rows, "tables": tables,
                        "selectors": {key: orientation(table) for key, table in tables.items()}})
    return sites, reports


def analyzer_bank(sign, angle, head, rng):
    """Compile a distinct contact for each sign before reading the same blank pointers."""
    native = pointer.pulses
    source = native.pulse("contact", 0, 1, angle % 2)
    dwell = Fraction(3 if sign > 0 else 1, 2)
    analyzer = native.pulse("contact", 0, 1, dwell)
    data = native.source_unitary((source, analyzer), 1)[:, 0]
    weights = (abs(head @ data) ** 2).sum(axis=1)
    ideal = (1 + sign * np.sin(np.pi * float(angle))) / 2
    if abs(weights.sum() - 1) > 3e-12 or abs(weights[0] - ideal) > 3e-12:
        raise ArithmeticError("native pointer differs from independent phase forecast")
    raw = grid_counts(rng, SHOTS, weights / weights.sum())
    return {"theta_over_pi": str(angle), "dwell_over_pi": str(dwell), "raw": raw.tolist(),
            "forecast": round(float(ideal), 14)}


def report():
    checks = material.certify()
    sites, runs = histories()
    archive = orientation(json.loads(ARCHIVE.read_bytes())["ordered_counts"])
    head = pointer.read_branches(pointer.address.schedule(), 0)
    rng = np.random.default_rng(129004)
    banks = {str(s): [analyzer_bank(s, theta, head, rng) for theta in ANGLES] for s in (1, -1)}
    radius = float(np.sqrt(np.log(2 * len(ANGLES) * 2 / 0.01) / (2 * SHOTS))) + 3e-12
    for bank in banks.values():
        for r in bank:
            if abs(r["raw"][0] / SHOTS - r["forecast"]) > radius:
                raise ArithmeticError("held contact forecast left the registered interval")
    for run in runs:
        decisions = run["selectors"]
        if decisions["backward"]["sign"] != -decisions["forward"]["sign"]:
            raise ArithmeticError("re-reading the history did not reverse its strict orientation")
        if any(decisions[x]["sign"] for x in ("pooled", "disconnected")):
            raise ArithmeticError("negative control issued a phase choice")
        # A dispatcher emits only a pulse command or refusal. Scoring angles
        # and terminal counts are unavailable to it; both pairings share banks.
        run["dispatch"] = {}
        for k, choice in decisions.items():
            s = choice["sign"]
            run["dispatch"][k] = {"declared_bank": str(s) if s else None,
                                  "opposite_bank": str(-s) if s else None}
    sources = [HERE / "PROTOCOL.md", Path(__file__), Path(material.__file__),
               Path(pointer.__file__), Path(pointer.pulses.__file__), Path(pointer.address.__file__),
               ROOT / "model/candidate_laws/interval_echo/candidate.py", ARCHIVE]
    return {"native_checks": checks, "sites": sites, "archive_selection": archive, "histories": runs,
            "shots_per_setting": SHOTS, "phase_seed": 129004, "simultaneous_99_radius": radius,
            "banks_by_analyzer_sign": banks,
            "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record = report()
    dest = HERE / "results.json"
    if args.write:
        dest.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    pointer.compare(json.loads(dest.read_bytes()), record)
    for run in record["histories"]:
        print(run["seed"], run["selectors"])
