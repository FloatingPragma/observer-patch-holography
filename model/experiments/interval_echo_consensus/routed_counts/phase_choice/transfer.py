"""Transmit a count decision to owner zero and drive a fresh phase measurement there."""

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_gauge.exterior_register import material as rail

from . import measure as trial

WORDS = {-1: 0, 0: 1, 1: 31}
COMMANDS = {v: k for k, v in WORDS.items()}
ANGLES = tuple(map(Fraction, ("1/6", "-1/6", "1/3", "-1/3")))
RECEIPT = trial.HERE / "transfer.json"


def source_word(code, station):
    geometry, base, _, _ = station
    sent = base.copy()
    for j, (a, b) in enumerate(geometry["dimers"]["bob"]):
        if code & (1 << j):
            sent[a], sent[b] = sent[b], sent[a]
    return sent


def stock_program(code, station, routing, chords):
    """Compile from the received public word and known stock, before actual execution."""
    planned = source_word(code, station)[station[2]].copy()
    target = trial.pointer.address.initial_word(1, 0)
    writes = []
    while (planned == 0).sum() > (target == 0).sum():
        a = int(np.flatnonzero(planned == 0)[0])
        b = int(np.flatnonzero(planned == 2)[0])
        chord = chords[len(writes)]
        writes.append({"source": (a, b), "chord": chord, "routing_pulses": 2 * len(routing.legs((a,b), chord))})
        planned[a] = planned[b] = 1
    if sorted(planned) != sorted(target):
        raise ValueError("incompatible launch inventory")
    swaps = []
    for i in range(len(planned)):
        if planned[i] == target[i]:
            continue
        j = next(j for j in range(i + 1, len(planned))
                 if planned[j] == target[i] and planned[j] != target[j])
        planned[i], planned[j] = planned[j], planned[i]
        swaps.append((i, j, len(routing.transpose(i, j))))
    return {"writes": writes, "endpoint_swaps": swaps}


def prepare_received(loads, program, routing):
    for write in program["writes"]:
        marked, pulses = trial.material.attempt(loads, routing, write["source"], write["chord"])
        if not marked or pulses != write["routing_pulses"]:
            raise ValueError("stock-conversion request failed")
    for i, j, _ in program["endpoint_swaps"]:
        for a, b in routing.transpose(i, j):
            loads[a], loads[b] = loads[b], loads[a]
    # Verification of the propagated word. No settings are selected here.
    if not np.array_equal(loads, trial.pointer.address.initial_word(1, 0)):
        raise ValueError("routed stock differs from actual phase launch")


def deliver(selected, cut, station, lookup):
    geometry, _, connected, isolated = station
    sent = source_word(WORDS[selected], station)
    tx = rail.read(sent, geometry["read_paths"]["bob"], lookup)
    received = sent[isolated if cut else connected].copy()
    rx = rail.read(received, geometry["read_paths"]["alice"], lookup)
    return tx, rx, received


def run():
    primary_path = trial.HERE / "results.json"
    primary = json.loads(primary_path.read_bytes())
    tables = [(f'{r["seed"]}:{name}', t) for r in primary["histories"]
              for name, t in r["tables"].items()]
    archive = json.loads(trial.ARCHIVE.read_bytes())["ordered_counts"]
    tables += [("archive:forward", archive), ("archive:backward", list(map(list, zip(*archive))))]
    station = rail.layout()
    lookup = rail.shared_report.certify()
    routing = trial.material.Routing(federation(0))
    calls = trial.pointer.address.schedule()
    active = {tuple(sorted((a, b))) for kind, a, b, _ in trial.pointer.address.pulse_word(calls)
              if kind == "contact"}
    active.add((0, 1))
    chords = [c for c in routing.chords if c not in active]
    head = trial.pointer.read_branches(calls, 0)
    rng = np.random.default_rng(129005)
    radius = float(np.sqrt(np.log(2 * 48 / 0.01) / (2 * trial.SHOTS))) + 3e-12
    receipts = []
    for name, t in tables:
        selection = trial.orientation(t)
        for cut in (False, True):
            tx, rx, received = deliver(selection["sign"], cut, station, lookup)
            code = rail.decode(rx)
            record = {"source": name, "selection": selection, "bridge_cut": cut,
                      "sent_clicks": tx, "received_clicks": rx, "received_word": code}
            sign = COMMANDS.get(code)
            if sign in (None, 0):
                record["no_read_per_angle"] = [trial.SHOTS] * len(ANGLES)
            else:
                record["stock_conversion"] = stock_program(code, station, routing, chords)
                prepare_received(received, record["stock_conversion"], routing)
                record["banks"] = [trial.analyzer_bank(sign, a, head, rng) for a in ANGLES]
                for bank in record["banks"]:
                    if abs(bank["raw"][0] / trial.SHOTS - bank["forecast"]) > radius:
                        raise ArithmeticError("fresh phase readings failed their prediction")
            receipts.append(record)
    return {"owners": station[0]["owners"], "message_contacts": station[0]["contact_count"],
            "message_bridge_contacts": station[0]["bridge_contacts"], "seed": 129005,
            "angles_over_pi": list(map(str, ANGLES)), "shots_per_angle": trial.SHOTS,
            "simultaneous_99_radius": radius, "records": receipts,
            "source_sha256": {str(p.relative_to(trial.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [primary_path, Path(__file__), Path(rail.__file__),
                                        Path(rail.shared_report.__file__), trial.HERE / "TRANSFER_PROTOCOL.md"]}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    trial.pointer.compare(json.loads(RECEIPT.read_bytes()), json.loads(json.dumps(result)))
    print("Verified delivered commands, stock conversion and all terminal readings.")
