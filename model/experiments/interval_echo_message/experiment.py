"""Bounded exposure acquisition. Default verification is read-only."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as reference
from model.candidate_laws.interval_echo import candidate as echo
from model.executable_model import federation

from .acquisition import evaluate, train
from .channel import HERE, SPEC, TOL, distance, encode, readout, require, sources
from .clock_join import probe
from .decoder import decode
from .headline import check


def panel(station):
    level = station["level"]
    board = federation(level)
    h = echo.generator(board)
    exchange = reference.generator(board, 0)
    residual = h - exchange - 6 * echo.KAPPA * sparse.eye(board.n)
    scalar_error = float(abs(residual).max())
    require(scalar_error <= TOL, "kappa removal changes the one-load generator beyond a scalar")
    launch = encode(board, station)
    calibrated = []
    for grid, s in enumerate(SPEC["exposures"]):
        p, norm = readout(h, launch, s)
        calibration, lookup = train(p, level, grid)
        calibrated.append((p, norm, calibration, lookup))
    # All receiver settings on this board are frozen before any evaluation outcome.
    windows = []
    controls_ok = True
    for grid, (s, entry) in enumerate(zip(SPEC["exposures"], calibrated, strict=True)):
        p, norm, calibration, lookup = entry
        zero, zero_norm = readout(exchange, launch, s)
        slow, slow_norm = readout(h - .5 * exchange, launch, s)
        equiv = float(np.max(abs(zero - p)))
        require(equiv <= TOL, "non-equivalent terminal laws at zero contrast")
        delta = float(np.max(abs(slow - p)))
        tv = distance(p)
        sham = np.repeat(p[:1], 2, axis=0)
        native = evaluate(p, lookup, level, grid, 0)
        no_encoding = evaluate(sham, lookup, level, grid, 1)
        require(distance(sham) == 0, "sham channel")
        for control in (native["permuted"], no_encoding):
            lo, hi = control["interval"]
            controls_ok &= lo <= .5 <= hi
        if s == SPEC["control_exposure"]:
            controls_ok &= delta >= SPEC["minimum_law_shift"]
        windows.append(dict(exposure=s, probabilities=p.tolist(), norm_residual=norm,
                            detection=p[:, :-1].sum(axis=1).tolist(), tv=tv, optimal=(1 + tv) / 2,
                            calibration_counts=calibration, lookup=list(lookup), native=native,
                            no_encoding=no_encoding, zero_kappa_probabilities=zero.tolist(),
                            zero_kappa_max_delta=equiv, zero_kappa_norm_residual=zero_norm,
                            half_exchange_probabilities=slow.tolist(), half_exchange_norm_residual=slow_norm,
                            half_exchange_max_delta=delta))
    planted = np.zeros_like(launch)
    planted[[0, 1], [0, 1]] = 1
    plant, _ = readout(h, planted, 0)
    # Known local port-to-bit mapping checks terminal sensitivity independently of training.
    plant_table = bytes([0, 1] + [0] * 11)
    plant_success = sum(plant[decode(y, plant_table), y] for y in range(13)) / 2
    controls_ok &= plant_success == 1 and distance(plant) == 1
    indexed = {r["exposure"]: r["native"]["interval"] for r in windows}
    resolved = indexed[SPEC["early"]][1] < SPEC["accuracy_target"] < indexed[SPEC["late"]][0]
    verdict = "CHANNEL_RESOLVED" if resolved else "BELOW_REGISTERED_SENSITIVITY"
    return dict(cells=board.cells, windows=windows, scalar_error=scalar_error,
                planted_probabilities=plant.tolist(), planted_success=float(plant_success),
                verdict=verdict if controls_ok else "NUMERICAL_OR_INSTRUMENT_FAILURE")


def acquire():
    sources()
    result = dict(protocol_sha256=hashlib.sha256((HERE / "protocol.json").read_bytes()).hexdigest(),
                  stations=[panel(station) for station in SPEC["stations"]])
    check(SPEC, result)
    return result


def compare(saved, actual):
    if type(saved) is not type(actual):
        raise ValueError("changed record type")
    if isinstance(saved, dict):
        require(saved.keys() == actual.keys(), "changed record fields")
        for key in saved:
            compare(saved[key], actual[key])
    elif isinstance(saved, list):
        for a, b in zip(saved, actual, strict=True):
            compare(a, b)
    elif type(saved) is float:
        require(abs(saved - actual) <= 1e-12, "changed numerical record")
    else:
        require(saved == actual, "changed record")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="create evidence at an unused path")
    args = parser.parse_args()
    result = acquire()
    result["clock_attempt"] = probe()
    if args.record:
        payload = json.dumps(result, allow_nan=False, separators=(",", ":"))
        # One public window per line keeps the complete vectors compact and inspectable.
        with args.record.open("x", encoding="utf8") as handle:
            handle.write(payload.replace('{"exposure":', '\n{"exposure":') + "\n")
    else:
        compare(json.loads((HERE / "receipt.json").read_text()), result)
    for site in result["stations"]:
        print(site["cells"], site["verdict"])
    print(result["clock_attempt"]["verdict"])
    if any(site["verdict"] != "CHANNEL_RESOLVED" for site in result["stations"]):
        raise SystemExit(1)
