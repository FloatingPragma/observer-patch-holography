"""Joint native cell-charge records, supplied Ramsey readout and matched marginals."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import CELL_SEAMS
from .. import measure as flight


HERE = Path(__file__).parent
ARCHIVE = HERE / "results.json"
TIME, EXTRA = 1_000_000, 4
PHASES = (np.pi/2, np.pi/3)
TOLERANCE = 2e-10


def pointer_decoder():
    """Complete native signed instruments and every allowed local word."""
    law = candidate.native_law
    echo, blank = law.echo(), np.eye(law.DIM, law.N_SRC)
    worst = 0.0
    for d in range(-2, 3):
        flag = np.repeat([False, *(v == d for _, v in law.LABELS)], law.N_SRC)
        actual = echo.T @ (flag[:, None] * (echo @ blank))
        expected = blank @ np.diag([a-b == d for a, b in law.SRC])
        worst = max(worst, float(np.max(abs(actual-expected))))
    if worst > 5e-14:
        raise ArithmeticError("Native signed pointer differs from its complete load effect")
    words = {(1,)*12}
    for excess in range(12):
        word = [1]*12
        word[excess] = 2
        words.add(tuple(word))
        word[excess] = 0
        words.add(tuple(word))
        for vacancy in range(12):
            if vacancy != excess:
                word = [1]*12
                word[excess], word[vacancy] = 2, 0
                words.add(tuple(word))
    decoder = {}
    for word in sorted(words):
        record = tuple(word[a]-word[b] for a, b in CELL_SEAMS)
        charge = sum(word)-12
        if record in decoder and decoder[record] != charge:
            raise ArithmeticError("The local signed records leave cell charge ambiguous")
        decoder[record] = charge
    if len(words) != 157 or len(decoder) != 157 or set(decoder.values()) != {-1, 0, 1}:
        raise ArithmeticError("The complete local charge domain differs from its declaration")
    encoded = json.dumps([[list(record), value] for record, value in sorted(decoder.items())]).encode()
    return {"local_words": len(words), "distinct_complete_records": len(decoder),
            "decoder_sha256": sha256(encoded).hexdigest(), "native_pointer_max_error": worst,
            "charge_outcomes": [-1, 0, 1],
            "domain": "One excess and one vacancy on the unit-filled board; no other local load patterns admitted."}


def hard_core_control(weights, x, y, n):
    """A positive maximum-entropy off-diagonal mixture with matching marginals."""
    row = np.bincount(x, weights=weights, minlength=n)
    column = np.bincount(y, weights=weights, minlength=n)
    mixture = row[:, None] * column[None, :]
    np.fill_diagonal(mixture, 0)
    for iteration in range(1, 2001):
        sums = mixture.sum(axis=1)
        mixture *= np.divide(row, sums, out=np.zeros_like(row), where=sums > 0)[:, None]
        sums = mixture.sum(axis=0)
        mixture *= np.divide(column, sums, out=np.zeros_like(column), where=sums > 0)[None, :]
        error = max(float(np.max(abs(mixture.sum(axis=1)-row))),
                    float(np.max(abs(mixture.sum(axis=0)-column))))
        if error < 2e-14:
            break
    else:
        raise ArithmeticError("Hard-core matched-marginal control did not converge")
    if mixture.min() < 0 or np.any(mixture.diagonal()) or abs(mixture.sum()-1) > 2e-12:
        raise ArithmeticError("Matched-marginal model control is not a valid hard-core mixture")
    return mixture[x, y], {"iterations": iteration, "maximum_marginal_residual": error,
                           "marginal_L1_residual": float(abs(mixture.sum(axis=1)-row).sum()+abs(mixture.sum(axis=0)-column).sum()),
                           "normalization_error": float(abs(mixture.sum()-1)),
                           "same_port_mass": 0,
                           "scope": "Prepared model control; no native preparation or reachability claim."}


def cell_records(weights, x, y, admission):
    records = []
    for owner in range(20):
        charge = (x//12 == owner).astype(int) - (y//12 == owner).astype(int)
        probabilities = np.bincount(charge+1, weights=weights, minlength=3)
        if probabilities.min() < 0 or abs(probabilities.sum()-1) > 2e-12:
            raise ArithmeticError("Cell-charge record loses an outcome")
        fields = []
        for phase in PHASES:
            cosine, sine = np.cos(phase), np.sin(phase)
            charged = probabilities[0]+probabilities[2]
            # Positive effect sums avoid cancellation at nearly dark ports.
            # Retain the numerical gain of the native input distribution.
            xminus = (1-cosine)*charged/2
            xplus = probabilities[1]+(1+cosine)*charged/2
            yplus = (1-sine)*probabilities[0]/2+probabilities[1]/2+(1+sine)*probabilities[2]/2
            yminus = (1+sine)*probabilities[0]/2+probabilities[1]/2+(1-sine)*probabilities[2]/2
            fields.append({"phase": phase, "X_plus_minus": [float(xplus), float(xminus)],
                           "Y_plus_minus": [float(yplus), float(yminus)],
                           "original_start_X_plus_minus_filter_miss": [float(admission*xplus), float(admission*xminus), float(1-admission)],
                           "original_start_Y_plus_minus_filter_miss": [float(admission*yplus), float(admission*yminus), float(1-admission)]})
        records.append({"owner": owner, "minus_zero_plus": probabilities.tolist(),
                        "original_start_minus_zero_plus_filter_miss": [*(admission*probabilities).tolist(), 1-admission],
                        "supplied_Ramsey": fields})
    different = float(weights[x//12 != y//12].sum())
    second_moment = sum(row["minus_zero_plus"][0]+row["minus_zero_plus"][2] for row in records)
    if abs(second_moment-2*different) > 2e-12:
        raise ArithmeticError("Native cell-charge second moment differs from the full pair census")
    return {"different_cell_probability": different, "summed_cell_charge_second_moment": second_moment,
            "cells": records}


def normalized_weights(state):
    raw = abs(state)**2
    normalization = float(raw.sum())
    if abs(normalization-1) > 5e-12:
        raise ArithmeticError("Native Ritz flight exceeds its normalization allowance")
    # This explicitly defines the normalized numerical Ritz state. Its change
    # from the raw law is retained and added to the physical-event allowance.
    return raw/normalization, abs(normalization-1)


def run():
    decoder = pointer_decoder()
    board = flight.federation(0)
    x, y, _, kinetic, h, cost = flight.defect.generators(board)
    original = json.loads((flight.HERE/"results.json").read_text())
    reduced = flight.PairReduction(h, cost)
    energies, basis, residual, _ = reduced.lowest(h, 24)
    if reduced.inertia(flight.CUT)["upper_lines_below"] != 24:
        raise ArithmeticError("The original moving pair's admitted band differs")
    projector_error = float(max(original["cluster_projector_error"], residual/(flight.CUT-energies[-1])))
    outputs = []
    for launch in original["flights"]:
        if launch["owner"] not in (0, 5, 10) or launch["admission"] < .001:
            continue
        index = flight.defect.ordered_index(*launch["ports"], board.n)
        overlap = basis[index]
        admission = float(overlap@overlap)
        coefficients = overlap/np.sqrt(admission)
        arrival = basis @ (np.exp(-1j*(energies-energies[0])*TIME)*coefficients)
        continued = basis @ (np.exp(-1j*(energies-energies[0])*(TIME+EXTRA))*coefficients)
        rival = expm_multiply((-1j*EXTRA)*kinetic, arrival)
        panels, norm_errors = {}, []
        for name, state in (("native_arrival", arrival), ("native_continuation", continued), ("exchange_only_continuation", rival)):
            weights, normalization_error = normalized_weights(state)
            panels[name] = cell_records(weights, x, y, admission)
            norm_errors.append(normalization_error)
            if name == "native_arrival":
                control, control_info = hard_core_control(weights, x, y, board.n)
                panels["matched_marginal_model"] = cell_records(control, x, y, admission)
        event_error = float(2*projector_error/np.sqrt(admission)+(TIME+EXTRA)*residual+1e-9+max(norm_errors))
        actual = panels["native_arrival"]["different_cell_probability"]
        changed = panels["exchange_only_continuation"]["different_cell_probability"]
        null = panels["matched_marginal_model"]["different_cell_probability"]
        outputs.append({"owner": launch["owner"], "initial_ports": launch["ports"],
                        "initial_filter_admission": admission, "filter_miss": 1-admission,
                        "native_event_probability_error": event_error,
                        "unconditional_event_error_upper": 2*projector_error+admission*event_error,
                        "explicit_numerical_normalization_L1_changes": norm_errors,
                        "matched_control": control_info, "panels": panels,
                        "gates": {"native_different_cells_below_001": bool(actual+event_error < .01),
                                  "exchange_only_different_cells_above_002": bool(changed-event_error > .02),
                                  "matched_model_effect_separation": null-actual-2*event_error,
                                  "signed_marginals_match_within_source_error": bool(control_info["maximum_marginal_residual"] < event_error)}})
    if len(outputs) != 3:
        raise ArithmeticError("The three enrolled moving-pair launches were not retained")
    return {"schema": "oph.interval_echo.cell_charge_probe.v1",
            "protocol_sha256": sha256((HERE/"PROTOCOL.md").read_bytes()).hexdigest(),
            "measure_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "source_receipt_sha256": sha256((flight.HERE/"results.json").read_bytes()).hexdigest(),
            "native_pointer_decoder": decoder, "native_flight": TIME, "extra_control_flight": EXTRA,
            "coupling_status": "The local signed pointers are native. The Ramsey ancilla phase theta*Q_C, drive, coupling strength and clock are supplied apparatus; electromagnetic identification is unestablished.",
            "scope": "Full final-address laws after original source admission; no arrival or same-owner conditioning. Ideal probabilities, no sampled trials or prediction credit.",
            "residual": residual, "projector_error": projector_error, "launches": outputs}


def compare(actual, archived):
    if type(actual) is not type(archived):
        raise ArithmeticError("Receipt type mismatch")
    if isinstance(actual, dict):
        if actual.keys() != archived.keys():
            raise ArithmeticError("Receipt fields differ")
        for key in actual:
            compare(actual[key], archived[key])
    elif isinstance(actual, list):
        if len(actual) != len(archived):
            raise ArithmeticError("Receipt outcome partition differs")
        for left, right in zip(actual, archived, strict=True):
            compare(left, right)
    elif isinstance(actual, float):
        if not np.isfinite(archived) or abs(actual-archived) > TOLERANCE:
            raise ArithmeticError("Receipt probabilities or bounds differ")
    elif actual != archived:
        raise ArithmeticError("Receipt metadata differs")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", type=Path, help="explicit exclusive first-record creation")
    args = cli.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.record:
        with args.record.open("x", encoding="utf-8") as destination:
            destination.write(payload)
    else:
        compare(result, json.loads(ARCHIVE.read_text(encoding="utf-8")))
    if any(not row["gates"]["native_different_cells_below_001"] or
           not row["gates"]["exchange_only_different_cells_above_002"] for row in result["launches"]):
        raise ArithmeticError("A recorded spatial-neutrality gate failed")
    for row in result["launches"]:
        print(row["owner"], {name: row["panels"][name]["different_cell_probability"] for name in row["panels"]})


if __name__ == "__main__":
    main()
