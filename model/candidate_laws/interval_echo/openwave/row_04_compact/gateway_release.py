"""Replay the bounded gateway source and the complete terminal records."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from model.executable_model import federation
from .experiment import gateway_table
from .gateway_acquisition import EXPOSURES, LAUNCHES, PULSES, gateway_acquisition, prepare_layout
from .gateway_bound import ATTEMPTS, ROUNDS, SCORE, board_charts, certify


HERE = Path(__file__).parent
ARMS = ("coherent", "dephased", "unfiltered_uniform", "stochastic_swap")


def experiment():
    charts = board_charts()
    source = gateway_acquisition()
    if source["rounds"] != ROUNDS or tuple(source["launches"]) != tuple(LAUNCHES) or len(LAUNCHES) != ATTEMPTS:
        raise ArithmeticError("Certificate and acquisition describe different preparation schedules")
    layouts = []
    for row in source["layouts"]:
        tips = tuple(row["gateways"])
        counts = [chart.get(tips, 0) for chart in charts.values()]
        settings = {}
        for pulse in PULSES:
            setting = row["settings"][str(pulse)]
            output = {"exhaustion": round(setting["exhaustion"], 10)}
            if pulse:
                bound = certify(tips, Fraction(pulse))
                for (low, high), measured in zip(bound["one_attempt_success_ppt"],
                                                 setting["success_by_sharp_launch"], strict=True):
                    if not low/SCORE-3e-12 <= measured <= high/SCORE+3e-12:
                        raise ArithmeticError("Sharp source acceptance disagrees with the integer enclosure")
                if (not .94 < bound["joint_home_floor_ppt"]/SCORE <= setting["all_time_joint_home_floor"]+3e-12
                        or setting["exhaustion"] > bound["exhaustion_ceiling_ppt"]/SCORE+3e-12):
                    raise ArithmeticError("Full entrance and exit law disagrees with the all-time certificate")
                output["joint_home_floor_ppt"] = bound["joint_home_floor_ppt"]
                output["exhaustion_ceiling_ppt"] = bound["exhaustion_ceiling_ppt"]
            settings[str(pulse)] = output
        layouts.append({"gateways": list(tips), "cells_per_board": counts, "pulses": settings,
                        "max_retry_contacts": row["settings"]["1.0"]["max_reset_contacts_over_47_retries"]})
    if sorted(tuple(row["gateways"]) for row in layouts) != sorted(set().union(*charts.values())):
        raise ArithmeticError("Incomplete gateway survey")
    readings = {}
    for pulse in PULSES:
        table = source["scored_appointments"][str(pulse)]
        if [row["t"] for row in table["readings"]] != list(EXPOSURES):
            raise ArithmeticError("Missing registered free-flight appointment")
        readings[str(pulse)] = {
            arm: [[round(float(v), 10) for v in
                   (*row[arm]["ports"], row[arm]["exterior_quiet"], row[arm]["source_exhaustion"])]
                  for row in table["readings"]] for arm in ARMS}
    return {
        "source_pin": "4fe10f9ac3e81423155fb3ac7baaf0fdb4c0b11c",
        "rounds": ROUNDS, "attempts": ATTEMPTS, "boards": list(charts), "layouts": layouts,
        "max_pulses": source["maximum_preparation_pulses"], "max_gateway_reads": source["maximum_gateway_reads"],
        "home": 13, "exposures": list(EXPOSURES),
        "terminal_order": "twelve local ports, accepted exterior quiet, source exhaustion",
        "terminal_by_pulse": readings,
    }


def verify(actual, stored):
    if type(actual) is not type(stored):
        raise ValueError("Changed gateway receipt field type")
    if isinstance(actual, dict):
        if actual.keys() != stored.keys():
            raise ValueError("Changed gateway receipt fields")
        for key, value in actual.items():
            verify(value, stored[key])
    elif isinstance(actual, list):
        for value, reference in zip(actual, stored, strict=True):
            verify(value, reference)
    elif isinstance(actual, float):
        # The portable archive rounds probability bins to ten decimals.
        if not (0 <= stored <= 1 and abs(actual-stored) <= 1.01e-10):
            raise ValueError("Changed gateway outcome probability")
    elif actual != stored:
        raise ValueError("Changed gateway protocol or integer certificate")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-clicks", type=int, nargs=3, metavar=("LEVEL", "CELL", "LAUNCH"),
                        help="print every first-click branch of the appointed sharp launch")
    parser.add_argument("--pulse", type=float, choices=PULSES)
    args = parser.parse_args()
    if args.first_clicks is None:
        if args.pulse is not None:
            parser.error("--pulse applies to --first-clicks")
        verify(experiment(), json.loads((HERE/"gateway_results.json").read_text(encoding="utf-8")))
        print("All gateway layouts, source budgets, continuous floors and terminal bins agree.")
        return
    level, home, port = args.first_clicks
    if level not in (0, 1, 2):
        parser.error("invalid board")
    charts = gateway_table(federation(level))
    if not 0 <= home < len(charts) or not 0 <= port < 12:
        parser.error("invalid home or launch")
    prep = prepare_layout(charts[home], args.pulse if args.pulse is not None else 1.0)
    print(json.dumps({"gateways": prep.gateways, "first_click": prep.first_click[:, :, port].tolist(),
                      "all_clear": float(prep.success_by_launch[port])}, allow_nan=False))


if __name__ == "__main__":
    main()
