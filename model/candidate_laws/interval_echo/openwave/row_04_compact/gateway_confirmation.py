"""Compare a registered, externally anchored release to an independent reproduction."""
import argparse
from base64 import b64encode
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from model.executable_model import federation
from tools.strict_json import loads
from .experiment import gateway_table
from .gateway_acquisition import prepare_layout, scored_appointments
from .gateway_bound import SCORE, certify
from .gateway_release import ARMS


HERE = Path(__file__).parent
PROTOCOL = HERE/"gateway_confirmation.json"
TASK = "openwave/xperiments/m10_oph/research/tasks/m10_4_task_details.md"
TASK_SHA256 = "af312744a6842aa093c8fdfc3961053bb6a5dd29f407f4b1ba11bff742cf3640"
PREREG_SHA256 = "2558d1c0f77bfdbd9dacf6d48cc92eac2317fdf720cae573359e3ff171d358be"


def registration_anchor(checkout):
    """Find the first exact task in a trusted, fetched official first-parent history."""
    protocol = PROTOCOL.read_bytes().replace(b"\r\n", b"\n")
    if sha256(protocol).hexdigest() != PREREG_SHA256:
        raise ValueError("Held gates differ from the registered protocol")
    base = json.loads(protocol)["prior_platform_revision"]

    def git(*arguments):
        return subprocess.run(["git", "-C", str(checkout), *arguments],
                              check=True, capture_output=True).stdout

    lineage = git("rev-list", "--first-parent", "refs/remotes/origin/main").decode().splitlines()
    if base not in lineage:
        raise ValueError("Official first-parent history lacks the registration baseline")
    if git("ls-tree", base, "--", TASK) and sha256(git("show", base+":"+TASK)).hexdigest() == TASK_SHA256:
        raise ValueError("Registration already present at the baseline")
    changes = git("log", "--first-parent", "--reverse", "--format=%H", base+"..refs/remotes/origin/main",
                  "--", TASK).decode().splitlines()
    for commit in changes:
        if not git("ls-tree", commit, "--", TASK):
            continue
        task = git("show", commit+":"+TASK)
        if sha256(task).hexdigest() == TASK_SHA256:
            if PREREG_SHA256.encode() not in task:
                raise ValueError("Task lacks the registered protocol digest")
            return commit
    return None


def verify_anchor(checkout, anchor):
    if type(anchor) is not str or not re.fullmatch(r"[0-9a-f]{40}", anchor):
        raise ValueError("Expected a lowercase 40-hex official commit")
    registered = registration_anchor(checkout)
    if registered is None:
        raise ValueError("The registered upstream task is absent")
    if registered != anchor:
        raise ValueError("Anchor differs from the first registered version")


def selector_digest(anchor):
    if not isinstance(anchor, str) or not re.fullmatch(r"[0-9a-f]{40}", anchor):
        raise ValueError("Expected the lowercase 40-hex official first-parent registration anchor")
    return sha256(("OPH-M10.4-gateway-v1:" + anchor).encode("ascii")).digest()


def selection(anchor):
    digest = selector_digest(anchor)
    return dict(home=digest[0] % 20, start=digest[1] % 12,
                direction=1 if digest[2] % 2 == 0 else -1, exposure=15+digest[3] % 24)


def acquisition(anchor):
    settings = selection(anchor)
    home = settings["home"]
    order = {k: settings[k] for k in ("start", "direction")}
    chart = tuple(gateway_table(federation(0))[home])
    protected = certify(chart, **order)
    preparation = prepare_layout(chart, **order)
    row, = scored_appointments(preparation, home=home, exposures=(settings["exposure"],))["readings"]
    return dict(anchor=anchor, selector_sha256=selector_digest(anchor).hex(),
                settings=settings, floor_ppt=protected["joint_home_floor_ppt"],
                bins={arm: [*row[arm]["ports"], row[arm]["exterior_quiet"], row[arm]["source_exhaustion"]]
                      for arm in ARMS})


def read_independent(contents):
    value = loads(contents.decode("utf-8"))
    json.dumps(value, allow_nan=False)
    return value


def terminal_law(bins, tolerance, *, native):
    if type(bins) is not dict or set(bins) != set(ARMS):
        raise ValueError("Missing a terminal arm")
    rows = {}
    for arm, payload in bins.items():
        if type(payload) is not list or len(payload) != 14 or any(
                type(x) not in (int, float) or not 0 <= x <= 1 for x in payload):
            raise ValueError("Expected a complete numeric terminal distribution")
        values = np.asarray(payload, dtype=float)
        if abs(values.sum()-1) > tolerance:
            if native:
                raise ArithmeticError("Native terminal law fails the normalization gate")
            raise ValueError("Invalid independent terminal law")
        rows[arm] = values
    return rows


def native_law(actual, gates):
    if type(actual["floor_ppt"]) is not int or not 0 <= actual["floor_ppt"] <= SCORE:
        raise ValueError("Invalid integer enclosure for the native release")
    return terminal_law(actual["bins"], gates["normalization_tolerance"], native=True)


def check(actual, independent, comparison=None):
    gates = json.loads(PROTOCOL.read_text(encoding="utf-8"))["gates"]
    native = native_law(actual, gates)
    if (type(independent) is not dict or set(independent) != {"anchor", "settings", "bins"}
            or independent["anchor"] != actual["anchor"] or type(independent["settings"]) is not dict
            or set(independent["settings"]) != set(actual["settings"])
            or any(type(v) is not int or v != actual["settings"][k] for k, v in independent["settings"].items())):
        raise ValueError("Independent calculation describes a different acquisition")
    peer = terminal_law(independent["bins"], gates["normalization_tolerance"], native=False)
    home = {}
    worst = {"absolute": -1., "arm": None, "outcome": None}
    for arm in ARMS:
        differences = abs(native[arm]-peer[arm])
        idx = int(np.argmax(differences))
        if differences[idx] > worst["absolute"]:
            worst = {"absolute": float(differences[idx]), "arm": arm, "outcome": idx}
        home[arm] = native[arm][:12].sum()
    if comparison is not None:
        comparison.update(worst)
    if worst["absolute"] > gates["independent_bin_tolerance"]:
        raise ArithmeticError(f"Independent terminal law differs for {worst['arm']}")
    if (actual["floor_ppt"] <= round(SCORE*gates["exact_joint_home_strictly_above"])
            or home["coherent"] <= gates["coherent_terminal_home_strictly_above"]
            or home["stochastic_swap"] >= gates["stochastic_swap_home_strictly_below"]
            or home["coherent"]-home["dephased"] < gates["minimum_coherent_minus_dephased"]):
        raise ArithmeticError("Held release fails its registered scientific gates")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("anchor", nargs="?", help="first official main commit containing the registered task")
    parser.add_argument("--openwave-checkout", type=Path, required=True,
                        help="trusted checkout with the fetched official refs/remotes/origin/main")
    parser.add_argument("--registration-status", action="store_true", help="inspect registration without selecting a flight")
    parser.add_argument("--compare", type=Path, help="independent calculation of all four arms")
    parser.add_argument("--record", type=Path, help="exclusive output file for the result")
    args = parser.parse_args()
    if args.registration_status:
        if args.anchor or args.compare or args.record:
            parser.error("--registration-status takes only --openwave-checkout")
        first = registration_anchor(args.openwave_checkout)
        print(json.dumps({"registered": first is not None, "anchor": first}))
        return
    if not (args.anchor and args.compare and args.record):
        parser.error("a held comparison requires ANCHOR, --compare and --record")
    verify_anchor(args.openwave_checkout, args.anchor)
    with args.record.open("x", encoding="utf-8") as destination:
        result = {"anchor": args.anchor, "selector_sha256": selector_digest(args.anchor).hex(),
                  "settings": selection(args.anchor), "comparison": {}}
        stage = "acquisition"
        try:
            native = acquisition(args.anchor)
            # A malformed nonfinite native output cannot be represented as JSON numbers.
            # Retain its diagnostic rendering as a string in a valid incomplete receipt.
            try:
                json.dumps(native, allow_nan=False)
            except ValueError:
                result["invalid_acquisition_text"] = json.dumps(native)
                raise
            result["acquisition"] = native
            stage = "native gates"
            native_law(native, json.loads(PROTOCOL.read_text(encoding="utf-8"))["gates"])
            stage = "independent receipt"
            contents = args.compare.read_bytes()
            result["independent_sha256"] = sha256(contents).hexdigest()
            result["independent_base64"] = b64encode(contents).decode("ascii")
            independent = read_independent(contents)
            stage = "comparison"
            check(result["acquisition"], independent, result["comparison"])
        except (ValueError, KeyError, TypeError, ArithmeticError, OSError, RecursionError) as failure:
            rejected = stage in ("native gates", "comparison") and isinstance(failure, ArithmeticError)
            result.update(passed=False if rejected else None, status="rejected" if rejected else "incomplete",
                          stage=stage, reason=str(failure))
            raise
        else:
            result.update(passed=True, status="passed")
        finally:
            json.dump(result, destination, indent=2, allow_nan=False)
            destination.write("\n")
    print("Held release passed with all independent bins verified.")


if __name__ == "__main__":
    main()
