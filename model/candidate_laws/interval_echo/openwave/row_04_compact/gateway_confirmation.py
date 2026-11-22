"""Compare a registered, externally anchored release to an independent reproduction."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from model.executable_model import federation
from .experiment import gateway_table
from .gateway_acquisition import prepare_layout, scored_appointments
from .gateway_bound import SCORE, certify
from .gateway_release import ARMS


HERE = Path(__file__).parent
PROTOCOL = HERE/"gateway_confirmation.json"
TASK = "openwave/xperiments/m10_oph/research/tasks/m10_4_task_details.md"
TASK_SHA256 = "f7dca7249ef997d399af6264f611035b2342939e4cec053a37a375b6ac6625b4"
PREREG_SHA256 = "2558d1c0f77bfdbd9dacf6d48cc92eac2317fdf720cae573359e3ff171d358be"


def verify_anchor(checkout, anchor):
    """Read the local copy of official origin/main. There is no network access."""
    selection(anchor)
    protocol = PROTOCOL.read_bytes().replace(b"\r\n", b"\n")
    if sha256(protocol).hexdigest() != PREREG_SHA256:
        raise ValueError("Held gates differ from the registered protocol")
    base = json.loads(protocol)["prior_platform_revision"]

    def git(*arguments):
        return subprocess.run(["git", "-C", str(checkout), *arguments],
                              check=True, capture_output=True).stdout

    lineage = git("rev-list", "--first-parent", "refs/remotes/origin/main").decode().splitlines()
    if base not in lineage or anchor not in lineage[:lineage.index(base)]:
        raise ValueError("Anchor absent from the official first-parent history after the baseline")
    changes = git("log", "--first-parent", "--reverse", "--format=%H", base+"..refs/remotes/origin/main",
                  "--", TASK).decode().splitlines()
    for commit in changes:
        if not git("ls-tree", commit, "--", TASK):
            continue
        task = git("show", commit+":"+TASK)
        if sha256(task).hexdigest() == TASK_SHA256:
            if PREREG_SHA256.encode() not in task or commit != anchor:
                raise ValueError("Anchor differs from the first registered version")
            return
    raise ValueError("The registered upstream task is absent")


def selection(anchor):
    if not isinstance(anchor, str) or not re.fullmatch(r"[0-9a-f]{40}", anchor):
        raise ValueError("Expected the lowercase 40-hex official first-parent registration anchor")
    digest = sha256(("OPH-M10.4-gateway-v1:" + anchor).encode("ascii")).digest()
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
    return dict(anchor=anchor, settings=settings, floor_ppt=protected["joint_home_floor_ppt"],
                bins={arm: [*row[arm]["ports"], row[arm]["exterior_quiet"], row[arm]["source_exhaustion"]]
                      for arm in ARMS})


def check(actual, independent):
    gates = json.loads(PROTOCOL.read_text(encoding="utf-8"))["gates"]
    if (set(independent) != {"anchor", "settings", "bins"} or independent["anchor"] != actual["anchor"]
            or set(independent["settings"]) != set(actual["settings"])
            or any(type(v) is not int or v != actual["settings"][k] for k, v in independent["settings"].items())
            or set(independent["bins"]) != set(ARMS)):
        raise ValueError("Independent calculation describes a different acquisition")
    home = {}
    for arm in ARMS:
        rows = []
        for payload in (actual["bins"][arm], independent["bins"][arm]):
            if type(payload) is not list or len(payload) != 14 or any(type(x) not in (int, float) for x in payload):
                raise ValueError("Expected a complete numeric terminal distribution")
            values = np.asarray(payload, dtype=float)
            if (not np.isfinite(values).all() or np.any(values < 0) or np.any(values > 1)
                    or abs(values.sum()-1) > gates["normalization_tolerance"]):
                raise ValueError("Invalid terminal law")
            rows.append(values)
        if np.max(abs(rows[0]-rows[1])) > gates["independent_bin_tolerance"]:
            raise ArithmeticError(f"Independent terminal law differs for {arm}")
        home[arm] = rows[0][:12].sum()
    if (actual["floor_ppt"] <= round(SCORE*gates["exact_joint_home_strictly_above"])
            or home["coherent"] <= gates["coherent_terminal_home_strictly_above"]
            or home["stochastic_swap"] >= gates["stochastic_swap_home_strictly_below"]
            or home["coherent"]-home["dephased"] < gates["minimum_coherent_minus_dephased"]):
        raise ArithmeticError("Held release fails its registered scientific gates")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("anchor", help="verified first official main commit containing the registered task")
    parser.add_argument("--openwave-checkout", type=Path, required=True,
                        help="trusted checkout with the fetched official refs/remotes/origin/main")
    parser.add_argument("--compare", type=Path, required=True, help="independent calculation of all four arms")
    args = parser.parse_args()
    verify_anchor(args.openwave_checkout, args.anchor)
    receipt = acquisition(args.anchor)
    check(receipt, json.loads(args.compare.read_text(encoding="utf-8")))
    print(json.dumps(receipt, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
