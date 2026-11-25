"""Replay one M10.5 appointment with the native and raw-edge implementations."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from model.executable_model import federation
from .. import gateway_independent as peer
from ..experiment import gateway_table
from ..gateway_acquisition import prepare_layout, scored_appointments
from ..gateway_bound import SCORE, certify
from ..gateway_release import ARMS

HERE = Path(__file__).resolve().parent
TERMS_HASH = "15a8374aa0b833fda546f4e5dfa67f5a1af5ddcfeba67a5b78035a6bd55c4c0d"
TASK_HASH = "fa7035ef0c707e2e23915e6350b62dc97eae907c3918d10831f55f8c22d9bb98"


def terms():
    raw = (HERE / "terms.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != TERMS_HASH:
        raise ValueError("Prospective terms differ from the candidate freeze")
    if hashlib.sha256((HERE / "TASK.md").read_bytes()).hexdigest() != TASK_HASH:
        raise ValueError("Task differs from the candidate freeze")
    return json.loads(raw)


def official_anchor(checkout, spec):
    """Use a trusted, freshly fetched OpenWave checkout and first-parent snapshots."""

    def git(*args):
        return subprocess.check_output(["git", "-C", str(checkout), *args])

    base, task = spec["platform_base"], spec["platform_task"]
    lineage = (
        git("rev-list", "--first-parent", "refs/remotes/origin/main")
        .decode()
        .splitlines()
    )
    if base not in lineage:
        raise ValueError("Fetch official main history through the frozen base")
    for commit in reversed(lineage[: lineage.index(base)]):
        entry = git("ls-tree", commit, "--", task)
        if (
            entry
            and hashlib.sha256(git("show", commit + ":" + task)).hexdigest()
            == TASK_HASH
        ):
            return commit
    return None


def selection(anchor, spec):
    digest = hashlib.sha256((spec["version"] + ":" + anchor).encode("ascii")).digest()
    return dict(
        home=digest[0] % 20,
        start=digest[1] % 12,
        direction=(-1 if digest[2] % 2 else 1),
        exposure=15 + digest[3] % 24,
    )


def native(settings):
    home = settings["home"]
    order = {k: settings[k] for k in ("start", "direction")}
    chart = tuple(gateway_table(federation(0))[home])
    prep = prepare_layout(chart, **order)
    bound = certify(chart, **order)["joint_home_floor_ppt"]
    (row,) = scored_appointments(prep, home=home, exposures=[settings["exposure"]])[
        "readings"
    ]
    bins = {
        arm: [
            *row[arm]["ports"],
            row[arm]["exterior_quiet"],
            row[arm]["source_exhaustion"],
        ]
        for arm in ARMS
    }
    return bound, bins


def compare(bins, independent, floor, gates):
    if set(bins) != set(ARMS) or set(independent) != set(ARMS):
        raise ValueError("Missing terminal arm")
    discrepancies = {}
    home = {}
    for arm in ARMS:
        raw, other = (
            np.asarray(source[arm], dtype=float) for source in (bins, independent)
        )
        for values in (raw, other):
            if (
                values.shape != (14,)
                or not np.isfinite(values).all()
                or min(values) < 0
                or max(values) > 1
                or abs(sum(values) - 1) > gates["normalization_tolerance"]
            ):
                raise ValueError("Invalid complete terminal law for " + arm)
        delta = abs(raw - other)
        discrepancies[arm] = {"bin": int(delta.argmax()), "error": float(delta.max())}
        home[arm] = float(raw[:12].sum())
    passed = (
        floor / SCORE > gates["exact_joint_home_strictly_above"]
        and home["coherent"] > gates["coherent_terminal_home_strictly_above"]
        and home["coherent"] - home["dephased"]
        >= gates["minimum_coherent_minus_dephased"]
        and home["stochastic_swap"] < gates["stochastic_swap_home_strictly_below"]
        and all(
            v["error"] <= gates["independent_bin_tolerance"]
            for v in discrepancies.values()
        )
    )
    return {"home": home, "worst_bin_by_arm": discrepancies, "passed": passed}


def replay(settings, record, spec):
    record["settings"] = settings
    floor, bins = native(settings)
    record.update(floor_ppt=floor, bins=bins)
    record["independent_bins"] = peer.probabilities(peer.BOARD, **settings)
    result = compare(bins, record["independent_bins"], floor, spec["gates"])
    record.update(result, status="passed" if result["passed"] else "rejected")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--pilot",
        action="store_true",
        help="historical cell 13, ascending, exposure 37",
    )
    group.add_argument("--status", type=Path, metavar="OPENWAVE")
    group.add_argument("--confirm", type=Path, metavar="OPENWAVE")
    parser.add_argument(
        "--record", type=Path, help="new output path for one acquisition"
    )
    args = parser.parse_args()
    spec = terms()
    anchor = None if args.pilot else official_anchor(args.status or args.confirm, spec)
    if args.status:
        if args.record:
            parser.error("--status only inspects the official history")
        print(
            json.dumps(
                {"registered": anchor is not None, "first_official_commit": anchor}
            )
        )
        return
    if args.confirm and anchor is None:
        parser.exit(
            2, "Official M10.5 task absent. No confirmatory acquisition selected.\n"
        )
    if not args.record:
        parser.error("a release requires --record")
    settings = (
        dict(home=13, start=0, direction=1, exposure=37)
        if args.pilot
        else selection(anchor, spec)
    )
    record = dict(
        kind="historical pilot" if args.pilot else "official replay",
        anchor=anchor,
        terms_sha256=TERMS_HASH,
        task_sha256=TASK_HASH,
        status="incomplete",
    )
    with args.record.open("x", encoding="utf-8") as stream:
        try:
            replay(settings, record, spec)
        except Exception as error:
            record["reason"] = str(error)
            raise
        finally:
            json.dump(record, stream, indent=2, allow_nan=False)
            stream.write("\n")
    print(json.dumps({k: record[k] for k in ("kind", "settings", "home", "passed")}))
    if not record["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
