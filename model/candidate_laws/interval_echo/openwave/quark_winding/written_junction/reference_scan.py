"""Field calibration with every reachable periodic bulk and compensating whole load."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import currents
from . import backgrounds, joint_field


def panel(row, geo, ext, rng):
    length = row["length"]
    g = backgrounds.bits(row["gauss"], length)
    data = backgrounds.block(g, np.array(geo["reference"]))
    words, heads, loads, _, tight = data[:5]
    if not np.any(tight):
        return dict(length=length, gauss=row["gauss"], tight=0, references=[])
    _, psi = backgrounds.terminal(data, ext)
    core_rows = [(i, list(joint_field.cores(loads[i]))) for i in np.flatnonzero(tight)]
    support = {
        int(j) for _, cores in core_rows for _, f, _ in cores for j in np.flatnonzero(f)
    }
    refs = []
    for phase in range(3):
        bulk = (np.arange(length) % 3 == phase).astype(int)
        for z in np.flatnonzero(bulk == 0):
            n = bulk.copy()
            n[z] += 1
            dist = min(
                min((j - int(z)) % length, (int(z) - j) % length) for j in support
            )
            for r in np.flatnonzero((loads == n).all(axis=1)):
                fractions = {}
                for sign in ("negative", "positive"):
                    good, total = 0, 0
                    for i, cores in core_rows:
                        subset = [(f, q) for tag, f, q in cores if tag == sign]
                        for f, q in subset:
                            wt = psi[i] ** 2 / len(subset)
                            delta = f @ (
                                heads[i] - heads[r] - np.roll(heads[i] - heads[r], 1)
                            )
                            assert abs(delta - (q - f[z])) < 1e-10
                            good += wt * (abs(delta - q) < 1e-10)
                            total += wt
                    fractions[sign] = good / total
                refs.append(
                    dict(
                        phase=phase,
                        site=int(z),
                        heads=int(words[r]),
                        distance=dist,
                        exact_fractions=fractions,
                    )
                )
    refs.sort(key=lambda s: (-s["distance"], s["phase"], s["site"], s["heads"]))
    result = dict(
        length=length, gauss=row["gauss"], tight=row["tight"], references=refs
    )
    if not refs:
        return result
    best = refs[0]
    full = joint_field.panel(row, best["site"], ext, geo, rng)
    result["selected"] = best
    # If an occupation has two consistent head strings, their difference is
    # constant and their cyclic divergences agree. Keep the designated head.
    for key in list(full["reads"]):
        if not key.startswith("phase_" + str(best["phase"]) + "/"):
            del full["reads"][key]
    result["selected_field"] = full["reads"]
    return result


def run(scans):
    rng = np.random.default_rng(12013003)
    return currents._portable(
        {
            k: [
                panel(r, scans["geometry"][str(r["length"])], k == "extension", rng)
                for r in v["selected"]
                if r["n"] > r["length"] // 3
            ]
            for k, v in scans["paired"].items()
        }
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    s = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    d = run(s)
    with args.output.open("x") as f:
        f.write(json.dumps(d, indent=2, allow_nan=False) + "\n")
