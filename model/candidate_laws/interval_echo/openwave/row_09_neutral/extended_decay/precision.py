"""Measure the longer signed delay aperture on two lengths and the crossing geometry."""

import argparse
import json
from .measure import HERE, run, compare_records

CASES = {
    "four": ((4,), 65536, 26092613, False),
    "three": ((3,), 65536, 26092614, False),
    "crossing": ((3,), 65536, 26092615, True),
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    for key, params in CASES.items():
        r = run(*params)
        dest = HERE / (key + ".json")
        if args.save:
            dest.write_text(json.dumps(r, indent=1) + "\n")
        else:
            compare_records(json.loads(dest.read_text()), r)
        print(
            key,
            r["dimers"],
            r["geometry"]["raw_seams_to_each_end"],
            r["receiver_formation"],
        )
        for panel in r["panels"]:
            print(
                panel["arm"],
                panel["false_loss_bound"],
                [v["loss_gain_other"][0] for v in panel["energy"]],
                [sum(w[1] for w in f["core_by_far_x"]) for f in panel["flight"]],
            )
