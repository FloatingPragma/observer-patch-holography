#!/usr/bin/env python3
"""Classify replay differences; does not change an original acceptance gate."""
import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np

from reproduce import BUNDLE, HERE, ORBITAL, save, sha
from independent import geometry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reproduction", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, default=HERE / "results/PROFILE_DIAGNOSTICS.json")
    args = parser.parse_args()
    results = {}
    for record in args.reproduction:
        receipt = json.loads(record.read_bytes())
        scratch = Path(receipt["scratch"])
        expected_path = scratch / "committed" / BUNDLE / ORBITAL / "RECEIPT.json"
        actual_path = scratch / "leaf_runs/orbital/RECEIPT.json"
        expected = json.loads(expected_path.read_bytes())
        actual = json.loads(actual_path.read_bytes())
        records_path = scratch / "committed" / BUNDLE / ORBITAL / "RECORDS.json.gz"
        tapes = json.loads(gzip.decompress(records_path.read_bytes()))
        _, lap = geometry()
        eigenvalues = (0., 5 - math.sqrt(5), 6., 5 + math.sqrt(5))
        projectors = []
        for value in eigenvalues:
            projection = np.eye(12)
            for other in eigenvalues:
                if other != value:
                    projection = projection @ (lap - other * np.eye(12)) / (value - other)
            projectors.append(projection)
        profiles, large = [], []
        for i, (a, b) in enumerate(zip(expected["population"], actual["population"])):
            for j, (x, y) in enumerate(zip(a["bands"], b["bands"])):
                u, v = np.asarray(x["ground"]["profile"]), np.asarray(y["ground"]["profile"])
                projector = projectors[j + 1]

                def identify_reference(entry):
                    spatial = np.repeat(np.array(entry["ground"]["profile"]), 3) / math.sqrt(3)
                    readouts = np.array(entry["actual_ground_mode_record_readouts"])
                    candidates = []
                    for port in range(12):
                        internal = projector[:, port] / math.sqrt(projector[port, port])
                        mode = np.kron(spatial, internal)
                        calculated = []
                        for channel in tapes["runs"][i]["channels"]:
                            y = math.sqrt(972 / (4 * 6318 * 32768)) * (
                                np.array(channel["occupation_counts"]) - 32768 / 972)
                            calculated.append(float(y @ mode))
                        candidates.append(float(np.max(abs(np.array(calculated) - readouts))))
                    best = int(np.argmin(candidates))
                    return dict(port=best, max_readout_error=candidates[best])

                profiles.append(dict(population=i, band=j,
                                     same_sign_error=float(np.max(abs(u - v))),
                                     opposite_sign_error=float(np.max(abs(u + v))),
                                     projector_error=float(np.max(abs(np.outer(u, u) - np.outer(v, v)))),
                                     energy_difference=abs(x["ground"]["slice_energy"] - y["ground"]["slice_energy"]),
                                     expected_internal_reference=identify_reference(x),
                                     actual_internal_reference=identify_reference(y),
                                     exact_projector_diagonal=entry_diagonal(j)))

        def walk(a, b, path=""):
            if isinstance(a, dict) and isinstance(b, dict):
                for key in a.keys() & b.keys():
                    walk(a[key], b[key], path + "/" + key)
            elif isinstance(a, list) and isinstance(b, list):
                for i, (x, y) in enumerate(zip(a, b)):
                    walk(x, y, path + "/" + str(i))
            elif isinstance(a, float) and isinstance(b, float) and abs(a - b) > 1e-10:
                large.append(dict(path=path, expected=a, actual=b,
                                  signed_difference=abs(a - b), magnitude_difference=abs(abs(a) - abs(b))))

        walk(expected, actual)
        results[record.parent.name] = dict(
            reproduction_sha256=sha(record), expected_receipt_sha256=sha(expected_path),
            actual_receipt_sha256=sha(actual_path), original_records_sha256=sha(records_path),
            profiles=profiles, large_differences=large,
            all_large_differences_are_sign_only=all(x["magnitude_difference"] < 1e-10 for x in large),
            note="Diagnostic threshold only. Ground signs and roundoff-sensitive argmax reference ports vary; original byte gate remains failed.")
    save(args.output, dict(script_sha256=sha(Path(__file__)), environments=results))
    print(json.dumps({name: dict(large_differences=len(r["large_differences"]),
                                all_sign_only=r["all_large_differences_are_sign_only"],
                                max_projector_error=max(x["projector_error"] for x in r["profiles"]))
                      for name, r in results.items()}, indent=2))


def entry_diagonal(band):
    # Vertex transitivity fixes every projector diagonal to rank / 12.
    return ("1/4", "5/12", "1/4")[band]


if __name__ == "__main__":
    main()
