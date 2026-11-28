"""Wavelength test for each physical five-connector record on the six-dimer strand."""

import argparse
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

from . import unsorted
from .measure import HERE, write_data
from .single_write import enroll


def run(counts_directory=None):
    level, route = enroll()
    gaps, weights, scales, names = unsorted.spectrum()
    if counts_directory is not None:
        counts_directory.mkdir(parents=True, exist_ok=False)
    rows = []
    for index, bits in enumerate(product((0, 1), repeat=len(route) - 1)):
        record = "".join(map(str, bits))
        selected = slice(index, index + 1)
        means, _ = unsorted.traces(gaps[selected], weights[selected])
        with TemporaryDirectory() as temporary:
            folder = Path(temporary)
            if counts_directory is not None:
                folder = counts_directory / record
                folder.mkdir()
            reading = unsorted.observe(
                32 + index, means, scales, names, folder, label=record
            )
        spatial = reading["spatial"]
        if not spatial["records_complete"]:
            raise ValueError(f"Incomplete or invalid spectrum for {record}: {spatial}")
        rows.append(
            dict(
                record=record,
                status=spatial["status"],
                profile_pass=spatial["profile_pass"],
                detail=spatial.get("detail"),
                odd_profile_q=spatial.get("candidate_q"),
                quadratic_cutoff=spatial.get("quadratic_cutoff"),
                launch_label_rejected=spatial["launch_label_rejected"],
                pair_sum_closure=reading["omitted_pair_closure"],
                count_sha256=reading["count_sha256"],
            )
        )
    passed = [row["record"] for row in rows if row["profile_pass"]]
    return dict(
        level=level,
        route=route,
        histories=len(rows),
        profile_pass_records=passed,
        native_weight_passing_profile=len(passed) / len(rows),
        rows=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--save-counts", type=Path, help="Create a directory for all 544 count arrays")
    args = parser.parse_args()
    result = run(args.save_counts)
    payload = write_data(result) + b"\n"
    receipt = HERE / "conditional_profiles.json"
    if args.record:
        with receipt.open("xb") as output:
            output.write(payload)
    else:
        if payload != receipt.read_bytes():
            raise SystemExit("Conditional profile replay differs from the receipt")
    print(f"{len(result['profile_pass_records'])}/{result['histories']} profiles pass")
