"""Execute the native producer and retain every explicit primitive-event row."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess

from prepare import prepare
from tape import encode, hash_file

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def native_path(path):
    path = Path(path).resolve()
    if os.name == "nt":
        return "/mnt/" + path.drive[0].lower() + path.as_posix()[2:]
    return str(path)


def command(binary, *args):
    prefix = ["wsl", "-d", "Ubuntu", "--"] if os.name == "nt" else []
    return prefix + [native_path(binary)] + list(args)


def run(q, variant, binary, output):
    output.mkdir(parents=True, exist_ok=True)
    stem = f"q{q}_{variant}"
    input_path = output / f"q{q}.input"
    info = prepare(q, input_path)
    log_path = output / f"{stem}.stderr"
    tape_path = output / f"{stem}.tape"
    with log_path.open("wb") as log, tape_path.open("wb") as tape:
        process = subprocess.Popen(command(binary, native_path(input_path), variant), stdout=subprocess.PIPE, stderr=log)
        try:
            result = encode(process.stdout, tape)
        except BaseException:
            process.kill()
            process.wait()
            raise
        if process.wait() != 0:
            raise RuntimeError(log_path.read_text(encoding="utf-8"))
    reported = json.loads(log_path.read_text(encoding="utf-8").splitlines()[-1])
    semantic_keys = ("events", "hops", "register_reads", "register_writes", "registers",
                     "max_abs_scaled_scalar", "max_read_depth", "total_read_depth", "tree_hops_per_round")
    costs = {key: reported[key] for key in semantic_keys}
    if result["events"] != costs["events"]:
        raise ValueError("producer event count differs from retained tape")
    receipt = {"schema": "oph.source_read_routing.run.v1", "variant": variant,
               "inputs": info, "costs": costs, "tape": dict(result,
                   path=tape_path.name, sha256=hash_file(tape_path), bytes=tape_path.stat().st_size),
               "producer_sha256": hash_file(HERE / "produce.cpp"),
               "implementation_resources": {key:value for key,value in reported.items() if key not in semantic_keys}}
    (output / f"{stem}.json").write_text(json.dumps(receipt, sort_keys=True, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(receipt, sort_keys=True), flush=True)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, required=True)
    parser.add_argument("--variant", choices=("baseline", "source", "branch", "scratch"), default="baseline")
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.q, args.variant, args.binary, args.output)
