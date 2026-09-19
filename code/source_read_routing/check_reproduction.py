"""Rebuild every small control through the public production workflow.

Archive replay alone cannot detect a broken producer, accounting command or
packer. This integration check starts from fresh inputs, executes all four
interventions, accounts a byte-identical second execution, packs every event,
and checks the result with both independent oracles and the retained tape hash.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import tempfile

from account import account
from check_control import check
from pack import pack
from run import run
from verify import load, require, verify

HERE = Path(__file__).resolve().parent


def reproduce(producer, verifier, work):
    raw, packed = work / "raw", work / "packed"
    for variant in ("baseline", "source", "branch", "scratch"):
        run(3, variant, producer, raw)
        receipt = raw / f"q3_{variant}.json"
        account(receipt, producer, work / "accounting")
        result = pack(receipt, packed, work / "packing")
        replay = verify(result, verifier, work / "replay")
        packet = load(result)
        oracle = check(packet, packed)
        reference = load(HERE / "controls" / result.name)
        require(packet["inputs"] == reference["inputs"], "regenerated input differs")
        require(packet["costs"] == reference["costs"], "regenerated census differs")
        require(replay["decoded_sha256"] == reference["tape"]["decoded_sha256"],
                "regenerated execution differs from the retained control")
        print(f"Reproduced {variant}: {oracle['events']} events, "
              f"{oracle['logical_pairs']} independently checked logical pairs.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--verifier", type=Path, required=True)
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    if args.work:
        reproduce(args.producer, args.verifier, args.work)
    else:
        with tempfile.TemporaryDirectory(prefix="oph-routing-reproduction-") as directory:
            reproduce(args.producer, args.verifier, Path(directory))
