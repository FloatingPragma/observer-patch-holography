"""Verify the exact formation census, or replay the inclusive pair flights."""

import argparse
from hashlib import sha1, sha256
import json
from math import isfinite
from pathlib import Path

from tools.strict_json import load
from . import flight, masks


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "results.json"
BINS = {"acquired", "random_swaps", "filed_input", "center", "minus", "plus", "independent_center_and_tip"}
SCALED = {"optimal_XZ", "source_exhaustion", "reset_unavailable"}


def sources():
    declared = [
        (flight.gateway, "44f7680df55b4e68b89c301896df3e10237f046d", "f10315f165ee362bd2d1cb771b53571d3bcfb5d7"),
        (flight.tidal, "d37d200a794b49446433ac5c7306f241fed18295", "1cd2d202ce7eb899441165f3e4c63508e8e18dd7"),
        (flight.transfer.donor, "44104ae607c93382f6115181337e26846ac868d9", "9c5eb99e088316c55aeb73b22cfc6a2945deea17"),
        (flight.transfer.LAW, "09f062f77f7f540f2783f128ae6314c5d8b818ea", "2c5aa95af9ccd8cb6214f20e961c38c2604093e6"),
    ]
    bindings = {}
    for module, pin, blob in declared:
        source = Path(module.__file__).read_bytes().replace(b"\r\n", b"\n")
        if sha1(b"blob %d\0" % len(source) + source).hexdigest() != blob:
            raise ValueError(f"Different source from the donor: {module.__name__}")
        bindings[module.__name__] = {"commit": pin, "git_blob": blob}
    board = flight.transfer.filed.BOARD
    if sha256(board.read_bytes()).hexdigest() != flight.transfer.filed.BOARD_HASH:
        raise ValueError("Wrong M10.4 board")
    return {"protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
            "board_sha256": flight.transfer.filed.BOARD_HASH, "source_blobs": bindings}


def run(full):
    result = {"schema": "oph.interval_echo.seam_origin.v1", "inputs": sources(), "census": masks.census()}
    if full:
        board = flight.federation(0)
        _, home = flight.open_pair.geometry()
        flight.transfer.check_meter()
        result.update(units=flight.UNITS, pair_outcomes=[*map(list, home), "quiet"],
                      historical=flight.historical(board, home), constituents=flight.constituents(board, home),
                      ensemble=flight.ensemble(board, home))
    return result


def compare(actual, saved, key=""):
    if type(actual) is not type(saved):
        raise ValueError(f"Different type at {key}")
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError(f"Different fields at {key}")
        for k, v in actual.items():
            compare(v, saved[k], k)
    elif isinstance(actual, list):
        if key in BINS:
            if (len(saved) != len(actual) or any(type(n) is not int or n < 0 for n in saved)
                    or abs(sum(saved) - flight.UNITS) > len(saved)/2
                    or any(abs(a-b) > 2 for a, b in zip(actual, saved, strict=True))):
                raise ValueError(f"Changed terminal law at {key}")
        else:
            for a, b in zip(actual, saved, strict=True):
                compare(a, b)
    elif key in SCALED and type(actual) is int:
        if not 0 <= saved <= flight.UNITS or abs(actual-saved) > 2:
            raise ValueError(f"Changed probability at {key}")
    elif type(actual) is float:
        if not isfinite(saved) or abs(actual-saved) > 2e-10:
            raise ValueError(f"Changed numerical result at {key}")
    elif actual != saved:
        raise ValueError(f"Changed exact record at {key}")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--record", type=Path, help="create a first full receipt; refuses an existing file")
    mode.add_argument("--full", action="store_true", help="verify every history and flight")
    args = cli.parse_args()
    full = args.full or args.record is not None
    if args.record and args.record.exists():
        cli.error("record path exists")
    result = run(full)
    if args.record:
        with args.record.open("x", encoding="utf-8", newline="\n") as out:
            json.dump(result, out, indent=1, allow_nan=False)
            out.write("\n")
    else:
        saved = load(ARCHIVE)
        compare(result, saved if full else {k: saved[k] for k in result})
    print("Filed local material probability:", result["census"]["filed_address"]["filed_local_probability"])
    if full:
        print("128-bank mean home and upper limit:", result["ensemble"]["mean_home"],
              result["ensemble"]["one_sided_upper_99"])


if __name__ == "__main__":
    main()
