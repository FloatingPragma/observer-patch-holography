"""Acquire the M10.3 held release using the public registration commit as the selector."""

import argparse
import base64
import hashlib
import json
import re
import zlib
from pathlib import Path

import numpy as np
from scipy import sparse

from row_04_check import MARGIN, SHOTS, evolve, local_law, screen


EXPOSURES = (3, 11, 29)
HERE = Path(__file__).resolve().parent


def acquire(commit, board):
    key = hashlib.sha256(b"M10.3:row4:" + bytes.fromhex(commit)).digest()
    eligible = [cell for cell in range(20) if cell not in (0, 13)]
    cell = eligible[int.from_bytes(key[:4], "big") % len(eligible)]
    edges, index, exchange, adjacent = screen(board)
    home = edges[(edges[:, 0] // 12 == cell) & (edges[:, 1] // 12 == cell)]
    gateways = np.unique(edges[edges[:, 0] // 12 != edges[:, 1] // 12])
    sources = home[np.isin(home, gateways).any(axis=1)]
    launch = sources[int.from_bytes(key[4:8], "big") % len(sources)]
    seed = int.from_bytes(key[8:16], "big")
    patch = index[home[:, 0], home[:, 1]]
    report = {"registration_commit": commit, "cell": cell, "launch": launch.tolist(),
              "outcomes": home.tolist() + ["quiet"], "shots": SHOTS,
              "simultaneous_radius": float(MARGIN), "arms": []}
    for arm, kappa in enumerate([8, 0]):
        h = exchange - 2 * kappa * sparse.diags(adjacent)
        records = []
        for j, state in enumerate(evolve(h, index[tuple(launch)], EXPOSURES)):
            probs = local_law(state, patch)
            rng = np.random.default_rng([seed, arm, j])
            draws = rng.choice(len(probs), size=SHOTS, p=probs / probs.sum()).astype("uint8")
            counts = np.bincount(draws, minlength=len(probs))
            records.append({"exposure": EXPOSURES[j], "counts": counts.tolist(),
                            "draws_zlib_base64": base64.b64encode(zlib.compress(draws.tobytes(), 9)).decode()})
        report["arms"].append({"kappa": kappa, "records": records})
    native_min = min(sum(r["counts"][:-1]) / SHOTS for r in report["arms"][0]["records"])
    final_control = sum(report["arms"][1]["records"][-1]["counts"][:-1]) / SHOTS
    report["passed"] = bool(native_min - MARGIN > 0.80 and final_control + MARGIN < 0.65)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("registration_commit", help="first accepted commit introducing the terms on official main")
    parser.add_argument("--board", type=Path, default=HERE / "FEDERATION_LEVEL0.npz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.registration_commit):
        parser.error("supply the full lowercase public registration commit")
    report = acquire(args.registration_commit, args.board)
    with args.output.open("x") as out:
        json.dump(report, out, indent=2)
        out.write("\n")
    print("PASS" if report["passed"] else "REJECT")
    raise SystemExit(0 if report["passed"] else 1)
