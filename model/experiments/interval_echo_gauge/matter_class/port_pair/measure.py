"""Routed Ramsey comparison of two labeled half-turns, retaining every start."""

import argparse
import json
from math import ceil, log, pi, sqrt
from pathlib import Path

import numpy as np

from .. import measure as matter
from ...scalar_grammar import measure as pulse
from ...scalar_grammar.shell_line.apparatus import transport


HERE = Path(__file__).resolve().parent
SEED = 26112531
SHOTS = 65536
# (port, signed exposure in half-turns)
MENUS = {
    "blank": (),
    "one": ((0, 1),),
    "same": ((0, 1), (0, 1)),
    "pair": ((0, 1), (1, 1)),
    "reversed_polarity": ((0, 1), (1, -1)),
    "cut": ((0, 1), (1, 1)),
}
SIGNED_TARGET = {
    "blank": 1,
    "one": 0,
    "same": -1,
    "pair": -1 / sqrt(5),
    "reversed_polarity": 1 / sqrt(5),
    "cut": 0,
}


def acquire():
    rng = np.random.default_rng(SEED)
    site = transport()
    words = [r[1] for r in site["address_clicks"]]
    (empty,) = {r[2] for r in site["address_clicks"]}
    words.append(empty)
    route = np.eye(3)[:, site["connected_columns"]]
    mounts, _ = matter.frame_program()
    banks = []
    for name, menu in MENUS.items():
        word, _ = matter.motion([(port, sign * pi) for port, sign in menu])
        source_gates = pulse.execute(word)
        for level in (0, 1):
            pair = ("pair", level, 2, 0.0)
            prepared = pulse.execute(pulse.prep(pair))[:, 0]
            received = route @ source_gates @ prepared
            for remounted in (False, True):
                # Physically turn the receiving frame and use its calibrated inverse reader.
                mounted = pulse.execute(mounts) @ received if remounted else received
                undo = pulse.inverse(mounts) if remounted else []
                final = pulse.execute(undo + pulse.reader(pair)) @ mounted
                probabilities = np.r_[abs(final) ** 2, 0.0]
                if name == "cut":
                    # Every source basis address gives the same empty raw word on the cut.
                    # With no local load all terminal gates leave the receiver in vacuum.
                    probabilities = np.array([0.0, 0.0, 0.0, 1.0])
                counts = pulse.grid_counts(rng, SHOTS, probabilities)
                banks.append(
                    dict(
                        menu=name,
                        level=level,
                        remounted=remounted,
                        counts=counts.tolist(),
                    )
                )
    return dict(
        stations=site["owners"],
        route_contacts=site["path_contacts"],
        detector_words=words,
        seed=SEED,
        starts_per_bank=SHOTS,
        banks=banks,
    )


def decode(result):
    """Return all-start signed means for the entire enrollment and their common width."""
    banks = result["banks"]
    n = result["starts_per_bank"]
    words = result["detector_words"]
    if n != SHOTS or words != ["000101", "010001", "010100", "010101"]:
        raise ValueError("Incomplete exposure or receiver calibration")
    expected = {
        (name, level, frame)
        for name in MENUS
        for level in (0, 1)
        for frame in (False, True)
    }
    contrasts = {}
    for bank in banks:
        key = bank["menu"], bank["level"], bank["remounted"]
        c = bank["counts"]
        if (
            key not in expected
            or key in contrasts
            or len(c) != 4
            or any(type(v) is not int or v < 0 for v in c)
            or sum(c) != n
        ):
            raise ValueError("Unknown, repeated or incomplete receiver bank")
        if key[0] == "cut" and c[3] != n:
            raise ValueError("Cut transfer produced an occupied remote reading")
        contrasts[key] = (c[key[1]] - c[2]) / n
    if set(contrasts) != expected:
        raise ValueError("Missing receiver bank")
    # Finite grid rounding of a two-output contrast contributes at most 8/2^40.
    radius = ceil(1e9 * (sqrt(2 * log(2 * len(expected) / 0.01) / n) + 8 / 2**40)) / 1e9
    return contrasts, radius


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    counts = acquire()
    path = HERE / "records.json"
    if args.write:
        path.write_text(json.dumps(counts, indent=1) + "\n")
    elif counts != json.loads(path.read_text()):
        raise SystemExit("Replayed receiver counts differ")
    contrasts, radius = decode(counts)
    for name, target in SIGNED_TARGET.items():
        observed = [value for (menu, _, _), value in contrasts.items() if menu == name]
        print(
            f"{name}: {min(observed):.6f} to {max(observed):.6f}, target {target:.6f}"
        )
        if max(abs(v - target) for v in observed) > radius:
            raise SystemExit(f"Signed comparison failed in {name}")
    print(f"Simultaneous 99% radius {radius:.9f}")


if __name__ == "__main__":
    main()
