"""Recompute the two-phase clock and spherical lens from their RER measurements."""

import argparse
from hashlib import sha256
import json
from math import isclose
from pathlib import Path

from .clock import run as clock
from .lens import run as lens

HERE = Path(__file__).resolve().parent


def verify(current, archived, key="result"):
    if (
        isinstance(current, dict)
        and isinstance(archived, dict)
        and current.keys() == archived.keys()
    ):
        for name in current:
            verify(current[name], archived[name], key + "." + name)
    elif (
        isinstance(current, list)
        and isinstance(archived, list)
        and len(current) == len(archived)
    ):
        for i, (a, b) in enumerate(zip(current, archived)):
            verify(a, b, key + "." + str(i))
    elif isinstance(current, float) and isinstance(archived, (int, float)):
        if not isclose(current, archived, rel_tol=1e-13, abs_tol=0.0):
            raise ArithmeticError(f"{key} differs: {current} vs {archived}")
    elif current != archived:
        raise ArithmeticError(f"{key} differs: {current} vs {archived}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="replace this packet's reference result"
    )
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("Run without -O to execute the record and interval checks")
    root = HERE.parents[6]
    fixture = json.loads((HERE / "native_fixture.json").read_bytes())
    for name, expected_hash in fixture["input_sha256"].items():
        if sha256((root / name).read_bytes()).hexdigest() != expected_hash:
            raise ValueError(f"Native source differs from the enrolled fixture: {name}")
    native = clock(fixture)
    shell = json.loads(
        (HERE.parents[1] / "dust_exterior" / "results.json").read_bytes()
    )
    publication = json.loads(
        (root / "model/experiments/interval_echo_clock/publication.json").read_bytes()
    )
    (stock,) = (
        s
        for s in publication["specimens"]
        if (s["cells"], s["width"], s["seed"]) == (80, 3, 31)
    )
    spherical = lens({**shell, "published_source": stock})
    data = {"native": native, "spherical": spherical}
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        verify(data, json.loads(path.read_text()))
    print(
        f"Native clock. {native['geometry_checks']} source configurations and "
        f"{native['native_columns_checked']} rational seam-column checks."
    )
    shift = native["probability_change"]["16"]["transmission"]
    print(f"Native transmission change at time 16: [{shift[0]:.17g}, {shift[1]:.17g}].")
    shot = spherical["held"]["shots"][0]["reading"]
    print(
        f"Held lens at 30 degrees: {shot[0]:.8f} rad, {shot[1]:.8f} excess time; "
        f"frequency ratio {spherical['held']['frequency_ratio']:.10f}."
    )


if __name__ == "__main__":
    main()
