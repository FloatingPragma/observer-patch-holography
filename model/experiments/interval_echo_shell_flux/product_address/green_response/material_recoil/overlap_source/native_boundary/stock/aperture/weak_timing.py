"""Distinguish temporal weak pulses and rare strong records at identical exposure."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigh

from . import contact_timing as clock

HERE = clock.HERE
weak = clock.weak
ARMS = ("recorded_short", "every_tick_short")


def decode(banks, old, tape, meters):
    result = clock.decode(banks, meters, tape, ARMS)
    q = Fraction(old["ticks"], 10000 * old["active_ticks"])
    # A random assembly receives the entire recorded history with probability q.
    # Condition on its public preparation flag and apply total probability.
    result["rare_recorded"] = {
        r: [float(q) * v for v in band]
        for r, band in clock.decode(old["source_banks"], meters, tape, ("recorded",))[
            "recorded"
        ].items()
    }
    return result


def run(old, tape, saved=None):
    h0, hw, meters = clock.apparatus(tape)
    if saved is not None:
        return decode(saved["source_banks"], old, tape, meters)
    bits = np.array(old["readings"]["2200"]) == 2
    n = len(bits)
    count = int(bits.sum())
    duty = Fraction(1, 10000)
    active_fraction = Fraction(count, n)
    narrow = duty / active_fraction
    assert n == old["ticks"] and count == old["active_ticks"] and 0 < narrow <= 1
    psi = weak.filtered(*eigh(h0))[0]
    empty = weak.unitary(h0, weak.STEP)
    blank = clock.time_record((empty,), psi, np.zeros(n, dtype=int))
    banks = {}
    rng = np.random.default_rng(clock.SEED + 2)
    banks["blank"] = clock.measure(blank, meters, rng)
    native = {}
    for name, flags, width in (
        (ARMS[0], bits, narrow),
        (ARMS[1], np.ones(n, dtype=int), duty),
    ):
        assert width * int(sum(flags)) / n == duty
        half = weak.unitary(h0, weak.STEP * (1 - float(width)) / 2)
        pulse = half @ weak.unitary(hw, weak.STEP * float(width)) @ half
        rho = clock.time_record((empty, pulse), psi, flags)
        banks[name] = clock.measure(rho, meters, rng)
        native[name] = clock.release_contrast(rho, blank, meters)
    native["rare_recorded"] = {
        r: float(narrow) * v
        for r, v in old["native_acceleration_change"]["recorded"].items()
    }
    decoded = decode(banks, old, tape, meters)
    for tag, row in native.items():
        for r, v in row.items():
            low, high = decoded[tag][r]
            assert low <= v <= high
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "WEAK_TIMING.md",
                HERE / "contact_timing.json",
                HERE / "mediated_release.json",
            )
        },
        exposure=str(duty),
        conditional_active_fraction=str(narrow),
        source_banks=banks,
        decoded=decoded,
        native_acceleration_change=native,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    modes = cli.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "weak_timing.json"
    old = json.loads((HERE / "contact_timing.json").read_bytes())
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    if args.counts_only:
        result = json.loads(path.read_bytes())
        assert result["decoded"] == run(old, tape, result)
    else:
        result = run(old, tape)
        if args.write:
            path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            clock.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), result)
    print(result["decoded"])
