"""Transfer a stationary-rate calibration to new, untrimmed native repair records."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigh

from . import contact_timing as clock

HERE = clock.HERE
weak = clock.weak
ENDS = (8192, 32768)
SEEDS = range(918501, 918505)


def age_banks(units, psi, flags):
    wave = psi.copy()
    block = np.empty((len(psi), 128), dtype=complex, order="F")
    total = np.zeros((len(psi), len(psi)), dtype=complex)
    samples = {}
    assert len(flags) == max(ENDS) and all(n % block.shape[1] == 0 for n in ENDS)
    for i, bit in enumerate(flags):
        wave = units[int(bit)] @ wave
        block[:, i % block.shape[1]] = wave
        if (i + 1) % block.shape[1] == 0:
            total += block @ block.conj().T
        if i + 1 in ENDS:
            rho = total / (i + 1)
            np.testing.assert_allclose(rho.trace(), 1, rtol=0, atol=1e-9)
            samples[i + 1] = rho
    return samples


def run(tape, old, saved=None):
    h0, hw, meters = clock.apparatus(tape)
    tags = ["regular", *(str(s) for s in SEEDS)]

    def infer(banks):
        return {str(n): clock.decode(banks[str(n)], meters, tape, tags) for n in ENDS}

    if saved is not None:
        return infer(saved["source_banks"])
    groups, edges = clock.stock.carrier()[2:4]
    tables = [np.cumsum(m, axis=1) for m in clock.stock.kernels()[:2]]
    cls = next(c for c in old["classes"].values() if "2200" in c["placements"])
    rate = Fraction(cls["seam_fractions"][1])
    assert rate == Fraction(424, 1055)
    target = Fraction(1, 10000)
    psi = weak.filtered(*eigh(h0))[0]
    off = weak.unitary(h0, weak.STEP)

    def pulse(width):
        half = weak.unitary(h0, weak.STEP * (1 - float(width)) / 2)
        return half @ weak.unitary(hw, weak.STEP * float(width)) @ half

    active = pulse(target / rate)
    regular = pulse(target)
    rng = np.random.default_rng(918510)
    blank = age_banks((off,), psi, np.zeros(max(ENDS), dtype=int))
    banks = {
        str(n): {"blank": clock.measure(rho, meters, rng)} for n, rho in blank.items()
    }
    native, observed = {}, {}
    for tag in tags:
        if tag == "regular":
            flags, width, on = np.ones(max(ENDS), dtype=bool), target, regular
        else:
            reads = clock.acquire("2200", int(tag), max(ENDS), groups, edges, tables)
            flags, width, on = np.array(reads) == 2, target / rate, active
        rho = age_banks((off, on), psi, flags)
        native[tag], observed[tag] = {}, {}
        # Encode every gate reading in chronological order, least significant bit first.
        if tag != "regular":
            observed[tag]["active_bits_hex"] = (
                np.packbits(flags, bitorder="little").tobytes().hex()
            )
        for n in ENDS:
            k = str(n)
            banks[k][tag] = clock.measure(rho[n], meters, rng)
            native[tag][k] = clock.release_contrast(rho[n], blank[n], meters)
            count = int(flags[:n].sum())
            observed[tag][k] = dict(active_ticks=count, exposure=str(width * count / n))
    decoded = infer(banks)
    for n, arms in decoded.items():
        for tag, entries in arms.items():
            for r, (lo, hi) in entries.items():
                assert lo <= native[tag][n][r] <= hi
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "REPAIR_EXPOSURE.md",
                HERE / "contact_timing.json",
                HERE / "mediated_release.json",
            )
        },
        service_rate=str(rate),
        granted_target=str(target),
        active_pulse_fraction=str(target / rate),
        observations=observed,
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
    path = HERE / "repair_exposure.json"
    old = json.loads((HERE / "contact_timing.json").read_bytes())
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    if args.counts_only:
        result = json.loads(path.read_bytes())
        assert result["decoded"] == run(tape, old, result)
    else:
        result = run(tape, old)
        if args.write:
            path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            clock.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), result)
    print({n: {k: v["6"] for k, v in a.items()} for n, a in result["decoded"].items()})
