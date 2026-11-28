"""Read the unattenuated native repair signal on the same long observation windows."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigh

from . import repair_exposure as repair

clock = repair.clock
HERE = repair.HERE


def run(tape, traces, saved=None):
    h0, hw, meters = clock.apparatus(tape)
    tags = [str(seed) for seed in repair.SEEDS]

    def infer(banks):
        return {n: clock.decode(v, meters, tape, tags) for n, v in banks.items()}

    if saved is not None:
        return infer(saved["source_banks"])
    psi = repair.weak.filtered(*eigh(h0))[0]
    units = [repair.weak.unitary(h, repair.weak.STEP) for h in (h0, hw)]
    rng = np.random.default_rng(918511)
    blank = repair.age_banks(units, psi, np.zeros(max(repair.ENDS), dtype=int))
    banks = {
        str(n): {"blank": clock.measure(rho, meters, rng)} for n, rho in blank.items()
    }
    direct = {}
    for tag in tags:
        original = traces["observations"][tag]
        flags = np.unpackbits(
            np.frombuffer(bytes.fromhex(original["active_bits_hex"]), dtype=np.uint8),
            bitorder="little",
        )
        assert all(
            int(flags[:n].sum()) == original[str(n)]["active_ticks"]
            for n in repair.ENDS
        )
        direct[tag] = {}
        for n, rho in repair.age_banks(units, psi, flags).items():
            banks[str(n)][tag] = clock.measure(rho, meters, rng)
            direct[tag][str(n)] = clock.release_contrast(rho, blank[n], meters)
    decoded = infer(banks)
    for n, arms in decoded.items():
        for tag, row in arms.items():
            for r, (lo, hi) in row.items():
                assert lo <= direct[tag][n][r] <= hi
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "FULL_REPAIR.md",
                HERE / "repair_exposure.json",
                HERE / "mediated_release.json",
            )
        },
        source_banks=banks,
        decoded=decoded,
        native_acceleration_change=direct,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    modes = cli.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "full_repair.json"
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    traces = json.loads((HERE / "repair_exposure.json").read_bytes())
    if args.counts_only:
        result = json.loads(path.read_bytes())
        assert run(tape, traces, result) == result["decoded"]
    else:
        result = run(tape, traces)
        if args.write:
            path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            clock.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), result)
    print(
        {
            n: {
                r: [min(v[r][0] for v in a.values()), max(v[r][1] for v in a.values())]
                for r in ("2", "3", "4", "6")
            }
            for n, a in result["decoded"].items()
        }
    )
