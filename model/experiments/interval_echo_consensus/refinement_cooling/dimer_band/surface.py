"""Owner-preserving repair using the same material reference as the dimer filter."""

import argparse
import json
from math import pi

import numpy as np
from scipy.linalg import expm

from .measure import HERE, gates, pointer

CAP = 4
SHOTS = 2048
OWNERS = 20


def branches(mode):
    h, decoder, native_filter = gates("native")
    energies, eigenvectors = np.linalg.eigh(h[:2, :2])
    signs = np.vstack((eigenvectors, np.zeros(2)))
    np.testing.assert_allclose(energies, [32, 33], atol=2e-13, rtol=0)
    decoded = decoder @ signs
    np.testing.assert_allclose(abs(decoded) ** 2, np.eye(3)[:, :2], atol=3e-13, rtol=0)
    proposal = np.eye(3) if mode == "parked_proposal" else expm(-1j * pi / 2 * h)
    test_filter = np.eye(3) if mode == "parked_filter" else native_filter

    # The B/vacant C pointer resolves both proposal rails. For a failed C/D
    # read the input was B, the state is now sharp C, and a full BC pulse
    # (two half contacts) restores B. Both pointers undo after projection.
    restored = native_filter @ native_filter @ np.eye(3)[:, 2]
    np.testing.assert_allclose(abs(restored) ** 2, [0, 1, 0], atol=3e-13, rtol=0)
    table = np.zeros((2, 3))
    for old in (0, 1):
        data = proposal @ decoded[:, old]
        for rail in (0, 1):
            sharp = np.eye(3)[:, rail] * data[rail]
            filtered = test_filter @ sharp
            table[old, 2] += abs(filtered[2]) ** 2
            filtered[2] = 0
            # The accepted dark component, after undoing the energy decoder,
            # is a native stationary mode. Project here for an algebra check.
            final = decoder.conj().T @ filtered
            read = abs((decoder @ final)[:2]) ** 2
            np.testing.assert_allclose(
                read, abs(signs.conj().T @ final) ** 2, atol=3e-13, rtol=0
            )
            np.testing.assert_allclose(read[1 - rail], 0, atol=1e-25, rtol=0)
            table[old, :2] += read
    np.testing.assert_allclose(table.sum(axis=1), 1, atol=3e-13, rtol=0)
    # These certified controls have dyadic probabilities. Remove only
    # machine roundoff before drawing the reproducible categorical counts.
    table = np.round(table, 13)
    return table / table.sum(axis=1, keepdims=True)


def collect():
    pointer()
    banks = {}
    for k, mode in enumerate(("native", "parked_filter", "parked_proposal")):
        table = branches(mode)
        rng = np.random.default_rng(38066 + k)
        owners = []
        for owner in range(OWNERS):
            rows = []
            for initial in (0, 1):
                left = SHOTS
                observed = []
                for attempt in range(CAP):
                    source = initial if attempt == 0 else 1
                    ground, excited, left = rng.multinomial(left, table[source])
                    observed.append([int(ground), int(excited)])
                rows.append(dict(old=initial, accepted=observed, timeout=int(left)))
            owners.append(dict(owner=owner, rows=rows))
        banks[mode] = owners
    return dict(owners=OWNERS, shots_per_owner_and_sign=SHOTS, cap=CAP, banks=banks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = collect()
    data = json.dumps(receipt, separators=(",", ":"), allow_nan=False) + "\n"
    target = HERE / "surface.json"
    if args.write:
        target.write_text(data)
    elif target.read_text() != data:
        raise ValueError("Thermal surface receipt mismatch")
    for mode, owners in receipt["banks"].items():
        rows = [r for o in owners for r in o["rows"]]
        total = np.sum([r["accepted"] for r in rows], axis=0).sum(axis=0)
        print(mode, *total, sum(r["timeout"] for r in rows))
