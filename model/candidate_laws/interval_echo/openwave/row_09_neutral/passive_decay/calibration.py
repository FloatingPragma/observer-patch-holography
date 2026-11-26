"""Replay the separate native clocks, signed private references and comparison law."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from .material import columns, hardware
from .measure import HERE, SHOTS, WAITS, hamiltonians
from .readers import clocks, charge


def run():
    board, dimers, _ = hardware()
    s, r, _, _ = hamiltonians(len(dimers))
    rng = np.random.default_rng(26092610)
    out = dict(
        shots=SHOTS,
        source_clocks=clocks(s, 32, rng),
        receiver_clocks=clocks(r, 33, rng),
        signed_receivers=charge(board, dimers),
    )
    # On binary blank dimers the changed sign law retains X/2 exchange.
    # All contrast-two forbidden words now cost 32 for both written signs.
    # Consequently the positive B pin costs 32 and both balanced weak links vanish.
    _, vs = eigh(s)
    _, vr = eigh(r)
    basis = np.kron(vs, vr)
    strong = s - np.diag([0.0, 1.0, 0.0, 1.0])
    soft = r + np.diag([0.0, 0.0, 0.0, 1.0])
    changed = np.kron(strong, np.eye(4)) + np.kron(np.eye(4), soft)
    e, v = eigh(changed)
    far_minus = np.kron(np.eye(8), np.array([[0.5, -0.5], [-0.5, 0.5]]))
    banks = []
    for idx in (0, 3):
        stock = basis[:, 4 * idx]
        outcomes = []
        for t in WAITS:
            psi = v @ (np.exp(-1j * e * t) * (v.T @ stock))
            p = abs(basis.T @ psi) ** 2
            far = np.vdot(psi, far_minus @ psi).real
            outcomes.append(
                dict(
                    wait=t,
                    joint=rng.multinomial(SHOTS, p / p.sum()).reshape(4, 4).tolist(),
                    far_minus=int(rng.binomial(SHOTS, far)),
                )
            )
        banks.append(dict(source=idx, settings=outcomes))
    out["changed_law"] = banks
    out["changed_max_column_residual"] = round(columns(comparison=True), 13)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    rec = run()
    if dest:
        dest.write_text(json.dumps(rec, indent=1) + "\n")
    else:
        assert rec == json.loads((HERE / "calibration.json").read_text())
        print("Verified all independent clocks, charge bins and changed-law counts.")
