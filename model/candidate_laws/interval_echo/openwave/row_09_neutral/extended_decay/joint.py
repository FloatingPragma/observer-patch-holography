"""Read the surviving core after a signed far transfer, including its final flight phase."""

import argparse
import json
import numpy as np
from scipy.linalg import eigh
from ..joined_decay.measure import core_pulse, delay_kernel, draw
from ..passive_decay.measure import SHOTS, compare_records, hamiltonians, physical
from ..passive_decay.material import columns
from .material import route
from .measure import HERE, rebuild

LABELS, PERIOD = 65536, 128


def acquire(h, initial, source, z, rng):
    e, v = eigh(h)
    dr = e[None, :] - e[:, None]
    probe = v.T @ z @ v
    traveling = np.exp(-16j * e) * (v.T @ initial)
    late = np.exp(-2j * np.pi * (LABELS - 1) / PERIOD * e)
    final = np.kron(source.T, np.eye(len(e) // 4)) @ v
    bins = np.arange(LABELS // (4 * PERIOD), 2 * LABELS // PERIOD + 1)
    outcome = np.zeros((2, 4))
    for sign in (1, -1):
        for i in range(0, len(bins), 16):
            filt = delay_kernel(dr, sign * bins[i : i + 16], labels=LABELS)
            vectors = (filt * probe) @ traveling
            cores = (vectors * late) @ final.T
            outcome[int(sign < 0)] += np.sum(
                abs(cores.reshape(-1, 4, len(e) // 4)) ** 2, axis=(0, 2)
            )
    rec = draw(rng, np.r_[outcome.ravel(), 1 - outcome.sum()])
    return dict(loss_core=rec[:4], gain_core=rec[4:8], other=rec[8])


def run():
    b, d, heads, _, _ = route((3,), True)
    s, *_ = hamiltonians(4)
    _, sv = eigh(s)
    native = rebuild(b, d, heads, False, True)
    cut = rebuild(b, d, heads, True, True)
    pulse = np.kron(core_pulse(sv), np.eye(8))
    z = np.diag(1 - 2 * (np.arange(32) % 2))
    rng = np.random.default_rng(26092616)
    stocks = []
    for arm, h, act in [
        ("joined_pulse", native, True),
        ("joined_idle", native, False),
        ("cut_pulse", cut, True),
    ]:
        _, v = eigh(h)
        psi = pulse @ v[:, 0] if act else v[:, 0]
        stocks.append((h, psi))
    primary = []
    for label, (h, psi) in zip(
        ("joined_pulse", "joined_idle", "cut_pulse"), stocks, strict=True
    ):
        primary.append(dict(arm=label, **acquire(h, psi, sv, z, rng)))
    # The comparison keeps the exact native stocks. Its three passive links
    # vanish, and B's pin loses its extra unit of energy.
    terms = physical(b, d, {pair: 32 for pair in heads})
    bare = native - np.diag(np.diag(native)) + np.diag(terms - terms[0])
    a, bits_b = (np.arange(32) >> 4) & 1, (np.arange(32) >> 3) & 1
    np.testing.assert_array_equal(np.diag(bare), 32 * a + 32 * bits_b - 32 * a * bits_b)
    x = np.kron(np.eye(16), np.array([[0.0, 1.0], [1.0, 0.0]]))
    minus = (np.eye(32) - x) / 2
    np.testing.assert_allclose(bare @ minus, minus @ bare, atol=1e-14)
    e, v = eigh(bare)
    changed = []
    for label, (_, psi) in zip(
        ("same_native_pulse", "same_native_idle"), stocks[:2], strict=True
    ):
        evolved = v @ (np.exp(-16j * e) * (v.T @ psi))
        p = float(np.vdot(evolved, minus @ evolved).real)
        result = acquire(bare, psi, sv, z, rng)
        changed.append(
            dict(
                arm=label,
                far_minus=int(rng.binomial(SHOTS, p)),
                loss=sum(result["loss_core"]),
                gain=sum(result["gain_core"]),
                other=result["other"],
            )
        )
    return dict(
        seed=26092616,
        shots=SHOTS,
        labels=LABELS,
        read_time=16,
        common_final_free_duration=float(2 * np.pi * (LABELS - 1) / PERIOD),
        joint=primary,
        changed=changed,
        changed_seam_residual=round(columns(comparison=True), 13),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    r = run()
    path = HERE / "joint.json"
    if args.save:
        path.write_text(json.dumps(r, indent=1) + "\n")
    else:
        compare_records(json.loads(path.read_text()), r)
    print(r)
