"""Half-flux period of the isolated ten-state wall cluster."""

import argparse
import json
from fractions import Fraction
from math import ceil
import numpy as np
from scipy.linalg import eigh
from . import measure as cur

HERE = cur.HERE
M = 2**23
DT = 1 / 128
SHOTS = 2**20
STEP = 2 * np.pi / (M * DT)
SEED = 122014


def enclosure(h, e, u):
    """Weyl bound from an integer Gram matrix, residual and polar factor."""
    s = 2**26
    A = np.rint(2 * h).astype("int64")
    np.testing.assert_allclose(2 * h, A, atol=1e-13, rtol=0)
    Q = np.rint(s * u).astype("int64")
    D = np.rint(s * e).astype("int64")
    limit = 2**63
    norm_a = int(abs(A).sum(axis=1).max())
    amp = int(abs(Q).max())
    if (
        max(s * norm_a * amp + 2 * amp * int(abs(D).max()), len(h) * amp**2 + s**2)
        >= limit
    ):
        raise ArithmeticError("Integer matrix bounds exceed 64 bits")
    gram = Q.T @ Q - s * s * np.eye(len(h), dtype="int64")
    residual = s * (A @ Q) - 2 * Q * D
    if len(h) * max(int(abs(gram).max()), int(abs(residual).max())) >= limit:
        raise ArithmeticError("Row norm may overflow")
    g = Fraction(int(abs(gram).sum(axis=1).max()), s * s)
    if g >= 1:
        raise ValueError("Singular frame")
    r = Fraction(max(int(abs(residual).sum(axis=i).max()) for i in [0, 1]), 2 * s * s)
    error = r + (Fraction(norm_a, 2) + Fraction(int(abs(D).max()), s)) * g
    return {
        "dyadic_scale": s,
        "gram_upper": str(g),
        "residual_upper": str(r),
        "energy_error_upper": str(error),
        "energy_error": float(error),
        "eigenvalue_numerators": D.tolist(),
    }


def resolution(sec, h):
    e, u = eigh(h)
    src = int(cur.belt.selected(sec)[0]) if len(h) == len(sec.states) else 0
    cert = enclosure(h, e, u)
    # Merge degenerate frequency lines before sampling; this is invariant
    # under rotation of the numerical basis in a degenerate eigenspace.
    vals, inv = np.unique(np.round(e, 10), return_inverse=True)
    w = np.round(np.bincount(inv, weights=abs(u[src]) ** 2), 12)
    return vals, w / w.sum(), e, cert


def qpe_samples(values, weights, rng):
    line = rng.choice(len(values), p=weights, size=SHOTS)
    u = values[line] / STEP
    nearest = np.rint(u).astype("int64")
    delta = u - nearest
    bins = nearest.copy()
    pending = np.arange(SHOTS)
    while len(pending):
        # Offset n has proposal q(0)=1/2, q(±n)=1/[2(n+1)(n+2)].
        # For the 2^23-bin clock 12q dominates the wrapped Fejer law.
        mag = (
            np.floor(1 / np.maximum(rng.random(len(pending)), 2**-53)).astype("int64")
            - 1
        )
        offset = mag * (2 * rng.integers(0, 2, len(pending)) - 1)
        x = delta[pending] - offset
        # Only evaluate the unique wrapped period's offsets.
        valid = (offset >= -M // 2) & (offset < M // 2)
        exact = np.zeros(len(pending))
        exact[valid] = (np.sinc(delta[pending][valid]) / np.sinc(x[valid] / M)) ** 2
        shifted = valid & (offset != 0)
        exact[shifted] *= (delta[pending][shifted] / x[shifted]) ** 2
        q = np.where(
            mag == 0, 0.5, 0.5 / (mag.astype(float) + 1) / (mag.astype(float) + 2)
        )
        accept = rng.random(len(pending)) * 12 * q < exact
        bins[pending[accept]] += offset[accept]
        pending = pending[~accept]
    bins = (bins + M // 2) % M - M // 2
    return np.unique(bins, return_counts=True)


def run():
    rng = np.random.default_rng(SEED)
    rows = []
    banks = {}
    for n in (3, 4):
        sec = cur.BinarySector(cur.belt.ring(10), n)
        sign = cur.counter_signs(sec, 0)
        for mode in ["native", "sharp", "exchange_only"]:
            h = sec.hamiltonian(0 if mode == "exchange_only" else cur.KAPPA)
            ids = cur.belt.selected(sec) if mode == "sharp" else np.arange(len(h))
            hx = h[np.ix_(ids, ids)]
            hz = (h * sign)[np.ix_(ids, ids)]
            resolutions = [resolution(sec, g) for g in [hx, hz]]
            lines = [r[2][-10:] for r in resolutions]
            distance = abs(lines[0][:, None] - lines[1]).min(axis=1)
            # The plus-counter spectrum nominations use all visible lines.
            vs, w, _, _ = resolutions[0]
            group_weight = np.array([w[np.argmin(abs(vs - x))] for x in lines[0]])
            distance[group_weight < 0.005] = -1
            picked = int(np.argmax(distance))
            det = float(distance[picked])
            target = float(lines[0][picked])
            radius = max(3 * STEP, det / 3)
            lo, hi = int(ceil((target - radius) / STEP)), int(
                np.floor((target + radius) / STEP)
            )
            bounds = []
            audits = []
            counts = []
            for bit, (vs, w, e, cert) in enumerate(resolutions):
                bins, freq = qpe_samples(vs, w, rng)
                banks[f"{n}_{mode}_{bit}"] = np.stack([bins, freq])
                found = int(freq[(bins >= lo) & (bins <= hi)].sum())
                counts.append(found)
                bounds.append(cur.cp(found, SHOTS, 0.01 / 12))
                audits.append(
                    dict(
                        energy=lines[bit].tolist(),
                        gap_below=float(e[-10] - e[-11]) if len(e) > 10 else None,
                        certificate=cert,
                    )
                )
            err = sum(x["certificate"]["energy_error"] + 1 / 2**27 for x in audits)
            rows.append(
                dict(
                    n=n,
                    mode=mode,
                    clock_window=[lo, hi],
                    chosen_energy=target,
                    min_distance=det,
                    certified_lower=max(0, det - err),
                    width=radius,
                    in_window=counts,
                    intervals=bounds,
                    branches=audits,
                )
            )
    return dict(shots=SHOTS, clock_bins=M, clock_step=DT, seed=SEED, rows=rows), banks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt, banks = run()
    dst = HERE / "band.json"
    if args.write:
        dst.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(HERE / "clock.npz", **banks)
    else:
        if json.loads(dst.read_bytes()) != receipt:
            raise ValueError("Band replay changed")
        with np.load(HERE / "clock.npz") as old:
            for key, v in banks.items():
                np.testing.assert_array_equal(v, old[key])
    for r in receipt["rows"]:
        print(r["n"], r["mode"], r["certified_lower"], r["in_window"])


if __name__ == "__main__":
    main()
