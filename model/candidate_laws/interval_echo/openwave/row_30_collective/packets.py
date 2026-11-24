"""Counterpropagating one-wall packets on the original three odd open strands."""

import argparse
import json
from pathlib import Path
import numpy as np
from .experiment import HERE, K, dual, verify

N_SHOTS = 2**20
SEED = 312143


def verify_conditioned_flip(site, n, angle):
    """Compile the one-wall rotation on 2/3 real adjacent dimers."""
    lo = max(0, site - 1)
    hi = min(n, site + 2)
    words = np.arange(1 << (hi - lo))
    target = site - lo
    occupied = (words >> target) & 1
    flip = np.eye(len(words))[:, words ^ (1 << target)]
    zleft = (
        (1 - 2 * ((words >> (target - 1)) & 1)) if site else np.ones(len(words), int)
    )
    zright = (
        (1 - 2 * ((words >> (target + 1)) & 1))
        if site < n - 1
        else -np.ones(len(words), int)
    )
    # The controlled-Z uses the occupied (logical bit one) ports. A physical
    # connector with the opposite convention is surrounded by native X flips.
    gate = np.ones(len(words))
    for j in range(lo, hi):
        if j != site:
            gate *= 1 - 2 * occupied * ((words >> (j - lo)) & 1)
    cz = np.diag(gate)
    unit = np.eye(len(words))
    sign = -1 if site == n - 1 else 1

    def xrot(a):
        return np.cos(a) * unit + 1j * np.sin(a) * flip

    actual = xrot(angle / 2) @ cz @ xrot(-sign * angle / 2) @ cz
    use = (1 - zleft * zright) / 2
    expected = np.diag(1 - use) + xrot(angle) @ np.diag(use)
    np.testing.assert_allclose(actual, expected, atol=1e-14)


def synthesize(target):
    # In the one-wall space an addressed dimer flip, conditioned on exactly
    # one neighbouring wall, is an isolated two-position rotation.
    a = np.zeros(len(target), dtype=complex)
    a[0] = 1
    for j in range(len(a) - 1):
        c = abs(target[j]) / np.linalg.norm(target[j:])
        s = np.sqrt(max(0, 1 - c * c))
        verify_conditioned_flip(j, len(a) - 1, np.arctan2(s, c))
        a[j], a[j + 1] = c * a[j], 1j * s * a[j]
    a *= np.exp(1j * (np.angle(target) - np.angle(a)))
    assert max(abs(a - target)) < 2e-12
    return a


def inspect(f, rng):
    # A native strip flip assigns +/- on the two wall locations and zero
    # otherwise. Use a relative quarter phase for its imaginary quadrature.
    xy = f[:-1].conj() * f[1:]
    weight = abs(f[:-1]) ** 2 + abs(f[1:]) ** 2
    totals = []
    for contrast in (2 * xy.real, 2 * xy.imag):
        p = np.column_stack(
            [(weight + contrast) / 2, (weight - contrast) / 2, 1 - weight]
        )
        reads = rng.multinomial(N_SHOTS, np.clip(p, 0, 1))
        totals.append(reads.tolist())
    re, im = (
        np.sum((np.asarray(r)[:, 0] - np.asarray(r)[:, 1]) / (2 * N_SHOTS))
        for r in totals
    )
    return float(np.arctan2(im, re)), totals


def occupations(f, u, e, v, t):
    co, si = np.cos(t * e), np.sin(t * e)
    aa, bb = (u * co) @ u.T, (v * co) @ v.T
    ab, ba = (u * si) @ v.T, -(v * si) @ u.T
    p = (aa + bb + 1j * (ba - ab)) / 2
    q = (aa - bb + 1j * (ba + ab)) / 2
    return (abs(q) ** 2).sum(axis=1) + abs(p @ f) ** 2 - abs(q @ f.conj()) ** 2


def run():
    boards = [
        b
        for b in json.loads((HERE / "results.json").read_text())["boards"]
        if b["seed"] == SEED
    ]
    meter = json.loads((HERE / "coupling_meter.json").read_text())["results"]
    m, c2 = (meter[name]["value"] for name in ("predicted_m", "predicted_c_squared"))
    out = []
    for b in boards:
        n = b["dimers"]
        assert n % 2
        rng = np.random.default_rng([833034, n])
        x = np.arange(n + 1)
        sigma = min(4.0, (n + 1) / 7)
        f = np.exp(-((x - n / 2) ** 2) / (4 * sigma**2) + 0.8j * x)
        f = synthesize(f / np.linalg.norm(f))
        k, reads = inspect(f, rng)
        other, other_reads = inspect(f.conj(), rng)
        grid = 2 * np.pi * np.fft.fftfreq(8192)
        # Fourier envelope uses the specified source width; carrier k is read.
        prob = abs(np.fft.fft(abs(f) * np.exp(1j * x * k), n=len(grid))) ** 2
        vel = c2 * np.sin(grid) / np.sqrt(m * m + c2 * (2 - 2 * np.cos(grid)))
        vpred = float(prob @ vel / prob.sum())
        times = np.arange(0, 3.01, 0.1)
        clock_shots = 2**20
        signs = []
        for norm in (np.cos(times), np.sin(times)):
            signs.append(rng.binomial(clock_shots, (1 + norm) / 2))
        hands = [2 * rr / clock_shots - 1 for rr in signs]
        observed = np.unwrap(np.arctan2(hands[1], hands[0]))
        arms = []
        for coupling in (K, 0):
            u, e, v = dual(n, coupling)
            tmat = (u - v) / 2
            # For one occupied bare orbital, <N_bog> = 1 + ||V||² -2||V^T f||².
            # Odd parity implies P(N_bog>=3) <= (<N_bog>-1)/2.
            excess = float((tmat * tmat).sum() - 2 * np.sum(abs(tmat.T @ f) ** 2))
            raw = []
            diff = []
            for idx in (0, 10, 20, 30):
                density = [occupations(a, u, e, v, times[idx]) for a in (f, f.conj())]
                rc = np.array(
                    [rng.binomial(N_SHOTS, np.clip(p, 0, 1)) for p in density]
                )
                # Symmetric source envelope and both charge signs give one
                # velocity per sign, with the same addresses and all outcomes.
                diff.append(float((x - n / 2) @ (rc[0] - rc[1]) / (2 * N_SHOTS)))
                raw.append(rc.tolist())
            t = observed[[0, 10, 20, 30]]
            speed = float(np.polyfit(t, diff, 1)[0])
            arms.append(
                dict(
                    coupling=coupling,
                    one_particle_lower_bound=max(0, 1 - excess / 2),
                    density_counts=raw,
                    half_centroid_differences=diff,
                    velocity=speed,
                )
            )
        out.append(
            dict(
                level=b["level"],
                dimers=n,
                cells=b["cells"],
                width=sigma,
                read_k=k,
                opposite_read_k=other,
                momentum_counts=[reads, other_reads],
                local_clock_counts=[s.tolist() for s in signs],
                read_times=t.tolist(),
                envelope_velocity_forecast=vpred,
                center_mode_velocity=float(
                    c2 * np.sin(k) / np.sqrt(m * m + c2 * (2 - 2 * np.cos(k)))
                ),
                native=arms[0],
                disabled=arms[1],
            )
        )
    return dict(source_seed=SEED, shots=N_SHOTS, boards=out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        verify(result, json.loads((HERE / "packets.json").read_text()))
    for b in result["boards"]:
        print(
            b["level"],
            b["dimers"],
            b["read_k"],
            b["native"]["velocity"],
            b["envelope_velocity_forecast"],
            b["disabled"]["velocity"],
            b["native"]["one_particle_lower_bound"],
        )
