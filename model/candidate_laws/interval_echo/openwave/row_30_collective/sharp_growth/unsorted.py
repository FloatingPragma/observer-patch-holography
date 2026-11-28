"""Apply the same phase and wavelength receivers directly to all sharp repairs."""

import argparse
import gzip
import hashlib
import json
from itertools import combinations, product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy.linalg import eigh

from .. import experiment
from ..confirmation import spatial
from .growth import SEED
from .measure import HERE, write_data
from .single_write import columns, enroll

N = 6
M = N + 1
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
LABELS = ["unsorted", "constant_0", "constant_1"]


def sources():
    spin, ops = columns(N)
    physical_pairs = list(combinations(range(2 * M), 2))
    addresses = list(combinations(range(M), 2))
    left, right = np.array(addresses).T
    aa = ops[[physical_pairs.index((2 * j, 2 * k)) for j, k in addresses]]
    profiles = np.sqrt(2 / (M + 1)) * np.sin(
        np.pi * np.outer(np.arange(1, M + 1), np.arange(1, 5)) / (M + 1)
    )

    def wedge(f, g):
        coefficients = f[left] * g[right] - f[right] * g[left]
        return coefficients @ aa.reshape(len(aa), -1), sum(abs(coefficients))

    settings, scales, names = [], [], []
    for i, j in PAIRS:
        f, g = profiles[:, i], profiles[:, j]
        src, src_scale = wedge(f, g)
        probes, factors = [src], [src_scale**2]
        for idx, ref in enumerate((f, g)):
            for k, basis in enumerate(np.eye(M)):
                op, scale = wedge(ref, basis)
                probes.append(op)
                factors.append(src_scale * scale)
        settings.append(np.reshape(probes, (-1, 2**N, 2**N)))
        scales.extend(factors)
        names.append(f"p{i + 1}_{j + 1}_")
    return spin, settings, np.asarray(scales), names


def spectrum():
    spin, settings, scales, names = sources()
    coords = np.arange(len(spin))
    facing = (spin[:, :-1] == -1) & (spin[:, 1:] == 1)
    gaps, weights = [], []
    for bits in product((0, 1), repeat=N - 1):
        h = np.diag((1 - facing) @ (4 * experiment.K + np.array(bits)))
        for i in range(N):
            h[coords, coords ^ (1 << i)] -= 0.5
        e, v = eigh(h)
        gaps.append((e[:, None] - e[None, :]).ravel())
        rows = []
        for ops in settings:
            probes = (v.T @ ops @ v).reshape(len(ops), -1)
            weight = probes[0].conj() * probes / len(spin)
            np.testing.assert_allclose(weight.imag, 0, atol=2e-14)
            rows.extend(weight.real)
        weights.append(np.array(rows).T)
    return np.array(gaps), np.array(weights), scales, names


def traces(gaps, weights, count=experiment.SAMPLES, dt=experiment.DT):
    """Fourier interpolation at uniform hand appointments with a bounded remainder."""
    t = np.arange(count) * dt
    size = 4 * count
    df = 2 * np.pi / (size * dt)
    g = gaps.ravel()
    peak = np.rint(g / df).astype(int)
    # The residual phase magnitude is at most pi/4 in the observation window.
    delta = (g - peak * df) * (count * dt)
    weight = weights.reshape(len(g), -1).T / len(gaps)
    bound = np.max(np.sum(np.abs(weight), axis=1)) * (np.pi / 4) ** 14 / 87178291200
    if bound > 1e-10:
        raise ArithmeticError("Uniform-time Fourier approximation exceeds tolerance")
    idx = peak % size
    out = np.zeros((weight.shape[0], count))
    powers, phase = np.ones(len(g)), np.ones(count, complex)
    for order in range(14):
        grid = np.stack(
            [np.bincount(idx, row * powers, minlength=size) for row in weight]
        )
        fft = np.fft.ifft(grid, axis=1)[:, :count] * size
        out += (fft * phase).real
        powers *= delta / (order + 1)
        phase *= 1j * t / (count * dt)
    # Audit against direct spectral sums at the beginning, interior and end.
    hold = [0, 11, count // 2, count - 1]
    direct = weight @ np.cos(g[:, None] * t[hold])
    np.testing.assert_allclose(out[:, hold], direct, atol=1e-10, rtol=0)
    return out, float(bound)


def observe(panel, out, scales, names, folder, *, label=None):
    times = experiment.DT * np.arange(experiment.SAMPLES)
    rng = np.random.default_rng([SEED, 50 + panel])
    hashes = {}

    def store(name, data):
        data = np.asarray(data, dtype="<u8")
        np.save(folder / (name + ".npy"), data, allow_pickle=False)
        hashes[name] = hashlib.sha256(data.tobytes()).hexdigest()

    raw_hand = []
    for name, hand in [("clock_cos", np.cos(times)), ("clock_sin", np.sin(times))]:
        cs = rng.binomial(spatial.CLOCK_SHOTS, (1 + hand) / 2)
        store(name, cs)
        raw_hand.append(2 * cs / spatial.CLOCK_SHOTS - 1)
    phase = np.unwrap(np.arctan2(raw_hand[1], raw_hand[0]))
    step = np.polyfit(np.arange(len(phase)), phase, 1)[0]
    sums = {}
    modes = [dict(spatial_launch=i + 1, k=None) for i in range(4)]
    for ix, (prefix, vv, norm) in enumerate(
        zip(
            names,
            out.reshape(5, 1 + 2 * M, -1),
            scales.reshape(5, 1 + 2 * M),
            strict=True,
        )
    ):
        p = (1 + vv / norm[:, None]) / 2
        assert -1e-12 <= p.min() and p.max() <= 1 + 1e-12
        cs = rng.binomial(spatial.SIGN_SHOTS, np.clip(p, 0, 1))
        store(prefix + "sum", cs[0])
        for j in range(2):
            store(prefix + f"scan{j}", cs[1 + j * M : 1 + (j + 1) * M].T)
        w = experiment.projection(2 * cs[0] / spatial.SIGN_SHOTS - 1, None, step)
        i, j = PAIRS[ix]
        sums[f"{i + 1},{j + 1}"] = w
    checked = spatial.score(
        dict(dimers=N, observer=dict(modes=modes, pair_sum_lines=sums)), folder
    )
    return dict(
        label=LABELS[panel] if label is None else label,
        count_sha256=hashes,
        pair_lines=sums,
        omega_first_from_three=(sums["1,2"] + sums["1,3"] - sums["2,3"]) / 2,
        omitted_pair_closure=abs(
            sums["2,4"] - (sums["1,4"] - sums["1,3"] + sums["2,3"])
        ),
        spatial=checked,
    )


def run(directory=None):
    level, route = enroll()
    gaps, weights, scales, names = spectrum()
    stored = json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))
    meters = [
        (e["results"]["J"]["value"], e["results"]["h"]["value"])
        for e in stored["meters"]
    ]
    meters.append(tuple(np.mean(meters, axis=0)))
    entries = []
    for i, idx in enumerate((slice(None), slice(0, 1), slice(-1, None))):
        out, err = traces(gaps[idx], weights[idx])
        with TemporaryDirectory() as tmp:
            folder = directory / LABELS[i] if directory else Path(tmp)
            folder.mkdir(exist_ok=True, parents=True)
            row = observe(i, out, scales, names, folder)
        row["fourier_uniform_error_bound"] = format(err, ".6g")
        above_nyquist = np.abs(gaps[idx]) >= np.pi / experiment.DT
        folded = (abs(weights[idx]) * above_nyquist[:, :, None]).sum(axis=(0, 1))
        row["largest_aliased_signed_mean_bound"] = format(
            float(max(folded / len(gaps[idx]) / scales)), ".6g"
        )
        record = row["spatial"]
        if record["profile_pass"]:
            k = np.asarray(record["even_trained_k"])
            seen = np.asarray(record["omega_from_pair_counts"])
            row["max_lattice_errors"] = [
                float(
                    np.max(
                        abs(
                            np.sqrt(4 * (j - h) ** 2 + 16 * j * h * np.sin(k / 2) ** 2)
                            - seen
                        )
                    )
                )
                for j, h in meters
            ]
        entries.append(row)
    return dict(level=level, route=route, meters=meters, panels=entries)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", action="store_true")
    cli.add_argument("--save-counts", type=Path)
    args = cli.parse_args()
    result = run(args.save_counts)
    rec = HERE / "unsorted.json"
    if args.record:
        with rec.open("xb") as f:
            f.write(write_data(result) + b"\n")
    else:
        assert rec.read_bytes() == write_data(result) + b"\n"
    for e in result["panels"]:
        print(
            e["label"],
            e["spatial"]["status"],
            e.get("max_lattice_errors"),
            e["omitted_pair_closure"],
        )
