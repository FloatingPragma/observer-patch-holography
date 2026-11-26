"""Free fork interference, reversible source preparation and terminal records."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ..measure import fixture, material
from ..scalar import wall_laplacian

HERE = Path(__file__).resolve().parent
SHOTS = 2**28
SCALE = 2**40
SEED = 830619
TIMES = np.r_[0.0, 0.0625, 0.25, 0.5, np.arange(1.0, 65.0)]
CONTRAST = (8, 4, 0)
SOURCES = ("B", "C", "plus", "minus")
WORDS = np.arange(128)
# Input B's inner gateway is the phase pivot. The specified paths use only
# addressed neighboring dimers. Each SWAP is three binary controlled flips.
PIVOT = 3
PATHS = ((3, 4), (3, 0, 5), (3, 0, 5, 6))


def cnot(indices, a, b):
    return indices ^ (((indices >> a) & 1) << b)


def swap(indices, a, b):
    return cnot(cnot(cnot(indices, a, b), b, a), a, b)


def encoder(launch):
    indices = WORDS.copy()
    for path in PATHS:
        before = list(zip(path[:-2], path[1:-1]))
        for a, b in before:
            indices = swap(indices, a, b)
        indices = cnot(indices, *path[-2:])
        for a, b in reversed(before):
            indices = swap(indices, a, b)
    return indices ^ (launch & (1 << PIVOT))


def sources(wall, orientation):
    launch = next(w for w in WORDS[wall == 4] if w % 2 == orientation)
    permutation = encoder(launch)
    base = launch & ~(1 << PIVOT)
    target = base ^ (1 << PIVOT)
    source = np.zeros((len(WORDS), 4))
    source[permutation[base]] = [1, 0, 2**-0.5, 2**-0.5]
    source[permutation[target]] = [0, 1, 2**-0.5, -(2**-0.5)]
    if permutation[base] != launch or wall[permutation[target]] != 16:
        raise ValueError("Source gate sequence failed its wall-word readback")
    # Inverse fanout, then X or Y analysis on the pivot.
    recovered = source[permutation]
    coherence = recovered[base].conj() * recovered[target]
    p = 0.5 + np.stack((coherence.real, coherence.imag), axis=-1)
    return source, np.stack((p, 1 - p), axis=-1)


def cdf40(probabilities):
    probabilities = np.clip(probabilities, 0, 1)
    total = probabilities.sum(axis=-1, keepdims=True)
    if np.max(abs(total - 1)) > 2e-11:
        raise ArithmeticError("Detector probabilities lost normalization")
    return np.rint(np.cumsum(probabilities / total, axis=-1) * SCALE).astype("int64")


def calculate():
    _, _, links = fixture(federation(0), 2)
    # Port 82 is the right end of A's gateway; 74 is the left end of its
    # private dimer. Their actual occupation word is read within cell six.
    remote = 2 * ((WORDS >> 1) & 1) + (1 - ((WORDS >> 2) & 1))
    full = np.empty((3, 2, 4, len(TIMES), 128), dtype=np.int64)
    local = np.empty((*full.shape[:-1], 4), dtype=np.int64)
    calibration = np.empty((2, 4, 2, 2), dtype=np.int64)
    for i, coupling in enumerate(CONTRAST):
        h, wall = material(links, coupling)
        energy, basis = np.linalg.eigh(h.toarray())
        phases = np.exp(-1j * TIMES[:, None] * energy)
        for r in (0, 1):
            states, measurement = sources(wall, r)
            calibration[r] = cdf40(measurement)
            for j, source in enumerate(states.T):
                amplitude = (phases * (basis.T @ source)) @ basis.T
                p = abs(amplitude) ** 2
                full[i, r, j] = cdf40(p)
                local[i, r, j] = cdf40(
                    np.stack([p[:, remote == w].sum(axis=1) for w in range(4)], axis=-1)
                )
    hand = 0.5 + 0.5 * np.stack((np.cos(TIMES), np.sin(TIMES)), axis=-1)
    return dict(
        full=full,
        local=local,
        source=calibration,
        hand=cdf40(np.stack((hand, 1 - hand), axis=-1)),
    )


def draw(laws):
    rng = np.random.default_rng(SEED)
    records = {}
    for key in ("full", "local", "source", "hand"):
        cdf = laws[key]
        p = np.diff(cdf, prepend=0, axis=-1) / SCALE
        records[key] = np.array(
            [rng.multinomial(SHOTS, row) for row in p.reshape(-1, p.shape[-1])],
            dtype=np.int32,
        ).reshape(p.shape)
    return records


def band(count):
    # Invert the two-sided Bernoulli Bernstein inequality with variance <= p.
    # The union includes <50,000 probabilities. Marginals of rounded
    # 128-category laws differ from their ideal marginals by <=128/SCALE.
    empirical = np.asarray(count) / SHOTS
    x = np.log(2 * 50000 / 0.01) / SHOTS
    a, b = np.sqrt(2 * x), 2 * x / 3
    low = (np.sqrt(np.maximum(0, empirical - b) + a * a / 4) - a / 2) ** 2
    high = (np.sqrt(empirical + b + a * a / 4) + a / 2) ** 2
    return np.maximum(0, low - 128 / SCALE), np.minimum(1, high + 128 / SCALE)


def compare_phases(counts):
    plus, minus = counts[2:4]
    lo_p, hi_p = band(plus)
    lo_m, hi_m = band(minus)
    lower, upper = lo_p - hi_m, hi_p - lo_m
    distance = np.maximum(0, np.maximum(lower, -upper))
    resolved = np.flatnonzero(distance)
    magnitude = abs(plus - minus)
    peak = (
        int(resolved[np.argmax(magnitude[resolved])])
        if resolved.size
        else int(np.argmax(magnitude))
    )
    return dict(
        earliest=int(resolved[0]) if resolved.size else None,
        peak=peak,
        plus=int(plus[peak]),
        minus=int(minus[peak]),
        interval=[float(lower[peak]), float(upper[peak])],
    )


def scalar_comparison(times, time_error, count):
    meters = json.loads((HERE.parents[1] / "coupling_meter.json").read_text())[
        "results"
    ]
    m, b = (meters[key]["value"] for key in ("predicted_m", "predicted_c_squared"))
    _, _, links = fixture(federation(0), 2)
    spatial, v = np.linalg.eigh(wall_laplacian(links))
    spectrum = np.sqrt(m * m + b * spatial) - m
    source = (np.eye(6)[2] + np.eye(6)[4]) / np.sqrt(2)
    time = float(times[1])
    amplitude = abs(v[0] @ (np.exp(-1j * time * spectrum) * (v.T @ source)))
    lmax = max(spatial)
    mlo = 0.99 * m
    variation = 0.01 * m * 1.01 * b * lmax / (2 * mlo * mlo) + 0.01 * b * lmax / (
        2 * mlo
    )
    error = (time + time_error) * variation + time_error * max(spectrum)
    forecast = [max(0.0, amplitude - error) ** 2, min(1.0, amplitude + error) ** 2]
    observed = list(band(count))
    return dict(
        forecast=[float(p) for p in forecast],
        count=int(count),
        observed=[float(p) for p in observed],
        compatible=bool(max(forecast[0], observed[0]) <= min(forecast[1], observed[1])),
    )


def read(records):
    # Infer each wall solely from the two physical load addresses.
    dimers, contacts, _ = fixture(federation(0), 2)
    occupied = {
        port: ((WORDS >> i) & 1) == bit
        for i, dimer in enumerate(dimers)
        for bit, port in enumerate(dimer)
    }
    wall = sum(
        (occupied[a] != occupied[b]).astype(int) << j
        for j, (a, b) in enumerate(contacts)
    )
    timer = 2 * records["hand"][:, :, 0].astype(float) / SHOTS - 1
    times = np.unwrap(np.arctan2(timer[:, 1], timer[:, 0]))
    eps = np.sqrt(2 * np.log(2 * 50000 / 0.01) / SHOTS) + 2 / SCALE
    timer_error = float(np.arcsin(np.sqrt(2) * eps))
    if np.max(np.diff(TIMES)) + 2 * timer_error >= np.pi:
        raise ValueError("Timer phase could wrap between appointments")
    result = dict(
        hand_readings=times.tolist(),
        hand_error=timer_error,
        source_quadratures={},
        settings={},
    )
    for r in (0, 1):
        calibrated = records["source"][r, :, :, 0]
        lo, hi = band(calibrated)
        result["source_quadratures"][str(r)] = dict(
            plus_counts=calibrated.tolist(),
            mean_lower=(2 * lo - 1).tolist(),
            mean_upper=(2 * hi - 1).tolist(),
        )
    for i, coupling in enumerate(CONTRAST):
        for r in (0, 1):
            full = records["full"][i, r]
            events = {
                name: full[:, :, (wall == mask)].sum(axis=-1)
                for name, mask in (("near", 1), ("far", 2))
            }
            events["remote"] = records["local"][i, r, :, :, 1:3].sum(axis=-1)
            row = {name: compare_phases(counts) for name, counts in events.items()}
            if i == 0:
                row["scalar_early"] = scalar_comparison(
                    times, timer_error, events["near"][2, 1]
                )
            result["settings"][f"k{coupling}_r{r}"] = row
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path)
    path = ap.parse_args().output
    laws = calculate()
    if path:
        path.mkdir(exist_ok=False)
        np.savez_compressed(path / "laws.npz", **laws)
        records = draw(laws)
        np.savez_compressed(path / "records.npz", **records)
        (path / "results.json").write_text(json.dumps(read(records), indent=2) + "\n")
        return
    with np.load(HERE / "laws.npz") as stored:
        for key, cdf in laws.items():
            # Compare computed CDFs to the portable acquisition quantiles.
            np.testing.assert_allclose(
                cdf / SCALE, stored[key] / SCALE, rtol=0, atol=8e-12
            )
        records = draw(stored)
    with np.load(HERE / "records.npz") as stored:
        for key, counts in records.items():
            np.testing.assert_array_equal(counts, stored[key])
    from ....electron_field.replay import compare

    compare(
        read(records), json.loads((HERE / "results.json").read_text()), "fork phase"
    )
    print("Verified every fork detector law, integer bank, clock and phase comparison.")


if __name__ == "__main__":
    main()
