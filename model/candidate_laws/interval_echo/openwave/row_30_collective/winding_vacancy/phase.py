"""Single-vacancy interferometry against a stationary state of the same unguarded guide."""

import argparse
import hashlib
import io
import json

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from .measure import HERE, N, generator, release

PROTOCOL = "ed16b145cf2679923ff5c0f3d7c8571f0dd4ae1077c176755a3a671716509527"
CELLS = 73
TIMES = np.arange(256) / 4
TRAIN = np.arange(5, 36)
HELD = np.arange(36, 68)


def clocks(energy, basis):
    x = np.arange(len(energy))
    r = (-1.0) ** x / np.sqrt(len(x))
    b = (-1.0) ** (x // 2) / np.sqrt(len(x))
    for state, line in [(r, 9), (b, 8)]:
        np.testing.assert_allclose(
            basis @ (energy * (basis.T @ state)), line * state, atol=2e-13
        )
    return r, b


def read(native, ref, signal, times, shots, rng):
    np.testing.assert_allclose(np.vdot(ref, signal), 0, atol=2e-14)
    data = []
    for phase in (1, -1, -1j, 1j):
        initial = (ref + phase * signal) / np.sqrt(2)
        p = release(*native, initial, times)
        data.append(np.array([rng.multinomial(shots, row / row.sum()) for row in p]))
    return np.array(data)


def quadratures(records):
    n = records[0, 0].sum()
    return ((records[0] - records[1]) + 1j * (records[2] - records[3])) / n


def slope(time, z):
    return time @ np.unwrap(np.angle(z)) / (time @ time)


def calibrate(raw, native, ref, b):
    rng = np.random.default_rng(302144)
    explicit = np.zeros(len(ref))
    explicit[[72, 74]] = np.array([1, -1]) / np.sqrt(2)
    times = np.arange(1, 9) / 32
    raw["centre_counts"] = read(native, ref, explicit, times, 2**24, rng)
    y = quadratures(raw["centre_counts"])
    centre = slope(times, y[:, 72] - y[:, 74])
    times = np.arange(1, 65) / 4
    raw["rest_counts"] = read(native, ref, b, times, 2**24, rng)
    y = quadratures(raw["rest_counts"])
    rest_line = slope(times, y[:, 2:] @ np.sign(ref[2:N] * b[2:N]))
    return dict(
        clock_minus_mean_work=centre,
        clock_minus_stationary_rest_line=rest_line,
        centred_rest=centre - rest_line,
    )


def marked_source(j):
    # All physical written contacts receive the same two-address, equal-phase launch.
    out = np.zeros(2 * CELLS)
    amplitude = np.sin(np.pi * j * np.arange(1, CELLS) / CELLS)
    out[1:-1] = np.repeat(amplitude, 2)
    return out / np.linalg.norm(out)


def spectral_lines(times, signals):
    u = np.linalg.svd(signals, full_matrices=False)[0][:, :2]
    step = np.linalg.lstsq(u[:-1], u[1:], rcond=None)[0]
    freq = np.sort(np.angle(np.linalg.eigvals(step)) % (2 * np.pi)) / (
        times[1] - times[0]
    )

    def residual(f):
        h = np.exp(1j * np.outer(times, f))
        error = h @ np.linalg.lstsq(h, signals, rcond=None)[0] - signals
        return np.r_[error.real.ravel(), error.imag.ravel()]

    f = least_squares(residual, freq, bounds=(0, 12), ftol=1e-12, gtol=1e-13).x
    h = np.exp(1j * np.outer(times, f))
    return f, np.linalg.lstsq(h, signals, rcond=None)[0]


def spatial_basis(cells, k):
    return np.column_stack([np.sin(k * cells), np.cos(k * cells)])


def spatial_line(amplitude):
    def fit(k):
        h = spatial_basis(TRAIN, k)
        coeff = np.linalg.lstsq(h, amplitude, rcond=None)[0]
        return np.linalg.norm(h @ coeff - amplitude)

    grid = np.linspace(1e-5, np.pi - 1e-5, 513)
    i = np.argmin([fit(k) for k in grid])
    lo, hi = grid[max(0, i - 1)], grid[min(i + 1, len(grid) - 1)]
    k = minimize_scalar(
        fit, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-13)
    ).x
    coeff = np.linalg.lstsq(spatial_basis(TRAIN, k), amplitude, rcond=None)[0]
    return k, spatial_basis(HELD, k) @ coeff


def decode(records, cal):
    y = quadratures(records)
    f, waves = spectral_lines(TIMES[::2], y[::2, 2 * TRAIN])
    predicted, lines = [], []
    for freq, amplitude in zip(f, waves, strict=True):
        k, profile = spatial_line(amplitude)
        predicted.append(profile)
        q = np.pi - k
        centre = cal["clock_minus_mean_work"]
        w = centre - freq
        m = cal["centred_rest"]
        lines.append(
            dict(
                clock_minus_line=freq,
                host_energy_from_clock=9 - freq,
                measured_k=k,
                measured_q=q,
                signed_centred_frequency=w,
                lattice_residual=abs(w) - np.sqrt(m**2 + 8 * np.sin(q / 2) ** 2),
                continuum_residual=abs(w) - np.sqrt(m**2 + 2 * q * q),
            )
        )
    forecast = np.exp(1j * np.outer(TIMES[1::2], f)) @ np.array(predicted)
    error = np.sqrt(np.mean(abs(forecast - y[1::2, 2 * HELD]) ** 2))
    return dict(
        held_spatial_rms=error, passed_spatial=bool(error < 0.00045), lines=lines
    )


def run():
    assert hashlib.sha256((HERE / "PHASE.md").read_bytes()).hexdigest() == PROTOCOL
    native = generator("01")
    ref, b = clocks(*native[1:])
    raw = {}
    cal = calibrate(raw, native, ref, b)
    rng = np.random.default_rng(303144)
    scans, rows = [], []
    for j in (69, 64, 57, 49, 36, 21):
        scan = read(native, ref, marked_source(j), TIMES, 2**19, rng)
        rows.append(dict(source_mode=j, **decode(scan, cal)))
        scans.append(scan)
    raw["line_counts"] = np.array(scans)
    lines = [row for scan in rows for row in scan["lines"]]
    q = np.array([line["measured_q"] for line in lines])
    w = np.abs([line["signed_centred_frequency"] for line in lines])
    h = np.column_stack((np.ones(len(q)), 1 - np.cos(q)))
    rival = np.linalg.lstsq(h, w, rcond=None)[0]
    residual = h @ rival - w
    max_lattice = max(abs(line["lattice_residual"]) for line in lines)
    max_continuum = max(
        abs(line["continuum_residual"]) for line in lines if line["measured_q"] < 0.45
    )
    archive = io.BytesIO()
    np.savez_compressed(archive, **raw)
    data = archive.getvalue()
    out = dict(
        protocol_sha256=PROTOCOL,
        counts_sha256=hashlib.sha256(data).hexdigest(),
        calibration=cal,
        rows=rows,
        fitted_hopping=rival.tolist(),
        max_hopping_error=max(abs(residual)),
        max_lattice_error=max_lattice,
        max_long_wave_continuum_error=max_continuum,
        passed=bool(
            max_lattice < 0.008
            and max_continuum < 0.008
            and all(row["passed_spatial"] for row in rows)
        ),
    )
    return json.dumps(out, indent=2, sort_keys=True) + "\n", data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    summary, data = run()
    for name, expected in (
        ("phase.json", summary.encode()),
        ("phase_counts.npz", data),
    ):
        path = HERE / name
        if args.write:
            path.write_bytes(expected)
        else:
            assert expected == path.read_bytes(), f"Replay mismatch for {name}"
    print("Written" if args.write else "Verified", "same-material phase records")


if __name__ == "__main__":
    main()
