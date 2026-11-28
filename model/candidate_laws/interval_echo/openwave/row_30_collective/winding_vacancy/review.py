"""Independent arithmetic checks using the archived records and the stated local equations."""

import argparse
from hashlib import sha256
from pathlib import Path
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent


def transfer_check():
    # Independent binary phase readback. The cosine transform uses measured recoil differences.
    sizes = 42
    values = {}
    with np.load(HERE / "harmonics_counts.npz", allow_pickle=False) as raw:
        angular = 2.0 ** np.arange(-5, 16)
        for field in ("phase", "foil"):
            lines = [0.0]
            for j in range(1, sizes + 1):
                bins = raw[f"{field}_{j}"]
                z = (2 * bins[0] / 2**20 - 1) + 1j * (2 * bins[1] / 2**20 - 1)
                freq = np.angle(z[0]) / angular[0]
                for t, datum in zip(angular[1:], z[1:], strict=True):
                    innovation = np.angle(datum * np.exp(-1j * t * freq))
                    freq += innovation / t
                lines.append(freq)
            even = np.r_[lines, lines[-2:0:-1]]
            values[field] = (2 * np.fft.rfft(even).real / len(even))[3:6].tolist()
    with np.load(HERE / "transfer_counts.npz", allow_pickle=False) as raw:
        order = [8.0, 9.0]
        for j in (69, 64, 57, 49, 36, 21):
            d = np.sqrt(16.25 + 4 * np.cos(np.pi * j / 73))
            order.extend((4.5 - d, 4.5 + d))
        t = np.r_[1 / 64, np.arange(128) / 8]
        for key in ("free", "check_1024"):
            bins = raw["gateway_" + key]
            y = (bins[..., 0] - bins[..., 1]).mean(axis=-1) / 2**20
            z = y[:, 0, :] + 1j * y[:, 1, :]
            err = np.angle(z * np.exp(1j * np.outer(np.array(order) - 4.5, t)))
            visible = abs(z) > 0.101
            values["visible_sample_phase_error_" + key] = float(
                np.max(abs(err[visible]))
            )
    return values


def run():
    size = 146
    d = np.full(size, 4.5)
    d[[0, -1]] = 8.5
    h = np.resize([-0.5, -4.0], size - 1)
    e, v = eigh_tridiagonal(d, h)
    j = np.arange(1, 73)
    k = j * np.pi / 73
    split = np.sqrt(16.25 + 4 * np.cos(k))
    independent = np.sort(np.r_[8, 9, 4.5 + split, 4.5 - split])
    roots = max(abs(e - independent))
    assert roots < 1e-12
    # Check the other producer's saved finite eigenenergies and sharp-charge records.
    with np.load(HERE / "counts.npz", allow_pickle=False) as raw:
        np.testing.assert_allclose(
            raw["alternating_matched_host_levels"], e, atol=1e-12
        )
        t = np.arange(1, 113) / 4
        psi = (np.exp(-1j * np.outer(t, e)) * v[72]) @ v.T
        probs = abs(psi[:, :144]) ** 2
        probs[:, :2] += abs(psi[:, 144:]) ** 2
        shots = raw["alternating_impulse_counts"]
        max_impulse = np.max(abs(shots / shots.sum(axis=1, keepdims=True) - probs))
    # Independent time and space extraction from the archived calibration/phase records.
    with np.load(HERE / "phase_counts.npz", allow_pickle=False) as raw:

        def z(x):
            return (x[0] - x[1] + 1j * (x[2] - x[3])) / x[0, 0].sum()

        t = np.arange(1, 65) / 4
        rest = z(raw["rest_counts"])
        x = np.arange(2, 144)
        sign = (-1.0) ** (x + x // 2)
        rest_clock = np.dot(t, np.unwrap(np.angle(rest[:, 2:] @ sign))) / np.dot(t, t)
        upper_rest = 9 - rest_clock
        upper = []
        for bank, mode in zip(
            raw["line_counts"], (69, 64, 57, 49, 36, 21), strict=True
        ):
            # Blindly deflate the strongest spectral peak before seeking the second.
            y = z(bank)[:, 4:140:2]
            t = np.arange(len(y)) * 0.25
            residual = y.copy()
            freqs = []
            for _ in range(2):
                nfft = 8192
                transform = np.fft.fft(residual, nfft, axis=0)
                loc = np.argmax(np.sum(abs(transform) ** 2, axis=1))
                f = 2 * np.pi * loc / (0.25 * nfft)

                def power(w):
                    return -(np.linalg.norm(np.exp(-1j * w * t) @ residual) ** 2)

                fine = minimize_scalar(
                    power, bounds=(f - 0.02, f + 0.02), method="bounded"
                ).x
                freqs.append(fine)
                design = np.exp(1j * t[:, None] * freqs)
                weights = np.linalg.lstsq(design, y, rcond=None)[0]
                residual = y - design @ weights
            # Coordinate refinement addresses leakage from the stronger resolved component.
            for _ in range(3):
                for i in range(2):
                    other = 1 - i
                    component = (
                        y - np.exp(1j * freqs[other] * t)[:, None] * weights[other]
                    )

                    def negative_power(w):
                        return -(np.linalg.norm(np.exp(-1j * w * t) @ component) ** 2)

                    freqs[i] = minimize_scalar(
                        negative_power,
                        bounds=(freqs[i] - 0.02, freqs[i] + 0.02),
                        method="bounded",
                    ).x
                    design = np.exp(1j * t[:, None] * freqs)
                    weights = np.linalg.lstsq(design, y, rcond=None)[0]
            upper_index = np.argmin(freqs)
            amps = weights[upper_index]
            # Translation phase from an independent spatial scan of the isolated temporal line.
            cells = np.arange(2, 70)

            def space_error(k):
                design = np.column_stack((np.sin(k * cells), np.cos(k * cells)))
                return np.linalg.norm(
                    amps - design @ np.linalg.lstsq(design, amps, rcond=None)[0]
                )

            grid = np.linspace(0.0001, np.pi - 0.0001, 257)
            g = np.argmin([space_error(k) for k in grid])
            fitted = minimize_scalar(
                space_error,
                bounds=(grid[max(g - 1, 0)], grid[min(g + 1, 256)]),
                method="bounded",
            ).x
            q = np.pi - fitted
            upper.append(dict(source_j=mode, q=q, frequency=9 - freqs[upper_index]))
    for row in upper:
        q = row["q"]
        predictions = dict(
            continuum=np.sqrt(upper_rest**2 + (32 / 7) * q * q),
            parabola=upper_rest + (2 / 7) * q * q,
            lattice=np.sqrt(upper_rest**2 + 4 * (32 / 7) * np.sin(q / 2) ** 2),
            hopping=upper_rest + (4 / 7) * (1 - np.cos(q)),
            native=4.5 + np.sqrt(12.25 + 8 * np.sin(q / 2) ** 2),
        )
        row["errors"] = {
            name: float(row["frequency"] - p) for name, p in predictions.items()
        }

    def extremes(rows):
        return {
            name: max(abs(r["errors"][name]) for r in rows)
            for name in rows[0]["errors"]
        }

    # Direct sharp-source position correlations from an independent real-space propagation.
    # Read every preparation and physical output. The final coordinate aliases port zero.
    with np.load(HERE / "density_counts.npz", allow_pickle=False) as raw:
        y = (
            raw["sharp"][
                :,
                :,
                :,
            ]
            / 8192
        )
        coord = np.arange(0, size, 2) % 144 // 2
        aidx = np.arange(0, size, 2)
        phi = np.exp(
            1j
            * (2 * (73 - np.array([69, 64, 57, 49, 42, 39])) * np.pi / 73)[
                :, None, None
            ]
            * (np.arange(72)[None, None, :] - coord[None, :, None])
        )
        out = np.einsum("tia,kia->tk", y[:, :, ::2], phi).real / len(coord)
        t = np.arange(256) / 4
        full = []
        for time in t:
            amp = (v[aidx] * np.exp(-1j * time * e)) @ v.T
            pr = abs(amp[:, ::2]) ** 2
            pr[:, 0] += pr[:, -1]
            full.append(np.einsum("ia,kia->k", pr[:, :72], phi).real / len(coord))
        sharp_direct = float(np.max(abs(np.array(full) - out)))
        # The constant guide energy is outside the actual finite one-vacancy spectrum.
        nearest = float(np.min(abs(e - 4.5)))
    # True E=4.5+sqrt(m²+2 q²-q^4/6+...), m=3.5.
    m = 3.5
    quartic = dict(
        native=-1 / (12 * m) - 1 / (2 * m**3),
        host_continuum=-((32 / 7) ** 2) / (8 * 8**3),
        host_lattice=-(32 / 7) / (24 * 8) - (32 / 7) ** 2 / (8 * 8**3),
        hopping=-(4 / 7) / 24,
    )
    return dict(
        algorithm="independent finite Jacobi, deflated periodogram and spatial scan; no producer imports",
        independent_transfers=transfer_check(),
        receipt_hashes={
            name: sha256((HERE / name).read_bytes()).hexdigest()
            for name in ("counts.npz", "phase_counts.npz", "density_counts.npz")
        },
        max_analytic_root_error=roots,
        max_sharp_impulse_error=max_impulse,
        max_direct_sharp_correlation_error=sharp_direct,
        nearest_stationary_energy_to_centre=nearest,
        rest_against_stationary_clock=rest_clock,
        upper_formation_rest=upper_rest,
        upper_kinetic_curvature=4 / 7,
        host_c_squared=32 / 7,
        rest_mass_times_front_speed_squared=(7 / 4) / 4,
        material_max_group_velocity=0.5,
        upper_host_continuum_speed=np.sqrt(32 / 7),
        quartic_coefficients=quartic,
        raw_upper=upper,
        max_long_errors=extremes([r for r in upper if r["q"] < 0.45]),
        max_all_errors=extremes(upper),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    output = json.dumps(run(), sort_keys=True, indent=2) + "\n"
    path = HERE / "review.json"
    if p.parse_args().write:
        path.write_text(output)
    else:
        assert path.read_text() == output, "Independent review mismatch"
    print("Independent winding arithmetic verified.")
