"""Proper-volume stops and independent field readings for a Weyl radiation source."""

from functools import lru_cache
import json
from math import pi
from pathlib import Path

import numpy as np
from scipy import linalg, optimize

from matched_einstein_limit.frame import FRAME
from ..measure import entropy, make_board, regional_covariances, sites_in_ball


@lru_cache(maxsize=8)
def quadrature(order):
    return np.polynomial.legendre.leggauss(order)


def source_word(side, length, cutoff, floor):
    """Flip both bands on the occupied shell. The total comoving work is an eigenvalue."""
    board = make_board(side, length, cutoff, floor, (0, 0, 0))
    p1 = 2 * pi * (side * np.fft.fftfreq(side) + 0.5) / length
    p = np.meshgrid(p1, p1, p1, indexing="ij", sparse=True)
    d = [np.sin(board.spacing * pj) / board.spacing for pj in p]
    energy = np.sqrt(sum(component**2 for component in d))
    flipped = sum(component**2 for component in p) < cutoff**2
    rho = float(2 * np.sum(energy[flipped]) / length**3)
    board.source_kernels = [np.fft.ifftn((1 - 2 * floor) * flipped * component / energy) for component in d]
    board.source_kernels.append(np.zeros_like(board.source_kernels[0]))
    board.energy = (1 - 2 * floor) * rho
    return board, int(flipped.sum()), rho


def geodesic_map(arc, mu, h, speed):
    """Endpoint and angular Gram Jacobian from unit spacelike launch directions."""
    mu = np.asarray(mu)
    boost = np.sqrt(1 + speed**2)
    v = speed * mu
    w = np.broadcast_to(arc, np.broadcast_shapes(np.shape(arc), np.shape(mu))).copy()
    for _ in range(7):
        x = h * w
        sine, cosine = np.sin(x), np.cos(x)
        sinc = np.sinc(2 * x / pi)
        proper = w * ((1 + v**2) / 2 + (1 - v**2) * sinc / 2 + v * sine * np.sinc(x / pi))
        q = cosine + v * sine
        w -= (proper - arc) / q**2
    x = h * w
    sine, cosine = np.sin(x), np.cos(x)
    q = cosine + v * sine
    check = w * ((1 + v**2) / 2 + (1 - v**2) * np.sinc(2 * x / pi) / 2 + v * sine * np.sinc(x / pi))
    if np.any(q <= 0) or np.any(w < 0) or not np.allclose(check, arc, rtol=2e-12, atol=1e-15):
        raise ValueError("rod exceeds the enrolled expanding chart or failed to reach its stop")
    # Integrate partial_v[cos(hz)+v sin(hz)]², then differentiate the rod stop.
    derivative_v = w * (sine * np.sinc(x / pi) + v * (1 - np.sinc(2 * x / pi)))
    w_mu = -speed * derivative_v / q**2
    time_mu = (-sine + v * cosine) * w_mu + speed * w * np.sinc(x / pi)
    sin_theta = np.sqrt(1 - mu**2)
    transverse_mu = sin_theta * w_mu - mu * w / sin_theta
    longitudinal_mu = boost * (mu * w_mu + w)
    norm_squared = transverse_mu**2 + longitudinal_mu**2 - time_mu**2
    jacobian = q**2 * w * sin_theta * np.sqrt(norm_squared)
    return w, jacobian


def surface(radius, h, speed, order):
    mu, weights = quadrature(order)
    _, jac = geodesic_map(np.asarray(radius)[..., None], mu, h, speed)
    return 2 * pi * (jac @ weights)


def volume(radius, h, speed):
    nodes, weights = quadrature(12)
    return radius * (weights @ surface(radius * (nodes + 1) / 2, h, speed, 12)) / 2


def proper_ball(radius, h, speed):
    target = 4 * pi * radius**3 / 3
    rod = optimize.brentq(lambda r: volume(r, h, speed) / target - 1,
                          radius * 0.8, radius * 1.2, xtol=1e-15, rtol=1e-14)
    area = surface(rod, h, speed, 18)
    return float(rod), float(area)


def conformal_ball(radius, h, speed):
    """A different causal diamond in a conformal-normal chart; integrate its own units."""
    mu, weights = quadrature(18)
    radial, wr = quadrature(14)

    def conformal_factor(r):
        y0 = r * speed * mu
        d = 1 + h * y0 - h**2 * r**2 / 4
        return (d + h * y0 - h**2 * r**2 / 2) / d**2

    def proper_volume(r):
        rr = r * (radial[:, None] + 1) / 2
        return pi * r * wr @ ((rr**2 * conformal_factor(rr)**3) @ weights)

    target = 4 * pi * radius**3 / 3
    rim = optimize.brentq(lambda r: proper_volume(r) / target - 1,
                          radius * 0.8, radius * 1.2, xtol=1e-15, rtol=1e-14)
    area = 2 * pi * rim**2 * weights @ conformal_factor(rim)**2
    return float(rim), float(area)


def rest_correlation_read(board, target_radius, rod, h):
    # In the homogeneous rest frame every geodesic rim event has the same conformal time.
    # The stationary field can therefore be read on that time plane inside its measured rim.
    waist, _ = geodesic_map(rod, np.array([0.0]), h, 0)
    sites_vacuum = sites_in_ball(board.spacing, target_radius, (0, 0, 0))
    sites_source = sites_in_ball(board.spacing, float(waist[0]), (0, 0, 0))
    ref, _ = regional_covariances(board, sites_vacuum)
    ground_on_source, difference = regional_covariances(board, sites_source)
    input_spectrum = linalg.eigvalsh(ref, check_finite=False)
    output_spectrum = linalg.eigvalsh(ground_on_source + difference, check_finite=False)
    core_radius = 0.985 * target_radius
    if core_radius >= min(target_radius, waist[0]):
        raise ValueError("enrolled inner aperture extends outside a proper ball")
    common_sites = sites_in_ball(board.spacing, core_radius, (0, 0, 0))
    common_ref, common_delta = regional_covariances(board, common_sites)
    common_change = entropy(linalg.eigvalsh(common_ref + common_delta)) - entropy(linalg.eigvalsh(common_ref))
    return {"comoving_waist": float(waist[0]), "vacuum_sites": len(sites_vacuum), "source_sites": len(sites_source),
            "same_sites": bool(np.array_equal(sites_source, sites_vacuum)),
            "all_reconstructed_sites_entropy_change": float(entropy(output_spectrum) - entropy(input_spectrum)),
            "common_core_radius": core_radius, "common_core_sites": len(common_sites),
            "registered_core_entropy_change": float(common_change)}


def clock_record(rate, separation=0.125, conformal_spacing=0.02):
    """Two null signals received by a comoving clock at a fixed optical address."""
    def proper(t):
        return t + rate * t**2 / 2

    times = np.array([-0.5, 0.5]) * conformal_spacing
    sent, received = proper(times), proper(times + separation)
    ratio = (received[1] - received[0]) / (sent[1] - sent[0])
    return {"optical_separation": separation, "sender_times": sent.tolist(), "receiver_times": received.tolist(),
            "clock_ratio": float(ratio), "inferred_expansion_rate": float((ratio - 1) / separation)}


def audit():
    side, length, cutoff, floor, stiffness = 48, 8.0, 1.8, 0.0005, 3.0
    _, calibration_count, calibration_density = source_word(side, length, 1.25, floor)
    calibration_clock = clock_record(np.sqrt(2 * calibration_density / stiffness))
    coupling = 3 * calibration_clock["inferred_expansion_rate"]**2 / calibration_density
    board, pairs, density = source_word(side, length, cutoff, floor)
    h = np.sqrt(2 * density / stiffness)
    geometry = []
    for radius in (0.12, 0.06, 0.03, 0.015, 0.0075):
        for frame in FRAME:
            speed = np.linalg.norm(np.asarray(frame[1:], dtype=float))
            rod, area = proper_ball(radius, h, speed)
            chart_radius, chart_area = conformal_ball(radius, h, speed)
            factor = 4 * pi * radius**4 / 15
            stiffness_reading = (4 * pi * radius**2 - area) / factor
            geometry.append({
                "radius": radius, "frame": list(map(float, frame)), "proper_radius": rod,
                "volume_at_fixed_rod_radius": float(volume(radius, h, speed)),
                "proper_area": area, "diamond_coordinate_radius": chart_radius, "diamond_area": chart_area,
                "area_coefficient": stiffness_reading, "stress_forecast": float(coupling * density * (1 + 4 * speed**2 / 3)),
                "diamond_minus_true_area_over_R4": (chart_area - area) / radius**4,
            })
    bulk = []
    for radius in (0.46, 0.55):
        rod, area = proper_ball(radius, h, 0)
        reading = rest_correlation_read(board, radius, rod, h)
        reading.update({"target_radius": radius, "proper_boundary_area": area,
                        "weighted_entropy_forecast": 8 * pi**2 * radius**4 * density / 15,
                        "registered_entropy_plus_area_change": reading["registered_core_entropy_change"] + (2 * pi / coupling) * (area - 4 * pi * radius**2)})
        bulk.append(reading)
    control_radius, control_speed = 0.015, np.sqrt(8)
    controls = {}
    for name, rate in [("declared", h), ("doubled_stiffness", h / np.sqrt(2)), ("flat_clock_and_rods", 0.)]:
        _, area = proper_ball(control_radius, rate, control_speed)
        controls[name] = (4 * pi * control_radius**2 - area) / (4 * pi * control_radius**4 / 15)
    return {"side": side, "length": length, "cutoff": cutoff, "readout_floor": floor, "flipped_band_pairs": pairs,
            "comoving_energy_density": density, "kinetic_stiffness": stiffness, "ordinary_coupling": coupling,
            "radiation_calibration": {"cutoff": 1.25, "flipped_band_pairs": calibration_count,
                                      "energy_density": calibration_density, "clock_record": calibration_clock},
            "held_clock": clock_record(h),
            "conformal_expansion_rate": h, "proper_balls": geometry, "finite_rest_entropy": bulk,
            "fixed_occupation_controls": {"radius": control_radius, "sinh_rapidity": control_speed, "area_coefficients": controls}}


if __name__ == "__main__":
    Path(__file__).with_name("results.json").write_text(json.dumps(audit(), indent=2) + "\n")
