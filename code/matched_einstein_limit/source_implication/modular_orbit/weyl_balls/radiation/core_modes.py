"""Eight physical interior mode channels on both homogeneous preparations and ten frames."""

import json
from math import pi
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.sparse.linalg import expm_multiply

from ..boosts import covariance_on_columns, frame_generator, hopping, smooth_step
from ..measure import entropy
from .measure import geodesic_map, proper_ball, source_word


def optical_clearance(rod, h, speed):
    """Sample the distances to the true rim, allowing for its optical retiming."""
    mu = np.cos(np.linspace(1e-6, pi - 1e-6, 241))
    w, _ = geodesic_map(rod, mu, h, speed)
    boost = np.sqrt(1 + speed**2)
    time = -h * w**2 * np.sinc(h * w / (2 * pi))**2 / 2 + speed * mu * w * np.sinc(h * w / pi)
    transverse = np.sqrt(1 - mu**2) * w
    longitudinal = boost * mu * w
    norm = transverse**2 + longitudinal**2 - time**2
    denominator = 1 - h * time - h**2 * norm / 4
    norm_time = (time + h * norm / 2) / denominator
    moved_time = boost * norm_time - speed * longitudinal / denominator
    moved_z = boost * longitudinal / denominator - speed * norm_time
    return float(np.min(np.sqrt(moved_z**2 + (transverse / denominator)**2) - np.abs(moved_time)))


def sample():
    side, length, cutoff, floor, radius, rapidity = 24, 8., 1.8, 0.0005, 0.51, 0.32
    board, pairs, rho = source_word(side, length, cutoff, floor)
    h = np.sqrt(2 * rho / 3)
    hamiltonian = hopping(side, board.spacing)
    positions = board.spacing * (np.indices((side,) * 3).reshape(3, -1).T - side // 2)
    squared = np.sum(positions**2, axis=1)
    core_radius = 0.86 * radius
    envelope = np.maximum(0, 1 - squared / core_radius**2) ** 3
    profiles = np.column_stack((envelope, envelope[:, None] * positions / radius))
    basis, normalizer = np.linalg.qr(profiles, mode="reduced")
    if np.min(np.abs(np.diag(normalizer))) < 1e-8:
        raise ValueError("insufficient resolution for the four aperture profiles")
    columns = np.kron(basis, np.eye(2))
    taper = np.prod(smooth_step((length / 2 - np.abs(positions)) / (length / 4)), axis=1)
    radial_position = np.repeat(squared * taper, 2)
    # On the central optical plane the temporal special conformal vector has
    # time component -|x|². Its adjoint mode inclusion has the indicated sign.
    conformal_generator = (hamiltonian.multiply(radial_position[:, None]) +
                           hamiltonian.multiply(radial_position[None, :])).tocsr() / 2
    rows = []
    axes = np.eye(3)
    directions = [np.zeros(3), *axes, *-axes, *((axes[i] + axes[j]) / np.sqrt(2)
                                                for i, j in [(0, 1), (1, 2), (2, 0)])]
    for direction in directions:
        tilt = rapidity if np.any(direction) else 0.0
        base = columns
        if tilt:
            boost = frame_generator(board, hamiltonian, direction)
            base = expm_multiply(-1j * tilt * boost, columns, traceA=0)
        moved = expm_multiply(-0.5j * h * conformal_generator, base, traceA=0)
        vacuum, _ = covariance_on_columns(board, base)
        on_curved_chart, difference = covariance_on_columns(board, moved)
        entropy_change = entropy(linalg.eigvalsh(on_curved_chart + difference)) - entropy(linalg.eigvalsh(vacuum))
        vacuum_shift = entropy(linalg.eigvalsh(on_curved_chart)) - entropy(linalg.eigvalsh(vacuum))
        speed = np.sinh(tilt)
        rod, area = proper_ball(radius, h, speed)
        rows.append({
            "direction": direction.tolist(), "rapidity": tilt, "true_proper_area": area,
            "sampled_optical_interior_clearance": optical_clearance(rod, h, speed),
            "registered_entropy_change": float(entropy_change), "vacuum_transport_change_at_finite_cutoff": float(vacuum_shift),
            "weighted_stress_forecast": float(8 * pi**2 * radius**4 * rho * (1 + 4 * speed**2 / 3) / 15),
            "registered_entropy_plus_area_change": float(entropy_change + pi * (area - 4 * pi * radius**2)),
        })
    return {"grid_side": side, "length": length, "momentum_cutoff": cutoff, "flipped_band_pairs": pairs,
            "kinetic_stiffness": 3., "rho": rho, "readout_floor": floor, "target_radius": radius,
            "packet_support_radius": core_radius, "mode_count": columns.shape[1], "frames": rows}


if __name__ == "__main__":
    Path(__file__).with_name("core_modes.json").write_text(json.dumps(sample(), indent=2) + "\n")
