"""Localized charged source prepared once and freely propagated by the measured hopping."""

import json
from math import pi
from pathlib import Path

import numpy as np
from scipy import linalg

from .measure import PAULI, entropy, make_board, regional_covariances, sites_in_ball


def packet_modes(board, launch):
    side, spacing = board.side, board.spacing
    momenta = 2 * pi * (side * np.fft.fftfreq(side) + 0.5) / (side * spacing)
    p = np.stack(np.meshgrid(momenta, momenta, momenta, indexing="ij"), axis=-1)
    d = np.sin(spacing * p) / spacing
    energy = np.linalg.norm(d, axis=-1)
    project_up = np.array([1, 0]) + np.einsum("...j,jb->...b", d / energy[..., None], PAULI[:, :, 0])
    dist_squared = np.sum((p - [0.4, 0.1, 2.5]) ** 2, axis=-1) / 1.9**2
    envelope = np.zeros_like(dist_squared)
    support = dist_squared < 1
    envelope[support] = np.exp(-1 / (1 - dist_squared[support]))
    # Physical grid indices are centered halfway through the antiperiodic torus.
    phase = np.exp(-1j * (p @ (np.array(launch) + side * spacing / 2)))
    profile = (envelope * phase)[..., None] * project_up
    return profile / np.linalg.norm(profile), energy


def audit():
    side, length, floor, fraction = 32, 8.0, 0.002, 0.65
    board = make_board(side, length, 2.5, floor)
    launch, times = (0.0, 0.0, -0.45), [0.0, 1.0]
    profile, energy = packet_modes(board, launch)
    twist = np.exp(1j * pi * np.indices((side,) * 3).sum(axis=0) / side)
    weight = (1 - 2 * floor) * fraction
    rows = []
    balls = [(0.64, (0., 0., -0.45)), (0.64, (0., 0., 0.45)), (0.80, (0., 0., 0.))]
    for radius, center in balls:
        sites = sites_in_ball(board.spacing, radius, center)
        cov, _ = regional_covariances(board, sites)
        eigen, basis = linalg.eigh(cov, check_finite=False)
        odds = np.log1p(-eigen) - np.log(eigen)
        index = tuple((sites + side // 2).T)
        ball_weight = (radius**2 - np.sum((board.spacing * sites - center) ** 2, axis=1)) / (2 * radius)
        for time in times:
            evolved = np.exp(-1j * time * energy)[..., None] * profile

            def observe(spectrum):
                field = np.fft.ifftn(spectrum, axes=(0, 1, 2), norm="ortho") * twist[..., None]
                return field[index]

            value = observe(evolved)
            acted = observe(energy[..., None] * evolved)
            vector = value.ravel()
            modal = weight * np.dot(odds, np.abs(basis.conj().T @ vector) ** 2)
            perturbed = cov + weight * np.outer(vector, vector.conj())
            shift = entropy(linalg.eigvalsh(perturbed, check_finite=False)) - entropy(eigen)
            proper_energy = 2 * pi * weight * np.einsum("i,ij,ij->", ball_weight, value.conj(), acted).real
            local_probability = weight * np.vdot(vector, vector).real
            rows.append({
                "radius": radius, "center": list(center), "time": time, "occupation_excess": float(local_probability),
                "weighted_hopping_energy": float(proper_energy), "modular_change": float(modal),
                "entropy_change": float(shift), "relative_entropy": float(modal - shift),
                "relative_entropy_upper_bound": float(local_probability**2 / (2 * floor * (1 - floor))),
            })
    return {"grid_side": side, "box_length": length, "floor": floor, "particle_fraction": fraction, "launch": launch, "readings": rows}


if __name__ == "__main__":
    Path(__file__).with_name("packet.json").write_text(json.dumps(audit(), indent=2) + "\n")
