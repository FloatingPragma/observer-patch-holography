"""Apply the same absolute gates and signed receiving contacts at two fixed repulsions."""

import numpy as np
from scipy.linalg import eigh

from . import ring


def eigensystem(counts, repulsion):
    h = ring.hamiltonian(4, *counts).toarray() / 16
    double_occupations = [
        sum(((word >> (2 * x)) & 3) == 3 for x in range(4))
        for word in ring.sector(4, *counts)[0]
    ]
    h += (repulsion - 1) * np.diag(double_occupations)
    return eigh(h)


def run(ring_record):
    neutral = (2, 2)
    neutral_cut = (
        ring_record["half_filled_ground"] + ring_record["half_filled_ground_gap"] / 2
    )
    neel = sum(1 << (2 * x + x % 2) for x in range(4))
    index = ring.sector(4, *neutral)[1][neel]
    report = []
    for repulsion in (1, 0.5):
        energies, frame = eigensystem(neutral, repulsion)
        ground_gate = frame[:, energies < neutral_cut]
        prepared = ground_gate @ ground_gate[index, :].conj()
        banks = []
        for original in ring_record["sign_banks"]:
            sign = original["sign"]
            window = original["two_energy_window"]
            counts = (3, 2) if sign > 0 else (2, 1)
            cut = (
                window["upper_gate_midpoint_in_corrected_units"]
                + ring_record["half_filled_ground"]
                + sign / 2
            )
            energies, frame = eigensystem(counts, repulsion)
            mask = energies < cut
            charged = frame[:, mask]

            def field(x):
                return ring.field(4, *neutral, x, int(sign < 0), sign > 0)

            injected = charged.T.conj() @ (field(0) @ prepared)
            sink = ground_gate.T.conj() @ (field(2).T @ charged)
            probabilities = []
            for t in (0, 16, window["half_beat_wait"]):
                drain = sink @ (np.exp(-1j * t * energies[mask]) * injected)
                probabilities.append(float(np.vdot(drain, drain).real / 4))
            if repulsion == 1:
                np.testing.assert_allclose(
                    probabilities,
                    [
                        window["three_waits"][label]["joint_spin_return"][2]
                        for label in ("parked", "wait_16", "half_beat")
                    ],
                    rtol=1e-10,
                    atol=1e-14,
                )
            banks.append(
                {
                    "sign": sign,
                    "charged_upper_gate_raw": float(cut),
                    "charged_gate_rank_at_fixed_spin_counts": int(mask.sum()),
                    "positive_spin_and_neutral_flag_per_original": probabilities,
                }
            )
        report.append(
            {
                "repulsion": repulsion,
                "neutral_gate_upper_raw": float(neutral_cut),
                "neutral_gate_rank_at_fixed_spin_counts": ground_gate.shape[1],
                "sign_banks": banks,
            }
        )
    return report
