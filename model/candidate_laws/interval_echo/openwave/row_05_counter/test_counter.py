"""Reconstruct contact energies, physical arm symmetries and long-time readouts."""

import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from .experiment import board_h

HERE = Path(__file__).parent


def test_two_charge_contact_from_full_native_seam():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    block = (exchange + 8 * contrast)[: law.N_SRC, : law.N_SRC]

    def cost(a, b):
        i = law.SRC.index((a, b))
        return block[i, i]

    for q in (-1, 1):
        full, left, right, links = board_h(q, 8, "all")
        cut, *_ = board_h(q, 8, "cut")
        connected = cost(0, 1 + q) - cost(0, 1) - cost(1, 1 + q) + cost(1, 1)
        assert np.isclose(connected, -17 if q == -1 else 15.5)
        site = np.array(
            [tuple(sorted((i, j))) in links for i in left for j in right], int
        )
        shift = full - cut - np.diag(site * connected)
        np.testing.assert_allclose(shift, np.eye(25) * shift[0, 0], atol=1.0e-12)


def test_single_link_reflection_exchanges_sources_and_commutes():
    for probe in (-1, 1):
        h, left, _, links = board_h(probe, 8, "one")
        b = next(x for x in links[0] if x in left)
        axis = left.index(b)
        refl = (2 * axis - np.arange(5)) % 5
        mirror = np.kron(np.eye(5)[refl], np.eye(5))
        np.testing.assert_allclose(h, mirror.T @ h @ mirror)
        modes = np.exp(2j * np.pi * np.arange(5) / 5) / np.sqrt(5)
        assert np.isclose(abs(np.vdot(modes.conj(), modes[refl])), 1)


def test_independent_full_unitary_power_and_archive_counts():
    # Direct coordinate-space pulses, without any momentum-block spectral projection.
    directory = json.loads((HERE / "beat.json").read_text())
    for study in directory["studies"]:
        probe = study["probe"]
        h, left, right, links = board_h(probe, 8, "all")
        zero, *_ = board_h(probe, 8, "cut")
        v = np.diag(h - zero)
        neigh = (np.arange(5)[:, None] - np.arange(5)[None, :]) % 5
        rim = np.eye(5) - ((neigh == 1) | (neigh == 4)) / 2
        cycle = np.diag(np.exp(-1j * study["duration"] * v)) @ np.kron(
            expm(-0.25j * rim), expm(-0.125j * rim)
        )
        dft = np.exp(2j * np.pi * np.outer(np.arange(5), np.arange(5)) / 5) / np.sqrt(5)
        for record in (study["held"][0], study["held"][-1]):
            u = np.linalg.matrix_power(cycle, record["periods"])
            for sign_read in record["arms"]:
                source = sign_read["sign"]
                psi = (
                    u @ np.kron(dft[:, source], (dft[:, 2] + dft[:, -2]) / np.sqrt(2))
                ).reshape(5, 5)
                spectral = dft.conj().T @ psi @ dft.conj()
                contrast = 2 * np.vdot(spectral[:, -2], spectral[:, 2])
                raw = np.array(sign_read["counts"])
                xy = (raw[:, :, 0] - raw[:, :, 1]).sum(axis=1) / directory["shots"]
                np.testing.assert_allclose(
                    [contrast.real, contrast.imag], xy, atol=0.003
                )
                assert raw.sum(axis=(1, 2)).tolist() == [directory["shots"]] * 2
        counts = np.array(study["calibration_counts"])
        assert np.all(counts.sum(axis=2) == directory["shots"])
        qp = (
            counts[:, :, 0] - counts[:, :, 1] + counts[:, :, 2] - counts[:, :, 3]
        ) / directory["shots"]
        freq = np.polyfit(
            directory["calibration_periods"],
            np.unwrap(np.angle(qp[:, 0] + 1j * qp[:, 1])),
            1,
        )[0]
        assert abs(freq - study["measured_frequency"]) < 1.0e-11 and freq > 0
        assert (
            max(abs(a["phase_error"]) for t in study["held"] for a in t["arms"]) < 0.018
        )


def test_charge_even_precession_and_unprotected_control():
    d = json.loads((HERE / "dispersive.json").read_text())
    for s in d["studies"]:
        rows = [a for t in s["times"] for a in t["signs"]]
        if s["active"]:
            assert max(abs(a["phase_error"]) for a in rows) > 0.1
            assert s["forecast_per_cycle"] < 0
        else:
            assert max(abs(a["phase_error"]) for a in rows) < 0.002
