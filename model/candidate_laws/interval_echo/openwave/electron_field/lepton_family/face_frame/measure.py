"""Read the exact port family band relative to an oriented material face."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.executable_model import cell, rotations
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet.measure import (
    face_edges,
    faces,
)
from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.family_projector.measure import (
    theorem_projectors,
)
from .field import rest_and_motion

HERE = Path(__file__).resolve().parent
OMEGA = np.exp(2j * np.pi / 3)
SHOTS = 2**22
SEED = 714026


def permutation(order):
    return np.eye(len(order))[:, list(order)]


def face_turn(face):
    a, b, c = face
    group = [p for p in rotations() if p[a] == b and p[b] == c and p[c] == a]
    (turn,) = group
    return turn


def character_frames(face, band):
    turn = permutation(face_turn(face))
    out = np.stack(
        [
            sum(
                OMEGA ** (-j * n) * np.linalg.matrix_power(turn, n) @ band
                for n in range(3)
            )
            / 3
            for j in range(3)
        ]
    )
    np.testing.assert_allclose(out.sum(axis=0), band, atol=1e-14)
    for p in out:
        np.testing.assert_allclose(p @ p, p, atol=1e-14)
        np.testing.assert_allclose(np.trace(p), 1, atol=1e-14)
    return out


def blind_spectroscopy(masses, weights, rng):
    """Two common clock quadratures, equal time spacing; fit free frequencies from counts."""
    step = 0.5
    times = step * np.arange(288)
    phases = np.exp(-1j * np.outer(times, masses)) @ weights
    prob = (1 + np.stack([phases.real, phases.imag], axis=1)) / 2
    counts = rng.binomial(SHOTS, np.clip(prob, 0, 1))
    data = (2 * counts / SHOTS - 1) @ np.array([1, 1j])
    train = data[:256]
    hankel = np.lib.stride_tricks.sliding_window_view(train, 32).T
    basis, sv, _ = np.linalg.svd(hankel, full_matrices=False)
    fits = []
    for rank in (1, 2, 3):
        u = basis[:, :rank]
        turn = np.linalg.lstsq(u[:-1], u[1:], rcond=None)[0]
        poles = np.linalg.eigvals(turn)
        frequencies = np.sort(-np.angle(poles) / step)
        design = np.exp(-1j * np.outer(times, frequencies))
        residues = np.linalg.lstsq(design[:256], train, rcond=None)[0]
        errors = abs(design[256:] @ residues - data[256:])
        fits.append(
            dict(
                number_of_lines=rank,
                energies=frequencies.tolist(),
                residue_magnitudes=abs(residues).tolist(),
                held_rms=float(np.sqrt(np.mean(errors**2))),
            )
        )
    # For both sources, cover both quadratures at indices 0..64 at once.
    # A single nonzero lower bound excludes every stationary two-tone series.
    delta = 2 * np.sqrt(np.log(2 * (2 * 2 * 65) / 0.01) / SHOTS)
    det_error = 3 * np.sqrt(3) * ((1 + delta) ** 3 - 1)
    minors = []
    for lag in range(1, 17):
        hankel3 = data[lag * (np.arange(3)[:, None] + np.arange(3))]
        estimated = float(abs(np.linalg.det(hankel3)))
        minors.append(
            dict(
                lag=lag,
                determinant_magnitude=estimated,
                simultaneous_lower=max(0.0, estimated - det_error),
            )
        )
    return dict(
        step=step,
        hankel_minors=minors,
        determinant_error_99=float(det_error),
        training_times=256,
        held_times=32,
        counts_xy=counts.tolist(),
        leading_hankel_singular_values=sv[:5].tolist(),
        fits=fits,
    )


def family_comparison(band, rng):
    """The three eigenvalues are tied to one covariant, retained *declared* mass term."""
    theta = 2 / 9
    roots = 1 + np.sqrt(2) * np.cos(theta + np.arange(3) * 2 * np.pi / 3)
    assert np.all(roots > 0)
    masses = roots**2
    frames = np.stack([character_frames(f, band) for f in faces()])
    mass = np.einsum("k,fkij->fij", masses, frames)
    for p in rotations():
        U = permutation(p)
        for i, face in enumerate(faces()):
            image = set(p[j] for j in face)
            j = next(j for j, f in enumerate(faces()) if set(f) == image)
            np.testing.assert_allclose(U @ mass[i] @ U.T, mass[j], atol=1e-13)
    np.testing.assert_allclose(
        frames.mean(axis=0), np.repeat(band[None] / 3, 3, axis=0), atol=2e-14
    )
    np.testing.assert_allclose(mass.mean(axis=0), 2 * band, atol=2e-14)
    # A second, independently drawn face after reading a line gives a flat
    # transition matrix. Retaining the same face gives perfect repeatability.
    old = frames[0]
    same = np.einsum("kij,lji->kl", old, old).real
    reset = np.einsum("kij,lji->kl", old, frames.mean(axis=0)).real
    np.testing.assert_allclose(same, np.eye(3), atol=1e-14)
    np.testing.assert_allclose(reset, np.full((3, 3), 1 / 3), atol=1e-14)

    # A path reference supplies balanced, coherent evolution-versus-idle
    # clock control; the readout uses no mass-dependent detector settings.
    times = np.arange(8) / 2
    source_ports = (0, 5)
    banks = []
    for src in source_ports:
        x = band[:, src] / np.linalg.norm(band[:, src])
        weights = np.einsum("i,fkij,j->k", x, frames, x).real / len(frames)
        np.testing.assert_allclose(weights, np.full(3, 1 / 3), atol=2e-14)
        series = []
        counts = []
        for t in times:
            phase = np.dot(weights, np.exp(-1j * t * masses))
            q = (1 + np.array([phase.real, phase.imag])) / 2
            draw = rng.binomial(SHOTS, np.clip(q, 0, 1))
            counts.append(draw.tolist())
            series.append(np.dot([1, 1j], 2 * draw / SHOTS - 1))
        # A three-tone recursion is fitted to six complex clock readings. Two
        # extra times are held. The exact moment rank is proved in the README.
        h0 = np.array([series[i : i + 3] for i in range(3)])
        h1 = np.array([series[i + 1 : i + 4] for i in range(3)])
        z = np.linalg.eigvals(np.linalg.solve(h0, h1))
        ehat = sorted((-np.angle(z) / 0.5).tolist())
        recurrence = np.linalg.solve(h0, np.array(series[3:6]))
        held_fit = list(series[:6])
        for _ in range(2):
            held_fit.append(np.dot(held_fit[-3:], recurrence))
        banks.append(
            dict(
                source_port=src,
                plus_counts_xy=counts,
                recovered_frequencies=ehat,
                held_complex_errors=[
                    float(abs(a - b))
                    for a, b in zip(held_fit[6:], series[6:], strict=True)
                ],
            )
        )
        banks[-1]["long_clock"] = blind_spectroscopy(masses, weights, rng)
    return dict(
        joint_field=rest_and_motion(frames[0], masses),
        theta=theta,
        unit="squared register amplitude with unit mass coefficient",
        energies=masses.tolist(),
        ordered_ratios=(np.sort(masses) / min(masses)).tolist(),
        koide=float(sum(masses) / sum(roots) ** 2),
        frame_count=len(frames),
        stabilizer_size=60 // len(frames),
        twirled_mass=2.0,
        same_frame_repeat=same.round(12).tolist(),
        independent_frame_repeat=reset.round(12).tolist(),
        shots_per_quadrature=SHOTS,
        sample_times=times.tolist(),
        banks=banks,
    )


def native(face, sign, kappa=8):
    """One unit of charge; three live even-parity direction heads and 27 blanks."""
    marked = face_edges(face)
    raw = blocks(kappa)
    h = np.zeros((96, 96), dtype=int)  # eight times the physical relative energy
    for u, v in cell().seams:
        hit = [j for j, (a, b) in enumerate(marked) if {a, b} == {u, v}]
        term = raw[(1, 1) if hit else "blank"]
        for head in range(8):
            coin = (
                (((head >> hit[0]) & 1) ^ int(marked[hit[0]][0] != u) ^ int(sign > 0))
                if hit
                else 0
            )
            filled = 9 * coin + LAW.SRC.index((1, 1))
            for site, other in ((u, v), (v, u)):
                pair = (1 + sign, 1) if site == u else (1, 1 + sign)
                i = 9 * coin + LAW.SRC.index(pair)
                col = 8 * site + head
                h[col, col] += term[i, i] - term[filled, filled]
                for k in np.flatnonzero(term[:, i]):
                    if k == i:
                        continue
                    assert LAW.SRC[k % 9] == pair[::-1]
                    new = head
                    if hit:
                        assert k // 9 == 1 - coin
                        new ^= 1 << hit[0]
                    h[8 * other + new, col] += term[k, i]
    np.testing.assert_array_equal(h, h.T)
    turn = permutation(face_turn(face))
    hr = permutation([((n << 1) & 7) | (n >> 2) for n in range(8)])
    rotation = np.kron(turn, hr)
    np.testing.assert_array_equal(rotation @ h @ rotation.T, h)
    return h, rotation


def native_assays(band, rng):
    face = faces()[0]
    eigprojectors = character_frames(face, band)
    states = []
    for p in eigprojectors:
        column = p[:, face[0]]
        states.append(column / np.linalg.norm(column))
    forms = []
    for kappa in (8, 4):
        h, turn = native(face, -1, kappa)
        np.testing.assert_array_equal(h, native(face, 1, kappa)[0])
        diag = h.diagonal()
        radius = np.sum(abs(h), axis=1) - abs(diag)
        energies, vectors = np.linalg.eigh(h / 8)
        # Collapse precisely repeated levels before comparing spectral measures.
        cuts = np.r_[0, np.where(np.diff(energies) > 1e-8)[0] + 1, len(energies)]
        rows = []
        for name, head in [
            ("ready", np.eye(8)[0]),
            ("plus_i", np.array([1j ** n.bit_count() for n in range(8)]) / np.sqrt(8)),
            (
                "minus_i",
                np.array([(-1j) ** n.bit_count() for n in range(8)]) / np.sqrt(8),
            ),
        ]:
            cols = np.stack([np.kron(v, head) for v in states], axis=1)
            np.testing.assert_allclose(
                turn @ cols, cols @ np.diag(OMEGA ** np.arange(3)), atol=1e-13
            )
            weights = abs(vectors.T @ cols) ** 2
            groups = np.array(
                [
                    weights[i:j].sum(axis=0)
                    for i, j in zip(cuts[:-1], cuts[1:], strict=True)
                ]
            )
            amplitudes = vectors @ (
                np.exp(-1j * energies)[:, None] * (vectors.T @ cols)
            )
            # Re-check the original head state and exact port band, without a family read.
            coherent = amplitudes.reshape(12, 8, 3)
            returned = band @ np.einsum("b,ibk->ik", head.conj(), coherent)
            prob = np.sum(abs(returned) ** 2, axis=0)
            counts = rng.binomial(SHOTS, np.clip(prob, 0, 1))
            leakage = (h / 8) @ cols - cols @ (cols.conj().T @ (h / 8) @ cols)
            rows.append(
                dict(
                    head=name,
                    mean_energy=(energies @ weights).tolist(),
                    spectral_weights=groups.tolist(),
                    conjugate_spectral_distance=float(
                        np.sum(abs(groups[:, 1] - groups[:, 2])) / 2
                    ),
                    leakage_squared=np.sum(abs(leakage) ** 2, axis=0).tolist(),
                    return_probability_time_one=prob.tolist(),
                    returned_out_of_starts=counts.tolist(),
                )
            )
        np.testing.assert_allclose(
            rows[0]["spectral_weights"],
            np.array(rows[0]["spectral_weights"])[:, [0, 2, 1]],
            atol=2e-12,
        )
        np.testing.assert_allclose(
            rows[1]["spectral_weights"],
            np.array(rows[2]["spectral_weights"])[:, [0, 2, 1]],
            atol=2e-12,
        )
        forms.append(
            dict(
                kappa=kappa,
                global_energy_enclosure=[
                    int(min(diag - radius)) / 8,
                    int(max(diag + radius)) / 8,
                ],
                energies=[
                    float(np.mean(energies[i:j]))
                    for i, j in zip(cuts[:-1], cuts[1:], strict=True)
                ],
                rows=rows,
            )
        )
    return dict(face=list(face), dimension=96, shots_per_read=SHOTS, forms=forms)


def run():
    projectors, link = theorem_projectors()
    band = projectors[1][0]
    rng = np.random.default_rng(SEED)
    return dict(
        family_theorem=link,
        seed=SEED,
        declared_field=family_comparison(band, rng),
        native_written_face=native_assays(band, rng),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Field ratios:", result["declared_field"]["ordered_ratios"])
    print(
        "Native energy bounds:",
        [r["global_energy_enclosure"] for r in result["native_written_face"]["forms"]],
    )
