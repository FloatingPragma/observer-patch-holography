"""Reconstruct the mirror-averaged material response from two measurable native source records."""

import argparse
from pathlib import Path
import json

import numpy as np
from ..matched_probe import measure as matched
from . import measure as reference
from model.executable_model import rotations
from ...closed_region.measure import check

HERE = Path(__file__).parent


def symmetries(board, seam):
    a, b = seam
    g = np.array(next(p for p in rotations() if 12 + p[a] == b))
    mirror = np.r_[12 + g, np.argsort(g)]
    for n in (1, 2):
        h, words, bits = matched.fp.generator(board, seam, n, 8, True)
        words_to_index = {word: i for i, word in enumerate(words)}
        qi = np.repeat(bits[:, :12].sum(axis=1), 4) - np.tile([0, 0, 1, 1], len(words))
        row, col = h.nonzero()
        np.testing.assert_array_equal(qi[row], qi[col])
        indices = [
            4 * words_to_index[tuple(sorted(mirror[list(word)]))] + (m + 2) % 4
            for word in words
            for m in range(4)
        ]
        assert (h[indices, :][:, indices] != h).nnz == 0


def run():
    board, dist, origins, blocks = matched.calibrated()
    eyes = np.eye(4)
    bright = (eyes[1] + eyes[2]) / np.sqrt(2)
    B = np.outer(bright, bright)
    C = np.diag([0.0, 1.0, 1.0, 0.0])
    symmetries(board, matched.fp.make_board()[1])
    _, terminal = matched.fp.restricted(8)
    assert set(terminal) == {(0, 1), (2, 0)}
    basis = [np.outer(eyes[i], eyes[i]) for i in (0, 1)] + [B]
    home = blocks[1, None]["number"]
    ports = sum(blocks[1, p]["number"] - home for p in origins) / len(origins)
    radial = np.einsum("txij,x->tij", ports, dist)
    cal = np.stack([matched.trace(r, radial) for r in basis])
    u = 2 * cal[1] - cal[2] - cal[0]
    v = 2 * (cal[2] - cal[1])
    recon = (
        cal[0, :, None, None] * np.eye(4) + u[:, None, None] * C + v[:, None, None] * B
    )
    error = float(np.max(abs(radial - recon)))
    assert error < 2e-10
    rng = np.random.default_rng(191204)
    cases = {
        name: matched.stocks()[name][1][1] for name in ("target", "phase_feedback")
    }
    cases["paired"] = matched.stocks()["target"][0][1]
    for rank in (1, 2, 3, 4):
        z = rng.normal(size=(4, rank)) + 1j * rng.normal(size=(4, rank))
        r = z @ z.conj().T
        cases[f"complex_rank_{rank}"] = r / np.trace(r)
    # On either root endpoint the two signs of memory parity have the same
    # compression. Each symmetric charge distribution depends only on v.
    _, seam, _, h, _, _, _ = reference.folded()
    energy, bath = np.linalg.eigh(h)
    beta = np.log(8 / 7)
    occupation = np.exp(-beta * (energy - energy[0]))
    bath = (bath * (occupation / occupation.sum())) @ bath.T

    def modular_law(r):
        embedded = np.zeros((8, 8), dtype=complex)
        words = [0, 4, 2, 6]
        embedded[np.ix_(words, words)] = r
        # K has real eigenprojectors. The imaginary antisymmetric part of the
        # input density has identically zero trace against each projector.
        spectrum = reference.modular_spectra(
            bath, sorted(seam), seam, [embedded.real] * 3
        )
        return np.array(spectrum["level_nats"]), np.array(spectrum["probabilities"][0])

    levels, dist_dark = modular_law(basis[0])
    _, dist_bright = modular_law(B)
    shots = 2**25
    # Each port and each launch uses fresh copies. Weighted Hoeffding for the
    # two half-weight launches and one full-weight blank gives squared weight 3/2.
    # Split one simultaneous error budget over standards, target ports and source reads.
    tests = (len(basis) + len(cases)) * len(home) + 3 * len(cases) + 2
    eps = np.sqrt(np.log(2 * tests / 0.01) / (2 * shots))
    radial_width = float(np.sqrt(1.5 * np.dot(dist, dist)) * eps)

    def port_counts(state, generator):
        counts = {}
        for launch in (None, *origins):
            prob = matched.trace(state, blocks[1, launch]["number"])
            counts[str(launch)] = generator.binomial(shots, prob.clip(0, 1))
        radial_read = sum(
            (counts[str(launch)] - counts["None"]) @ dist for launch in origins
        ) / (len(origins) * shots)
        return counts, radial_read

    charge_rng = np.random.default_rng(191206)

    def charge_counts(prob):
        prob = np.maximum(prob, 0)
        return charge_rng.multinomial(shots, prob / prob.sum())

    standard_charge_counts = [charge_counts(prob) for prob in (dist_dark, dist_bright)]
    kd, kb = [c @ levels / shots for c in standard_charge_counts]
    assert abs(kb - kd) > 1
    calibration_rng = np.random.default_rng(191205)
    standard_rows = [port_counts(state, calibration_rng) for state in basis]
    measured_cal = np.stack([entry[1] for entry in standard_rows])
    measured_u = 2 * measured_cal[1] - measured_cal[2] - measured_cal[0]
    measured_v = 2 * (measured_cal[2] - measured_cal[1])
    assert np.max(abs(measured_cal - cal)) < radial_width
    # All states in these mixtures have maximally mixed coin and endpoint.
    # The physical triangle is generated by dark unmatched, dark matched and bright.
    matched_dark = (eyes[1] - eyes[2]) / np.sqrt(2)
    vertex = [np.diag([0.5, 0, 0, 0.5]), np.outer(matched_dark, matched_dark), B]
    same_charge = []
    for coincidence in (0.25, 0.50, 0.75, 1.0):
        echo = 0.25
        mixture = (
            (1 - coincidence) * vertex[0]
            + (coincidence - echo) * vertex[1]
            + echo * vertex[2]
        )
        tensor = mixture.reshape(2, 2, 2, 2)
        for first, second in ((0, 2), (1, 3)):
            np.testing.assert_allclose(
                tensor.trace(axis1=first, axis2=second), np.eye(2) / 2
            )
        np.testing.assert_allclose(
            matched.trace(mixture, np.stack([C, B])), [coincidence, echo]
        )
        _, distribution = modular_law(mixture)
        np.testing.assert_allclose(
            distribution, (1 - echo) * dist_dark + echo * dist_bright, atol=2e-12
        )
        delta = matched.trace(mixture - cases["paired"], radial)
        same_charge.append(
            dict(coincidence=coincidence, unwrite=echo, radial_change=delta.tolist())
        )
    rows = {}
    for name, r in cases.items():
        local = matched.trace(r, np.stack([C, B]))
        predicted = cal[0] + u * local[0] + v * local[1]
        measured_local = rng.binomial(shots, local.clip(0, 1)) / shots
        inferred = (
            measured_cal[0]
            + measured_u * measured_local[0]
            + measured_v * measured_local[1]
        )
        exact = matched.trace(r, radial)
        np.testing.assert_allclose(predicted, exact, rtol=0, atol=2e-10)
        _, distribution = modular_law(r)
        np.testing.assert_allclose(
            distribution,
            dist_dark * (1 - local[1]) + dist_bright * local[1],
            rtol=0,
            atol=2e-12,
        )
        charge = charge_counts(distribution)
        v_from_charge = (charge @ levels / shots - kd) / (kb - kd)
        charge_prediction = (
            measured_cal[0]
            + measured_u * measured_local[0]
            + measured_v * v_from_charge
        )
        # The true mean interpolates its two calibration means with positive
        # weights. Target and interpolated calibration each cost one mean width.
        v_error = 2 * float(np.ptp(levels)) * eps / abs(kb - kd)
        charge_width = (
            4 * radial_width + eps * abs(measured_u) + v_error * abs(measured_v)
        )
        raw, obs = port_counts(r, rng)
        assert np.all(abs(obs - charge_prediction) < charge_width)
        # For physical source frequencies 0 <= v <= u <= 1, the sum of
        # absolute interpolation weights is at most three. Add the independent
        # target response uncertainty and the two source-frequency errors.
        inferred_width = 4 * radial_width + eps * (abs(measured_u) + abs(measured_v))
        assert np.all(abs(obs - inferred) < inferred_width)
        assert np.max(abs(obs - predicted)) < radial_width
        rows[name] = dict(
            local_reads=local.tolist(),
            measured_local=measured_local.tolist(),
            inferred_from_charge_and_coincidence=charge_prediction.tolist(),
            modular_counts=charge.tolist(),
            unwrite_from_modular_mean=v_from_charge,
            inferred_from_source_counts=inferred.tolist(),
            prediction=predicted.tolist(),
            observed=obs.tolist(),
            source_matrix_real=r.real.tolist(),
            source_matrix_imag=r.imag.tolist(),
            counts={k: x.tolist() for k, x in raw.items()},
        )
    return dict(
        times=[*matched.fp.TIMES, "infinite_time_mean"],
        trials_per_port=shots,
        sampling_99_half_width=radial_width,
        fixed_charge_family=same_charge,
        fully_sampled_half_width=inferred_width.tolist(),
        modular_calibration=dict(
            nats=levels.tolist(),
            dark_counts=standard_charge_counts[0].tolist(),
            bright_counts=standard_charge_counts[1].tolist(),
            dark_mean=kd,
            bright_mean=kb,
            unwrite_99_half_width=v_error,
            prediction_99_half_width=charge_width.tolist(),
        ),
        source_99_half_width=eps,
        standard_order=["unmatched", "matched", "bright"],
        standard_counts=[
            {k: v.tolist() for k, v in entry[0].items()} for entry in standard_rows
        ],
        measured_standards=measured_cal.tolist(),
        calibration=dict(
            baseline=cal[0].tolist(), coincidence=u.tolist(), unwrite=v.tolist()
        ),
        maximum_matrix_residual=error,
        held=rows,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    result = run()
    path = HERE / "two_reads.json"
    if p.parse_args().write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), result)
