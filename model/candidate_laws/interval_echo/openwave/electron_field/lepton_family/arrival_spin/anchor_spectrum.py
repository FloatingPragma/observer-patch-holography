"""One-anchor comparison from an unlabelled vacancy and simulated Ramsey records."""

import argparse
import base64
import json
from pathlib import Path
import sys
import zlib

import mpmath as mp
import numpy as np
from scipy.linalg import expm, hankel

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (  # noqa: E402
    rest_partner as rest,
)

SHOTS = 2**22
STEP = 4
SURVEY = STEP * np.arange(128)
BLOCK = STEP * np.arange(64)
DELAYS = 2 ** np.arange(9, 27)


def instrument(effect):
    board, route = rest.guide.longest_unpinned_route()
    patch = rest.Patch(board, route, min(rest.geometries(board, route)))
    columns, pair = patch.local_space()
    source, parked = columns[:, 0], columns @ pair[:, 0]
    reference = np.zeros_like(source)
    reference[: patch.n] = np.isin(patch.ports[: patch.n], route)
    reference /= np.linalg.norm(reference)
    masses, characters = np.linalg.eigh(effect)
    sharp_address_weights = abs(characters[0]) ** 2
    first_fails = np.zeros(32)
    weights = []
    other = 0.0
    for m, address_weight in zip(masses, sharp_address_weights):

        def generator(enabled):
            return rest.guide.generator(patch.b * enabled[:, None], np.array([[m]]))

        decode = generator(np.isin(patch.tags, ["cost", "ring1"]))
        ring = generator(np.isin(patch.tags, ["cost", "ring0", "ring1"]))
        bridge = generator(patch.tags == "bridge")
        accepted = expm(-1j * np.pi / 6 * decode) @ source
        cycle = expm(-1j * np.pi / np.sqrt(8) * bridge) @ expm(
            -1j * np.pi / (2 * np.sqrt(10)) * ring
        )
        for i in range(len(first_fails)):
            accepted = cycle @ accepted
            reject = patch.bridge_flux != 0
            first_fails[i] += address_weight * np.linalg.norm(accepted[reject]) ** 2
            accepted[reject] = 0
        np.testing.assert_allclose(reference @ accepted, 0, atol=1e-13)
        dark_weight = abs(parked @ accepted) ** 2
        weights.append(address_weight * dark_weight)
        other += (
            address_weight
            * np.linalg.norm(accepted - parked * (parked @ accepted)) ** 2
        )
        h = generator(patch.tags != "bridge")
        np.testing.assert_allclose(h @ reference, m * reference, atol=2e-13)
        np.testing.assert_allclose(h @ parked, -m * parked, atol=2e-13)
        # Signed vacancy stays at the hub even after restoring the bridge and all matter hops.
        np.testing.assert_allclose(parked[patch.ports != patch.hub], 0, atol=1e-14)
        free = generator(np.ones(len(patch.b), bool))
        np.testing.assert_allclose(free @ parked, -m * parked, atol=2e-13)
        np.testing.assert_allclose(parked @ source, 0, atol=1e-14)
        shifted = h + np.eye(len(h)) / 8
        np.testing.assert_allclose(
            [reference @ shifted @ reference, parked @ shifted @ parked],
            [1 / 8 + m, 1 / 8 - m],
            atol=2e-13,
        )
    np.testing.assert_allclose(first_fails.sum() + sum(weights) + other, 1, atol=2e-13)
    # Canonicalize floating roundoff before seeded finite-count sampling.
    return (
        np.round(2 * masses, 15),
        np.round(weights, 14),
        float(np.round(other, 14)),
        first_fails,
    )


def acquire(times, omega, weights, other, rng, filter_failure=()):
    """Each setting retains +/- clicks, accepted off-mode and preparation failures."""
    enrolled = weights.sum()
    z = np.exp(1j * np.outer(times, omega)) @ (weights / enrolled)
    reads = np.stack([z.real, z.imag], axis=-1)
    active = rng.binomial(SHOTS, enrolled, size=reads.shape)
    plus = rng.binomial(active, np.clip(np.round((1 + reads) / 2, 12), 0, 1))
    off = rng.binomial(SHOTS - active, np.clip(other / (1 - enrolled), 0, 1))
    outcomes = [plus, active - plus, off]
    remaining_shots = SHOTS - active - off
    remaining_probability = 1 - enrolled - other
    for probability in filter_failure:
        failed = rng.binomial(remaining_shots, probability / remaining_probability)
        outcomes.append(failed)
        remaining_shots -= failed
        remaining_probability -= probability
    return np.stack([*outcomes, remaining_shots], axis=-1)


def fringe(counts):
    plus, minus = counts[:, :, 0], counts[:, :, 1]
    xy = (plus - minus) / (plus + minus)
    return xy[:, 0] + 1j * xy[:, 1]


def discover(counts):
    z = fringe(counts)
    matrix = hankel(z[:48], z[47:])
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    # Simultaneous Hoeffding entry bound (two quadratures, union bound), then Frobenius.
    n_min = counts[:, :, :2].sum(axis=2).min()
    cutoff = np.sqrt(matrix.size * 4 * np.log(4 * len(z) / 0.001) / n_min)
    rank = int(np.count_nonzero(singular > cutoff))
    if rank == 0:
        raise ValueError("No resolved Ramsey frequencies")
    basis = u[:, :rank]
    shift = np.linalg.lstsq(basis[:-1], basis[1:], rcond=None)[0]
    omega = np.sort(np.angle(np.linalg.eigvals(shift))) / STEP
    return omega, {
        "threshold": float(cutoff),
        "leading_singular_values": singular[:8].tolist(),
    }


def amplitudes(times, z, omega):
    return np.linalg.lstsq(np.exp(1j * np.outer(times, omega)), z, rcond=None)[0]


def infer(survey, blocks):
    omega, diagnosis = discover(survey)
    for delay, counts in zip(DELAYS, blocks, strict=True):
        phase = np.angle(amplitudes(BLOCK, fringe(counts), omega))
        cycles = np.rint((delay * omega - phase) / (2 * np.pi))
        omega = (phase + 2 * np.pi * cycles) / delay
    coefficient = amplitudes(SURVEY, fringe(survey), omega)
    return omega, coefficient, diagnosis


def pack(array):
    raw = np.asarray(array, dtype="<u4")
    return {
        "shape": list(raw.shape),
        "zlib_base64_uint32_le": base64.b64encode(
            zlib.compress(raw.tobytes(), level=9)
        ).decode("ascii"),
    }


def unpack(record):
    raw = zlib.decompress(base64.b64decode(record["zlib_base64_uint32_le"]))
    return np.frombuffer(raw, dtype="<u4").reshape(record["shape"]).astype(np.int64)


def selection(omega, weights, fitted, other, rng):
    # The successful output of a balanced Ramsey read is routed back into the parked mode.
    kept = weights.copy()
    filters = []
    failures = []
    for reject in fitted[1:]:
        t = np.pi / (reject - fitted[0])
        phase = -fitted[0] * t
        probability = np.cos((omega * t + phase) / 2) ** 2
        failures.append(kept @ (1 - probability))
        kept = kept * probability
        filters.append({"delay": float(t), "phase": float(phase)})
    counts = acquire(SURVEY, omega, kept, other, rng, failures)
    resolved, diagnosis = discover(counts)
    # One further long-delay validation uses only independently sampled selected records.
    blocks = [
        acquire(delay + BLOCK, omega, kept, other, rng, failures) for delay in DELAYS
    ]
    refined, _, _ = infer(counts, blocks)
    return {
        "settings_from_discovery_counts": filters,
        "acceptance_per_original_enrollment": float(kept.sum()),
        "modeled_heavy_contamination_given_acceptance": float(
            kept[1:].sum() / kept.sum()
        ),
        "filter_failure_probabilities": list(map(float, failures)),
        "fresh_line_count": len(resolved),
        "fresh_frequencies": refined.tolist(),
        "rank_diagnostic": diagnosis,
        "survey_counts": pack(counts),
        "refinement_counts": pack(blocks),
    }


def panel(instrument_data, seed, repeats):
    omega, weights, other, failures = instrument_data
    acquisitions = []
    main = {}
    for repeat in range(repeats):
        rng = np.random.default_rng(seed + repeat)
        survey = acquire(SURVEY, omega, weights, other, rng)
        blocks = [
            acquire(delay + BLOCK, omega, weights, other, rng) for delay in DELAYS
        ]
        fitted, coeff, diag = infer(survey, blocks)
        acquisitions.append(fitted.tolist())
        if repeat:
            continue
        held_times = np.array([71, 901, 7123, 41297, 600003, 2525121, 51867890])
        held = acquire(held_times, omega, weights, other, rng)
        predicted = np.exp(1j * np.outer(held_times, fitted)) @ coeff
        err = fringe(held) - predicted
        n_min = held[:, :, :2].sum(axis=2).min()
        main = {
            "line_count": len(fitted),
            "survey_counts": pack(survey),
            "refinement_counts": pack(blocks),
            "rank_diagnostic": diag,
            "heldout_times": held_times.tolist(),
            "heldout_counts": pack(held),
            "largest_heldout_quadrature_error": float(
                max(abs(err.real).max(), abs(err.imag).max())
            ),
            "heldout_tolerance": float(6 / np.sqrt(n_min)),
            "geometric_mode_probability": float(weights.sum()),
            "accepted_other_probability": other,
            "first_bridge_failure_probabilities": failures.tolist(),
            "first_failure_calibration_counts_then_accepted": rng.multinomial(
                SHOTS, [*failures, 1 - failures.sum()]
            ).tolist(),
        }
        if len(fitted) > 1:
            main["light_filter"] = selection(omega, weights, fitted, other, rng)
    main["frequencies_by_acquisition"] = acquisitions
    return main


def scale_audit(read_ratios):
    # Fourier eigenvalues of the circulant, evaluated without building its matrix.
    with mp.workdps(60):
        theta = mp.mpf(2) / 9

        def ratios(angle):
            levels = sorted(
                (1 + mp.sqrt(2) * mp.cos(angle + j * 2 * mp.pi / 3)) ** 2
                for j in range(3)
            )
            return [v / levels[0] for v in levels]

        expected = ratios(theta)
        electron, muon = mp.mpf("0.51099895069"), mp.mpf("105.6583755")
        measured_ratio = mp.mpf("206.7682827")
        adjusted = mp.findroot(lambda a: ratios(a)[1] - measured_ratio, theta)
        required_electron = muon / expected[1]
        return {
            "ratio_source": "https://physics.nist.gov/cuu/pdf/RevModPhys.97.025002.pdf",
            "codata_eq_187_ratio": float(measured_ratio),
            "ratio_standard_uncertainty": 4.6e-6,
            "closed_form_ratios": [mp.nstr(v, 48) for v in expected],
            "muon_difference_ppm": float(1e6 * (expected[1] / measured_ratio - 1)),
            "difference_in_quoted_ratio_uncertainties": float(
                (expected[1] - measured_ratio) / mp.mpf("0.0000046")
            ),
            "range_of_eight_count_ratio_estimates": [
                float(min(read_ratios)),
                float(max(read_ratios)),
            ],
            "electron_keV_implied_by_muon": float(1000 * required_electron),
            "implied_minus_measured_electron_eV": float(
                1e6 * (required_electron - electron)
            ),
            "muon_excess_keV_exact_e_anchor": float(
                1000 * (electron * expected[1] - muon)
            ),
            "muon_excess_keV_literal_511_anchor": float(
                1000 * (mp.mpf("0.511") * expected[1] - muon)
            ),
            "phase_fitted_using_CODATA_ratio": float(adjusted),
            "fitted_minus_declared_phase_radians": float(adjusted - theta),
            "conditional_tau_MeV_with_fitted_phase": float(
                electron * ratios(adjusted)[2]
            ),
        }


def run():
    effect, _ = rest.guide.register_effect()
    setup = instrument(effect)
    main = panel(setup, 74021, 8)
    scalar = np.eye(len(effect)) * np.trace(effect).real / len(effect)
    control = panel(instrument(scalar), 74121, 1)
    fits = np.array(main["frequencies_by_acquisition"])
    ratios = fits / fits[:, :1]
    # Laboratory inputs appear only in this terminal comparison.
    measured = np.array([0.51099895069, 105.6583755, 1776.93])
    ideal = setup[0] / setup[0][0] * measured[0]
    return {
        "shots_per_setting": SHOTS,
        "survey_step": STEP,
        "refinement_delays": DELAYS.tolist(),
        "seed": 74021,
        "control_seed": 74121,
        "scale_audit": scale_audit(ratios[:, 1]),
        "discovery": main,
        "scalar_control": control,
        "comparison": {
            "unit": "MeV",
            "pdg": "https://pdg.lbl.gov/2026/tables/rpp2026-sum-leptons.pdf",
            "measured": measured.tolist(),
            "reported_uncertainty": [1.6e-10, 2.3e-6, 0.09],
            "measured_electron_anchor_first_acquisition": (
                ratios[0] * measured[0]
            ).tolist(),
            "repeat_standard_deviation": (
                ratios.std(axis=0, ddof=1) * measured[0]
            ).tolist(),
            "literal_511_keV_anchor_first_acquisition": (ratios[0] * 0.511).tolist(),
            "exact_generator_electron_anchor": ideal.tolist(),
            "exact_generator_minus_measurement": (ideal - measured).tolist(),
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as stream:
        json.dump(run(), stream, indent=2, allow_nan=False)
        stream.write("\n")
