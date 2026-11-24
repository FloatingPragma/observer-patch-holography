"""Read inertia at the far end of a guide after heralding one local vacancy."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (  # noqa: E402
    anchor_spectrum,
    recoil,
    rest_partner,
)

guide = rest_partner.guide
SHOTS = 2**18
SURVEY = np.arange(1024) * 0.5
BLOCK = np.arange(256) * 0.5
DELAYS = 2 ** np.arange(8, 22)
SEED = 178037


def port_read(route, ends):
    """Read the endpoint by subtracting its background from the actual written heads."""
    head = guide.head_paths(route)[0]
    background = head["reference_head_counts_minus_charge"]
    edges = [tuple(sorted(edge)) for edge in zip(route, route[1:])]
    masks = []
    for end in ends:
        site = route[end]
        charges = []
        for word in head["written_words_in_port_order"]:
            total = sum(
                edge[0 if word & (1 << j) else 1] == site
                for j, edge in enumerate(edges)
            )
            charges.append(total - background[end])
        mask = np.r_[charges, charges[:-1]] == -1
        masks.append(mask)
    return masks


def channels(bonds, effect, route, detector=-1):
    source_mask, distal = port_read(route, (0, detector))
    n = len(route)
    source_mask[n:] = False  # local vertex-grade read after the charge click
    (src,) = np.flatnonzero(source_mask)
    (dst,) = np.flatnonzero(distal[:n])  # resolve vertex grade at an interior detector
    e, chars = np.linalg.eigh(effect)
    np.testing.assert_equal(src, 0)
    np.testing.assert_equal(dst, detector % n)
    # Independent full-register short flight, including coherent cyclic-address amplitudes.
    h = guide.generator(bonds, effect)
    full = expm(-1j * h)[:, 3].reshape(-1, 3)
    entries, spectra = [], []
    for m, char in zip(e, chars.T, strict=True):
        single = guide.generator(bonds, np.array([[m]]))
        levels, v = np.linalg.eigh(single)
        entry = (np.exp(-1j * levels) @ (v[src] * v[1])) * char[0].conjugate()
        np.testing.assert_allclose(entry, char.conj() @ full[src], atol=2e-13)
        entries.append(abs(entry) ** 2)
        spectra.append((levels.round(12), v[src] * v[dst]))
    return e, np.array(entries), spectra


def select(entries, masses, fitted, target):
    # Balanced local vertex/incidence split, free hold on grades +F and -F,
    # recombine with a phase, then keep the vertex. Both grades have port zero's head.
    keep = entries.copy()
    failures, settings = [], []
    if target is not None:
        for j, freq in enumerate(fitted):
            if j == target:
                continue
            delay = np.pi / abs(freq - fitted[target])
            phase = -delay * fitted[target]
            transmit = np.cos((2 * masses * delay + phase) / 2) ** 2
            failures.append(float(keep @ (1 - transmit)))
            keep *= transmit
            settings.append([float(delay), float(phase)])
    return keep, [float(1 - sum(entries)), *failures], settings


def arrivals(times, spectra):
    return np.array(
        [np.exp(-1j * np.outer(times, levels)) @ coeff for levels, coeff in spectra]
    )


def acquire(times, flight, masses, spectra, keep, failures, rng):
    # The two histories start with the same charge at the same source port. Only
    # the late arm accrues its held vertex phase before the common final flight.
    early = arrivals(times + flight, spectra)
    late = np.exp(-1j * np.outer(masses, times)) * arrivals([flight], spectra)
    probability = keep @ ((abs(early) ** 2 + abs(late) ** 2) / 2)
    interference = keep @ (early.conj() * late)
    counts = []
    for q in (interference.real, interference.imag):
        bins = np.column_stack(
            [
                (probability + q) / 2,
                (probability - q) / 2,
                sum(keep) - probability,
                *[np.full(len(times), f) for f in failures],
            ]
        )
        np.testing.assert_allclose(bins.sum(axis=1), 1, atol=1e-13)
        bins = np.maximum(0, bins.round(11))
        counts.append(np.array([rng.multinomial(SHOTS, b / b.sum()) for b in bins]))
    return np.stack(counts, axis=1)


def fringe(counts):
    signed = (counts[:, :, 0] - counts[:, :, 1]) / SHOTS
    return signed[:, 0] + 1j * signed[:, 1]


def amplitude(t, z, freq):
    window = np.hanning(len(t))
    return np.exp(-1j * np.outer(freq, t)) @ (window * z) / window.sum()


def read(survey, refined):
    z = fringe(survey)
    nfft = 16384
    power = abs(np.fft.fft(z * np.hanning(len(z)), nfft)[: nfft // 2])
    # Ignore the tapered DC main lobe. Only peaks 10% as high as the strongest
    # positive line enter the reader. The threshold receives no family parameters.
    peaks = find_peaks(power, height=power[80:].max() / 10)[0]
    indices = peaks[peaks >= 80][:3]
    if len(indices) < 3:
        raise ValueError("Fewer than three spatial recoil lines in the arrival counts")
    freq = 2 * np.pi * indices / (0.5 * nfft)
    base = fringe(survey[: len(BLOCK)])
    for delay, block in zip(DELAYS, refined, strict=True):
        first = amplitude(BLOCK, base, freq)
        second = amplitude(BLOCK, fringe(block), freq)
        phi = np.angle(second * first.conj())
        freq = (phi + 2 * np.pi * np.rint((delay * freq - phi) / (2 * np.pi))) / delay
    return freq


def measurement(flight, masses, spectra, keep, failures, rng):
    survey = acquire(SURVEY, flight, masses, spectra, keep, failures, rng)
    blocks = [
        acquire(delay + BLOCK, flight, masses, spectra, keep, failures, rng)
        for delay in DELAYS
    ]
    freq = read(survey, blocks)
    # At zero delay both histories freely traverse the same guide. Either pointer
    # basis measures the unconditional incidence of a distal signed-charge read.
    shot = survey[0, 0]
    return freq, {
        "common_flight": flight,
        "lowest_three_positive_frequencies": freq.tolist(),
        "distal_given_source_and_tag_at_zero_delay": float(
            sum(shot[:2]) / sum(shot[:3])
        ),
        "survey_counts": anchor_spectrum.pack(survey),
        "long_wait_counts": anchor_spectrum.pack(blocks),
    }


def verify_local_circuit(bonds, effect, route, settings, keep, masses, spectra):
    # Evolve the actual coherent register. The two-level splitter acts only on
    # the vertex and incidence grade at the source, which share a signed head word.
    h = guide.generator(bonds, effect)
    nf = len(effect)
    sample = expm(-1j * h)[:nf, nf]
    h_local = guide.generator(np.zeros((1, 1)), effect)
    half = np.kron(np.array([[1, 1], [1, -1]]) / np.sqrt(2), np.eye(nf))
    for t, phase in settings:
        psi = half @ np.r_[sample, np.zeros(nf)]
        psi = expm(-1j * t * h_local) @ psi
        psi[nf:] *= np.exp(1j * phase)
        sample = (half @ psi)[:nf]
    np.testing.assert_allclose(np.vdot(sample, sample).real, sum(keep), atol=1e-13)
    delay, flight = 3.5, 6
    dst = slice(nf * (len(route) - 1), nf * len(route))
    early = expm(-1j * (flight + delay) * h)[dst, :nf] @ sample
    late = expm(-1j * flight * h)[dst, :nf] @ (expm(-1j * delay * effect) @ sample)
    first = arrivals([flight + delay], spectra)[:, 0]
    second = arrivals([flight], spectra)[:, 0] * np.exp(-1j * delay * masses)
    expected = [
        keep @ abs(first) ** 2,
        keep @ abs(second) ** 2,
        keep @ (first.conj() * second),
    ]
    np.testing.assert_allclose(
        [np.vdot(early, early), np.vdot(late, late), np.vdot(early, late)],
        expected,
        atol=3e-12,
    )


def spatial_read(panels, keep, failures, masses, frequencies, rng):
    times = SURVEY[:512]
    counts = [
        acquire(times, 6, masses, spectra, keep, failures, rng) for spectra in panels
    ]
    scores = spatial_scores(counts, frequencies, times)
    return {
        "ranked_best_modes": np.argmin(scores, axis=1).tolist(),
        "relative_residuals_by_mode": scores.tolist(),
        "survey_counts_site_time_basis_bin": anchor_spectrum.pack(counts),
    }


def spatial_scores(counts, frequencies, times):
    n = len(counts)
    coefficients = [amplitude(times, fringe(c), np.r_[0, frequencies]) for c in counts]
    # The local late-arm amplitude multiplies every frequency at a detector. Its
    # phase and magnitude cancel by referencing the local DC line in the regression.
    z = np.array(coefficients)
    waves = np.cos(np.pi * np.outer(np.arange(n) + 0.5, np.arange(n)) / n)
    modes = z[:, :1] * waves
    scores = []
    for line in z[:, 1:].T:
        constants = modes.conj().T @ line / np.sum(abs(modes) ** 2, axis=0)
        scores.append(
            np.linalg.norm(line[:, None] - modes * constants, axis=0)
            / np.linalg.norm(line)
        )
    return np.array(scores)


def run():
    board, route = guide.longest_unpinned_route()
    bonds = guide.path_incidence(board, route)
    effect, _ = guide.register_effect()
    prior = json.loads((HERE / "anchor_spectrum.json").read_text())
    fitted = prior["discovery"]["frequencies_by_acquisition"][0]
    rng = np.random.default_rng(SEED)
    studies = {}
    scalar = np.eye(3) * effect.trace().real / 3
    for name, register in (("family", effect), ("scalar", scalar)):
        masses, entries, spectra = channels(bonds, register, route)
        # Close the first contact before the source herald, then independently
        # close only the middle contact during the two release histories.
        cut = bonds.copy()
        cut[0] = 0
        _, source_cut, _ = channels(cut, register, route)
        middle = bonds.copy()
        middle[(len(route) - 1) // 2] = 0
        _, _, blocked = channels(middle, register, route)
        np.testing.assert_allclose(source_cut, 0, atol=1e-25)
        runs = []
        for target in range(len(fitted)) if name == "family" else (None,):
            keep, failures, settings = select(entries, masses, fitted, target)
            verify_local_circuit(
                bonds, register, route, settings, keep, masses, spectra
            )
            first, survey = measurement(6, masses, spectra, keep, failures, rng)
            mass, stiff, third = recoil.fit_first_two(
                len(route), first[0:1], first[1:2]
            )
            held, fresh = measurement(8, masses, spectra, keep, failures, rng)
            cut_counts = acquire(SURVEY[:32], 6, masses, blocked, keep, failures, rng)
            np.testing.assert_array_equal(cut_counts[:, :, :2], 0)
            runs.append(
                {
                    "tag": target,
                    "filter_settings_delay_phase": settings,
                    "source_click_probability": float(sum(entries)),
                    "prep_failure_probabilities": failures,
                    "joint_tag_probability": float(keep.sum()),
                    "first_flight": survey,
                    "dispersion_mass": float(mass[0]),
                    "stiffness": float(stiff[0]),
                    "third_frequency_prediction": float(third[0]),
                    "fresh_flight": fresh,
                    "third_frequency_residual": float(held[2] - third[0]),
                    "passes_third_frequency": bool(abs(held[2] - third[0]) < 1e-5),
                    "release_cut_counts": anchor_spectrum.pack(cut_counts),
                }
            )
        studies[name] = {"source_cut_probability": float(sum(source_cut)), "runs": runs}
    independent = np.random.default_rng(SEED + 1)
    panels = [channels(bonds, effect, route, site)[2] for site in range(len(route))]
    masses, entries, _ = channels(bonds, effect, route)
    for target, row in enumerate(studies["family"]["runs"]):
        keep, failures, _ = select(entries, masses, fitted, target)
        frequencies = row["first_flight"]["lowest_three_positive_frequencies"]
        row["spatial_read"] = spatial_read(
            panels, keep, failures, masses, frequencies, independent
        )
    # Subject the scalar control to the identical three filters, including the
    # low-yield middle setting. Compare the observed frequencies under each tag.
    masses, entries, spectra = channels(bonds, scalar, route)
    scalar_tags = []
    for target in range(len(fitted)):
        keep, failures, _ = select(entries, masses, fitted, target)
        frequencies, raw = measurement(6, masses, spectra, keep, failures, independent)
        mass, _, _ = recoil.fit_first_two(
            len(route), frequencies[0:1], frequencies[1:2]
        )
        scalar_tags.append(
            {
                "tag": target,
                "joint_tag_probability": float(sum(keep)),
                "dispersion_mass": float(mass[0]),
                "flight": raw,
            }
        )
    studies["scalar"]["same_three_tag_settings"] = scalar_tags
    return {
        "route": route,
        "distinct_cells": len({p // 12 for p in route}),
        "shots_per_basis": SHOTS,
        "seed": SEED,
        "spatial_and_scalar_tag_seed": SEED + 1,
        "source_wait": 1,
        "survey_times": [0, 0.5, len(SURVEY)],
        "block_times": [0, 0.5, len(BLOCK)],
        "long_wait_delays": DELAYS.tolist(),
        "bins": [
            "distal signed charge, pointer +",
            "distal signed charge, pointer -",
            "away from distal port",
            "missed local source",
            "first filter reject if applied",
            "second filter reject if applied",
        ],
        "studies": studies,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.output.open("x") as stream:
        json.dump(run(), stream, indent=2, allow_nan=False)
        stream.write("\n")
