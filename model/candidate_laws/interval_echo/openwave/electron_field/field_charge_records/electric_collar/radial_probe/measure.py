"""Read a reciprocal radial probe on the same late triangle/reference arrivals."""

import argparse
from hashlib import sha256
import json
from math import inf, log, nextafter, sqrt
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

from ....fan_path.fixture import FanPath, archived_fixture
from ....fan_path.local_launch.face_current.measure import LABEL, READS, pair_amplitudes
from .apparatus import RadialCollar


HERE = Path(__file__).resolve().parent
TIMES = (1, 2, 4)
TRIALS = 2**40
SEED = 10451046
Y = np.array([[0, -1j], [1j, 0]])


def inf_norm(matrix):
    return float(np.max(np.sum(abs(matrix), axis=1)))


def spectral_aperture(h):
    """Enclose evolution from a checked symmetric dense diagonalization."""
    n = len(h)
    twice = 2 * h
    np.testing.assert_array_equal(twice, np.rint(twice))
    np.testing.assert_array_equal(h, h.T)
    assert n <= 512 and np.max(abs(h)) < 64
    values, basis = eigh(h)
    absolute = abs(basis)
    columns = absolute.sum(axis=0)
    # Nine dot/sum roundoff factors cover Gram, reconstruction, and outward
    # norms, including the two diagonal scalings and cancellation.
    guard = 16 * n * 2.0**-53
    gram = (1 + guard) * (
        inf_norm(basis.T @ basis - np.eye(n))
        + guard * (1 + np.max(absolute.T @ absolute.sum(axis=1)))
    )
    reconstruction = (1 + guard) * (
        inf_norm((basis * values) @ basis.T - h)
        + guard * (inf_norm(h) + np.max(absolute @ (abs(values) * columns)))
    )
    assert gram < 1e-8 and max(abs(values)) < 128
    operator_error = reconstruction + 3 * max(abs(values)) * gram
    wave_error = nextafter(float(3 * gram + max(TIMES) * operator_error + 2e-12), inf)
    return values, basis, wave_error


def radial_effects(collar, h):
    labels = connected_components(h, directed=False)[1]
    starts = np.arange(collar.starts).reshape(2, collar.backgrounds, 2)
    assert np.array_equal(labels[starts[:, :, 0]], labels[starts[:, :, 1]])
    destinations = np.array([p for _, p, _ in collar.basis])
    effects = np.zeros((len(TIMES), 2, collar.backgrounds, 2, 2), dtype=complex)
    worst = 0
    for component in np.unique(labels[starts]):
        group = np.flatnonzero(labels == component)
        entries = np.argwhere(labels[starts[:, :, 0]] == component)
        initial = starts[entries[:, 0], entries[:, 1]].ravel()
        columns = np.searchsorted(group, initial)
        energies, basis, error = spectral_aperture(h[group][:, group].toarray())
        worst = max(worst, error)
        for row, time in enumerate(TIMES):
            phases = np.exp(-1j * time * energies)
            flight = basis @ (phases[:, None] * basis.T[:, columns])
            for k, (direction, bg) in enumerate(entries):
                end = destinations[group] == collar.radial[1 - direction]
                transfers = flight[end, 2 * k : 2 * k + 2]
                effects[row, direction, bg] = transfers.conj().T @ transfers
    return effects, worst


def late_arrival(path, label, source_receipt):
    observed = source_receipt["receiver_reference_pauli_moments"]
    wave = expm_multiply(-1j * 30.8125 * path.h, path.source())
    amps, cell, hub = pair_amplitudes(path, wave, label)
    samples = [amps @ (np.array([1, -sign * 1j]) / sqrt(2)) for sign in (1, -1)]
    density = np.array([sample.T @ sample.conj() / 2 for sample in samples])
    admissions = density.trace(axis1=1, axis2=2).real
    currents = np.trace(density @ (Y / 2), axis1=1, axis2=2).real
    np.testing.assert_allclose(
        admissions,
        [(observed[0][0] + s * observed[0][2]) / 2 for s in (1, -1)],
        atol=2e-13,
    )
    np.testing.assert_allclose(
        currents,
        [(observed[2][0] + s * observed[2][2]) / 4 for s in (1, -1)],
        atol=2e-13,
    )
    assert currents[0] > 0 > currents[1]
    phase = np.diag([1, -1])
    changed = phase @ density @ phase
    np.testing.assert_array_equal(
        changed.diagonal(axis1=1, axis2=2), density.diagonal(axis1=1, axis2=2)
    )
    np.testing.assert_array_equal(
        np.trace(changed @ (Y / 2), axis1=1, axis2=2).real, -currents
    )
    incoming, target = path.bridges[-1]
    for i in range(path.blocks[3, 0].rows.start, path.blocks[3, 0].rows.stop):
        x, word = path.basis[i]
        j = path.ids[incoming, word]
        assert (
            x == target
            and path.bridge_term[i, j] == -4
            and path.bridge_term[i, i] == path.bridge_term[j, j] == 4
        )
    refusal = [1 - cell, cell - hub, hub - float(sum(admissions))]
    return np.stack([density, changed]), admissions, currents / admissions, refusal


def sample_counts(probabilities, admissions, refusal, rng):
    counts = np.empty((*probabilities.shape[:-1], 7), dtype=np.int64)
    for index in np.ndindex(probabilities.shape[:-1]):
        plus, minus = probabilities[index]
        bins = np.array(
            [*refusal, admissions[0] - plus, plus, admissions[1] - minus, minus]
        )
        assert min(bins) >= -2e-14 and abs(sum(bins) - 1) < 2e-12
        bins = np.maximum(bins, 0)
        counts[index] = rng.multinomial(TRIALS, bins / sum(bins))
    np.testing.assert_array_equal(counts.sum(axis=-1), TRIALS)
    return counts


def radial_diagnostics(collar, raw, counts, response_error):
    """Compare same-field directions and the independently opened reciprocal lane."""
    difference = raw["mobile"] - raw["neutral"]
    outward, inward = difference[:, :, 0], difference[:, :, 1]
    directions = raw["mobile"], raw["neutral"]
    sampled = (
        np.array(counts["mobile"])[..., [4, 6]]
        - np.array(counts["neutral"])[..., [4, 6]]
    ) / TRIALS
    tested_bins = 4 * 4 * 2 * len(TIMES) * 2 * collar.backgrounds * 3
    bin_error = sqrt(log(2 * tested_bins / 0.01) / (2 * TRIALS))
    failures = 4 * bin_error + response_error

    # Start both directions on a locally open rail by reversing only the
    # probe head in the inward preparation. All other written heads stay fixed.
    bit = collar.probe_mask >> collar.inside_bits
    outward_open = [
        b
        for b in range(collar.backgrounds)
        if bool(b & bit) == (collar.probe_sign == -1)
    ]
    opened = []
    for b in outward_open:
        partner = b ^ bit
        for direction, background in ((0, b), (1, partner)):
            for w in range(2):
                i = (direction * collar.backgrounds + background) * 2 + w
                assert collar.terms["probe_hop"][i, i] == 4
        opened.append((b, partner))

    comparisons = {}
    for label, pairs in (
        (
            "same_exterior_trial_both_directions_observed",
            [(b, b) for b in range(collar.backgrounds)],
        ),
        ("reciprocal_both_lanes_initially_open", opened),
    ):
        best = None
        for outward_word, inward_word in pairs:
            diff = outward[:, :, outward_word] - inward[:, :, inward_word]
            for wait, ref in np.ndindex(diff.shape[1:]):
                scores = diff[:, wait, ref]
                if scores[0] * scores[1] >= 0 or min(abs(scores)) <= response_error:
                    continue
                if (
                    label.startswith("same")
                    and min(
                        np.min(arm[:, wait, direction, word, ref])
                        for arm in directions
                        for direction, word in ((0, outward_word), (1, inward_word))
                    )
                    <= response_error
                ):
                    continue
                swing = abs(scores[0] - scores[1])
                if best is None or swing > best["absolute_phase_swing"]:
                    estimates = (
                        sampled[:, wait, 0, outward_word, ref]
                        - sampled[:, wait, 1, inward_word, ref]
                    )
                    best = {
                        "absolute_phase_swing": float(swing),
                        "outward_and_inward_prepared_words": [
                            outward_word,
                            inward_word,
                        ],
                        "reference_y": (1, -1)[ref],
                        "collar_wait": TIMES[wait],
                        "original_and_reversed_phase_scores_per_launch": scores.tolist(),
                        "original_and_reversed_simulated_scores_per_launch": estimates.tolist(),
                        "simultaneous_99pct_interval_per_phase_score": [
                            [float(x - failures), float(x + failures)]
                            for x in estimates
                        ],
                    }
        assert best is not None
        comparisons[label] = best
    best = None
    pooled = (outward - inward).sum(axis=-1)
    for wait, background in np.ndindex(pooled.shape[1:]):
        scores = pooled[:, wait, background]
        if scores[0] * scores[1] >= 0 or min(abs(scores)) <= 2 * response_error:
            continue
        if (
            min(np.min(arm[:, wait, :, background].sum(axis=-1)) for arm in directions)
            <= 2 * response_error
        ):
            continue
        swing = abs(scores[0] - scores[1])
        if best is None or swing > best["absolute_phase_swing"]:
            estimates = (
                sampled[:, wait, 0, background] - sampled[:, wait, 1, background]
            ).sum(axis=-1)
            radius = 4 * bin_error + 2 * response_error
            best = {
                "absolute_phase_swing": float(swing),
                "outward_and_inward_prepared_words": [background, background],
                "collar_wait": TIMES[wait],
                "original_and_reversed_phase_scores_per_launch": scores.tolist(),
                "original_and_reversed_simulated_scores_per_launch": estimates.tolist(),
                "simultaneous_99pct_interval_per_phase_score": [
                    [float(x - radius), float(x + radius)] for x in estimates
                ],
            }
    assert best is not None
    comparisons["same_exterior_reference_traced"] = best
    return comparisons, bin_error


def static_clock(collar, scalar, acceptance, rng):
    """Record the two independently placed arms of the exact static echo."""
    wait = 1 / 8
    # A calibrated controlled stationary/neutral echo has this same scalar
    # phase for every prepared exterior word with either probe position.
    phase = np.exp(-1j * scalar * wait)
    p = [(1 + phase.real) / 2, (1 + phase.imag) / 2]
    return {
        "clock_wait": wait,
        "common_relative_phase_radians": -scalar * wait,
        "contrast_order": ["real", "imaginary"],
        "outcome_order": ["pair refused", "+", "-"],
        "same_field_placements": [
            {
                "probe_port": port,
                "plus_probability_given_pair": p,
                "raw_simulated_counts": [
                    rng.multinomial(
                        TRIALS, [1 - acceptance, acceptance * x, acceptance * (1 - x)]
                    ).tolist()
                    for x in p
                ],
            }
            for port in collar.radial
        ],
    }


def report(path, sign, words, arrival, incoming_error, rng):
    collar = RadialCollar(path, sign, words)
    hamiltonians, scalar = collar.generators()
    zero_static = collar.static_difference(hamiltonians, scalar)
    density, admissions, currents, refusal = arrival
    raw = {}
    errors = {}
    archived = {}
    for law, h in hamiltonians.items():
        effects, error = radial_effects(collar, h)
        probabilities = np.einsum("tdbij,prji->ptdbr", effects, density).real
        pad = nextafter(
            float(incoming_error + (2 + error) * error * max(admissions) + 2e-12), inf
        )
        errors[law] = pad
        raw[law] = probabilities
        archived[law] = {
            "transfer_probabilities": probabilities.tolist(),
            "absolute_error_per_launch_upper": pad,
            "counts": sample_counts(probabilities, admissions, refusal, rng).tolist(),
        }
    np.testing.assert_allclose(
        raw["mobile_cut"], raw["neutral_cut"], atol=2e-12, rtol=0
    )
    directions = raw["mobile"] - raw["neutral"]
    radial = directions[:, :, 0] - directions[:, :, 1]
    response_error = 2 * (errors["mobile"] + errors["neutral"])
    diagnostics, sample_bin_error = radial_diagnostics(
        collar,
        raw,
        {law: bank["counts"] for law, bank in archived.items()},
        response_error,
    )
    phase_error = 2 * response_error
    difference = radial[0] - radial[1]
    idx = np.unravel_index(np.argmax(abs(difference)), difference.shape)
    phase_signal = float(difference[idx])
    t, b, y = map(int, idx)
    return (
        collar,
        hamiltonians,
        radial,
        response_error,
        {
            "source_sign": path.sign,
            "probe_sign": sign,
            "source_word_order": list(words),
            "phase_specific_head_current_given_face_and_reference": [
                currents.tolist(),
                (-currents).tolist(),
            ],
            "reference_admission_per_launch": admissions.tolist(),
            "refusal_probabilities": refusal,
            "sharp_chart_records": len(collar.basis),
            "clamped_source_scalar_energy": scalar,
            "exact_connected_near_far_energy_difference": zero_static,
            "static_echo": static_clock(collar, scalar, float(sum(admissions)), rng),
            "laws": archived,
            "simultaneous_99pct_error_per_simulated_dynamic_bin": sample_bin_error,
            "radial_mobility_score_per_launch": radial.tolist(),
            "radial_mobility_score_error_upper": response_error,
            "resolved_two_direction_phase_examples": diagnostics,
            "largest_phase_comparison": {
                "current_phase_difference_per_original_launch": phase_signal,
                "absolute_difference_error_upper": phase_error,
                "exterior_word": b,
                "collar_wait": TIMES[t],
                "reference_y": (1, -1)[y],
                "original_phase_radial_per_launch": float(radial[0, t, b, y]),
                "opposite_phase_radial_per_launch": float(radial[1, t, b, y]),
            },
        },
    )


def verify_conjugacy(a, b):
    left, matrices, _, _, _ = a
    right, opposite, _, _, _ = b
    mask = (1 << len(left.edges)) - 1
    order = [right.ids[x, p, w ^ mask] for x, p, w in left.basis]
    assert len(set(order)) == len(left.basis) == len(right.basis)
    for law in matrices:
        assert (matrices[law] != opposite[law][order][:, order]).nnz == 0


def run():
    receipt = json.loads(READS.read_text())
    late = next(
        row
        for case in receipt["cases"]
        if case["fan_count"] == 4
        for stock in case["stock_arms"]
        if stock["preparation"] == "clock_bypassed"
        for row in stock["flights"]
        if row["material_bridge_rate"] == 4 and row["appointment"] == 30.8125
    )
    original = next(
        row for row in late["pair_readings"] if row["port_face_and_words"] == LABEL
    )
    rng = np.random.default_rng(SEED)
    with archived_fixture() as archive:
        p, m = [FanPath(archive, 4, s, asymmetric=True) for s in (1, -1)]
        mask = (1 << len(p.field.edges)) - 1
        order = [m.ids[x, tuple(w ^ mask for w in words)] for x, words in p.basis]
        assert (p.h != m.h[order][:, order]).nnz == 0
        np.testing.assert_array_equal(m.source()[order], -p.source())
        comparisons = {}
        for path in (p, m):
            words = LABEL[-2:] if path.sign == 1 else [w ^ mask for w in LABEL[-2:]]
            signals = late_arrival(path, [*LABEL[:2], *words], original)
            for test_sign in (1, -1):
                comparisons[path.sign, test_sign] = report(
                    path,
                    test_sign,
                    words,
                    signals,
                    late["observable_absolute_error_upper_per_pair"] + 2e-12,
                    rng,
                )
        for sign in (1, -1):
            verify_conjugacy(comparisons[1, sign], comparisons[-1, -sign])
        collar = comparisons[1, 1][0]
        odd = sum(q * p * case[2] for (q, p), case in comparisons.items()) / 4
        odd_error = sum(case[3] for case in comparisons.values()) / 4
        return {
            "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
            "inherited_word_reads_sha256": sha256(READS.read_bytes()).hexdigest(),
            "physical_geometry": {
                "faces": [list(face) for face in collar.faces],
                "near_and_far_probe": list(collar.radial),
                "intact_hops": collar.ranges,
                "probe_degrees": collar.probe_degrees,
                "eligible_shortest_choices": collar.eligible,
                "written_exterior_seam_bit_order": [
                    list(edge) for edge in collar.added
                ],
                "isolated_active_lane_transfer_probability_each_direction": [
                    float(np.sin(4 * t) ** 2) for t in TIMES
                ],
                "isolated_blocked_lane_transfer_probability": 0,
            },
            "local_waits": list(TIMES),
            "phase_order": ["original", "word_b_phase_reversed"],
            "direction_order": ["outward", "inward"],
            "reference_y_order": [1, -1],
            "probability_array_axes": [
                "phase",
                "wait",
                "direction",
                "background",
                "reference_y",
            ],
            "radial_score_array_axes": ["phase", "wait", "background", "reference_y"],
            "raw_count_bins": [
                "outside receiver",
                "receiver off hub",
                "different word",
                "Y+ stay",
                "Y+ transfer",
                "Y- stay",
                "Y- transfer",
            ],
            "simulated_starts_per_bank": TRIALS,
            "simulated_seed": SEED,
            "signed_odd_radial_score_per_launch": odd.tolist(),
            "signed_odd_radial_score_error_upper": odd_error,
            "cases": [item[-1] for item in comparisons.values()],
        }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    data = run()
    encoded = json.dumps(data, separators=(",", ":"), allow_nan=False)
    if args.output:
        with args.output.open("x") as handle:
            handle.write(encoded + "\n")
    else:
        compare = run_path(str(HERE.parents[2] / "replay.py"))["compare"]
        compare(
            json.loads(encoded),
            json.loads((HERE / "results.json").read_text()),
            "radial_probe",
        )
        print("Verified both radial probe directions and all matched controls.")


if __name__ == "__main__":
    main()
