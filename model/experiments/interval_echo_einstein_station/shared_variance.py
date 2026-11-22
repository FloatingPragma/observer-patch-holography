"""Two written collar pairs with a common or independent physical preparation tag.

All arms have identical individual collar-pair states and the same distribution
of both stored pair energies. Free flights encounter their joint state. The
ordinary reference removes only the one stationary unit on the shared port.
"""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.experiments.contact import radar
from . import cap, refinement, volume


HERE = Path(__file__).resolve().parent
TIPS = (0, 5, 7, 8)
TAG = (6, 10)
ARMS = ("reversed_second_tag", "independent_tags", "shared_tag")


def source_words(source):
    """Joint source probabilities from two independent native writes given their tags."""
    units = np.asarray(source["source_units"], dtype=int).reshape(3, 9)
    lower, middle, upper = units
    np.testing.assert_array_equal(lower + upper, 2 * middle)
    # A common fair tag selects two independent writes given that value.
    pairs = np.stack(
        [
            np.outer(lower, upper) + np.outer(upper, lower),
            2 * np.outer(middle, middle),
            np.outer(lower, lower) + np.outer(upper, upper),
        ]
    )
    energies = np.asarray(source["source_word_energies"], dtype=int)
    directed = (upper - middle).reshape(3, 3)
    np.testing.assert_array_equal(pairs[2] - pairs[1], 2 * np.outer(directed, directed))
    # The probability denominator is 2592, so the total variation from
    # independent tags is 72/(2*2592) = 1/72. In each word, a mixed cubic
    # is the first energy monomial that can see the directed change.
    np.testing.assert_equal(np.abs(pairs[2] - pairs[1]).sum(), 72)
    for i in range(3):
        for j in range(3 - i):
            np.testing.assert_equal((energies**i) @ directed @ (energies**j), 0)
    cubic = int(energies @ directed @ (energies**2))
    if cubic == 0 or cubic != -int((energies**2) @ directed @ energies):
        raise ValueError("Missing opposite cubic moments on the oriented energy words")
    sums = np.add.outer(energies, energies).ravel()
    for arm in pairs:
        for axis in (0, 1):
            np.testing.assert_array_equal(arm.sum(axis=axis), 72 * middle)
        # Even the two separate pair-energy readings made together stay fixed.
        for a, b in product(np.unique(sums), repeat=2):
            selection = (sums == a)[:, None] & (sums == b)[None, :]
            np.testing.assert_array_equal(
                arm[selection].sum(), pairs[1][selection].sum()
            )
    return pairs / 2592, middle / 36


def shared_station(level=0, extra_sites=()):
    """Build the existing native radar with two written pairs around port 1."""
    _, _, hand_yes, source = volume.radar_source_and_hand(federation(0))
    base, _, placement = refinement.station(level, source, TAG)
    board = federation(level)
    remote = placement["serviced_bridge"][1]
    sites = [4, 9, remote, remote + 5]
    # The tagged words are prepared before the separate material hand is armed.
    # The b->9 sitter and disjoint tip spares use blank swaps, even where
    # a spare touches a previously written (parked) memory.
    written = {(0, 1), (1, 5), (1, 7), (1, 8)}
    for vertices, spares in [((0, 1, 5), (11, 9, 4)), ((7, 1, 8), (0, 9, 3))]:
        if set(spares) & (set(vertices) | set(TAG)) or len(set(spares)) != 3:
            raise ValueError("Preparation would overwrite its source or protected tag")
        if any(
            p not in board.neighbours[v] or tuple(sorted((v, p))) in written
            for v, p in zip(vertices, spares, strict=True)
        ):
            raise ValueError(
                "The ordered-word preparation needs three blank service rails"
            )
    # Existing station parks the unused seams at port 1; reinstating written
    # contacts (1,7) and (1,8) gives precisely the same blank-exchange removal.
    if any(i not in board.neighbours[1] or board.degree[i] != 5 for i in TIPS):
        raise ValueError("The four degree-five source endpoints must contact port 1")
    if any(y not in board.neighbours[x] for x, y in (TIPS[:2], TIPS[2:])):
        raise ValueError("Each source triangle needs its own blank probe splitter")
    hand = source["hand_seam"]
    for origin in (*sites, *extra_sites):
        if origin in {1, *TIPS, *hand, *TAG} or not 0 <= origin < base.shape[0]:
            raise ValueError("The detector must use an unmarked physical port")
        spares = set(range(12)) - {1, *TIPS, *hand, *TAG, origin}
        if not any(
            p != q
            for p in spares & set(board.neighbours[hand[0]])
            for q in spares & set(board.neighbours[hand[1]])
        ):
            raise ValueError("No free hand-metering rails for this launch")
    return base, sites, hand_yes, source


def word_generators(base, energies):
    for indices in product(range(3), repeat=4):
        shifts = np.zeros(base.shape[0])
        shifts[list(TIPS)] = -np.asarray(energies)[list(indices)]
        yield base + sparse.diags(shifts, format="csr")


def connected_modular_charge(weights):
    """Read the joint collar minus the two marginal modular generators."""
    shots, seed = 2**34, 2026112262
    rng = np.random.default_rng(seed)
    held = np.array(
        [rng.multinomial(shots, w.ravel()).reshape(9, 9) for w in weights[1:]]
    )
    held_difference = (held[1] - held[0]) / shots
    change = weights[2] - weights[1]
    # With the marginals fixed, all additive errors in their log scores cancel
    # in the population contrast. Bound the remaining joint-log calibration
    # error by total variation and a simultaneous bound on the reference cells.
    np.testing.assert_allclose(change.sum(axis=0), 0, atol=1e-16)
    np.testing.assert_allclose(change.sum(axis=1), 0, atol=1e-16)
    ref_eps = np.sqrt(np.log(2 * weights.size / 0.005) / (2 * shots))
    held_eps = np.sqrt(np.log(2 * 6 / 0.005) / (2 * shots))
    rows = []
    for ref in weights:
        counts = rng.multinomial(shots, ref.ravel()).reshape(9, 9)
        measured_ref = counts / shots
        if measured_ref.min() <= ref_eps:
            raise ValueError(
                "Joint log calibration requires resolved positive word weights"
            )

        def score(p):
            return np.log(np.outer(p.sum(axis=1), p.sum(axis=0)) / p)

        measured_score = score(measured_ref)
        calibration_radius = abs(change).sum() * -np.log1p(
            -ref_eps / measured_ref.min()
        )
        radius = calibration_radius + 2 * held_eps * np.ptp(measured_score)
        rows.append(
            dict(
                reference_counts=counts.tolist(),
                native_shared_minus_independent_nats=float(np.sum(change * score(ref))),
                measured_shared_minus_independent_nats=float(
                    np.sum(held_difference * measured_score)
                ),
                simultaneous_99pct_radius_nats=float(radius),
            )
        )
    return dict(
        shots_per_setting=shots,
        seed=seed,
        held_counts=held.tolist(),
        reference_order=["reversed_second_tag", "independent_tags", "shared_tag"],
        rows=rows,
        scope="Sharp commuting collar-word records under three independently prepared faithful references. Identification of their connected log charge with the physical boost remainder requires a reference and a stress conversion on the same patch.",
    )


def launch_acceleration(base, source, weights):
    """Launch two unmeasured tracers and read their relative angular coordinate."""
    n = base.shape[0]
    directions = cap.chord_coordinates()
    coordinate = np.zeros(n)
    coordinate[:12] = directions @ directions[0]
    edges = ((5, 4), (11, 2))
    starts = np.zeros((n, len(edges)))
    for k, (i, j) in enumerate(edges):
        if not np.isclose(base[i, j], -0.5) or federation(0).degree[j] != 5:
            raise ValueError(
                "The tracer launch requires a native split and a degree phase"
            )
        # A native half splitter gives (1,i)/sqrt(2); a degree phase -pi/2
        # on the second rail prepares the real, initially zero-current launch.
        starts[[i, j], k] = 1 / np.sqrt(2)

    def measurements(h):
        velocity = h @ starts
        accel = 2 * np.sum(
            coordinate[:, None] * (velocity**2 - starts * (h @ velocity)), axis=0
        )
        # Three snapshots on fresh preparations; all outcomes outside the patch
        # enter the outside count and contribute zero to this patch coordinate.
        flight = expm_multiply(-1j * h, starts, start=0, stop=1 / 32, num=3)
        patch = abs(flight[:, :12, :]) ** 2
        outside = np.maximum(0, 1 - patch.sum(axis=1, keepdims=True))
        return accel, np.concatenate([patch, outside], axis=1)

    blank, blank_reads = measurements(base)
    acc, reads = zip(
        *(
            measurements(h)
            for h in word_generators(base, source["source_word_energies"])
        ),
        strict=True,
    )
    mixing = weights.reshape(-1, 81)
    expected = np.vstack([blank, mixing @ np.array(acc)])
    probabilities = np.concatenate(
        [blank_reads[None, ...], np.einsum("aw,wtps->atps", mixing, reads)], axis=0
    )
    np.testing.assert_allclose(
        expected[1:], np.broadcast_to(expected[2], (3, 2)), atol=3e-12, rtol=0
    )
    analytic = (
        0.5 * np.mean(source["source_word_energies"]) * (coordinate[5] - coordinate[4])
    )
    np.testing.assert_allclose(expected[2] - blank, [analytic, 0], rtol=0, atol=4e-12)
    shots, seed = 2**40, 2026112217
    rng = np.random.default_rng(seed)
    times = np.array([0.0, 1 / 64, 1 / 32])
    gap = source["hand_gap"]
    p_clock = (
        1 + np.column_stack([np.cos(gap * times[1:]), np.sin(gap * times[1:])])
    ) / 2
    clock_reads = rng.binomial(shots, p_clock)
    xy = 2 * clock_reads / shots - 1
    hand_time = np.r_[0, np.arctan2(xy[:, 1], xy[:, 0]) / gap]

    # Extrapolate even position means through fourth order; report the finite
    # stencil bias against the independently computed exact launch derivative.
    def stencil(t):
        return np.linalg.solve(np.array([np.ones(3), t**2, t**4]), [0, 2, 0])

    scored = np.r_[coordinate[:12], 0]
    rows = []
    for label, exact, menu in zip(
        ("vacant", *ARMS), expected, probabilities, strict=True
    ):
        counts = np.array(
            [
                [
                    rng.multinomial(shots, radar.acquisition_weights(p) / radar.LATTICE)
                    for p in snapshot.T
                ]
                for snapshot in menu
            ]
        )
        means = counts @ scored / shots
        coefficients = stencil(hand_time)
        measured = coefficients @ (means[:, 0] - means[:, 1])
        mean_variances = (counts @ scored**2 / shots - means**2) / (shots - 1)
        conditional_se = np.sqrt(coefficients**2 @ mean_variances.sum(axis=1))
        ideal = np.einsum("tps,p->ts", menu, scored)
        finite = stencil(times) @ (ideal[:, 0] - ideal[:, 1])
        rows.append(
            dict(
                arm=label,
                tracer_accelerations=exact.tolist(),
                native_relative_acceleration=float(exact[0] - exact[1]),
                native_finite_stencil=float(finite),
                measured_finite_stencil=float(measured),
                port_sampling_se_given_hand_times=float(conditional_se),
                patch_and_outside_counts=counts.tolist(),
            )
        )
    return dict(
        tracer_launch_rails=edges,
        angular_coordinate="Home slow-triplet dot product with port 0, zero score on the separately recorded outside event.",
        shots_per_setting=shots,
        seed=seed,
        native_exposures=times.tolist(),
        independent_hand_times=hand_time.tolist(),
        yes_counts_cos_sin=clock_reads.tolist(),
        ordinary_relative_acceleration=analytic,
        common_history_exact_acceleration_change=0,
        rows=rows,
        scope="Held angular projections of separately released at-rest-in-mean tracers. A physical geodesic-deviation or curvature interpretation requires additional metric and freely falling tracer calibration.",
    )


def source_read(weights, probabilities, shots, rng):
    """Two separated material probes, one per retained two-seam word."""
    conditional = np.column_stack((probabilities, 1 - probabilities))
    joint = conditional.T @ weights @ conditional
    centre = joint[0].sum()
    scores = np.outer([1 - centre, -centre], [1 - centre, -centre])
    counts = rng.multinomial(shots, joint.ravel()).reshape(2, 2)
    return dict(
        complete_counts=counts.tolist(),
        single_probe_return=float(centre),
        native_centred_product=float(np.sum(joint * scores)),
        measured_centred_product=float(np.sum(counts * scores) / shots),
        simultaneous_radius=float(np.ptp(scores) * np.sqrt(np.log(800) / (2 * shots))),
    )


def sample_outward(rate, shots, rng):
    # Compact binomial draws in independent batches retain exact quiet and
    # arrival totals. Python's integer sum is unbounded; JSON strings avoid
    # a loss of precision on readers using IEEE doubles.
    total = 0
    for _ in range(shots // (2**42 * 4096)):
        total += int(rng.binomial(2**42, rate, size=4096).sum(dtype=np.int64))
    return total


def waveform_peaks(base, sites, weights, waveforms):
    """Use the same whole-window maximum for every arm and both directions."""
    rays = [ray for i, j in refinement.PAIRS for ray in ((i, j), (j, i))]
    mixed = np.einsum("aw,wmijt->amijt", weights.reshape(len(weights), -1), waveforms)
    onset = [float(abs(base[sites[j], sites[i]]) ** 2) for i, j in rays]
    if min(onset) <= 0:
        raise ValueError("The quadratic onset formula here requires a direct bond")
    result = dict(
        directed_ports=[[sites[i], sites[j]] for i, j in rays],
        exposure_per_bin=radar.DT,
        number_of_bins=radar.TICKS,
        terminal_short_time_probability_over_time_squared=onset,
        scope="Population forecasts. First-hit probabilities include earlier quiet appointments. Terminal probabilities use separate freely evolving launches at each appointment. Peak locations resolve only this supplied grid.",
    )
    for method, data in zip(
        ("first_hit", "unmonitored_terminal"),
        mixed.transpose(1, 0, 2, 3, 4),
        strict=True,
    ):
        curves = np.stack([data[:, i, j] for i, j in rays], axis=1)
        top = np.sort(curves, axis=-1)[..., -2:]
        result[method] = dict(
            peak_bins_by_arm=(1 + curves.argmax(axis=-1)).tolist(),
            minimum_peak_margin_over_arms=(top[..., 1] - top[..., 0])
            .min(axis=0)
            .tolist(),
            maximum_abs_shared_minus_independent_probability_by_ray=np.max(
                abs(curves[2] - curves[1]), axis=-1
            ).tolist(),
        )
    return result


def run():
    base, sites, hand_yes, source = shared_station()
    weights, reference = source_words(source)
    full_strength = weights
    # Halve each conditional signed shift with independent local randomization.
    # Averaging over the common tag leaves exactly a quarter of the joint contrast.
    attenuated = (3 * weights[1] + weights[2]) / 4
    native_words = np.asarray(source["source_units"]).reshape(3, 9) / 36
    halves = (native_words[[0, 2]] + native_words[1]) / 2
    np.testing.assert_allclose(attenuated, halves.T @ halves / 2, atol=1e-16, rtol=0)
    weights = np.concatenate((weights, attenuated[None]))
    gaps = np.asarray(source["source_word_energies"])
    score = np.log(1 / (9 * reference))
    score -= reference @ score
    amplitudes = [float(score @ w @ score) for w in weights]
    unit = np.log(3) / 12
    np.testing.assert_allclose(
        amplitudes, [-(unit**2), 0, unit**2, unit**2 / 4], atol=2e-16, rtol=0
    )
    # The degree-five quadrature and half-splitter on each triangle are the
    # native operations certified by ordered_collar_control. Their two probe
    # paths have disjoint tips and see commuting phases at the fixed centre.
    gap_difference = np.subtract.outer(gaps, gaps).ravel()
    probabilities = (1 + np.sin(gap_difference * np.pi / 3)) / 2
    enrolled = json.loads((HERE / "backreaction.json").read_text())[
        "correlated_contacts"
    ]["ordered_three_levels"]
    for arm in enrolled["arms"]:
        native = np.array(arm["source_units_of_1_over_36"]).ravel() @ probabilities / 36
        np.testing.assert_allclose(
            native, arm["readers"]["quadrature"]["probabilities"][0], atol=2e-13, rtol=0
        )
    rng = np.random.default_rng(2026112255)
    probes = [source_read(w, probabilities, 2**22, rng) for w in weights]
    np.testing.assert_allclose(
        [a["native_centred_product"] for a in probes],
        [-1 / 432, 0, 1 / 432, 1 / 1728],
        atol=1e-15,
    )
    joint, enrolled_outward, waveforms = [], [], []
    max_leg_disagreement, tail = 0.0, 0.0
    for h in word_generators(base, gaps):
        legs, error = refinement.projector_check(h, sites)
        renewal, terminal = refinement.sparse_hit_laws(h, sites, return_intensity=True)
        waveforms.append((legs[..., :-1], terminal))
        max_leg_disagreement = max(
            max_leg_disagreement, float(np.max(abs(legs - renewal)))
        )
        tail = max(tail, error)
        joint.append(radar.round_trip_law(legs))
        enrolled_outward.append(renewal[0, 1, :-1].sum())
    if max_leg_disagreement > 1e-11:
        raise ValueError("Independently computed first-hit probabilities disagree")
    by_word = np.asarray(joint)
    fields = np.einsum(
        "sij,ijxyab->sxyab", weights, by_word.reshape(9, 9, *by_word.shape[1:])
    )
    clock = (
        np.unwrap(np.arctan2(2 * hand_yes[1] - 1, 2 * hand_yes[0] - 1))
        / source["hand_gap"]
    )
    half_hand = clock[np.add.outer(np.arange(radar.TICKS), np.arange(radar.TICKS))] / 2
    vacant = refinement.summarize(by_word[0], half_hand)
    forecasts = [refinement.summarize(law, half_hand) for law in fields]
    # Select the outbound event for acquisition after exploratory forecasts of
    # this dyad. Count fresh simulations of every launch, including timeouts.
    shots = 2**68
    counts, rates = [], []
    acquisition_rates = weights.reshape(-1, 81) @ enrolled_outward
    for field, acquisition_rate in zip(fields, acquisition_rates, strict=True):
        rate = float(1 - field[0, 1, -1, -1])
        rates.append(rate)
        bins = radar.acquisition_weights([acquisition_rate, 1 - acquisition_rate])
        hits = sample_outward(bins[0] / radar.LATTICE, shots, rng)
        counts.append({"arrived": str(hits), "timeout": str(shots - hits)})
    np.testing.assert_allclose(acquisition_rates, rates, rtol=0, atol=1e-13)
    native_shifts = [r - rates[1] for r in rates]
    measured_shifts = [
        (int(row["arrived"]) - int(counts[1]["arrived"])) / shots for row in counts
    ]
    # Bound exact-arithmetic Chebyshev truncation and keep interval conditional
    # on the high-precision numerical probability table and supplied preparation.
    log_terms = np.log(8 / 0.01)
    radius = 2 * np.sqrt(log_terms / (2 * shots)) + 1 / radar.LATTICE
    return dict(
        cells=20,
        source_tip_pairs=[TIPS[:2], TIPS[2:]],
        stationary_ordinary_port=1,
        tag_rails=list(TAG),
        hand_seam=source["hand_seam"],
        radar_dyads=[sites[:2], sites[2:]],
        arms=[*ARMS, "half_strength_each_write"],
        outward_probability_lattice=radar.LATTICE,
        word_law="Two conditionally independent categorical preparations, each using the three-level native writer. The second preparation uses the opposite, a fresh independent, or the same fair tag. No tag is read by the flying detector.",
        supplied="Conditional preparations, physical tag protection, stationary central load, detector appointments and reference. Protected source words survive both flight legs.",
        between_word_score_covariance_nats_squared=amplitudes,
        mutual_information_nats=[
            float(np.sum(w * np.log(w / weights[1]))) for w in weights
        ],
        shared_small_strength_quartic_coefficient=float(
            0.5 * np.sum((weights[2] - weights[1]) ** 2 / weights[1])
        ),
        full_shared_vs_fresh_total_variation=1 / 72,
        ramsey_shots=2**22,
        seed=2026112255,
        two_contact_probe_measurements=probes,
        vacant_clock_law=vacant,
        occupied_clock_laws=forecasts,
        waveform_peaks=waveform_peaks(base, sites, weights, waveforms),
        maximum_leg_disagreement=max_leg_disagreement,
        chebyshev_half_step_truncation=tail,
        outward_enrollments_per_arm=str(shots),
        outward_batch_size=2**42,
        outward_sampling="Exploratory 4-to-9 outbound event selected before independent seeded acquisition. These are synthetic binomial samples of the renewed first-hit table. The simultaneous radius is conditional on this event, ideal binomial sampling and the represented probability table.",
        outward_counts=counts,
        outward_rates=rates,
        maximum_outward_renewal_disagreement=float(
            np.max(np.abs(acquisition_rates - rates))
        ),
        outward_changes_vs_independent=native_shifts,
        measured_outward_changes_vs_independent=measured_shifts,
        simultaneous_99pct_outward_change_radius=float(radius),
        launch_acceleration=launch_acceleration(base, source, full_strength),
        connected_modular_charge=connected_modular_charge(full_strength),
        scope="Native conditional duration forecasts and separate simulated tag-blind outward detector counts. The central load and two-word energy distribution are held fixed across the four occupied arms. Flight scattering and launch acceleration are distinct finite observables; neither has an independently established curvature calibration.",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "shared_variance.json"
    result = json.loads(json.dumps(run(), allow_nan=False))
    if parser.parse_args().write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        expected = json.loads(path.read_text())
        volume.cap.compare(result, expected)
        np.testing.assert_allclose(
            result["outward_changes_vs_independent"],
            expected["outward_changes_vs_independent"],
            rtol=0,
            atol=1e-13,
        )
        for method in ("first_hit", "unmonitored_terminal"):
            key = "maximum_abs_shared_minus_independent_probability_by_ray"
            np.testing.assert_allclose(
                result["waveform_peaks"][method][key],
                expected["waveform_peaks"][method][key],
                rtol=0,
                atol=1e-13,
            )
        measured, saved = (
            r["occupied_clock_laws"][2][0]["conditional_half_hand"]
            - r["occupied_clock_laws"][1][0]["conditional_half_hand"]
            for r in (result, expected)
        )
        np.testing.assert_allclose(measured, saved, rtol=0, atol=1e-13)
    print(
        result["outward_changes_vs_independent"],
        result["measured_outward_changes_vs_independent"],
        result["simultaneous_99pct_outward_change_radius"],
    )
