"""Replay neutral-reference phase and held arrival at the archived 14–24 gateway."""

import argparse
from base64 import b64decode, b64encode
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from zlib import compress, decompress

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
FIELD = HERE.parent
sys.path.insert(0, str(FIELD))

DT = 1 / 32
SOURCE_SIZE = 4096
PAIR_SIZE = 3072
CLOCK_SIZE = 768
PHASE_SHOTS = 2**26
ARRIVAL_SHOTS = 2**22
RESOLUTION = 0.002
COVERAGE = 0.99
ARRIVAL_OUTCOMES = 7
GROUND_BRACKET = (Fraction(3547, 1000), Fraction(3549, 1000))
HIGHER_FLOOR = Fraction(93, 20)
LABELS = ("vacuum", "first_outer_odd", "last_outer_odd")
SEED = 142406
COLLECTIVE = "experiments/interval_echo_collective_spin/"
ADJUST = {"vacuum": 0, "first_outer_odd": 1, "last_outer_odd": 1}


class Flight:
    def __init__(self, h):
        self.levels, self.vectors = eigh(h)

    def evolve(self, t, source):
        return self.vectors @ (
            np.exp(-1j * t * self.levels) * (self.vectors.T @ source)
        )

    def trace_against_launch(self, t, source):
        weights = abs(self.vectors.T @ source) ** 2
        return np.exp(-1j * np.outer(t, self.levels)) @ weights

    def release_filter(self, time, frequency, source):
        factors = np.mean(np.exp(-1j * np.outer(time, self.levels - frequency)), axis=0)
        return self.vectors @ (factors * (self.vectors.T @ source))


def phase_bank(contrast, rng):
    quadratures = np.column_stack((contrast.real, contrast.imag))
    counts = rng.binomial(PHASE_SHOTS, np.clip((1 + quadratures) / 2, 0, 1))
    measured = (2 * counts[:, 0] / PHASE_SHOTS - 1) + 1j * (
        2 * counts[:, 1] / PHASE_SHOTS - 1
    )
    return measured, {
        "settings": len(contrast),
        "plus_counts": b64encode(compress(counts.astype("<u4").tobytes())).decode(),
    }


def peaks(contrast):
    time = np.arange(len(contrast)) * DT
    window = np.hanning(len(time))
    fringe = contrast * window / window.sum()
    step = 2 * np.pi / (len(time) * DT)
    scan = abs(np.fft.ifft(fringe)) * len(time)
    candidates, _ = find_peaks(scan[: len(scan) // 2], height=0.002)
    fitting = []
    for index in candidates:
        opt = minimize_scalar(
            lambda frequency: -abs(fringe @ np.exp(1j * time * frequency)),
            bounds=((index - 0.75) * step, (index + 0.75) * step),
            method="bounded",
            options={"xatol": 1e-9},
        )
        assert opt.success
        fitting.append(
            {
                "fourier_bin": int(index),
                "angular_frequency": float(opt.x),
                "amplitude": float(-opt.fun),
            }
        )
    return fitting


def reads(contrast, rng):
    measured, record = phase_bank(contrast, rng)
    lines = peaks(measured)
    return record | {"observed_lines": lines}


def first_prominent(record):
    return next(line for line in record["observed_lines"] if line["amplitude"] > 0.01)


def infer_isolated_frequency(
    bank, bracket, excited_floor, spectral_ceiling, eligible_bins
):
    """Use two disjoint weighted phase reads with a simultaneous count bound."""
    counts = np.frombuffer(
        decompress(b64decode(bank["plus_counts"])), dtype="<u4"
    ).reshape(bank["settings"], 2)
    contrast = 2 * counts / PHASE_SHOTS - 1
    contrast = contrast[:, 0] + 1j * contrast[:, 1]
    n = len(contrast) // 2
    frequency = 2 * np.pi * first_prominent(bank)["fourier_bin"] / (2 * n * DT)
    wait = n * DT
    window = np.hanning(n)
    window /= window.sum()
    # One Hoeffding bound for each real/imaginary component, half, sign and bin.
    count_bound = np.sqrt(
        4 * (window @ window) * np.log(8 * eligible_bins / (1 - COVERAGE)) / PHASE_SHOTS
    )
    leakage = 0.0
    if excited_floor is not None:
        min_angle = (excited_floor - frequency) * DT
        assert 0 < min_angle < np.pi
        hann_angle = 2 * np.pi / (n - 1)
        assert (spectral_ceiling - frequency) * DT + hann_angle < np.pi
        denominators = 2 * np.sin(
            (min_angle + np.array([0, hann_angle, -hann_angle])) / 2
        )
        assert min(denominators) > 0
        leakage = 4 * (1 - np.cos(hann_angle)) / ((n - 1) * np.prod(denominators))
    demodulated = contrast * np.exp(1j * frequency * DT * np.arange(2 * n))
    first, last = demodulated.reshape(2, n) @ window
    eps = count_bound + leakage
    assert eps < min(abs(first), abs(last))
    angle_error = sum(np.arcsin(eps / abs(z)) for z in (first, last))
    assert angle_error < np.pi / 4
    # The bracket certifies no branch cut and positive symmetric Hann response.
    assert max(abs(frequency - edge) for edge in bracket) * wait + angle_error < np.pi
    assert max(abs(frequency - edge) for edge in bracket) * (n - 1) * DT / 2 < np.pi / 2
    value = frequency - np.angle(last * first.conjugate()) / wait
    uncertainty = angle_error / wait
    return {
        "frequency": float(value),
        "confidence_interval": [float(value - uncertainty), float(value + uncertainty)],
        "two_half_count_error_bound": float(count_bound),
        "off_line_leakage_bound": float(leakage),
        "minimum_measured_half_window_contrast": float(min(abs(first), abs(last))),
    }


def two_copy_insertion(fans, local_field, neutral_sectors):
    """Verify the declared feeder's load and written-head changes, input by input."""
    from threshold import Instrument, local_columns

    free = local_columns(Instrument())
    states = [(0, (1, 1)), (1, (1, 1)), (1, (2, 0)), (0, (0, 2))]
    for a, b in ((states[0], states[2]), (states[1], states[3])):
        assert a[0] - a[1][0] == b[0] - b[1][0] and sum(a[1]) == sum(b[1])
        assert all((coin, tuple(load)) != b for coin, load, _ in free[a])
        assert all((coin, tuple(load)) != a for coin, load, _ in free[b])
    neutral_dimension = 1
    for sign, space, sector, feeder in zip(
        (-1, 1), fans, neutral_sectors, (1, 0), strict=True
    ):
        neutral_dimension *= space.portal[0].size
        bridge_head = 6 if sign == -1 else 0
        counts = [
            np.bincount(
                [a if bit & word else b for (a, b), bit in local_field.bits.items()],
                minlength=6,
            )
            for word in local_field.words[sector]
        ]
        for i, left in enumerate(counts):
            for j, right in enumerate(counts):
                expected = np.concatenate((left, right))
                expected[bridge_head] += 1
                position, written = space.basis[space.portal[0][i, j]]
                observed = np.bincount(space.heads(written), minlength=12)
                assert position == 0
                np.testing.assert_array_equal(observed, expected)
                expected[0] += feeder
                observed[0] += 1 - feeder
                assert observed[0] - expected[0] == sign
                np.testing.assert_array_equal(observed[1:], expected[1:])
    # A declared 4(I-X) preparation contact gives a pure swap at pi/8.
    levels, basis = eigh(np.array([[4, -4], [-4, 4]]))
    pulse = (basis * np.exp(-1j * np.pi * levels / 8)) @ basis.T
    np.testing.assert_allclose(pulse, [[0, 1], [1, 0]], rtol=0, atol=1e-14)
    return neutral_dimension


def multinomial_marginal_intervals(counts, settings):
    """Simultaneous two-sided exact binomial intervals for the listed categories."""
    tail = (1 - COVERAGE) / (2 * settings * ARRIVAL_OUTCOMES)
    return [
        [
            float(beta.ppf(tail, k, ARRIVAL_SHOTS - k + 1)) if k else 0.0,
            float(beta.ppf(1 - tail, k + 1, ARRIVAL_SHOTS - k))
            if k < ARRIVAL_SHOTS
            else 1.0,
        ]
        for k in counts
    ]


def arrivals(space, flight, source, code, vacuum, frequency, rng, archive_times):
    joint = np.zeros((len(source), len(vacuum) * code.shape[1]))
    joint[space.portal[1].ravel()] = np.kron(np.eye(len(vacuum)), code)
    shared_vacuum = np.zeros(len(source))
    shared_vacuum[space.portal[1].ravel()] = np.kron(vacuum, vacuum)
    np.testing.assert_array_equal(joint.T @ shared_vacuum, 0)
    right_port = space.portal[1].ravel()
    right = np.array([port >= 6 for port, _ in space.basis])
    initial = np.array([port == 0 for port, _ in space.basis])
    assert len(right_port) == 16 and np.count_nonzero(initial) == 16
    filter_times = np.arange(SOURCE_SIZE) * DT
    filtered = flight.release_filter(filter_times, frequency, source)
    result = []
    outcome_order = [
        "filter_declined",
        "source_gateway_charge",
        "other_left_charge",
        "right_nongateway_charge",
        "right_gateway_both_field_vacua",
        "right_gateway_local_doublet",
        "right_gateway_other_field",
    ]
    assert len(outcome_order) == ARRIVAL_OUTCOMES
    settings = 2 * len(LABELS) * (1 + len(archive_times))
    for name, launch in (("unfiltered", source), ("clock_filtered", filtered)):
        for time in [0, *archive_times]:
            state = flight.evolve(time, launch)
            total = np.vdot(state, state).real
            left_gateway = np.vdot(state[initial], state[initial]).real
            arriving = np.vdot(state[right], state[right]).real
            gateway = np.vdot(state[right_port], state[right_port]).real
            orbital = np.sum(abs(joint.T @ state) ** 2)
            no_field = abs(shared_vacuum @ state) ** 2
            bins = np.array(
                [
                    1 - total,
                    left_gateway,
                    total - left_gateway - arriving,
                    arriving - gateway,
                    no_field,
                    orbital,
                    gateway - no_field - orbital,
                ]
            )
            assert min(bins) > -1e-12
            bins = np.clip(bins, 0, 1)
            bins /= bins.sum()
            counts = rng.multinomial(ARRIVAL_SHOTS, bins).tolist()
            result.append(
                {
                    "read": name,
                    "delay_in_cadence_steps": round(time / DT),
                    "probabilities": bins.tolist(),
                    "counts": counts,
                    "simultaneous_marginal_intervals": multinomial_marginal_intervals(
                        counts, settings
                    ),
                }
            )
    return {"outcome_order": outcome_order, "settings": result}


def run_fixture(field, original, degrees, ports, archive_times):
    from binding_audit.measure import count_below
    from search import JoinedSpace, vacuum_code

    complement = tuple(
        sum(i in edge for edge in field.edges) - h for i, h in enumerate(original)
    )
    vacancy = JoinedSpace(field, (complement, complement), -1)
    excess = JoinedSpace(field, (original, original), 1)
    guides = (-1, 1)
    h, h_plus = [space.hamiltonian(degrees) for space in (vacancy, excess)]
    assert h.shape == (160, 160)
    mask = (1 << len(vacancy.edges)) - 1
    mapping = [excess.ids[port, word ^ mask] for port, word in vacancy.basis]
    np.testing.assert_array_equal(h, h_plus[np.ix_(mapping, mapping)])
    neutral_insertions = two_copy_insertion(
        (vacancy, excess), field, (complement, original)
    )
    vacuum, odd = vacuum_code(field, complement)
    vacuum_plus, _ = vacuum_code(field, original)
    neutral_h = field.ring_laplacian(field.words[complement])
    clock = Flight(neutral_h)
    diagonal_positive = h_plus.copy()
    diagonal_positive[diagonal_positive == -4] = 0
    diagonal_control = Flight(diagonal_positive)
    moving = Flight(h)
    moving_plus = Flight(h_plus)
    abs_during_flight = max(abs(h).sum(axis=1).max(), abs(h_plus).sum(axis=1).max())
    # Both signed spectra are nonnegative, and the joint band lies below Nyquist.
    assert 2 * abs_during_flight + 1 < np.pi / DT
    inertia = [count_below(h, cut) for cut in (*GROUND_BRACKET, HIGHER_FLOOR)]
    assert inertia == [0, 1, 1]
    control_source = np.zeros(len(h_plus))
    control_source[excess.portal[0].ravel()] = np.kron(vacuum_plus, vacuum_plus)
    # Closing all compensator matter contacts leaves its launch stationary.
    parked_plus = np.kron(
        field.ring_laplacian(field.words[original]), np.eye(4)
    ) + np.kron(np.eye(4), field.ring_laplacian(field.words[original]))
    np.testing.assert_array_equal(
        parked_plus @ control_source[excess.portal[0].ravel()], 0
    )
    rng = np.random.default_rng(SEED)
    clock_t = DT * np.arange(CLOCK_SIZE)
    # The two read outcomes are the X and Y projections of (vacuum + odd)/sqrt(2).
    clock_samples = []
    for i in range(2):
        start = (vacuum + odd[:, i]) / np.sqrt(2)
        phases = np.array(
            [
                2 * np.conj(vacuum @ state) * (odd[:, i] @ state)
                for state in (clock.evolve(t, start) for t in clock_t)
            ]
        )
        clock_samples.append(reads(phases, rng))
    source_t = DT * np.arange(SOURCE_SIZE)
    hold_t = DT * np.arange(PAIR_SIZE)
    comp_moving = moving_plus.trace_against_launch(hold_t, control_source)
    comp_control = diagonal_control.trace_against_launch(hold_t, control_source)
    source_data = {}
    held_arrivals = {}
    sources = {}
    for label, local in zip(LABELS, [vacuum, *odd.T], strict=True):
        charge = np.zeros(len(h))
        charge[vacancy.portal[0].ravel()] = np.kron(local, vacuum)
        sources[label] = charge
        neutral_ref = np.array(
            [np.conj(local @ clock.evolve(t, local)) for t in source_t]
        )
        source_data[label] = reads(
            moving.trace_against_launch(source_t, charge) * neutral_ref, rng
        )
    measured_single = first_prominent(source_data["vacuum"])["angular_frequency"]
    assert min(GROUND_BRACKET) < measured_single < max(GROUND_BRACKET)
    held_single_amplitude = moving.trace_against_launch(hold_t, sources["vacuum"])
    held_pair = reads(held_single_amplitude * comp_moving, rng)
    held_diagonal = reads(held_single_amplitude * comp_control, rng)
    measured_pair = first_prominent(held_pair)["angular_frequency"]
    measured_diagonal = first_prominent(held_diagonal)["angular_frequency"]
    clock_first, clock_last = [
        first_prominent(c)["angular_frequency"] for c in clock_samples
    ]
    observed_ratio = measured_pair / clock_last
    predicted_ratio = 2 * measured_single / clock_first
    control_ratio = measured_diagonal / clock_last
    assert abs(observed_ratio - predicted_ratio) < RESOLUTION
    assert abs(control_ratio - predicted_ratio) > RESOLUTION
    eligible_bins = (2 * CLOCK_SIZE + SOURCE_SIZE + PAIR_SIZE) // 2
    clock_cert, opposite_cert = [
        infer_isolated_frequency(bank, (0.99, 1.01), None, 2, eligible_bins)
        for bank in clock_samples
    ]
    single_cert = infer_isolated_frequency(
        source_data["vacuum"],
        GROUND_BRACKET,
        HIGHER_FLOOR,
        abs_during_flight,
        eligible_bins,
    )
    pair_cert = infer_isolated_frequency(
        held_pair,
        [2 * bound for bound in GROUND_BRACKET],
        HIGHER_FLOOR + min(GROUND_BRACKET),
        2 * abs_during_flight,
        eligible_bins,
    )
    a, b = single_cert["confidence_interval"]
    c, d = clock_cert["confidence_interval"]
    e, f = pair_cert["confidence_interval"]
    g, h_clock = opposite_cert["confidence_interval"]
    difference_interval = [e / h_clock - 2 * b / c, f / g - 2 * a / d]
    assert max(map(abs, difference_interval)) < RESOLUTION
    for label, state in sources.items():
        held_arrivals[label] = arrivals(
            vacancy,
            moving,
            state,
            odd,
            vacuum,
            measured_single,
            rng,
            archive_times,
        )
    # Add a half unit to each defect. Terms are scalar within each
    # charged block, and the neutral arm and both ring clocks are unchanged.
    local_shift = 0.5
    shifted = Flight(h + local_shift * np.eye(len(h)))
    for label, state in sources.items():
        for time in archive_times:
            np.testing.assert_allclose(
                shifted.evolve(time, state),
                np.exp(-1j * local_shift * time) * moving.evolve(time, state),
                atol=5e-13,
                rtol=0,
            )
    shifted_pair = reads(
        held_single_amplitude * comp_moving * np.exp(-2j * local_shift * hold_t),
        rng,
    )
    shifted_pair_frequency = first_prominent(shifted_pair)["angular_frequency"]
    assert abs(shifted_pair_frequency - measured_pair - 2 * local_shift) < RESOLUTION
    return {
        "free_ports": ports,
        "degrees_in_port_order": degrees,
        "signed_loads": guides,
        "neutral_background_heads": [complement, original],
        "free_configurations_per_signed_copy": len(h),
        "neutral_basis_at_insertion": neutral_insertions,
        "feeder_gate_ports_by_copy": {
            "vacancy_copy": ports[0][0],
            "excess_copy": ports[0][0],
        },
        "native_neutral_to_adjacent_pair_contact_columns": 0,
        "free_pair_configurations": len(h) ** 2,
        "free_pair_frequency_row_sum_bound": 2 * abs_during_flight,
        "exact_released_vacancy_inertia": {
            "cuts": [str(cut) for cut in (*GROUND_BRACKET, HIGHER_FLOOR)],
            "below_cut": inertia,
        },
        "source_phase_banks_with_parked_positive": source_data,
        "independent_neutral_clock_phase_banks": clock_samples,
        "held_pair_phase_bank_with_moving_positive": held_pair,
        "diagonal_only_positive_control_bank": held_diagonal,
        "pair_comparison": {
            "single_vacancy_frequency_in_first_clock_units": measured_single
            / clock_first,
            "predicted_pair_frequency_in_second_clock_units": predicted_ratio,
            "held_pair_frequency_in_second_clock_units": observed_ratio,
            "diagonal_only_control_frequency_in_second_clock_units": control_ratio,
            "held_minus_prediction": observed_ratio - predicted_ratio,
            "control_minus_prediction": control_ratio - predicted_ratio,
            "resolution": RESOLUTION,
        },
        "simultaneous_half_window_frequency_bounds": {
            "first_neutral_clock": clock_cert,
            "second_neutral_clock": opposite_cert,
            "parked_compensator_vacancy": single_cert,
            "moving_compensator_pair": pair_cert,
            "pair_clock_ratio_minus_prediction": difference_interval,
            "joint_coverage_at_least": COVERAGE,
        },
        "defect_count_energy_shift": {
            "shift_per_defect": local_shift,
            "same_fixed_defect_flight_probabilities_for_all_delays": True,
            "observed_moving_pair_frequency_shift": shifted_pair_frequency
            - measured_pair,
            "shifted_phase_bank": shifted_pair,
        },
        "held_arrival_banks": held_arrivals,
        "joint_arrival_interval_coverage_at_least": COVERAGE,
    }


def run():
    manifest = dict(
        line.split("  ", maxsplit=1)[::-1]
        for line in (FIELD / "records.sha256").read_text().splitlines()
    )
    model_root = FIELD.parents[4]
    sys.path.insert(0, str(model_root))
    with TemporaryDirectory() as work, ZipFile(FIELD / "records.zip") as archive:
        root = Path(work)
        assert sorted(archive.namelist()) == sorted(manifest)
        assert all(
            sha256(archive.read(name)).hexdigest() == digest
            for name, digest in manifest.items()
        )
        (root / model_root.name).symlink_to(model_root, target_is_directory=True)
        archive.extractall(root)
        package = root / COLLECTIVE
        sys.path.insert(0, str(package))
        from search import fan

        record = json.loads((package / "results.json").read_text())
        choice = record["headline"]["best_worst_qubit_success"]
        release = choice["release"]
        assert release["sign"] == 1
        background = tuple(release["head_sectors"][0])
        assert tuple(release["head_sectors"][1]) == background
        ports = [choice["actual_ports"][side] for side in ("left", "right")]
        assert [side[0] for side in ports] == [14, 24]
        times = [release["best_receiver_only_time"], release["best_worst_time"]]
        assert times == [27.0625, 30.8125]
        result = run_fixture(
            fan(4),
            background,
            [6, *choice["degrees"][0], 6, *choice["degrees"][1]],
            ports,
            times,
        )
    return {
        "archived_sources_sha256": {
            file: manifest[file]
            for file in [COLLECTIVE + "search.py", COLLECTIVE + "results.json"]
        },
        "cadence_scheduled_native_units": DT,
        "phase_trials_per_quadrature": PHASE_SHOTS,
        "arrival_trials_per_setting": ARRIVAL_SHOTS,
        "seed": SEED,
        "packed_phase_count_encoding": "base64 zlib little-endian uint32, real then imaginary plus count for each setting",
        "fixture": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as out:
            out.write(result)
    else:
        from replay import compare

        compare(
            json.loads(result),
            json.loads((HERE / "records.json").read_text()),
            "fan_recoil",
        )
        print("Reproduced paired-clock phase and held arrival counts.")


if __name__ == "__main__":
    main()
