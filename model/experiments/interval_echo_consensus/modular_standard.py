"""Compile an experimental seam law from the stationary records of a finite comparison service.

Run from the repository root with python -m model.experiments.interval_echo_consensus.modular_standard.
Use --write to regenerate its compact results. The active candidate remains independently specified.
"""

from argparse import ArgumentParser
from collections import Counter
from decimal import Decimal, localcontext
import base64
import hashlib
import json
import zlib
from math import pi
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo.openwave.submission import row_04_check as pair
from model.executable_model import cell, federation
from . import offer_service, prospective_offer_capture
from .standard_observers import observer_tests

HERE = Path(__file__).resolve().parent
OPENWAVE = HERE.parents[1] / "candidate_laws/interval_echo/openwave"
GATE_TOLERANCE = 1e-11


def decimal_log(ratio):
    return (Decimal(ratio.numerator) / Decimal(ratio.denominator)).ln()


def energies(service):
    """Derive energies using only stationary mode weights and the exchange doublet as unit."""
    zero, one_even, one_odd, two_even, two_odd = service["modes"]
    with localcontext() as context:
        context.prec = 70
        beta = decimal_log(one_even / one_odd)
        first = decimal_log(zero / one_even) / beta
        second = decimal_log(zero / two_even) / beta
        if abs(decimal_log(two_even / two_odd) / beta - 1) > Decimal("1e-60"):
            raise ValueError("Exchange splitting differs between contrast magnitudes")
        return float(beta), float(first), float(second)


def quadratic_energy_residual(service):
    """Compute E_two − 4 E_one from the exact held odds before rounding either energy."""
    _, one_even, one_odd, _, _ = service["modes"]
    with localcontext() as context:
        context.prec = 70
        return float(-decimal_log(service["held_ratios"][0]) / decimal_log(one_even / one_odd))


def calibrate_pulse(gain, seed, witnesses=None):
    """Use only a one-port before/after offer comparison to set an unknown dimer drive.

    Each shot has a preparation, an offer symbol and, where valid, an output.
    The paired indicator variance accounts for using the same initial preparation.
    """
    generator = np.random.default_rng(seed)
    exchange, contrast = law.local_term_parts()
    start = law.SRC.index((1, 0))
    end = law.SRC.index((0, 1))
    idle_down_up = np.array([0, 0, -1, 1])
    labels = np.arange(3)[:, None]
    after = labels + idle_down_up[None, :]
    valid = (after >= 0) & (after <= 2)
    before_middle = np.broadcast_to(labels == 1, after.shape)
    after_middle = valid & (after == 1)

    def middle_probability(command):
        # The material device contains a fixed strength unknown to the feedback step.
        gate = expm(-1j * gain * command * (exchange + 1.5 * contrast))
        return abs(gate[end, start]) ** 2

    def read(command, shots):
        middle = middle_probability(command)
        # Sample the independent apparatus decisions separately. This keeps
        # exactly fair coins at 1/2 across native floating-point implementations.
        middle_count = generator.binomial(shots, middle)
        left_count = generator.binomial(shots - middle_count, 0.5)
        load_counts = (left_count, middle_count, shots - middle_count - left_count)
        draws = np.array([generator.multinomial(n, [0.25]*4) for n in load_counts])
        counts = draws / shots
        initial = float(counts[before_middle].sum())
        eligible = float(counts[valid].sum())
        final = float(counts[after_middle].sum() / eligible)
        influence = (
            (after_middle.astype(float) - final * valid) / eligible
            - before_middle
            + initial
        )
        se = float(np.sqrt(np.sum(counts * influence**2) / shots))
        return dict(
            command=command,
            before_middle=initial,
            after_middle_given_valid=final,
            difference=final - initial,
            standard_error=se,
            words=draws.tolist(),
        )

    lower, upper = 0.15, 2.25
    training = []
    for _ in range(11):
        command = (lower + upper) / 2
        item = read(command, 2**19)
        training.append(item)
        if item["difference"] > 0:
            lower = command
        else:
            upper = command
    held = read((lower + upper) / 2, 2**21)
    # Propagate the achieved command using the known simulated device only
    # after feedback. No magnitude-two events were drawn in this calibration.
    middle = middle_probability(held["command"])
    endpoint = (1 - middle) / 2
    low, centre = (2*endpoint + middle)/4, (2*endpoint + 2*middle)/4
    z, u, v = 2*low**2 + centre**2, 2*low*centre, low**2
    m = standard_count() if witnesses is None else witnesses
    predicted = (v*z**3 / u**4)**m
    return dict(
        device_gain=gain,
        seed=seed,
        training=training,
        held=held,
        residual_within_readout_precision=abs(held["difference"]) <= 5 * held["standard_error"],
        model_prediction_at_achieved_command=dict(
            generator_exposure=gain*held["command"], middle_probability=middle,
            valid_offers=1, witnesses=m, magnitude_two_even_ratio=predicted,
            magnitude_two_odd_ratio=predicted,
        ),
    )


def standard_count():
    reference = cell()
    counts = {
        len(offer_service.incident_witnesses(reference, e)) for e in reference.seams
    }
    if len(counts) != 1:
        raise ValueError(
            "The shielded primitive supplies several inequivalent reference counts"
        )
    return counts.pop()


def selected_offers(witnesses, reader):
    """Use the existing phase-sensitive pi-contact tolerance and keep the adjacent failure."""
    for offers in range(33):
        s = offer_service.finite_service(offers, witnesses, reader=reader)
        error = abs(np.expm1(-1j * pi * quadratic_energy_residual(s)))
        if error < GATE_TOLERANCE:
            return offers
    raise ValueError("No finite source horizon satisfies the required phase tolerance")


def extend_collar_energy(splittings):
    """The common-model blank, coin, load and echo extension of a measured five-mode spectrum."""
    phase_generator = np.zeros((1 + law.N_LAB, 1 + law.N_LAB))
    for coin in (0, 1):
        for magnitude, (even_energy, odd_energy) in enumerate(splittings, 1):
            plus = 1 + law.LABELS.index((coin, magnitude))
            minus = 1 + law.LABELS.index((coin, -magnitude))
            mean = (even_energy + odd_energy) / 2
            difference = (even_energy - odd_energy) / 2
            phase_generator[np.ix_([plus, minus], [plus, minus])] = [
                [mean, difference],
                [difference, mean],
            ]
    echo = law.echo()
    return echo @ np.kron(phase_generator, np.eye(law.N_SRC)) @ echo.T


def stationary_operator(witnesses=None):
    """Return the experimentally proposed local generator, compiled from a calibrated finite source."""
    m = standard_count() if witnesses is None else witnesses
    source = offer_service.exact_stationary_source(offers=1, witnesses=m)
    zero, one_even, one_odd, two_even, two_odd = source["mode_probabilities"]
    beta = np.log(one_even / one_odd)
    levels = [
        (np.log(zero / pos) / beta, np.log(zero / neg) / beta)
        for pos, neg in ((one_even, one_odd), (two_even, two_odd))
    ]
    return extend_collar_energy(levels), source


def compile_native(service):
    """Use the logarithm of the measured fixed-state eigenvalues on the ten written labels."""
    zero, one_even, one_odd, two_even, two_odd = service["modes"]
    beta, first, second = energies(service)
    splittings = [
        (float(decimal_log(zero / pos)) / beta, float(decimal_log(zero / neg)) / beta)
        for pos, neg in ((one_even, one_odd), (two_even, two_odd))
    ]
    operator = extend_collar_energy(splittings)
    injection = np.eye(law.DIM)[:, : law.N_SRC]
    interchange = np.eye(law.N_SRC)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    potential = np.array([0.0, first, second])
    reduced = (np.eye(law.N_SRC) - interchange) / 2 + np.diag(
        [potential[abs(a - b)] for a, b in law.SRC]
    )
    np.testing.assert_allclose(
        operator @ injection, injection @ reduced, atol=2e-13, rtol=0
    )
    return operator


def local_spectrum(operator):
    """Two complementary Ramsey outcomes on empty load ports recover all four energy gaps."""
    source_index = law.SRC.index((0, 0))

    def label(d):
        return np.eye(law.DIM)[
            :, (1 + law.LABELS.index((0, d))) * law.N_SRC + source_index
        ]

    reference = label(0)
    exposure = pi / 128
    unitary = expm(-1j * exposure * operator)
    response = []
    for magnitude in (1, 2):
        for parity in (1, -1):
            mode = (label(magnitude) + parity * label(-magnitude)) / np.sqrt(2)
            initial = (reference + mode) / np.sqrt(2)
            quadrature = (reference - 1j * mode) / np.sqrt(2)
            evolved = unitary @ initial
            x = abs(np.vdot(initial, evolved)) ** 2
            y = abs(np.vdot(quadrature, evolved)) ** 2
            gap = np.arctan2(2 * y - 1, 2 * x - 1) / exposure
            response.append(
                dict(
                    magnitude=magnitude,
                    parity=parity,
                    x_probability=x,
                    y_probability=y,
                    energy=gap,
                )
            )
    return response


def pair_release(energy, home, times, original, *, launch_index=0, seed=917421, record_draw_hash=False):
    board = federation(0)
    sector = PairSector(board)
    edges = sorted(
        tuple(sorted(e)) for e in board.seams if all(p // 12 == home for p in e)
    )
    gateways = {v for e in board.external for v in e}
    # Cell 13 reproduces the published registration. Cell 7 uses a second internal launch.
    candidates = [e for e in edges if any(p in gateways for p in e) == original]
    source = candidates[launch_index % len(candidates)]
    patch = [sector.index(*e) for e in edges]
    initial = sector.index(*source)
    report = []
    for active in (True, False):
        h = (
            sector.exchange()
            if not active
            else sector.exchange()
            - 2 * energy * sparse.diags(sector.adjacent.astype(float))
        )
        state = np.zeros(sector.P, dtype=complex)
        state[initial] = 1
        previous = 0
        readings = []
        arm = int(not active) + (0 if original else 2)
        for index, t in enumerate(times):
            state = expm_multiply(-1j * (t - previous) * h, state)
            p = pair.local_law(state, patch)
            draw = np.random.default_rng([seed, arm, index]).choice(
                len(p), size=pair.SHOTS, p=p / p.sum()
            )
            counts = np.bincount(draw, minlength=len(p))
            readings.append(
                dict(
                    exposure=t,
                    probability=float(1 - p[-1]),
                    counts=counts.tolist(),
                    outcomes=p.tolist(),
                    norm_error=abs(np.vdot(state, state).real - 1),
                )
            )
            if record_draw_hash:
                readings[-1]["draw_bytes_sha256"] = hashlib.sha256(draw.astype(np.uint8).tobytes()).hexdigest()
            previous = t
        report.append(readings)
    native, control = report
    lower = min(1 - row["counts"][-1] / pair.SHOTS - pair.MARGIN for row in native)
    upper = 1 - control[-1]["counts"][-1] / pair.SHOTS + pair.MARGIN
    return dict(
        cells=board.cells,
        sector_size=sector.P,
        source=source,
        addresses=edges,
        quiet_index=len(edges),
        native=native,
        exchange_only=control,
        native_lower=lower,
        control_upper=upper,
        passes=bool(lower > 0.8 and upper < 0.65),
    )


def registered_pair_release(energy):
    """Reproduce the externally enrolled launch with the standard's binary-sector energy.

    Sampling uses the same registered streams for direct comparison of every
    encoded outcome, rather than claiming statistically independent repeats.
    """
    commit = "bb697927ff37446130d4c6e8c68126b6a3eac00d"
    times = (3, 11, 29)
    key = hashlib.sha256(b"M10.3:row4:" + bytes.fromhex(commit)).digest()
    eligible = [cell_id for cell_id in range(20) if cell_id not in (0, 13)]
    home = eligible[int.from_bytes(key[:4], "big") % len(eligible)]
    launch_index = int.from_bytes(key[4:8], "big")
    seed = int.from_bytes(key[8:16], "big")
    calculated = pair_release(energy, home, times, True, launch_index=launch_index,
                              seed=seed, record_draw_hash=True)
    file = HERE.parents[2] / "openwave/interval_echo/row_04/registered_release.json"
    payload = file.read_bytes()
    archived = json.loads(payload)
    if (archived["registration_commit"] != commit or archived["cell"] != home
            or archived["launch"] != list(calculated["source"]) or archived["shots"] != pair.SHOTS
            or not np.isclose(archived["simultaneous_radius"], pair.MARGIN, rtol=2e-15, atol=0)
            or archived["outcomes"] != [list(edge) for edge in calculated["addresses"]] + ["quiet"]
            or [arm["kappa"] for arm in archived["arms"]] != [8,0]):
        raise ValueError("The deposited Row 4 acquisition differs from its externally registered design")
    cells_and_quiet = len(calculated["addresses"]) + 1
    for arm, native in zip(archived["arms"], (calculated["native"], calculated["exchange_only"]), strict=True):
        for deposited, derived in zip(arm["records"], native, strict=True):
            raw = zlib.decompress(base64.b64decode(deposited["draws_zlib_base64"], validate=True))
            observed = np.frombuffer(raw, dtype=np.uint8)
            if (observed.size != pair.SHOTS or observed.max() >= cells_and_quiet
                    or deposited["exposure"] != derived["exposure"]
                    or deposited["counts"] != derived["counts"]
                    or np.bincount(observed,minlength=cells_and_quiet).tolist() != deposited["counts"]
                    or hashlib.sha256(raw).hexdigest() != derived["draw_bytes_sha256"]):
                raise ValueError("The derived pair replay changed a registered draw or outcome count")
    if not archived["passed"] or not calculated["passes"]:
        raise ValueError("The registered finite patch-persistence gates failed")
    return dict(registration_commit=commit, deposit_sha256=hashlib.sha256(payload).hexdigest(),
                home_cell=home, seed=seed, shots_per_setting=pair.SHOTS, simultaneous_radius=float(pair.MARGIN),
                checked_count_bins=2*len(times)*cells_and_quiet, checked_draws=2*len(times)*pair.SHOTS,
                scope="Freely released two-load pair on the registered finite home patch",
                native_lower_threshold=0.80, final_control_upper_threshold=0.65, **calculated)


def wave_forecast(e1):
    data = json.loads((OPENWAVE / "row_30_selective/channel_review.json").read_text())
    readings = []
    for row in data["frequencies_from_seam"]:
        q = row["q"]
        forecast = abs(e1) / 2 * np.sqrt(3 - np.sqrt(6 + 2 * np.cos(q)))
        readings.append(
            dict(
                q=q,
                observed=row["observed"],
                forecast=forecast,
                percent_error=100 * (forecast / row["observed"] - 1),
            )
        )
    return readings


def compact_arm(witnesses, offers, reader, word=(2, 1, 1)):
    s = offer_service.finite_service(
        offers, witnesses, reader=reader, word=word, cap=2048
    )
    beta, first, second = energies(s)
    spectrum = local_spectrum(compile_native(s))
    np.testing.assert_allclose(
        [r["energy"] for r in spectrum],
        [first, first + 1, second, second + 1],
        atol=2e-12,
        rtol=0,
    )
    zero, x, y, u, v = s["modes"]
    predicted_even = x**4 / zero**3
    predicted_odd = predicted_even * y / x
    return dict(
        witnesses=witnesses,
        offers=offers,
        reader=reader,
        offer_word=word,
        beta=beta,
        one_energy=first,
        two_energy=second,
        held_even_ratio=float(u / predicted_even),
        held_odd_ratio=float(v / predicted_odd),
        held_even_relative_error=float(u / predicted_even - 1),
        held_odd_relative_error=float(v / predicted_odd - 1),
        second_energy_minus_quadratic=quadratic_energy_residual(s),
        accepted_mode_probabilities=list(map(float, s["modes"])),
        spectrum=spectrum,
    )


def run():
    count = standard_count()
    depth = selected_offers(count, "common")
    fair_depth = selected_offers(count, "fair")
    arms = {}
    for name, m, n, reader, word in (
        ("standard", count, depth, "common", (2, 1, 1)),
        ("previous_horizon", count, depth - 1, "common", (2, 1, 1)),
        ("removed_contact", count - 1, depth, "common", (2, 1, 1)),
        ("fair_parity", count, fair_depth, "fair", (2, 1, 1)),
        ("biased_offers", count, depth, "common", (1, 1, 2)),
        ("reversed_chart", count, depth, "common_reversed_chart", (2, 1, 1)),
        ("sign_blind_odd", count, depth, "common_sign_blind", (2, 1, 1)),
    ):
        arms[name] = compact_arm(m, n, reader, word)
    standard = arms["standard"]
    fair = arms["fair_parity"]
    if (
        standard["one_energy"] != count
        or arms["removed_contact"]["one_energy"] != count - 1
    ):
        raise ValueError(
            "Removing one actual witness did not lower the common standard by one"
        )
    if (
        abs(standard["held_even_relative_error"]) > GATE_TOLERANCE
        or arms["biased_offers"]["held_even_ratio"] < 2
    ):
        raise ValueError("The symmetry prediction or its biased-offer control failed")
    s = offer_service.finite_service(depth, count, reader="common", cap=2048)
    a = standard["one_energy"]
    b = standard["two_energy"]
    witness_screens = {}
    for level in (0, 1, 2):
        board = federation(level)
        witnesses = Counter(
            len(offer_service.incident_witnesses(board, e)) for e in board.seams
        )
        witness_screens[str(board.cells)] = {
            str(k): v for k, v in sorted(witnesses.items())
        }
    operator, exact = stationary_operator()
    masked, unused = stationary_operator(count - 1)
    exchange, agreement = law.local_term_parts()
    exact["full_operator_max_error"] = float(
        np.max(abs(operator - exchange - count * agreement))
    )
    exact["masked_full_operator_max_error"] = float(
        np.max(abs(masked - exchange - (count - 1) * agreement))
    )
    np.testing.assert_allclose(
        operator, exchange + count * agreement, atol=2e-13, rtol=0
    )
    np.testing.assert_allclose(
        masked, exchange + (count - 1) * agreement, atol=2e-13, rtol=0
    )
    exact["readouts"] = local_spectrum(operator)
    exact["removed_witness_readouts"] = local_spectrum(masked)
    exact["compiled_states"] = operator.shape[0]
    physical_layouts = offer_service.replay()
    pulse_calibration = [
        calibrate_pulse(gain, 35819 + i, count) for i, gain in enumerate((0.85, 1.12))
    ]
    if not all(case["residual_within_readout_precision"] for case in pulse_calibration):
        raise ValueError("One-port pulse calibration failed the held stationarity read")
    direct = prospective_offer_capture.measure()
    if not direct["enrolled_test"]["within_3_standard_errors"]:
        raise ValueError(
            "The rotated three-contact prospective reading missed its registered gate"
        )
    independent_read = {
        "training_pair_attempts": direct["training"]["pair_attempts"],
        "held_pair_attempts": direct["held"]["pair_attempts"],
        "held_cohorts": direct["held"]["buffered_cohorts"],
        "held_counts": direct["held"]["accepted_counts"],
        "test": direct["enrolled_test"],
    }
    reversible = {}
    for cap in (None, 512):
        comparison = offer_service.reversible_diagnostics(depth, count, cap=cap)
        reversible["until_decided" if cap is None else str(cap)] = {
            key: comparison[key]
            for key in (
                "detailed_balance_error",
                "gap",
                "tv_1pct_spectral_ceiling",
                "tv_1pct_all_point_starts_first_epoch",
                "mean_raw_pair_attempts_per_epoch",
                "mean_raw_pair_attempts_per_change",
                "mean_raw_pair_attempts_per_entry_into_magnitude_two",
                "mean_raw_pair_attempts_per_entry_into_odd_magnitude_two",
            )
        }
    replacements = {
        key: float(s[key])
        for key in (
            "overflow",
            "disagreement",
            "parity_reject",
            "pointer_flip",
            "trace_distance_1pct_epochs",
        )
    }
    return dict(
        source_arms=arms,
        stationary_pulse_standard=exact,
        pulse_calibration=pulse_calibration,
        board_witness_counts=witness_screens,
        physical_layouts=physical_layouts,
        epoch_outcomes=replacements,
        reversible_comparison=reversible,
        prospective_three_contact_acquisition=independent_read,
        buffer_cost={
            k: float(s[k])
            for k in (
                "mean_pairs_per_epoch",
                "mean_pairs_per_flip",
                "mean_pairs_per_either_held",
                "mean_pairs_per_held_odd",
            )
        },
        registered_row_04_pair=registered_pair_release(a),
        standard_pair=[
            pair_release(a, 13, (1, 7, 37), True),
            pair_release(a, 7, (2, 11, 29), False),
        ],
        observers=observer_tests(a, b),
        fair_observers=observer_tests(fair["one_energy"], fair["two_energy"]),
        waves={
            key: wave_forecast(arms[key]["one_energy"])
            for key in ("standard", "fair_parity", "removed_contact")
        },
    )


# These fields come directly from finite rational odds. Retain the represented
# double so a unit prediction is distinguishable from a resolved small miss.
EXACT_RATIO_FIELDS = {"held_even_ratio", "held_odd_ratio"}
RELATIVE_FIELDS = {
    "held_even_relative_error", "held_odd_relative_error", "second_energy_minus_quadratic",
    "accepted_mode_probabilities", "signed_probabilities", "mode_probabilities",
    "endpoint_probabilities", "first_failure_one_chain", "pair_survival", "per_window_survival",
}
BRANCH_FIELDS = {"one_batch_branches", "buffered_epoch_branches", "epoch_outcomes"}


def stable(v, field=None):
    if isinstance(v, dict):
        return {str(k): stable(x, str(k)) for k, x in v.items()}
    if isinstance(v, (tuple, list)):
        return [stable(x, field) for x in v]
    if isinstance(v, (float, np.floating)):
        return float(v) if field in EXACT_RATIO_FIELDS else float(format(v, ".13g"))
    if isinstance(v, np.generic):
        return v.item()
    return v


def check(expected, actual, path=()):
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError("Receipt fields changed")
        for key, value in expected.items():
            check(value, actual[key], (*path, key))
    elif isinstance(expected, list):
        for old, new in zip(expected, actual, strict=True):
            check(old, new, path)
    elif isinstance(expected, float):
        field = path[-1] if path else None
        if field in EXACT_RATIO_FIELDS:
            if actual != expected:
                raise ValueError(f"Exact finite ratio changed at {'.'.join(path)}")
        else:
            # Stochastic rational-source probabilities and stable analytical
            # residuals retain relative precision down to arbitrarily small
            # rates. BLAS-derived norms near zero retain an absolute tolerance.
            strict = field in RELATIVE_FIELDS or any(k in BRANCH_FIELDS for k in path)
            np.testing.assert_allclose(actual, expected, rtol=3e-9, atol=0 if strict else 5e-10,
                                       err_msg=".".join(path))
    elif expected != actual:
        raise ValueError((expected, actual))


def main(argv=None):
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    result = stable(run())
    receipt = HERE / "modular_standard.json"
    if args.write:
        receipt.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        check(json.loads(receipt.read_text()), result)
    s = result["source_arms"]["standard"]
    print(
        f"{s['witnesses']} witnesses, {s['offers']} offers, first and second energies {s['one_energy']}, {s['two_energy']}."
    )
    print("Pair gates:", [r["passes"] for r in result["standard_pair"]])


if __name__ == "__main__":
    main()
