"""Read spatial phase and a receiver-local phase delay on one guided charge."""

import argparse
from hashlib import sha256
import json
from math import asin, atan2, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm
from scipy.sparse import diags

from ... import measure as source
from ...apparatus import LAW, blocks
from ...arm_wavelength import interval
from ..guides import detour
from ..impulse import certify
from ..measure import clock


HERE = Path(__file__).resolve().parent
SHOTS = 2**56
SEED = 513029
LONG = 1
SHORT = 1 / 16
TAIL = 0.01 / 4096
NUMERICAL_MARGIN = 2e-15
CLOCK_DELAYS = 2**36
CLOCK_TICK = 1 / 32


class UnresolvedRecord(ValueError):
    """A valid bank whose measured signal has insufficient resolution."""


def bank(probabilities, rng, admitted, shots=SHOTS):
    """Keep separately the detector complement and the two upstream refusals."""
    values = np.asarray(probabilities, dtype=float) / 3
    remainder = admitted / 3 - values.sum()
    p = np.asarray([*values, remainder, (1 - admitted) / 3, 2 / 3])
    if p.min() < -5e-16 or abs(p.sum() - 1) > 5e-16:
        raise ValueError("Incomplete outcome bank")
    p = p.clip(0, 1)
    return rng.multinomial(shots, p / p.sum()).tolist()


def bounds(count):
    lo, hi = interval(count, SHOTS, TAIL)
    return max(0.0, lo - NUMERICAL_MARGIN), min(1.0, hi + NUMERICAL_MARGIN)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def minus(a, b):
    return a[0] - b[1], a[1] - b[0]


def ratio(a, b):
    if b[0] <= 0:
        raise UnresolvedRecord("Unresolved positive denominator")
    endpoints = [x / y for x in a for y in b]
    return min(endpoints), max(endpoints)


def hop_bank(rng):
    """A separate blank 99-column transfer on the time of the long delay."""
    evolution = expm(-1j * LONG * blocks()["blank"] / 8)
    origin = LAW.SRC.index((0, 1))
    arrival = LAW.SRC.index((1, 0))
    probability = abs(evolution[arrival, origin]) ** 2
    return [int(rng.binomial(SHOTS, probability)), SHOTS]


def read_hop(record):
    count, starts = record
    if starts != SHOTS:
        raise ValueError("Hop reference has the wrong denominator")
    return [2 * asin(sqrt(p)) for p in bounds(count)]


def verify_parity_instrument():
    """A quarter reference phase then the blank half swap measures both parities."""
    H = blocks()["blank"] / 8
    basis = [LAW.SRC.index(pair) for pair in ((0, 1), (1, 0))]
    U = expm(-0.5j * pi * H[:9, :9])
    V = U[np.ix_(basis, basis)] @ np.diag([1, -1j])
    for observed, intended in zip(
        V, np.array([[1, 1], [1, -1]]) / sqrt(2), strict=True
    ):
        np.testing.assert_allclose(
            np.outer(observed.conj(), observed),
            np.outer(intended, intended),
            atol=2e-14,
        )
    return certify(np.pi / 2)


def spatial_sites(legs, ports):
    length = len(legs[0]) + 1
    stations = range(length // 3, (2 * length) // 3 + 1)
    index = {p: j for j, p in enumerate(ports)}
    rows = np.array([[index[p] for p in leg] for leg in legs]).T
    return stations, rows


def spatial_banks(wave, legs, ports, rng):
    stations, rows = spatial_sites(legs, ports)
    wave = wave.reshape(-1, 64)
    admitted = np.linalg.norm(wave) ** 2
    centers = np.linalg.norm(wave[rows[np.array(stations) - 1]]) ** 2
    parity_banks = []
    for side in (-1, 1):
        for parity in (0, 1):
            cohort = np.array([j for j in stations if j % 2 == parity])
            a = wave[rows[cohort - 1]]
            b = wave[rows[cohort + side - 1]]
            parity_banks.append(
                bank(
                    [np.linalg.norm(a + b) ** 2 / 2, np.linalg.norm(a - b) ** 2 / 2],
                    rng,
                    admitted,
                )
            )
    return {
        "middle_count": bank([centers], rng, admitted),
        "contact_parity_counts": parity_banks,
        "ordinals": list(stations),
        "physical_station_ports": [[leg[j - 1] for leg in legs] for j in stations],
    }


def read_spatial(records):
    middle = records["middle_count"]
    banks = records["contact_parity_counts"]
    if (
        len(middle) != 4
        or len(banks) != 4
        or any(len(r) != 5 for r in banks)
        or any(sum(r) != SHOTS for r in [middle, *banks])
    ):
        raise ValueError("Incomplete spatial count bank")
    contrasts = [
        minus(bounds(n_plus), bounds(n_minus)) for n_plus, n_minus, *_ in banks
    ]
    signed = np.sum(contrasts, axis=0)
    cosine = ratio(signed, [4 * p for p in bounds(middle[0])])
    if min(cosine) <= -1 or max(cosine) >= 1:
        raise UnresolvedRecord("Unresolved spatial phase")
    return {
        "cosine": list(cosine),
        "radians_per_step": np.arccos(cosine[::-1]).tolist(),
    }


def temporal_banks(coefficients, modes, energy, rng):
    """Use a supplied path control to delay the actual charge in one branch."""
    admitted = np.linalg.norm(coefficients) ** 2
    arrival = modes[64:128]
    reference = arrival @ coefficients
    records = []
    for time in (SHORT, LONG):
        delayed = arrival @ (coefficients * np.exp(-1j * energy * time))
        for phase in (1, 1j):
            records.append(
                bank(
                    [
                        np.linalg.norm(reference + phase * delayed) ** 2 / 4,
                        np.linalg.norm(reference - phase * delayed) ** 2 / 4,
                    ],
                    rng,
                    admitted,
                )
            )
    return records


def phase_interval(real, imaginary):
    midpoint = atan2(sum(imaginary), sum(real))
    corners = []
    for x in real:
        for y in imaginary:
            corner = atan2(y, x)
            corners.append(midpoint + (corner - midpoint + pi) % (2 * pi) - pi)
    if real[0] <= 0 <= real[1] and imaginary[0] <= 0 <= imaginary[1]:
        raise UnresolvedRecord("Phase is unresolved at the origin")
    return min(corners), max(corners)


def read_temporal(records, hop):
    if len(records) != 4 or any(len(r) != 5 or sum(r) != SHOTS for r in records):
        raise ValueError("Incomplete temporal count bank")
    pieces = [minus(bounds(n_plus), bounds(n_minus)) for n_plus, n_minus, *_ in records]
    short = phase_interval(*pieces[:2])
    long = phase_interval(*pieces[2:])
    unwrapped = [
        tuple(x + 2 * pi * n for x in long)
        for n in range(-16, 17)
        if max(long[0] + 2 * pi * n, short[0] * LONG / SHORT)
        <= min(long[1] + 2 * pi * n, short[1] * LONG / SHORT)
    ]
    if len(unwrapped) != 1:
        raise UnresolvedRecord("Temporal phase winding is unresolved")
    return list(ratio(unwrapped[0], hop))


def blocks_at_middle(h, legs, ports, hopping):
    length = len(legs[0]) + 1
    checked = 0
    index = {p: i for i, p in enumerate(ports)}
    for leg in legs:
        for pos in range(length // 3, (2 * length) // 3 + 1):
            left, site, right = map(index.__getitem__, leg[pos - 2 : pos + 1])
            for word in range(64):
                row = h.getrow(64 * site + word)
                if dict(zip(row.indices, row.data, strict=True)) != {
                    64 * site + word: 17.0,
                    64 * left + word: -hopping,
                    64 * right + word: -hopping,
                }:
                    raise ValueError("Unexpected native bulk operator")
                checked += 1
    return checked


def source_lines(ev, input_weights, length, half_bandwidth):
    zero = ev.searchsorted(17 - half_bandwidth)
    second = ev.searchsorted(17 - half_bandwidth * np.cos(1.5 * pi / length))
    rows = [zero, second]
    for fraction in (0.425, 0.625, 0.82):
        n = round(fraction * length)
        lo, hi = 17 - half_bandwidth * np.cos(
            np.array([n - 0.5, n + 0.5]) * pi / length
        )
        allowed = slice(ev.searchsorted(lo), ev.searchsorted(hi))
        rows.append(allowed.start + np.argmax(abs(input_weights[allowed])))
    return rows


def read_longitudinal(records, reference):
    spatial = read_spatial(records["spatial"])
    energy = read_temporal(records["temporal"], reference)
    return {
        "spatial": spatial,
        "energy": energy,
        "inferred_V": list(add(energy, spatial["cosine"])),
    }


def secants(lines):
    def slope(a, b):
        return ratio(
            minus(b["energy"], a["energy"]),
            minus(b["spatial"]["radians_per_step"], a["spatial"]["radians_per_step"]),
        )

    first, middle, last = lines
    left = slope(first, middle)
    right = slope(middle, last)
    return {
        "left": list(left),
        "right": list(right),
        "lower_slope_difference": left[0] - right[1],
        "positive_massive_convexity_rejected": bool(left[0] > right[1]),
    }


def read_preparation(
    settings, index, energies, input_weights, modes, legs, ports, reference, rng
):
    coeffs = clock(energies, settings) * input_weights
    wave = modes @ coeffs
    admitted = np.linalg.norm(coeffs) ** 2
    record = {
        "prepared_center": float(settings["center"]),
        "source_admission_before_detector": float(admitted / 3),
        "finite_clock_main_line_fraction": float(abs(coeffs[index]) ** 2 / admitted),
        "spatial": spatial_banks(wave, legs, ports, rng),
        "temporal": temporal_banks(coeffs, modes, energies, rng),
    }
    try:
        record["measured"] = read_longitudinal(record, reference)
    except UnresolvedRecord as exc:
        record["measured"] = None
        record["unresolved_reason"] = str(exc)
    return record


def panel(count, scale, instrument_ports, reference, rng, forecast_offset=None):
    replacement, ports, legs, h, actions, q, small = detour(count)
    if instrument_ports.intersection(
        p for added in replacement.values() for p in added
    ):
        raise ValueError("The detector's parked stock overlaps the detour")
    if scale != 1:
        h = (scale * h + (1 - scale) * diags(h.diagonal())).tocsr()
        for permutation in actions.values():
            if (h[permutation, :][:, permutation] != h).nnz:
                raise ValueError("The hopping comparison breaks arm symmetry")
        small = (q.T @ h @ q).toarray()
    n = len(legs[0]) + 1
    eigenvalues, vectors = eigh(small, driver="evd")
    input_weights = source.BASE[:, 0] @ (q[source.ONE] @ vectors)
    modes = q @ vectors
    readings = []
    for index in source_lines(eigenvalues, input_weights, n, scale):
        settings = dict(
            center=eigenvalues[index], tick=CLOCK_TICK, number_delays=CLOCK_DELAYS
        )
        readings.append(
            read_preparation(
                settings,
                index,
                eigenvalues,
                input_weights,
                modes,
                legs,
                ports,
                reference,
                rng,
            )
        )
    prepared = [line["measured"] for line in readings]
    if forecast_offset is None:
        if any(row is None for row in prepared):
            raise ValueError("Unresolved native panel")
        forecast_offset = prepared[0]["inferred_V"]
        held = prepared[1:]
    else:
        held = prepared
    differences = [
        list(minus(row["inferred_V"], forecast_offset)) if row else None for row in held
    ]
    upper = prepared[-3:]
    return {
        "length": n,
        "source_and_receiving_hubs": ports[:2],
        "off_diagonal_flight_scale": scale,
        "bulk_word_rows_checked": blocks_at_middle(h, legs, ports, scale / 2),
        "readings": readings,
        "three_upper_harmonic_secants": secants(upper) if all(upper) else None,
        "native_forecast_offset": forecast_offset,
        "inferred_band_edge_from_first_line": list(minus(forecast_offset, (1, 1)))
        if scale == 1
        else None,
        "held_energy_minus_native_cosine_forecast": differences,
    }


def synthetic_convex_control():
    k = (1.3, 1.9, 2.6)
    points = [
        {
            "spatial": {"radians_per_step": [x, x]},
            "energy": [sqrt(16**2 + 4 * x**2)] * 2,
        }
        for x in k
    ]
    return secants(points)


def audit_source_doublet(energies, amplitudes, modes, actions, f):
    """Compress each measured effect to the original two logical input columns."""
    cycle = actions[source.CY]
    first = np.zeros(len(modes))
    first[source.ONE] = source.BASE[:, 0]
    other = (first[cycle] - first[cycle[cycle]]) / sqrt(3)
    sign = other[source.ONE] @ source.BASE[:, 1]
    np.testing.assert_allclose(abs(sign), 1, atol=1e-14)

    def two_columns(wave):
        return np.stack(
            [wave, sign * (wave[cycle] - wave[cycle[cycle]]) / sqrt(3)], axis=1
        )

    prepared = two_columns(modes @ amplitudes)
    stations, rows = spatial_sites(f["legs"], f["ports"])
    waves = prepared.reshape(-1, 64, 2)
    effects = [prepared, waves[rows[np.array(stations) - 1]]]
    for side in (-1, 1):
        for parity in (0, 1):
            cohort = np.array([j for j in stations if j % 2 == parity])
            left = waves[rows[cohort - 1]]
            right = waves[rows[cohort + side - 1]]
            effects += [(left + right) / sqrt(2), (left - right) / sqrt(2)]
    reference = prepared[64:128]
    for time in (SHORT, LONG):
        delayed = two_columns(modes @ (amplitudes * np.exp(-1j * energies * time)))[
            64:128
        ]
        for phase in (1, 1j):
            effects += [
                (reference + phase * delayed) / 2,
                (reference - phase * delayed) / 2,
            ]
    residual = 0.0
    for effect in effects:
        output = effect.reshape(-1, 2)
        gram = output.conj().T @ output
        average = np.trace(gram).real / 2
        residual = max(residual, np.linalg.norm(gram - average * np.eye(2), ord=2))
    if residual > 1e-14:
        raise ValueError("An enrolled detector resolves the logical input")
    return {
        "compressed_effect_count": len(effects),
        "maximum_non_scalar_effect_operator_norm_before_source_aperture": float(
            residual
        ),
    }


def archived_arrival_line(reference, forecast_offset, rng):
    """Apply the independent record readers at the exact deposited source clock."""
    payload = (HERE.parents[1] / "results.json").read_bytes()
    archive = json.loads(payload)
    f = source.fixture()
    energies, _, input_weights, q, vectors, _, actions = source.build()
    modes = q @ vectors
    record = read_preparation(
        archive["clock"],
        archive["selected_standard_index"],
        energies,
        input_weights,
        modes,
        f["legs"],
        f["ports"],
        reference,
        rng,
    )
    if record["measured"] is None:
        raise ValueError("The archived arrival line is unresolved")
    difference = minus(record["measured"]["inferred_V"], forecast_offset)
    record["held_energy_minus_native_cosine_forecast"] = list(difference)
    record["same_clock_tick"] = archive["clock"]["tick"]
    record["same_clock_delays"] = archive["clock"]["number_delays"]
    record["source_and_receiving_hubs"] = f["ports"][:2]
    record["archived_source_sha256"] = sha256(payload).hexdigest()
    record["standard_input_effect_audit"] = audit_source_doublet(
        energies, clock(energies, archive["clock"]) * input_weights, modes, actions, f
    )
    return record


def run():
    rng = np.random.default_rng(SEED)
    hop = hop_bank(rng)
    reference = read_hop(hop)
    instrument = verify_parity_instrument()
    ports = {
        instrument["private_port"],
        instrument["staging_port"],
        *instrument["compensating_dimer"],
    }
    runs = [panel(count, 1, ports, reference, rng) for count in range(3)]
    changed = panel(0, 0.5, ports, reference, rng, runs[0]["native_forecast_offset"])
    bridge = archived_arrival_line(reference, runs[0]["native_forecast_offset"], rng)
    math_control = synthetic_convex_control()
    if not all(
        run["three_upper_harmonic_secants"]["positive_massive_convexity_rejected"]
        and all(
            lo <= 0 <= hi for lo, hi in run["held_energy_minus_native_cosine_forecast"]
        )
        for run in runs
    ):
        raise ValueError("Native record challenge differs")
    if not any(
        not interval[0] <= 0 <= interval[1]
        for interval in changed["held_energy_minus_native_cosine_forecast"]
        if interval is not None
    ):
        raise ValueError("Hopping change does not challenge the measured coefficient")
    if math_control["positive_massive_convexity_rejected"]:
        raise ValueError("The convex synthetic comparison must remain admissible")
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "threshold_receipt_sha256": sha256(
            (HERE.parent / "threshold.json").read_bytes()
        ).hexdigest(),
        "starts_per_record": SHOTS,
        "seed": SEED,
        "confidence_tail_per_binomial": TAIL,
        "absolute_probability_margin": NUMERICAL_MARGIN,
        "short_and_long_delay_settings": [SHORT, LONG],
        "blank_hop_counts": hop,
        "blank_hop_phase_standard": reference,
        "phase_instrument": instrument,
        "native": runs,
        "changed_hopping": changed,
        "archived_arrival_line": bridge,
        "mathematical_convex_sensitivity": math_control,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = run()
    if args.record:
        with args.record.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        recorded = json.loads((HERE / "records.json").read_text())
        if recorded != result:
            raise SystemExit("Guided recoil replay differs from its frozen records.")
    for row in result["native"]:
        s = row["three_upper_harmonic_secants"]
        print(
            row["length"],
            row["inferred_band_edge_from_first_line"],
            "upper secant decrease at least",
            s["lower_slope_difference"],
        )


if __name__ == "__main__":
    main()
