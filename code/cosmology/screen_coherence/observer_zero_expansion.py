"""Exact finite-window radar expansion and unit conversion for the level-13 echo bank."""

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ECHO_FILE = "observer_zero_capacity.json"
ECHO_SHA256 = "b53beb926ac942f6d45e8da606d06726b76908c05991ed1361907273119aaa00"
C = 299_792_458
# CODATA 2022 central value and standard uncertainty, in seconds.
PLANCK_TIME_S = 5.391247e-44
PLANCK_TIME_UNCERTAINTY_S = 0.000060e-44


def read_bank(path):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != ECHO_SHA256:
        raise ValueError(
            "Expected the published observer-zero complete level-13 echo bank"
        )
    bank = json.loads(raw)
    labels = [tuple(label) for label in bank["reflectors"]]
    radar, send = {}, {}
    for label, waits in zip(labels, zip(*bank["echo_waits"], strict=True), strict=True):
        echo_ranges = {
            tick: Fraction(wait, 2)
            for tick, wait in zip(bank["origin_emission_ticks"], waits, strict=True)
        }
        send[label] = echo_ranges
        radar[label] = {s + r: r for s, r in echo_ranges.items()}
    return {"radar": radar, "send": send}


def slopes(ranges, step):
    left, middle, right = ranges
    speed = (right - left) / (2 * step)
    acceleration = (right - 2 * middle + left) / step**2
    return speed, acceleration, speed / middle, acceleration / middle


def window(tracks, center, step, tau, cohort=None):
    """Use only sources observed at all three times, with fixed equal port weights."""
    times = [center - step, center, center + step]
    if step <= 0 or tau <= 0:
        raise ValueError(
            "Both the half window and the seconds per tick must be positive"
        )
    candidates = tracks if cohort is None else cohort
    paired = [p for p in candidates if all(t in tracks[p] for t in times)]
    result = dict(
        times_ticks=times,
        times_s=[t * tau for t in times],
        candidates=len(candidates),
        beacons=len(paired),
        remote_observers=len({owner for owner, _ in paired}),
        sources=paired,
    )
    if not paired:
        return result | {"measurement": None}
    ranges = [
        sum((tracks[p][t] for p in paired), Fraction(0)) / len(paired) for t in times
    ]
    v, acceleration, H, A = slopes(ranges, step)
    changes = [p for p in paired if len({tracks[p][t] for t in times}) > 1]
    result["measurement"] = dict(
        changing_sources=changes,
        mean_range_relay_units=ranges,
        mean_range_m=[x * C * tau for x in ranges],
        speed_over_c=v,
        speed_m_per_s=v * C,
        H_per_tick=H,
        H_per_s=H / tau,
        range_acceleration_relay_units_per_tick2=acceleration,
        range_acceleration_m_per_s2=acceleration * C / tau,
        relative_acceleration_per_tick2=A,
        relative_acceleration_per_s2=A / tau**2,
    )
    return result


def overview(tracks, tau):
    radar, send = tracks["radar"], tracks["send"]
    common = sorted(set.intersection(*(set(track) for track in radar.values())))
    if not common:
        raise ValueError("No common exact reflection times")
    first, last = common[0], common[-1]
    step = 4
    early = window(radar, first - step, step, tau)
    settled = window(radar, first + step, step, tau)
    observed = early["measurement"]
    last_echo = max(s + 2 * r for track in send.values() for s, r in track.items())
    # Two markers from the fixed run's signed local history and its manifest.
    markers = [
        (
            "last_local_repair",
            10,
            "https://dovybzrnzxfo1.cloudfront.net/analysis/observer_zero_information.json",
        ),
        ("global_settlement", 96, "https://dovybzrnzxfo1.cloudfront.net/manifest.json"),
        (
            "last_received_echo_in_complete_replay",
            last_echo,
            "https://dovybzrnzxfo1.cloudfront.net/analysis/" + ECHO_FILE,
        ),
    ]
    steady = all(
        track[t] == track[first]
        for track in radar.values()
        for t in track
        if t >= first
    )
    # The fixed sample is complete at every send time. Keep the one-tick
    # derivatives at every interior sample, including negative and zero values.
    sends = sorted(set.intersection(*(set(track) for track in send.values())))
    means = [sum(track[t] for track in send.values()) / len(send) for t in sends]
    if any(b != a + 1 for a, b in zip(sends, sends[1:])):
        raise ValueError("The send-time series must have consecutive clock ticks")
    series = []
    for i, (t, r) in enumerate(zip(sends, means, strict=True)):
        row = {"send_tick": t, "mean_range_relay_units": r}
        if 0 < i < len(sends) - 1:
            v, acc, H, A = slopes(means[i - 1 : i + 2], 1)
            row.update(
                speed_over_c=v,
                acceleration_relay_units_per_tick2=acc,
                H_per_tick=H,
                relative_acceleration_per_tick2=A,
            )
        series.append(row)
    return dict(
        integer_radar_ticks_with_full_coverage=[
            t for t in common if t.denominator == 1
        ],
        first_common_radar_tick=first,
        last_common_radar_tick=last,
        all_individual_ranges_constant_from_first_common_time=steady,
        elapsed_from_preparation=[
            dict(
                event=event,
                tick=t,
                seconds_at_supplied_tick_duration=t * tau,
                seconds_if_one_planck_time_per_tick=float(t) * PLANCK_TIME_S,
                source=url,
            )
            for event, t, url in markers
        ],
        service_cadence="One local reference flip per scheduled round is supplied in the clocked replay. Global settlement is a manifest marker. Later echoes use the stationary continuation.",
        early_common_time_window=early,
        same_early_beacons_later=window(
            radar, first + step, step, tau, early["sources"]
        ),
        full_common_time_window=settled,
        all_beacons_at_send_ticks_0_4_8=window(send, 4, step, tau),
        full_cohort_send_history=series,
        one_planck_time_hypothesis=dict(
            status="Conditional unit assignment using the laboratory Planck time",
            codata="https://physics.nist.gov/cuu/pdf/RevModPhys.97.025002.pdf",
            seconds_per_tick=PLANCK_TIME_S,
            standard_uncertainty_s=PLANCK_TIME_UNCERTAINTY_S,
            metres_per_radar_unit=C * PLANCK_TIME_S,
            early_window_H_per_s=float(observed["H_per_tick"]) / PLANCK_TIME_S,
            early_window_relative_acceleration_per_s2=float(
                observed["relative_acceleration_per_tick2"]
            )
            / PLANCK_TIME_S**2,
            early_window_range_acceleration_m_per_s2=float(
                observed["range_acceleration_relay_units_per_tick2"]
            )
            * C
            / PLANCK_TIME_S,
            local_clock="One ideal blank dimer flip per radar tick at the supplied service cadence. The isolated two-state generator has eigenvalues 8 and 9 with flip at pi and return at 2*pi generator-time units.",
            clock_receipts="https://dovybzrnzxfo1.cloudfront.net/analysis/observer_zero_information.json",
            gravity_test="For the same two energy eigenstates, compare w=DeltaE*tau/hbar with g=G*DeltaM/(c^3*tau) from independent weak-field gravity readings. With DeltaM=DeltaE/c^2 and the radar signal identified with light, g/w=(t_P/tau)^2.",
            clock_gap_radians_per_tick=math.pi,
            required_gravitational_length_in_radar_units_if_one_planck_tick=math.pi,
            observed_gravitational_length_in_radar_units=None,
            measured_planck_times_per_tick=None,
        ),
    )


def report(tracks, tau, coordinate="radar", center=None, step=4):
    if tau <= 0:
        raise ValueError("Seconds per tick must be positive")
    result = dict(
        source="https://dovybzrnzxfo1.cloudfront.net/analysis/" + ECHO_FILE,
        source_sha256=ECHO_SHA256,
        receiver=[0, 0],
        c_m_per_s=C,
        seconds_per_tick=tau,
        metres_per_relay_range_unit=C * tau,
        clock_attachment="The caller supplies seconds per observer tick. The published example assigns one second per tick. This history has no independent laboratory clock comparison.",
        ruler="One radar range unit is the distance assigned to one tick of signal flight, with c fixed to the specified value. This is material relay radar.",
        cohort="One vote per terminally reachable remote beacon port. A common-time window uses the intersection of source labels across the three exact radar times. Sources with missing readings at any of the three times are excluded throughout that window.",
        estimators="Mean range R. Symmetric speed v=(R_plus-R_minus)/(2h), acceleration g=(R_plus-2R_mid+R_minus)/h^2, H=v/R_mid, A=g/R_mid. These are finite differences. The midpoint time is send + (receive-send)/2. The send series instead groups equal launch times.",
    )
    if center is None:
        result.update(overview(tracks, tau))
    else:
        result.update(
            time_coordinate=coordinate,
            window=window(tracks[coordinate], center, step, tau),
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bank", type=Path, default=HERE / ECHO_FILE)
    parser.add_argument(
        "--seconds-per-tick",
        type=Fraction,
        default=Fraction(1),
        help="Unit convention or independent clock calibration, e.g. 1/1000",
    )
    parser.add_argument(
        "--center", type=Fraction, help="Report a window centered on this tick"
    )
    parser.add_argument(
        "--half-window",
        type=Fraction,
        help="Half width of the custom window in ticks (default 4)",
    )
    parser.add_argument(
        "--coordinate",
        choices=("radar", "send"),
        help="Clock coordinate for a custom window (default radar)",
    )
    parser.add_argument(
        "--output", type=Path, help="Write JSON to this path (default stdout)"
    )
    args = parser.parse_args()
    if args.center is None and (
        args.half_window is not None or args.coordinate is not None
    ):
        parser.error("Specify --center when selecting a custom window or coordinate")
    try:
        result = report(
            read_bank(args.bank),
            args.seconds_per_tick,
            args.coordinate or "radar",
            args.center,
            args.half_window if args.half_window is not None else 4,
        )
    except ValueError as exc:
        parser.error(str(exc))
    payload = json.dumps(result, indent=1, default=str) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload)


if __name__ == "__main__":
    main()
