"""Test which vacuum densities observer zero's repeated radial clock readings distinguish."""

from fractions import Fraction
from hashlib import sha256
import json
from math import pi
from pathlib import Path


SURVEY = Path(__file__).with_name("level202_distance_redshift.json")
CALIBRATION_TICK = 16


def clock_stamps(row):
    return (
        row["reflection_at"],
        row["echo_at"],
        row["beacon_sent"],
        row["beacon_received"],
    )


def static_stamps(send, flight, period):
    """Immediate reflection with synchronized proper clocks and equal radial light-flight times."""
    reflection = send + flight
    emissions = [reflection - period / 2, reflection + period / 2]
    return reflection, reflection + flight, emissions, [t + flight for t in emissions]


def readout(raw):
    survey = json.loads(raw)
    period = Fraction(survey["four_tick_period"])
    panels = {
        tick: {tuple(row["source"]): row for row in survey["panels"][str(tick)]}
        for tick in survey["query_epochs"]
    }
    ranges = {
        source: Fraction(row["echo_at"] - row["query_at"], 2)
        for source, row in panels[CALIBRATION_TICK].items()
        if row["echo_at"] is not None
        and row["beacon_received"] is not None
        and None not in row["beacon_received"]
    }
    comparisons = {}
    for name, radius in (("Minkowski", None), ("Einstein_static_dust", 100)):
        if radius is None:
            curvature = Fraction(0)
            flight = ranges
        else:
            curvature = Fraction(1, radius**2)
            # Radial angular coordinates on the round three-sphere, fixed once at tick 16.
            angles = {source: distance / radius for source, distance in ranges.items()}
            if any(angle >= pi for angle in angles.values()):
                raise ValueError("The specified radial rays must be shorter than the antipodal arc")
            flight = {source: radius * angle for source, angle in angles.items()}
        matches = {
            str(tick): sum(
                static_stamps(tick, duration, period) == clock_stamps(rows[source])
                for source, duration in flight.items()
            )
            for tick, rows in panels.items()
        }
        # Static Friedmann and acceleration equations for dust give Lambda = K
        # and 8 pi G rho_dust / c^2 = 2 K for spatial sectional curvature K.
        comparisons[name] = dict(
            curvature_radius_in_light_ticks=radius,
            Lambda_times_light_tick_squared=str(curvature),
            eight_pi_G_dust_mass_density_times_tick_squared=str(2 * curvature),
            exact_timestamp_matches_by_query_tick=matches,
        )
    last_tick = max(panels)
    fitting = [
        name for name, model in comparisons.items()
        if ranges and model["exact_timestamp_matches_by_query_tick"][str(last_tick)] == len(ranges)
    ]
    return dict(
        source=SURVEY.name,
        source_sha256=sha256(raw).hexdigest(),
        level=survey["level"],
        receiver=survey["receiver"],
        enrolled=survey["enrolled"],
        range_calibration_tick=CALIBRATION_TICK,
        beacons_with_range_and_period_at_calibration=len(ranges),
        range_in_light_ticks=(
            [float(min(ranges.values())), float(max(ranges.values()))] if ranges else None
        ),
        comparison_scope="Only the calibrated beacons, using their reflection, echo, emission and reception stamps at every survey epoch. Beacon visibility and early routing have separate material dynamics.",
        static_geometries=comparisons,
        last_query_tick=last_tick,
        static_geometries_matching_all_calibrated_beacons_at_last_query=fitting,
        observer_vacuum_energy_density_J_m3=None,
        interpretation="The radial measurements test the optical path to each beacon. Inferring vacuum energy on the Einstein branch requires independent gravitational curvature and matter-stress readings or a measured de Sitter horizon with a physical calibration.",
    )


if __name__ == "__main__":
    print(json.dumps(readout(SURVEY.read_bytes()), indent=2, allow_nan=False))
