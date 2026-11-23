"""Transfer an archived SPARC acceleration scale to external gas-rich UDGs."""

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
CALIBRATION = HERE.parent / "receipts" / "held_galaxies.json"
GM_SUN = 1.3271244e20  # IAU 2015 Resolution B3 nominal solar mass parameter, SI
KPC = 3.085677581491367e19


def fit_metadata():
    transfer = json.loads(CALIBRATION.read_text())["compact_exterior_profile"][
        "amplitude_transfer_from_disjoint_galaxies"
    ]
    return transfer["fitted_a0_m_s2"], transfer[
        "training_galaxy_bootstrap_a0_95pct_m_s2"
    ]


def asymptote(mass, a0):
    return (mass * GM_SUN * a0) ** 0.25 / 1000


def effective_scale(mass, speed):
    return (speed * 1000) ** 4 / (GM_SUN * mass)


def required_tilt(inclination, observed, predicted):
    return math.degrees(
        math.asin(observed * math.sin(math.radians(inclination)) / predicted)
    )


def sample(rows, scale, lower_scale):
    result = []
    for agc, dist, tilt, log_mass, mass_error, speed, err_low, err_high in rows:
        mass = 10**log_mass
        predicted = asymptote(mass, scale)
        bounds = []
        for span in (1, 2):
            low_mass = 10 ** (log_mass - span * mass_error)
            max_observed = speed + span * err_high
            lower_prediction = asymptote(low_mass, lower_scale)
            bounds.append(
                dict(
                    width=span,
                    lower_prediction_kms=lower_prediction,
                    upper_observation_kms=max_observed,
                    minimum_gap_kms=lower_prediction - max_observed,
                )
            )
        result.append(
            dict(
                agc=agc,
                observed_kms=speed,
                predicted_kms=predicted,
                ratio=predicted / speed,
                effective_a0_m_s2=effective_scale(mass, speed),
                needed_tilt_deg=required_tilt(tilt, speed, predicted),
                bounds=bounds,
            )
        )
    return result


def load(path):
    return json.loads(path.read_text())


def resolved(rows, scale, scale_floor):
    forecasts = []
    for ring in rows:
        radius = ring["radius_kpc"] * KPC
        vbar = ring["baryonic_kms"]
        vcirc = ring["circular_kms"]
        err = ring["plotted_error_kms"]
        gbar = (vbar * 1000) ** 2 / radius
        gobs = (vcirc * 1000) ** 2 / radius

        def prediction(a0):
            return math.sqrt(radius * (gbar + math.sqrt(a0 * gbar))) / 1000

        forecasts.append(
            dict(
                radius_kpc=ring["radius_kpc"],
                gbar_m_s2=gbar,
                selected=gbar < 1.2e-11,
                baryonic_kms=vbar,
                observed_kms=vcirc,
                plotted_error_kms=err,
                predicted_kms=prediction(scale),
                prediction_at_scale_floor_kms=prediction(scale_floor),
                standardized_residual=(prediction(scale) - vcirc) / err,
                scale_floor_residual=(prediction(scale_floor) - vcirc) / err,
                signed_inverse_a0_m_s2=(gobs - gbar) * abs(gobs - gbar) / gbar,
            )
        )
    return forecasts


def report():
    data = load(HERE / "observations.json")
    scale, (lower, upper) = fit_metadata()
    rows = sample(data["table"], scale, lower)
    result = dict(
        acceleration_scale=scale,
        training_interval=[lower, upper],
        calibration_sha256=hashlib.sha256(CALIBRATION.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(
            (HERE / "observations.json").read_bytes()
        ).hexdigest(),
        six_galaxies=rows,
        median_predicted_to_measured=median(r["ratio"] for r in rows),
        misses=[
            sum(r["bounds"][i]["minimum_gap_kms"] > 0 for r in rows) for i in (0, 1)
        ],
    )
    figures = load(HERE / "rings.json")["panels"]
    result["ring_inputs_sha256"] = hashlib.sha256(
        (HERE / "rings.json").read_bytes()
    ).hexdigest()
    result["resolved_rings"] = {
        fig: resolved(data["rings"], scale, lower) for fig, data in figures.items()
    }
    primary = result["resolved_rings"]["figure8"]
    mass, mass_error = 1.47e9, 0.14e9
    terminal = primary[-1]
    result["followup"] = dict(
        mass_solar=mass,
        mass_error_solar=mass_error,
        observed_tilt_deg=31,
        error_tilt_deg=2,
        observed_distance_mpc=78.7,
        error_distance_mpc=1.5,
        asymptotic_prediction_kms=asymptote(mass, scale),
        required_tilt_deg=required_tilt(
            31, terminal["observed_kms"], asymptote(mass, scale)
        ),
        max_abs_baryonic_two_curve_discrepancy_kms=max(
            abs(r["baryonic_kms"] - r["baryonic_from_second_curve_kms"])
            for r in figures["figure8"]["rings"]
        ),
        sign_of_residual="positive means predicted speed is higher than observed",
        data_geometry="Independent optical 31-degree inclination; gas rotation and kinematic-model uncertainties remain correlated across rings.",
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="regenerate the archived transfer"
    )
    args = parser.parse_args()
    packet = report()
    out = HERE / "results.json"
    if args.write:
        out.write_text(json.dumps(packet, indent=2, allow_nan=False) + "\n")
    elif packet != load(out):
        raise SystemExit("External transfer differs from the archive")
    print(packet["misses"], packet["median_predicted_to_measured"])
