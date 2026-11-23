"""Reconstruct the homogeneous Einstein–scalar equations from clocked probe records."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# The following inputs belong to the forward source. Only its probe tapes reach the decoder.
COUPLING, VACUUM, MASS_SQUARED = 1.6, 0.12, 0.7
LAUNCHES = ((0.32, 0.42, 0.52, 0.80),
            (0.27, 0.58, -0.65, 0.73),
            (0.50, 0.31, 0.73, -0.88))
TICK = 2.0 ** -40
SPACINGS = (0.05, 0.025, 0.0125)


def evolve(launch, coupling):
    hx, hy, field, velocity = launch
    density = (velocity**2 + MASS_SQUARED * field**2) / 2
    hz = (VACUUM + coupling*density - hx*hy)/(hx+hy)
    initial = np.array([1., 1., 1., hx, hy, hz, field, velocity, 0., 0., 0.])

    return propagate(initial, 0., 1.8, MASS_SQUARED, coupling, VACUUM)


def propagate(initial, start, end, mass_squared, coupling, vacuum):
    def derivative(_, state):
        arms, expansion = state[:3], state[3:6]
        field, velocity = state[6:8]
        total = expansion.sum()
        potential = mass_squared * field**2 / 2
        return np.r_[arms * expansion, -total*expansion + vacuum + coupling*potential,
                     velocity, -total*velocity - mass_squared*field, 1/arms]

    solution = solve_ivp(derivative, (start, end), initial, dense_output=True,
                         method="DOP853", atol=1e-13, rtol=1e-13)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.sol


def records(solution, appointment, spacing):
    # Two legs to a comoving reflector traverse twice its coordinate separation.
    # Find a recorded trip centered on each appointment for equal-time comparison.
    tape = []
    for offset in (-1, 0, 1):
        time = appointment + offset*spacing
        trips = []
        for axis in range(3):
            def crossing_error(half_trip):
                return (solution(time+half_trip)[8+axis] - solution(time-half_trip)[8+axis]
                        - 2*spacing)
            half_trip = brentq(crossing_error, 0., 4*spacing, xtol=1e-14)
            trips.append([round((time-half_trip)/TICK), round((time+half_trip)/TICK)])
        tape.append({"proper_clock_ticks": round(time/TICK), "radar_send_return_ticks": trips,
                     "scalar_meter_ticks": round(solution(time)[6]/TICK)})
    return {"comoving_arm_spacing": spacing, "tick_unit": TICK, "readings": tape}


def local_state(tape):
    rows = tape["readings"]
    tick = tape["tick_unit"]
    times = np.array([row["proper_clock_ticks"] for row in rows])*tick
    field = np.array([row["scalar_meter_ticks"] for row in rows])*tick
    trips = np.array([row["radar_send_return_ticks"] for row in rows])
    arms = (trips[:,:,1]-trips[:,:,0])*tick/(2*tape["comoving_arm_spacing"])
    step = (times[2]-times[0])/2
    h = (arms[2]-arms[0])/(2*step*arms[1])
    dd_arm = (arms[2]-2*arms[1]+arms[0])/(step**2*arms[1])
    v = (field[2]-field[0])/(2*step)
    acceleration = (field[2]-2*field[1]+field[0])/step**2
    return times[1], arms[1], h, dd_arm, field[1], v, acceleration


def reconstruct(tape):
    _, _, h, dd_arm, field, v, acceleration = local_state(tape)
    scalar_kinetic = acceleration + h.sum()*v
    g00 = float(h[0]*h[1] + h[0]*h[2] + h[1]*h[2])
    spatial = np.array([-(dd_arm[j]+dd_arm[k]+h[j]*h[k])
                        for j,k in ((1,2),(0,2),(0,1))])
    return g00, spatial, field, v, scalar_kinetic


def calibrate(tape):
    g00, spatial, field, velocity, kinetic = reconstruct(tape)
    mass_squared = -kinetic / field
    coupling = (g00 + spatial.mean()) / velocity**2
    vacuum = g00 - coupling * (velocity**2 + mass_squared*field**2)/2
    return {"mass_squared": float(mass_squared), "coupling": float(coupling), "vacuum": float(vacuum)}


def held_errors(tapes, calibrated):
    einstein_errors, field_errors = [], []
    for tape in tapes:
        g00, spatial, field, velocity, kinetic = reconstruct(tape)
        mass = calibrated["mass_squared"]
        density, pressure = [(velocity**2 + sign*mass*field**2)/2 for sign in (1, -1)]
        predicted_00 = calibrated["vacuum"] + calibrated["coupling"]*density
        predicted_space = calibrated["coupling"]*pressure - calibrated["vacuum"]
        einstein_errors.extend(np.abs(np.r_[g00-predicted_00, spatial-predicted_space]).tolist())
        field_errors.append(float(abs(kinetic + mass*field)))
    return {"max_einstein_component_error": max(einstein_errors),
            "max_scalar_equation_error": max(field_errors)}


def forecast(tapes, calibrated):
    time, arms, rates, _, field, velocity, _ = local_state(tapes[0])
    last = tapes[-1]["readings"][1]["proper_clock_ticks"] * tapes[-1]["tick_unit"]
    initial = np.r_[arms, rates, field, velocity, 0., 0., 0.]
    predicted = propagate(initial, time, last, **calibrated)
    errors = []
    for tape in tapes[1:]:
        time, observed, _, _, field, _, _ = local_state(tape)
        state = predicted(time)
        errors.append({"arm_errors": abs(state[:3]-observed).tolist(),
                       "scalar_error": float(abs(state[6]-field))})
    return errors


def run():
    reference = [evolve(launch, COUPLING) for launch in LAUNCHES]
    test_field = [evolve(launch, 0.) for launch in LAUNCHES[1:]]
    rows = []
    for spacing in SPACINGS:
        calibration_tape = records(reference[0], 0.32, spacing)
        calibrated = calibrate(calibration_tape)
        arms = {}
        for name, sources in (("coupled", reference[1:]), ("test_field_on_vacuum_geometry", test_field)):
            tapes = [records(source, time, spacing) for source in sources for time in (0.4, 0.72, 1.15)]
            predictions = [error for start in (0, 3)
                           for error in forecast(tapes[start:start+3], calibrated)]
            arms[name] = {"held_tapes": tapes, "record_forecasts": predictions,
                          "max_forecast_error": max(max(*row["arm_errors"], row["scalar_error"])
                                                    for row in predictions),
                          **held_errors(tapes, calibrated)}
        rows.append({"spacing": spacing, "calibration_tape": calibration_tape,
                     "calibrated": calibrated, "arms": arms})
    for metric in ("max_einstein_component_error", "max_scalar_equation_error", "max_forecast_error"):
        errors = [row["arms"]["coupled"][metric] for row in rows]
        assert errors[0] > errors[1] > errors[2] and errors[2] < 0.005, errors
    assert rows[-1]["arms"]["test_field_on_vacuum_geometry"]["max_einstein_component_error"] > 0.05
    return {"source": "controlled homogeneous Bianchi I Einstein–scalar action",
            "forward_coefficients": [COUPLING, VACUUM, MASS_SQUARED], "refinements": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="save the integer probe records to a fresh path")
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    for row in result["refinements"]:
        print(row["spacing"], row["calibrated"],
              {name: {key: value for key, value in arm.items() if key.startswith("max_")}
               for name, arm in row["arms"].items()})


if __name__ == "__main__":
    main()
