"""Weigh the occupied excitations and test their energy as a gravitating source."""

import argparse
from hashlib import sha256
import json

from mpmath import iv
import numpy as np

from . import thermal_mass as maps
from .occupied_work import experiment

flight = experiment.flight
HERE = maps.HERE
LAB = flight.HERE


def ordinary(twice, rng):
    h = twice / 2
    n = len(h)
    e, basis = np.linalg.eigh(h)
    wave = basis[:, 0]
    lower = np.floor(e[0] * 1e6) / 1e6 - 1e-7
    shifted = h - lower * np.eye(n)
    chol = np.linalg.cholesky(shifted)
    # A reconstructed positive Gram matrix bounds all negative residual
    # eigenvalues. Include dot-product and subtraction roundoff in the row norm.
    gamma = 4 * n * np.finfo(float).eps
    residual = abs(shifted - chol @ chol.T) + gamma * (
        abs(shifted) + abs(chol) @ abs(chol).T
    )
    lo = lower - np.linalg.norm(residual, ord=np.inf)
    ray = wave @ h @ wave / (wave @ wave)
    upper = ray + gamma * (1 + abs(wave) @ abs(h) @ abs(wave))
    assert 0 < lo <= e[0] + 1e-12 <= upper
    tau = 1 / 128
    ages = 2**40
    assert np.abs(h).sum(1).max() * tau < np.pi / 2
    amplitudes = basis.T @ wave
    prob = (1 + (amplitudes**2) @ np.sin(tau * e)) / 2
    positive = int(rng.binomial(ages, prob))
    eps = np.sqrt(flight.RISK / (2 * ages))
    # The trial-wave Ramsey phase differs from the ground-frequency phase
    # by at most t times its spectral residual and ground bracket width.
    error = tau * (np.linalg.norm(h @ wave - ray * wave) + upper - lo) + 2e-10
    sine = 2 * positive / ages - 1 + np.array([-1, 1]) * (2 * eps + error)
    measured = np.arcsin(sine) / tau
    assert measured[0] < lo < upper < measured[1]
    return dict(
        reference="all-one chart and both written contacts in inactive sectors",
        spectroscopic_interval=measured.tolist(),
        ground_energy_certificate=[float(lo), float(upper)],
        trial_rayleigh=float(ray),
        phase_dwell=tau,
        shots=ages,
        positive_flags=positive,
    )


def run():
    rng = np.random.default_rng(78419)
    board, _, edges, _, _ = flight.chart()
    bath = json.loads((LAB / "occupied_equilibrium.json").read_bytes())["results"]
    surveyed = json.loads((HERE / "thermal_mass_maps.json").read_bytes())
    inputs, _, meta, _ = maps.rp.galaxy.transfer.budget.load_inputs()
    count = sum(len(d["radius_kpc"]) for d in inputs.values())
    tests = {}
    for load, tag in ((0, "vacancy"), (2, "excess")):
        twice, _ = flight.apparatus(board, edges, load, -1, 1)
        standard = ordinary(twice, rng)
        slope = iv.mpf(bath[tag]["internal_energy_per_joint_increment"])
        scale = flight.lab.interval_output(
            slope / iv.mpf(standard["spectroscopic_interval"])
        )
        ceiling = max(map(abs, scale)) / 4
        forecasts = {}
        for orientation, sc in (
            ("forward", scale),
            ("reverse", [-scale[1], -scale[0]]),
        ):
            full = maps.forecast(surveyed, sc, inputs)
            forecasts[orientation] = {}
            for tracer, caps in full.items():
                for cap, menus in caps.items():
                    for menu, r in menus.items():
                        forecasts[orientation][f"{tracer}/{cap}/{menu}"] = dict(
                            within_two=r["within_two"],
                            median_abs_kms=r["median_abs_kms"],
                            maximum_speed_change_kms=r[
                                "maximum_baryon_speed_change_kms"
                            ],
                            by_galaxy={
                                k: v["within_two"] for k, v in r["galaxies"].items()
                            },
                        )
        tests[tag] = dict(
            ordinary=standard,
            energy_per_joint_change=flight.lab.interval_output(slope),
            energy_per_ordinary_mass=scale,
            any_independent_population_ceiling=ceiling,
            map_forecasts=forecasts,
        )
    upper = max(row["any_independent_population_ceiling"] for row in tests.values())
    prior = json.loads((HERE / "outer_disk.json").read_bytes())
    mass = (
        prior["all_atomic_He_msun"]
        + prior["all_molecular_He_msun"]
        + 1e9 * meta["NGC2841"]["L36_gsun"]
    )
    need = prior["strongest"]["required_extra_mass_floor_msun"]
    cmb = json.loads((HERE / "thermal_mass.json").read_bytes())[
        "planck_cold_per_baryon"
    ]
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "OCCUPIED_MASS.md",
                HERE / "thermal_mass_maps.json",
                HERE / "thermal_mass.json",
                HERE / "outer_disk.json",
                LAB / "occupied_equilibrium.json",
            )
        },
        galaxy_count=len(inputs),
        ring_count=count,
        results=tests,
        combined_mass_ratio_upper=upper,
        cosmic_ratio=cmb,
        minimum_cosmic_factor=cmb / upper,
        ngc2841=dict(
            inventory_msun=mass,
            granted_extra_msun=mass * upper,
            required_outer_shell_msun=need,
            minimum_factor=need / (mass * upper),
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    path = HERE / "occupied_mass.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    print(
        out["combined_mass_ratio_upper"], out["minimum_cosmic_factor"], out["ngc2841"]
    )
