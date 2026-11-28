"""Test the measured single-seam clock response against resolved galaxy orbits."""

import argparse
import json
from hashlib import sha256

import numpy as np
from . import local_well as local

HERE = local.HERE


def run():
    labfile = local.lab.HERE / "mediated_release.json"
    medium = local.lab.HERE / "medium_field.json"
    field = json.loads(medium.read_bytes())["from_counts"]["shell_strength"]
    equilibrium = local.lab.HERE / "undisturbed_medium.json"
    station_field = json.loads(equilibrium.read_bytes())[
        "stationary_shell_field_from_counts"
    ]
    prior = HERE / "source_profile.json"
    measured = json.loads(labfile.read_bytes())["decoded"][
        "shell_clock_slowdown_from_counts"
    ]
    inputs = json.loads(prior.read_bytes())["resolved_maps"]["inputs"]
    charts = local.profiles.mapped_annuli(local.profiles.maps.enroll()[0])
    meta, _ = local.profiles.galaxies.catalogue()
    shape = (np.arange(len(measured)), np.asarray(measured))
    files = [
        labfile,
        medium,
        equilibrium,
        local.lab.HERE / "UNDISTURBED_MEDIUM.md",
        local.lab.HERE / "CLOCK_TO_FLIGHT.md",
        prior,
        local.lab.HERE / "MEDIATED_RELEASE.md",
        HERE / "LOCAL_WELL_PROTOCOL.md",
        local.profiles.maps.HERE / "leroy_table7.dat",
        local.profiles.maps.HERE / "leroy_sample.dat",
        local.profiles.galaxies.sparc.DATA / "table1.dat",
    ]
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in files},
        ring_window_fraction=local.WINDOW,
        measured_shell_slowdown=measured,
        held=local.evaluate(inputs, charts, meta, shape),
        stationary=local.evaluate(
            inputs,
            charts,
            meta,
            (np.arange(len(station_field)), np.asarray(station_field)),
        ),
        flow=local.evaluate(
            inputs, charts, meta, (np.arange(len(field)), np.asarray(field))
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "mediated_clock.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    for name in ("held", "flow", "stationary"):
        print(name, {k: v for k, v in result[name].items() if k != "galaxies"})
