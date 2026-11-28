"""Compare the balanced linear material-release channel to SPARC mass models."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog

from . import moving_budget
from cosmology.rar_deep_regime.individual_trail_budget import measure as galaxies


def signed_columns(rows):
    v = np.array([[r[k] for k in ("vgas_kms", "vdisk_kms", "vbul_kms")] for r in rows])
    return v * abs(v)


def zero_excess_fit(meta, rows):
    """Grant an entire galaxy common distance, inclination and mass conversions.

    k=(D/D0)(sin(i)/sin(i0))² converts a squared modeled speed to the
    catalogue's nominal inclination and distance. Bounds k*l <= x <= k*u
    remain linear even for a negative squared gas contribution.
    """
    columns = signed_columns(rows)
    v, err = [np.array([r[k] for r in rows]) for k in ("vobs_kms", "e_vobs_kms")]
    low, high = np.maximum(0.0, v - 2 * err) ** 2, (v + 2 * err) ** 2
    baseline = columns @ [1.0, 0.5, 0.7]
    inc, di = meta["inclination_deg"], meta["e_inclination_deg"]
    distance = meta["e_dist_mpc"] / meta["dist_mpc"]
    angle = np.radians(np.clip([inc - 2 * di, inc + 2 * di], 0, 90))
    kmin, kmax = (
        np.maximum(0, [1 - 2 * distance, 1 + 2 * distance])
        * (np.sin(angle) / np.sin(np.radians(inc))) ** 2
    )
    amin, amax = np.array([0.7, 0.3, 0.4]), np.array([1.5, 0.8, 1.0])
    strongest = np.maximum(columns * amin, columns * amax).sum(axis=1)
    weakest = np.minimum(columns * amin, columns * amax).sum(axis=1)
    low_env = np.minimum(kmin * weakest, kmax * weakest)
    high_env = np.maximum(kmin * strongest, kmax * strongest)
    within = (low <= baseline) & (baseline <= high)
    possible = (low <= high_env) & (low_env <= high)
    # x=(scaled gas, scaled disk, scaled bulge, k). One vector per disk.
    orbit = np.column_stack((columns, np.zeros(len(columns))))
    matrix = np.vstack(
        (
            orbit,
            -orbit,
            np.column_stack((np.eye(3), -amax)),
            np.column_stack((-np.eye(3), amin)),
        )
    )
    cap = np.concatenate((high, -low, np.zeros(6)))
    fit = linprog(
        np.zeros(4),
        A_ub=matrix,
        b_ub=cap,
        bounds=[(0, None)] * 3 + [(kmin, kmax)],
        method="highs",
    )
    assert fit.status in (0, 2)
    if fit.success:
        assert np.all(matrix @ fit.x <= cap + 1e-6) and possible.all()
    return dict(
        rings=len(rows),
        within_two=int(within.sum()),
        individually_possible=int(possible.sum()),
        nominal_whole=bool(within.all()),
        shared_parameters_whole=bool(fit.success),
        baryons_too_slow_even_at_upper_limit=int((low > high_env).sum()),
    )


def run():
    native = moving_budget.LAB / "material_release.json"
    lab = json.loads(native.read_bytes())
    for order in lab["schedules"].values():
        for paths in order["contrasts"].values():
            for t in (1, 2):
                for j in range(2):
                    assert (
                        Q(paths[f"{t},0"]["current_per_unit"][j])
                        + Q(paths[f"0,{t}"]["current_per_unit"][j])
                        == 0
                    )
    meta, groups = galaxies.catalogue()
    panels = {}
    for name, rows in groups.items():
        good = [
            r
            for r in rows
            if r["rad_kpc"] > 0 and r["vobs_kms"] > 0 and r["e_vobs_kms"] > 0
        ]
        if good:
            panels[name] = zero_excess_fit(meta[name], good)
    enrolled, *_ = moving_budget.inventory.rp.galaxy.transfer.budget.load_inputs()

    def summarize(names):
        return {
            key: sum(panels[name][key] for name in names)
            for key in next(iter(panels.values()))
        } | {"galaxies": len(names)}

    # Previously enrolled disks can use a proper subset of the catalogue
    # radii. Keep the exact original enrollment for the nominal comparison.
    previous = {}
    for name, tab in enrolled.items():
        v, b, e = [
            np.asarray(tab[k]) for k in ("observed_kms", "baryons_kms", "error_kms")
        ]
        previous[name] = dict(rings=len(v), within_two=int((abs(b - v) <= 2 * e).sum()))
    qualified = [
        n
        for n in panels
        if meta[n]["quality"] in (1, 2) and meta[n]["inclination_deg"] >= 30
    ]
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                native,
                galaxies.sparc.DATA / "table1.dat",
                galaxies.sparc.DATA / "table2.dat",
                moving_budget.HERE / "source_profile.json",
            )
        },
        orientation_balanced_linear_excess=0,
        catalogue_rows_total=sum(map(len, groups.values())),
        all_catalogue=summarize(panels),
        quality_1_2_inclination_at_least_30=summarize(qualified),
        previously_enrolled=previous,
        galaxies=panels,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    path = moving_budget.HERE / "release_rotation.json"
    if args.write:
        path.write_text(output)
    else:
        assert path.read_text() == output
    print(json.loads(output)["all_catalogue"])
