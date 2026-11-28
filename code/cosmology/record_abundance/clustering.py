"""Ordinary occupation versus adjacent disagreement in a stationary native load process."""

import argparse
from hashlib import sha256
from io import BytesIO
import json
from math import pi
import re
import tarfile
from urllib.request import urlopen

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from . import stationary as stock

HERE = stock.measure.HERE
RECEIPT = HERE / "clustering.json"
COPIES, SEED = 8192, 24031128
APPOINTMENTS = (0, 11, 64)
TABLE = (
    (10.1, 0.91, 0.68),
    (27.7, 7.3, 2.5),
    (75.9, 48.3, 3.9),
    (208, 252, 15),
    (568, 679, 82),
    (1554, 956, 180),
)
LE = (
    "0711.3908v2",
    "Joint_Analysis_A1689.tex",
    "f8ff2e8010de97281ebada5fb157abe0ed7488fe992d2dfe1f31aa899afdbf0a",
)
GO = (
    "1309.3565v1",
    "ms.tex",
    "340533654ea063293528070903ba3524701de8b0cbff97f41cd4b34a50c2958e",
)


def source(spec):
    ref, name, digest = spec
    with urlopen("https://arxiv.org/src/" + ref, timeout=30) as response:
        with tarfile.open(fileobj=BytesIO(response.read())) as archive:
            raw = archive.extractfile(name).read()
    if sha256(raw).hexdigest() != digest:
        raise ValueError("Original TeX has changed")
    return raw.decode("latin1")


def verify_public_tables():
    text = source(LE)
    table = text[text.index(r"$r$ [h$^{-1}$ kpc] & $\rho_{\rm total}$") :].split(
        r"\end{tabular}", 1
    )[0]
    rows = re.findall(
        r"^\s*([\d.]+)\s*&[^&\n]+&\s*\$([\d.]+)\s*\\pm\s*([\d.]+)\$", table, re.M
    )
    np.testing.assert_array_equal(np.array(rows, float), TABLE)
    if "$z=0.183$" not in text or not re.search(
        r"f_\{bary\}=0\.136\s*\\pm\s*0\.005", source(GO)
    ):
        raise ValueError("Cluster redshift or baryon-mass inference differs")


def cluster_caps(p):
    # Each of three incident edges contributes at most half a unit to its owner.
    dark_cap = 1 / (2 * p * (1 - p))
    total_cap = max(dark_cap, 1 / p)
    z = 0.183
    critical_today = 2.77536627e11  # h² Msun Mpc^-3 for H0=100h km/s/Mpc.
    # Use Omega_m0 <= 1, so this reference is at least the actual matter mean.
    max_mean = critical_today * (1 + z) ** 3
    rows = []
    for r, mass, error in TABLE:
        size = 4 * pi / 3 * (r / 1000) ** 3
        ratio = 1e12 * mass / (size * max_mean)
        lower = 1e12 * (mass - error) / (size * max_mean)
        rows.append(
            {
                "radius_kpc_h_inverse": r,
                "mass_1e12_msun_h_inverse": mass,
                "quoted_mass_error": error,
                "density_over_maximum_cosmic_mean": ratio,
                "using_lower_mass": lower,
                "lower_over_total_model_cap": lower / total_cap,
            }
        )
    # The weighted baryon fraction is for clusters with both stellar and gas inventories.
    fraction, error = 0.136, 0.005
    nonbaryonic_over_critical = 500 * (1 - fraction - error)
    return {
        "dark_density_cap": dark_cap,
        "occupied_mass_cap": 1 / p,
        "total_density_cap": total_cap,
        "a1689_redshift": z,
        "fixed_apertures": rows,
        "massive_cluster_baryon_fraction": fraction,
        "massive_cluster_baryon_fraction_error": error,
        "r500_nonbaryonic_density_over_critical_at_upper_baryon_fraction": nonbaryonic_over_critical,
        "maximum_local_fraction_of_cluster_nonbaryonic_mass": dark_cap
        / nonbaryonic_over_critical,
    }


def sample_panel(p, settings, owners, distance_bins, rng):
    n = len(distance_bins)
    a, b = owners.T
    # Value 2 tracks the only difference between the planted and empty-root
    # stocks. The common exchange words move the same tagged owner in both.
    encoded = (rng.random((COPIES, n)) < p).astype(np.uint8)
    encoded[:, 0] = 2
    tracer = np.eye(1, n)[0]
    groups = np.eye(4)[distance_bins]
    incidence = (groups[a] + groups[b]) / 2
    responses = {}
    for t in range(APPOINTMENTS[-1] + 1):
        if t:
            x, y = settings[(t - 1) % len(settings)]
            swap = rng.integers(0, 2, (COPIES, len(x)), dtype=np.uint8).astype(bool)
            stock.screen.paired_moves(encoded, x, y, swap)
            stock.screen.apply_mean(tracer, x, y, 0.5)
        if t not in APPOINTMENTS:
            continue
        assert (np.count_nonzero(encoded == 2, axis=1) == 1).all()
        predicted = (1 - 2 * p) * (tracer[a] + tracer[b]) @ incidence
        rows = []
        for panel in np.array_split(encoded, 64):
            plus, minus = panel > 0, panel == 1
            delta = (plus[:, a] != plus[:, b]).astype(np.int8) - (
                minus[:, a] != minus[:, b]
            )
            rows.append(delta @ incidence)
        sample = np.concatenate(rows)
        np.testing.assert_allclose(predicted.sum(), 3 * (1 - 2 * p), atol=1e-13, rtol=0)
        responses[str(t)] = {
            "tracer_probability": (tracer @ groups).tolist(),
            "predicted_charge_change_in_seam_units": predicted.tolist(),
            "sample_mean": sample.mean(axis=0).tolist(),
            "sample_se": (sample.std(axis=0, ddof=1) / np.sqrt(COPIES)).tolist(),
            "total_sample_mean": float(sample.sum(axis=1).mean()),
            "total_sample_se": float(sample.sum(axis=1).std(ddof=1) / np.sqrt(COPIES)),
        }
    return {
        "occupation": p,
        "long_wavelength_fractional_response": (1 - 2 * p) / (1 - p),
        "appointments": responses,
    }


def run():
    board = stock.measure.federation(2)
    owners = np.asarray(board.external) // 12
    settings = [tuple(g[:2]) for g in stock.screen.observations.schedules(board)[0]]
    graph = sparse.csr_matrix(
        (np.ones(len(owners)), owners.T), shape=(board.cells,) * 2
    )
    distance = shortest_path(graph, directed=False, unweighted=True, indices=0)
    bins = np.digitize(distance, [0.5, 1.5, 3.5])
    ps = (stock.screen.DENSITY, 0.5, 0.84)
    rng = np.random.SeedSequence(SEED).spawn(len(ps))
    return {
        "protocol_sha256": sha256((HERE / "CLUSTERING.md").read_bytes()).hexdigest(),
        "samples_per_occupation": COPIES,
        "seed": SEED,
        "cells": board.cells,
        "distance_bin_sizes": np.bincount(bins).tolist(),
        "primary_teX_sha256": {LE[0]: LE[2], GO[0]: GO[2]},
        "equilibrium_perturbations": [
            sample_panel(p, settings, owners, bins, np.random.default_rng(seed))
            for p, seed in zip(ps, rng, strict=True)
        ],
        "observational_comparison": cluster_caps(ps[0]),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="verify the numbers against the two versioned source manuscripts",
    )
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.fetch:
        verify_public_tables()
    result = run()
    if args.write:
        if not args.fetch:
            parser.error(
                "The original observational sources are required to create the receipt"
            )
        with RECEIPT.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        stock.check(json.loads(RECEIPT.read_text()), result)
    print(json.dumps(result["observational_comparison"], indent=2))
