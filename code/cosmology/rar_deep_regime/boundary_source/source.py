"""Construct outer apertures and the source directly from surface densities."""

import csv
import itertools
from collections import defaultdict
from pathlib import Path

import numpy as np

from cosmology.rar_deep_regime.rar_deep_regime_diagnostic import DATA, read_table1

HERE = Path(__file__).resolve().parent
THRESHOLD_MAX = 10_000.0  # solar masses per pc²
AREA_UNITS = 1e6  # convert kpc² to pc²


def table(filename):
    with (HERE / filename).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def catalog_id(name):
    if name == "Ho II":
        return "UGC04305"
    if name.startswith("NGC"):
        return f"NGC{int(name[3:]):04d}"
    return name.replace(" ", "")


def edges(centres):
    return np.r_[
        0,
        0.5 * (centres[:-1] + centres[1:]),
        centres[-1] + (centres[-1] - centres[-2]) / 2,
    ]


def filled(radius, rows, column):
    values = np.array([float(row[column]) if row[column] else np.nan for row in rows])
    measured = np.isfinite(values)
    if column == "SigmaH2":
        return np.nan_to_num(values), int(measured.sum())
    if not np.any(measured):
        raise ValueError(f"Unmeasured {column}")
    return np.interp(radius, radius[measured], values[measured]), int(measured.sum())


def rebin(old_edges, surface, width):
    """Exact area-weighted field on radial intervals, including the partial edge."""
    new = np.r_[np.arange(0, old_edges[-1] - 1e-10, width), old_edges[-1]]
    lo = np.maximum(new[:-1, None], old_edges[None, :-1])
    hi = np.minimum(new[1:, None], old_edges[None, 1:])
    weights = np.maximum(0, hi**2 - lo**2) / np.diff(new**2)[:, None]
    return new, weights @ surface


def shell_factor(radius, shells, strength):
    """Integrate (r u(r))² exactly over each linear segment."""
    w = strength / sum(strength)
    nodes = np.r_[0, shells[shells < radius], radius]
    values = np.minimum(1, nodes[:, None] / shells[None, :]) @ w
    return float(
        np.sum(
            np.diff(nodes)
            * (values[:-1] ** 2 + values[:-1] * values[1:] + values[1:] ** 2)
        )
        / (3 * radius)
    )


def variation(ring_edges, densities):
    sigma = np.clip(densities, 0, THRESHOLD_MAX)
    return 2 * np.pi * AREA_UNITS * ring_edges[1:] * abs(np.diff(np.r_[sigma, 0]))


def profile(
    meta,
    leroy,
    rows,
    orbit,
    width=0.8,
    star=1.0,
    atomic=1.0,
    molecular=1.0,
    distance=None,
):
    distance = meta["dist_mpc"] if distance is None else distance
    scale = distance / float(leroy["Dist"])
    radii = np.array([float(r["r"]) for r in rows]) * scale
    if np.any(np.diff(radii) <= 0):
        raise ValueError("Radii must increase")
    fields = [filled(radii, rows, s) for s in ("Sigma*", "SigmaHI", "SigmaH2")]
    surface = np.column_stack([f[0] for f in fields]) * [star, atomic, molecular]
    ring_edges, coarse = rebin(edges(radii), surface, width)
    if np.min(coarse) < 0:
        raise ValueError("Negative surface density")
    sigma = coarse.sum(axis=1)
    masses = np.pi * AREA_UNITS * np.diff(ring_edges**2) @ coarse
    strength = variation(ring_edges, sigma)
    total = float(sum(strength))
    rd = float(leroy["l*"]) * scale
    return {
        "mass_msun": float(sum(masses)),
        "stellar_msun": masses[0],
        "atomic_with_he_msun": masses[1],
        "molecular_with_he_msun": masses[2],
        "n_input_rings": len(rows),
        "n_measured_star_HI_H2": [f[1] for f in fields],
        "max_sigma_msun_pc2": max(sigma),
        "edge_kpc": ring_edges[-1],
        "rd_kpc": rd,
        "S_msun_kpc": total,
        "F": shell_factor(orbit, ring_edges[1:], strength),
        "terminal_jump_fraction": strength[-1] / total,
        "atomic_fraction_of_sparc": masses[1]
        / (1.36 * meta["MHI_gsun"] * 1e9 * (distance / meta["dist_mpc"]) ** 2),
    }


def enroll():
    """Read only ID and radius from the orbital file in this stage."""
    catalogue = {r["name"]: r for r in read_table1()}
    last = defaultdict(float)
    for line in (DATA / "table2.dat").read_text().splitlines():
        name = line[0:11].strip()
        last[name] = max(last[name], float(line[19:25]))
    rows = defaultdict(list)
    for r in table("leroy_table7.dat"):
        rows[r["Name"]].append(r)
    kept, dropped = {}, {}
    for lrow in table("leroy_sample.dat"):
        key = catalog_id(lrow["Name"])
        if key not in catalogue:
            dropped[key] = "no SPARC counterpart"
            continue
        meta = catalogue[key]
        if meta["quality"] > 2 or meta["inclination_deg"] < 30:
            dropped[key] = "quality or inclination"
            continue
        src = profile(meta, lrow, rows[lrow["Name"]], last[key])
        if last[key] < max(src["edge_kpc"], 3 * src["rd_kpc"]):
            dropped[key] = "orbit interior to aperture or 3 Rd"
            continue
        kept[key] = (meta, lrow, rows[lrow["Name"]], last[key], src)
    return kept, dropped


def assignments(sources):
    order = sorted(sources, key=lambda k: (sources[k]["mass_msun"], k))
    candidates = []
    for a, b in itertools.combinations(order, 2):
        if sources[b]["mass_msun"] > 2 * sources[a]["mass_msun"]:
            continue
        intensity = {
            k: sources[k]["S_msun_kpc"] / np.sqrt(sources[k]["mass_msun"])
            for k in (a, b)
        }
        hi, lo = sorted([a, b], key=intensity.get, reverse=True)
        candidates.append((intensity[hi] / intensity[lo], hi, lo))
    _, hi, lo = max(candidates)
    return order[::2], order[1::2], (hi, lo)
