"""Spatial interferometry of the two independently specified fractional walls."""

import argparse
import hashlib
import json
from math import asin, pi, sin
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from ...quark_winding.belt_walls import belt
from ..polar_recoil.measure import angle_error, line_error, unwrap


HERE = Path(__file__).resolve().parent
SIZE = 10
SHOTS = 2**20
WAITS = [0, *(2.0**j for j in range(-4, 9))]
HARMONICS = (1, -1, 2, -2, 3, -3, 4, -4)
SETTINGS = [(n, k, m) for n in (3, 4) for k in (8, 0) for m in HARMONICS]
FLOOR = .95


def wall_indices(sector):
    """Address each sharp wall by the start of the anomalous gap."""
    addresses = sorted((belt.wall(sector.states[i], SIZE)[0], i) for i in belt.selected(sector))
    if [x for x, _ in addresses] != list(range(SIZE)):
        raise ValueError("The source does not give one sharp word at every address")
    return [int(i) for _, i in addresses]


def characters(dimension, indices, harmonic):
    result = np.zeros((dimension, 2), dtype=complex)
    result[indices, 0] = 1 / np.sqrt(SIZE)
    result[indices, 1] = np.exp(2j*pi*harmonic*np.arange(SIZE)/SIZE) / np.sqrt(SIZE)
    return result


def dominant(energies, columns):
    boundaries = np.r_[0, np.where(np.diff(energies) > 1e-9)[0] + 1]
    weights = np.add.reduceat(abs(columns)**2, boundaries, axis=0)
    return weights.max(axis=0)


def acquire():
    rng = np.random.default_rng(301310)
    banks = []
    poles, ports = next(iter(belt.belts().items()))
    for loads in (3, 4):
        sector = BinarySector(belt.ring(SIZE), loads)
        indices = wall_indices(sector)
        for kappa in (8, 0):
            energies, vectors = np.linalg.eigh(sector.hamiltonian(kappa))
            evolution = np.exp(-1j*np.outer(energies, WAITS))
            for harmonic in HARMONICS:
                columns = vectors.T @ characters(len(energies), indices, harmonic)
                weights = dominant(energies, columns)
                prepared = columns.sum(axis=1) / np.sqrt(2)
                waves = vectors @ (evolution * prepared[:, None])
                probabilities = abs(waves)**2
                if np.max(abs(probabilities.sum(axis=0) - 1)) > 1e-11:
                    raise ArithmeticError("Incomplete free release")
                outside = np.ones(len(energies), dtype=bool)
                outside[indices] = False
                outcomes = np.vstack((probabilities[indices], probabilities[outside].sum(axis=0)))
                counts = []
                for panel in outcomes.T:
                    p = panel.round(12)
                    counts.append(rng.multinomial(SHOTS, p / p.sum()).tolist())
                banks.append(dict(loads=loads, kappa=kappa, harmonic=harmonic,
                                  weights=weights.round(11).tolist(), counts=counts))
    return dict(protocol_sha256=hashlib.sha256((HERE/"PROTOCOL.md").read_bytes()).hexdigest(),
                poles=list(poles), belt_ports=list(ports), shots=SHOTS, waits=WAITS, rows=banks)


def recoil(row):
    counts = np.asarray(row["counts"])
    if (counts.shape != (len(WAITS), SIZE+1) or counts.dtype.kind not in "iu"
            or np.any(counts < 0) or np.any(counts.sum(axis=1) != SHOTS)):
        raise ValueError("Incomplete wall or complement counts")
    weights = np.asarray(row["weights"])
    if weights.shape != (2,) or not np.all((weights >= 0) & (weights <= 1)):
        raise ValueError("Invalid source audit")
    modes = np.fft.fft(counts[:, :SIZE] / SHOTS, axis=1)
    strengths = abs(modes[0, 1:SIZE//2+1])
    order = np.argsort(strengths)
    k = int(order[-1] + 1)
    if k != abs(row["harmonic"]) or strengths.max() - strengths[order[-2]] < .01:
        return dict(status="unresolved wavelength")
    if min(weights) < FLOOR:
        return dict(status="unresolved spectral line")
    signal = modes[:, k].conj() if row["harmonic"] > 0 else modes[:, k]
    try:
        errors = []
        for z in signal[1:]:
            # Rounding eleven probabilities at 12 places and renormalizing
            # perturbs their spatial Fourier moment by less than 2e-11.
            radius = abs(z)*sin(angle_error(z, SHOTS)) + 2e-11
            errors.append((z, 2*line_error(FLOOR) + asin(radius/abs(z))))
        interval = unwrap((-16, 16), WAITS[1:], errors)
    except ValueError:
        return dict(status="unresolved temporal phase")
    return dict(status="resolved", interval=interval, harmonic=k)


def forecast(anchor, harmonic, kind):
    n = abs(harmonic)
    ratio = n if kind == "continuum" else sin(pi*n/SIZE) / sin(pi/SIZE)
    extrema = [d*r for d in anchor for r in (ratio, ratio**2)]
    return [min(extrema), max(extrema)]


def read(raw):
    poles, ports = next(iter(belt.belts().items()))
    if (raw["waits"] != WAITS or raw["shots"] != SHOTS or raw["poles"] != list(poles)
            or raw["belt_ports"] != list(ports)):
        raise ValueError("Unrecognized enrollment")
    by_key = {}
    for row in raw["rows"]:
        key = (row["loads"], row["kappa"], row["harmonic"])
        if key in by_key:
            raise ValueError("Repeated launch")
        by_key[key] = recoil(row)
    if set(by_key) != set(SETTINGS):
        raise ValueError("Incomplete launch menu")
    results = []
    for loads, kappa, harmonic in SETTINGS:
        result = by_key[loads, kappa, harmonic]
        anchor = by_key[loads, kappa, 1].get("interval")
        item = dict(loads=loads, kappa=kappa, launch=harmonic, **result)
        if harmonic != 1:
            for kind in ("continuum", "lattice"):
                limits = result.get("interval")
                if anchor is None or anchor[0] <= 0 <= anchor[1] or limits is None:
                    item[kind] = dict(status="unresolved")
                    continue
                lo, hi = forecast(anchor, harmonic, kind)
                gap = max(0, limits[0] - hi, lo - limits[1])
                item[kind] = dict(status="excluded" if gap else "compatible",
                                  envelope=[lo, hi], gap=gap)
        results.append(item)
    return results


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire", type=Path, help="write an unused simulation record")
    cli.add_argument("--counts-only", action="store_true", help="decode the deposited records")
    args = cli.parse_args()
    if args.counts_only and args.acquire:
        cli.error("Choose one mode")
    files = HERE/"counts.json"
    data = json.loads(files.read_text()) if args.counts_only else acquire()
    if data["protocol_sha256"] != hashlib.sha256((HERE/"PROTOCOL.md").read_bytes()).hexdigest():
        raise ValueError("Changed protocol")
    results = read(data)
    if args.acquire:
        with args.acquire.open("x") as dest:
            json.dump(data, dest, separators=(",", ":"), allow_nan=False)
            dest.write("\n")
    elif not args.counts_only and data != json.loads(files.read_text()):
        raise ArithmeticError("Acquisition differs from archive")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
