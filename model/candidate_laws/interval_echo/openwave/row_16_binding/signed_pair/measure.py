"""Held signed binding and quadrupole records on a fixed native family."""

import argparse
import gzip
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.executable_model import carrier
from ..apparatus import NATIVE_LINE, time_port
from .material import competing_orbit, Assembly, axes, bound, decompositions, vertices
from .aperture import flights, WAITS
from .probe import SHOTS, Quadrupole, calibrate, draw

HERE = Path(__file__).resolve().parent


def cap_bank(rng, gain):
    sec = BinarySector(carrier.cell(), 6)
    h = sec.hamiltonian(8, False)
    e, v = eigh(h)
    coords = vertices()
    frame = axes(coords[0])
    indices = [sec.index[tuple(sorted([i, *sec.car.neighbours[i]]))] for i in range(12)]
    orbit = np.zeros((len(e), 3))
    orbit[indices] = coords @ frame.T / 2  # sum_i r_i r_i^T=4I
    assert np.allclose(orbit.T @ orbit, np.eye(3))
    low = e < 144
    exact = v[:, low] @ (v[:, low].T @ orbit)
    singular = np.linalg.svd(exact, compute_uv=False)
    energy = np.sum(exact * (h @ exact), axis=0) / np.sum(exact**2, axis=0)
    out = v @ (time_port(e, NATIVE_LINE)[:, None] * (v.T @ orbit))
    words = np.zeros((len(e), 12), int)
    for i, w in enumerate(sec.states):
        words[i, list(w)] = 1
    reader = Quadrupole(words, words, frame, -1, gain)
    p0 = abs(out[:, 0]) ** 2 / 12
    p1 = abs((out[:, 1] + 1j * out[:, 2]) / np.sqrt(2)) ** 2 / 12
    sharp = np.zeros(len(e))
    sharp[indices[0]] = 1
    banks = {
        k: reader.bank(p, rng)
        for k, p in [
            ("axial_0", p0),
            ("circular_plus", p1),
            ("isotropic", (p0 + 2 * p1) / 3),
            ("sharp_cap", sharp),
        ]
    }
    return {
        "vector_singular_values": singular.tolist(),
        "vector_energies": energy.tolist(),
        "projected_energy_residual": float(np.linalg.norm(h @ exact - exact * energy)),
        "above_aperture_weight": float(
            np.max(np.sum(abs(v[:, ~low].T @ out) ** 2, axis=0))
        ),
        "probability_per_source": float(p0.sum()),
        "banks": banks,
    }


def run():
    rng = np.random.default_rng(161401)
    gain, calibration = calibrate(rng)
    gain_cal = (calibration["counts"][2][0] - calibration["counts"][2][1]) / SHOTS
    inventory = {}
    for a, b in product(range(5), range(4)):
        if a + b:
            inventory[f"{a},{b}"] = bound(a, b)
    threshold = inventory["2,1"][0] + inventory["2,2"][0]
    banks = {}
    for k in (8, 0):
        assembly = Assembly(4, 3, k)
        seeds = assembly.sources()
        masks, like, opposite = decompositions(assembly.w)
        reader = Quadrupole(assembly.w + 1, abs(assembly.w), assembly.axes, 0, gain)
        data, degree = flights(assembly.h, assembly.bracket, seeds, 209 if k else 0)
        sample = {}
        for j, name in enumerate(("sharp", "residue_0", "residue_1")):
            src = seeds[:, j]
            acceptance = float(np.vdot(src, src).real)
            if acceptance < 1e-12:
                sample[name] = {"source_probability": acceptance}
                continue
            unit = src / np.sqrt(acceptance)
            mean = np.vdot(unit, assembly.h @ unit).real
            variance = np.linalg.norm(assembly.h @ unit - mean * unit) ** 2
            row = {
                "source_probability": acceptance,
                "mean_energy": float(mean),
                "variance": float(variance),
                "contact_bank": assembly.contact_bank(unit, rng, SHOTS),
                "above_threshold_bound": float(
                    variance / (variance + (threshold - mean) ** 2)
                )
                if mean < threshold
                else 1,
                "flight": [],
            }
            for i, t in enumerate(WAITS):
                for select in (0, 1):
                    p = abs(data[2 * i + select, :, j]) ** 2
                    assert p.sum() <= 1 + 1e-10
                    topo = (masks.T @ (2 ** np.arange(4))).astype(int)
                    x = np.bincount(topo, weights=p, minlength=16)
                    row["flight"].append(
                        {
                            "wait": t,
                            "energy_port": bool(select),
                            "shape_census": draw(rng, [*x, 1 - p.sum()]),
                            "mean_equal_edges": float(p @ like / p.sum()),
                            "mean_opposite_edges": float(p @ opposite / p.sum()),
                            "quadrupole": reader.bank(p, rng),
                        }
                    )
            sample[name] = row
        banks[str(k)] = {
            "joined_dimension": len(assembly.w),
            "sharp": assembly.w[assembly.i].tolist(),
            "orbit": assembly.orbit,
            "degree": degree,
            "bracket": assembly.bracket,
            "sources": sample,
        }
    separated = {}
    for name, plus, minus in [
        ("positive", 2, 1),
        ("neutral", 2, 2),
        ("four_excess", 4, 0),
        ("three_vacancy", 0, 3),
    ]:
        a = Assembly(plus, minus)
        src = a.sources()[:, 0]
        q = plus - minus
        phase = 260 * np.pi / 31 * q
        # The same source pointer visits the entire cell before its Y analysis.
        address_sign = np.sign(np.sin(260 * np.pi / 31))
        ans = address_sign * np.sin(phase)
        separated[name] = {
            "sector": (plus, minus),
            "word": a.w[a.i].tolist(),
            "pointer_y_counts": draw(rng, [(1 + ans) / 2, (1 - ans) / 2]),
            "contact_bank": a.contact_bank(src, rng, SHOTS),
        }
    return {
        "shots": SHOTS,
        "calibration": calibration,
        "unit_gain": gain,
        "observed_unit_gain": gain_cal,
        "inventories": inventory,
        "independent": separated,
        "competing_rotational_orbit": competing_orbit(),
        "joined": banks,
        "cap": cap_bank(rng, gain),
    }


def normalize(result):
    if isinstance(result, dict):
        return {k: normalize(v) for k, v in result.items()}
    if isinstance(result, (list, tuple)):
        return [normalize(v) for v in result]
    if isinstance(result, float):
        return round(result, 10)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    obj = normalize(run())
    if args.output:
        args.output.write_bytes(
            gzip.compress(
                (json.dumps(obj, separators=(",", ":")) + "\n").encode(), mtime=0
            )
        )
    else:
        assert (
            json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes())) == obj
        )
        print(
            "Replayed the signed stocks, energy and shape banks, and orbital-vector control."
        )
