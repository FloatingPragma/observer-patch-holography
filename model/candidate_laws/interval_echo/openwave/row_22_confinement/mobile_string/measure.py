"""Manufacture seeded readout counts, or decode the stored receiver banks."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import material as mat

HERE = Path(__file__).parent
SHOTS = 1 << 22
LAUNCHES = 1 << 24
SEED = 2026112622
TICK = np.arange(257) / 8
CAL_HOLDS = np.array([np.pi / 128, np.pi / 4, np.pi / 32])
FREE = {"native": (mat.G, mat.H), "degree_off": (0, mat.H), "exchange_off": (mat.G, 0)}


def sample_ports(rng, signal):
    return rng.binomial(SHOTS, np.clip((1 + signal) / 2, 0, 1)).astype(np.uint32)


def quadratures(theta):
    return np.stack((np.cos(theta), np.sin(theta)))


def acquire():
    rng = np.random.default_rng(SEED)
    archive = {}
    calibration = mat.calibration_energies()
    np.testing.assert_allclose(
        calibration, [4 * mat.J, 2 * mat.H, mat.G], atol=1e-12, rtol=0
    )
    for name, (g, h) in FREE.items():
        law = mat.free_laws(g, h)
        archive["word_law_" + name] = law.astype(np.float32)
        for lev in range(3):
            counts = np.empty_like(law, dtype=np.uint32)
            for i, launches in enumerate(law):
                for k, p in enumerate(launches):
                    counts[i, k] = rng.multinomial(LAUNCHES, p / p.sum())
            archive[f"word_{lev}_{name}"] = counts
    for lev in range(3):
        archive[f"ports_{lev}"] = mat.route(lev)
        archive[f"clock_{lev}"] = sample_ports(rng, quadratures(calibration[1] * TICK))
        for name, g in (("native", mat.G), ("degree_off", 0)):
            energies = calibration.copy()
            if name == "degree_off":
                energies[-1] = 0
            archive[f"cal_{lev}_{name}"] = sample_ports(
                rng, quadratures(energies * CAL_HOLDS)
            )
            bank, grounds, weights = [], [], []
            for length in mat.LENGTHS:
                ground, panels = mat.static_laws(length, g)
                bank.append([sample_ports(rng, x) for x, _ in panels])
                grounds.append(ground)
                weights.append([w for _, w in panels])
            archive[f"static_{lev}_{name}"] = bank
            archive["ground_" + name] = grounds
            archive["weight_" + name] = weights
    np.savez_compressed(HERE / "receivers.npz", **archive)


def port_phase(counts):
    return np.arctan2(
        2 * counts[1].astype(float) / SHOTS - 1, 2 * counts[0].astype(float) / SHOTS - 1
    )


def ramsey(counts, clock):
    angles = port_phase(counts)
    rate = angles[0] / mat.HOLDS[0]
    for phase, t in zip(angles[1:], mat.HOLDS[1:]):
        rate = (phase + 2 * np.pi * round((t * rate - phase) / (2 * np.pi))) / t
    visibility = np.hypot(*(2 * counts.astype(float) / SHOTS - 1)).min()
    return {
        "line": round(rate / clock, 7),
        "min_visibility": round(float(visibility), 7),
    }


def decode():
    with np.load(HERE / "receivers.npz") as arc:
        levels = {}
        masks = mat.event_masks()
        for lev in range(3):
            rate = float(
                np.polyfit(TICK, np.unwrap(port_phase(arc[f"clock_{lev}"])), 1)[0]
            )
            panel = {"physical_ports": arc[f"ports_{lev}"].tolist(), "clock_rate": rate}
            levels[lev] = panel
            for name in FREE:
                data = panel[name] = {}
                if name != "exchange_off":
                    cal = port_phase(arc[f"cal_{lev}_{name}"]) / (CAL_HOLDS * rate)
                    data["contact_4J_swap_2h_private_g"] = cal.tolist()
                    grounds = arc["ground_" + name]
                    weights = arc["weight_" + name]
                    static = data["static"] = {}
                    for i, (length, plates) in enumerate(
                        zip(mat.LENGTHS, arc[f"static_{lev}_{name}"])
                    ):
                        reads = [ramsey(counts, rate) for counts in plates]
                        for j, row in enumerate(reads):
                            row["ground_product_weight"] = float(weights[i, j])
                        static[length] = {
                            "calculated_ground_work": float(grounds[i]),
                            "sharp_intact": reads[0],
                            "sharp_screened": reads[1],
                            "calibrated_zero_swap_prediction": min(
                                cal[0] + cal[2] * length, 2 * cal[0] + 2 * cal[2]
                            ),
                        }
                    slope = (
                        static[3]["sharp_intact"]["line"]
                        - static[2]["sharp_intact"]["line"]
                    )
                    data["slope_2_to_3"] = slope
                    data["held_intact_line_minus_extrapolated"] = {
                        q: static[q]["sharp_intact"]["line"]
                        - (static[2]["sharp_intact"]["line"] + (q - 2) * slope)
                        for q in (4, 5)
                    }
                laws = arc["word_law_" + name]
                words = arc[f"word_{lev}_{name}"]
                event_table = data["release"] = []
                for i, (law, bank) in enumerate(zip(laws, words)):
                    row = []
                    for k, (launch, counts) in enumerate(zip(law, bank)):
                        assert int(counts.sum()) == LAUNCHES
                        row.append(
                            [
                                {
                                    "count": int(counts[m].sum()),
                                    "probability": float(
                                        launch[m].sum(dtype=np.float64)
                                    ),
                                }
                                for m in masks[k]
                            ]
                        )
                    event_table.append(row)
        return {
            "protocol": "a77043d82",
            "seed": SEED,
            "ramsey_trials_per_quadrature": SHOTS,
            "word_launches_per_preparation": LAUNCHES,
            "free_holds": mat.RELEASE.tolist(),
            "record_order": "upper and reversed launches; signed-endpoint survival, extra interior pair with surviving endpoints, enclosed opposite four-site island with surviving endpoints",
            "origin": "seeded binomial quadratures and multinomial word counts from calculated exact finite native-sector Born laws; spectra and eigenweights independently labelled",
            "levels": levels,
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquire", action="store_true")
    if parser.parse_args().acquire:
        acquire()
    result = decode()
    (HERE / "results.json").write_text(json.dumps(result, indent=1) + "\n")
