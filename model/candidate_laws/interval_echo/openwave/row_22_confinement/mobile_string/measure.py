"""Manufacture seeded readout counts, or decode the stored receiver banks."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import material as mat

HERE = Path(__file__).parent
SHOTS = 1 << 22
# Component Hoeffding radius 5/sqrt(N), enclosing complex-disk radius below.
VISIBILITY_RADIUS = 5 * np.sqrt(2 / SHOTS)
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
    """Unwrap in the first hold's principal interval, with a resolved-turn check."""
    angles = port_phase(counts)
    visibility = np.hypot(*(2 * counts.astype(float) / SHOTS - 1))
    weak = visibility <= VISIBILITY_RADIUS
    angular_error = np.full_like(visibility, np.pi)
    angular_error[~weak] = np.arcsin(VISIBILITY_RADIUS / visibility[~weak])
    ambiguous_holds = []
    if abs(angles[0]) + angular_error[0] >= np.pi:
        ambiguous_holds.append(float(mat.HOLDS[0]))
    rate = angles[0] / mat.HOLDS[0]
    for i, t in enumerate(mat.HOLDS[1:], 1):
        predicted = t * rate
        phase = angles[i] + 2 * np.pi * round((predicted - angles[i]) / (2 * np.pi))
        uncertainty = angular_error[i] + t / mat.HOLDS[i - 1] * angular_error[i - 1]
        if abs(phase - predicted) + uncertainty >= np.pi:
            ambiguous_holds.append(float(t))
        rate = phase / t
    attempted = round(float(rate / clock), 7)
    row = {"min_visibility": round(float(visibility.min()), 7)}
    if weak.any() or ambiguous_holds:
        return {**row, "status": "unresolved", "line": None, "attempted_line": attempted,
                "low_visibility_holds": mat.HOLDS[weak].tolist(), "ambiguous_turn_holds": ambiguous_holds}
    return {**row, "status": "resolved", "line": attempted}


def intact_trend(static):
    first = static[2]["sharp_intact"]["line"]
    second = static[3]["sharp_intact"]["line"]
    slope = second - first if first is not None and second is not None else None
    held = {}
    for length in (4, 5):
        line = static[length]["sharp_intact"]["line"]
        held[length] = line - (first + (length - 2) * slope) if line is not None and slope is not None else None
    return slope, held


def check_analysis(current, stored):
    """Match every outcome and integer exactly, allowing numerical roundoff in computed floats."""
    if type(current) is not type(stored):
        raise ValueError("stored analysis changed a value type")
    if isinstance(current, dict):
        if current.keys() != stored.keys():
            raise ValueError("stored analysis changed its fields")
        for key, value in current.items():
            check_analysis(value, stored[key])
    elif isinstance(current, list):
        for left, right in zip(current, stored, strict=True):
            check_analysis(left, right)
    elif isinstance(current, float):
        np.testing.assert_allclose(current, stored, rtol=0, atol=1e-10, equal_nan=False)
    elif current != stored:
        raise ValueError(f"stored analysis differs: {stored!r}, expected {current!r}")


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
                    data["slope_2_to_3"], data["held_intact_line_minus_extrapolated"] = intact_trend(static)
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
            "ramsey_visibility_radius": float(VISIBILITY_RADIUS),
            "ramsey_analysis": "A line requires every quadrature-vector disk to exclude zero and every nearest-turn choice to survive the angular error bounds. An unresolved attempt is retained only as a diagnostic.",
            "word_launches_per_preparation": LAUNCHES,
            "free_holds": mat.RELEASE.tolist(),
            "record_order": "upper and reversed launches; signed-endpoint survival, extra interior pair with surviving endpoints, enclosed opposite four-site island with surviving endpoints",
            "origin": "seeded binomial quadratures and multinomial word counts from calculated exact finite native-sector Born laws; spectra and eigenweights independently labelled",
            "levels": levels,
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--acquire", action="store_true")
    mode.add_argument("--check", action="store_true", help="verify the current analysis of the deposited bank")
    args = parser.parse_args()
    if args.acquire:
        acquire()
    result = decode()
    current = HERE / "results.json"
    if args.check:
        check_analysis(json.loads(json.dumps(result)), json.loads(current.read_text()))
    else:
        current.write_text(json.dumps(result, indent=1) + "\n")
    panels = [rows[preparation] for level in result["levels"].values()
              for law in FREE if "static" in level[law] for rows in level[law]["static"].values()
              for preparation in ("sharp_intact", "sharp_screened")]
    resolved = sum(row["status"] == "resolved" for row in panels)
    print(f"Ramsey analysis: {resolved} resolved and {len(panels) - resolved} unresolved preparations.")
