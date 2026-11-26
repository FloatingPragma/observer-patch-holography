"""Acquire, retain and decode the independent material and voltage receiver banks."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from . import readers, ring

HERE = Path(__file__).parent
SHOTS = 2**22
SEED = 2601120
GRAPHS = ((0, 14), (1, 14), (2, 14), (0, 16))


def calibration(rng):
    bins = {}
    for n in (0, 1, 2):
        phases = np.array(
            [readers.phase(n, t) * readers.phase(0, t).conjugate() for t in ring.ANGLES]
        )
        p = (1 + np.stack((phases.real, phases.imag), axis=1)) / 2
        bins[n] = {
            "bright": rng.binomial(SHOTS, p.clip(0, 1)).tolist(),
            "probability": p.tolist(),
        }
    term = ring.candidate.native_law.seam_term(ring.K)
    swap = expm(-1j * np.pi * term)
    for n in range(3):
        i, j = [ring.candidate.native_law.SRC.index(w) for w in ((0, n), (n, 0))]
        assert abs(swap[j, i]) > 1 - 1e-12
    return bins


def inferred_phase(rows):
    axis = 2 * np.array(rows["bright"]) / SHOTS - 1
    angles = np.arctan2(axis[:, 1], axis[:, 0])
    return float((angles[1] - angles[2]) / (np.pi / 2))


def read_motion(meter, standards):
    base = inferred_phase(standards[1] if 1 in standards else standards["1"])
    inferred = {}
    for key, scans in meter.items():
        inferred[key] = [
            {
                "charge_per_wall_step": inferred_phase(s["ports"]["wall"]) / (2 * base),
                "ordinary_per_site_step": inferred_phase(s["ports"]["standard"]) / base,
                "parked_per_wall_step": inferred_phase(s["ports"]["parked"])
                / (2 * base),
            }
            for s in scans
        ]
    return inferred


def acquire():
    rng = np.random.default_rng(SEED)
    bank = {}
    graphs = []
    for level, n in GRAPHS:
        station = readers.stations(level, n)
        cal = calibration(rng)
        meter = readers.interferometry(n, rng, SHOTS)
        board = dict(
            level=level,
            n=n,
            receivers=station,
            standards=cal,
            held_motion=meter,
            inferred_motion=read_motion(meter, cal),
            release={},
        )
        graphs.append(board)
        base, sites = ring.basis(n)
        bank[f"words_{n}"] = base
        for amount in (1, 2):
            for sign in (1, -1):
                gate, two, tapered, flat = ring.regional(n, sign, amount)
                table = {}
                board["release"][f"{amount}_{sign}"] = table
                for code, k in (("native", ring.K), ("contrast_off", 0)):
                    law = ring.released(n, amount, k)
                    bank[f"law_{n}_{amount}_{code}"] = law.astype(np.float32)
                    counts = np.stack(
                        [rng.multinomial(SHOTS, p / p.sum()) for p in law]
                    )
                    bank[f"counts_{level}_{n}_{amount}_{sign}_{code}"] = counts.astype(
                        np.uint32
                    )
                    table[code] = {
                        "joint_gate_counts": counts[:, gate].sum(axis=1).tolist(),
                        "two_wall_counts": counts[:, two].sum(axis=1).tolist(),
                        "joint_gate_probabilities": (law @ gate).tolist(),
                        "two_wall_probabilities": (law @ two).tolist(),
                    }
                p = ring.released(n, amount, ring.K)[4]
                table["two_voltage_probes"] = readers.receiver_phases(
                    n, amount, sign, p, gate, tapered, flat, rng, SHOTS
                )
        # Separate sharp preparations of each translated alternating bulk phase.
        ref = {}
        for amount in (1, 2):
            for parity in (0, 1):
                sharp = amount * ((np.arange(n) + parity) % 2)
                fields = []
                for arm in ((-5, -4, -3), (2, 3, 4)):
                    shape = np.zeros(n)
                    for j in arm:
                        shape[j % n] += 0.5
                        shape[(j + 1) % n] += 0.5
                    source = readers.phase_word(sharp, shape, 2 * np.pi)
                    zero = readers.phase_word(np.zeros(n, int), shape, 2 * np.pi)
                    flat = readers.phase_word(sharp, (shape > 0).astype(int), 2 * np.pi)
                    flat0 = readers.phase_word(
                        np.zeros(n, int), (shape > 0).astype(int), 2 * np.pi
                    )
                    p = (
                        1
                        + np.real([source * zero.conjugate(), flat * flat0.conjugate()])
                    ) / 2
                    fields.append(rng.binomial(SHOTS, p.clip(0, 1)).tolist())
                ref[f"{amount}_{parity}"] = fields
        board["bulk_reference_full_turn_bright_taper_and_flat"] = ref
    manifest = {
        "seed": SEED,
        "shots": SHOTS,
        "holds": ring.HOLDS.tolist(),
        "gradient_angles": ring.ANGLES.tolist(),
        "regional_angles": (np.pi * np.array([0, 0.5, 1, 2])).tolist(),
        "word_layout": "positive occupancies indexed by sorted high-count words; negative occupancies are their bitwise complements within n sites; occupied sites carry the enrolled amount",
        "regional_bins": "accept=no,yes, then left and right bright,dark; both translated sharp bulk preparations give the same relative phase",
        "origin": "Seeded binomial phase and multinomial joint counts drawn from calculated finite native Born probabilities. Sharp digit codes come from two restored routed native receiver heads.",
        "boards": graphs,
    }
    return manifest, bank


def same_record(actual, stored):
    """Counts and labels replay exactly; calculated phases may differ by roundoff."""
    if isinstance(stored, dict):
        return (
            isinstance(actual, dict)
            and actual.keys() == stored.keys()
            and all(same_record(actual[k], v) for k, v in stored.items())
        )
    if isinstance(stored, list):
        return (
            isinstance(actual, list)
            and len(actual) == len(stored)
            and all(same_record(a, b) for a, b in zip(actual, stored))
        )
    if isinstance(stored, float):
        return bool(np.isclose(actual, stored, atol=3e-11, rtol=3e-11))
    return actual == stored


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    writing = parser.parse_args().write
    report, bank = acquire()
    raw = json.dumps(report, indent=1) + "\n"
    if writing:
        (HERE / "records.json").write_text(raw)
        np.savez_compressed(HERE / "words.npz", **bank)
    else:
        with np.load(HERE / "words.npz") as saved:
            if set(bank) != set(saved.files) or any(
                (
                    not np.allclose(bank[k], saved[k], rtol=3e-7, atol=3e-12)
                    if k.startswith("law_")
                    else not np.array_equal(bank[k], saved[k])
                )
                for k in bank
            ):
                raise ArithmeticError("Raw word replay changed")
        if not same_record(
            json.loads(raw), json.loads((HERE / "records.json").read_text())
        ):
            raise ArithmeticError("Voltage or decoded record changed")
    for b in report["boards"]:
        print(
            b["level"],
            b["n"],
            b["release"]["1_1"]["native"]["joint_gate_counts"][4],
            b["inferred_motion"]["1_1"],
        )
