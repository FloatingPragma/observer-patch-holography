"""Compare two sender resets in the half-filled field and its held matter excitation."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30.experiment import select_route
from model.executable_model.carrier import federation
from ..measure import ring_modes, tangents
from ..selection import measure as selection

HERE = Path(__file__).resolve().parent
DISTANCES = (1, 4, 8, 16, 51, 55)
SHOTS = 2**18
SEED = 246055


def reset(c, bit):
    out = c.copy()
    out[0, :] = 0
    out[:, 0] = 0
    out[0, 0] = bit
    return out


def exact_words(interaction, close=False, native_diagonal=False):
    """Full, four-site binary word evolution from two independently specified backgrounds."""
    h = np.zeros((16, 16))
    bonds = [(0, 1), (1, 2), (2, 3)] + ([(3, 0)] if close else [])
    for word in range(16):
        for a, b in bonds:
            mask = (1 << a) | (1 << b)
            n = (word & mask).bit_count()
            if n == 1:
                h[word ^ mask, word] = -0.5
                if native_diagonal:
                    h[word, word] += candidate.KAPPA + 0.5
            elif n == 2:
                h[word, word] += interaction
    u = expm(-3j * h)
    receiver_occupied = (np.arange(16) & (1 << 2)) != 0
    rows = []
    for reference in (0, 1 << 1):
        p = [np.sum(abs(u[receiver_occupied, reference | bit])**2) for bit in (0, 1)]
        rows.append({"sender_empty": p[0], "sender_full": p[1], "excess": p[1] - p[0]})
    return rows


def run(archive=None):
    board = federation(2)
    route, outside = select_route(board)
    sender, spare = route[0], outside[0]
    if spare in route or spare not in board.neighbours[sender]:
        raise ArithmeticError("The sender reset lacks an adjacent independent rail")
    k, modes, filled = ring_modes(len(route))
    source_receipt = json.loads((selection.HERE / "results.json").read_text())
    intervals = selection.read_calibration(source_receipt["counts_symmetric_antisymmetric_quiet"])
    lower = {i for i, (_, hi) in enumerate(intervals) if hi < 0}
    cooled = selection.trajectory("inverted", -np.cos(k), selection.work_choices(intervals), lower, intervals)
    if cooled["unfinished"] or cooled["final_upper_count"] or lower != set(np.flatnonzero(filled)):
        raise ArithmeticError("The enrolled cooling trajectory did not prepare the matter reference")
    basis = modes[:, sorted(lower)]
    cold = basis @ basis.conj().T
    variations = tangents(k, modes, filled)
    hole, particle = variations["coherent_pair"]
    excited = cold + np.outer(particle, particle.conj()) - np.outer(hole, hole.conj())
    h = modes[:, np.flatnonzero(filled)[-1]]
    plus = variations["signed_rotation"][1]
    rotation = cold + np.outer(plus, plus.conj()) - np.outer(h, h.conj())
    sources = {"cold": cold, "held_replacement": excited, "held_rotation_pi_over_four": rotation}
    pairs = {name: [reset(c, bit) for bit in (0, 1)] for name, c in sources.items()}
    reference = []
    readings = {name: [] for name in sources}
    for t in range(81):
        phase = np.exp(1j * t * np.cos(k))
        rows = (modes[list(DISTANCES)] * phase) @ modes.conj().T
        forecast = abs(np.exp(1j * np.outer(DISTANCES, k)) @ phase / len(route))**2
        reference.append(forecast)
        for name, preps in pairs.items():
            values = np.stack([np.einsum("ij,ij->i", rows @ c, rows.conj()).real for c in preps])
            readings[name].append(values)
    reference = np.array(reference)
    readings = {name: np.array(rows) for name, rows in readings.items()}
    errors = {name: np.max(abs(rows[:, 1] - rows[:, 0] - reference)) for name, rows in readings.items()}
    if max(errors.values()) > 2e-12:
        raise ArithmeticError("Independent reset arms disagree with the signal forecast")
    peaks = []
    if archive is None:
        rng = np.random.default_rng(SEED)
        counts = {name: rng.binomial(SHOTS, probs).astype(np.int32) for name, probs in readings.items()}
    else:
        counts = {}
        for name, probs in readings.items():
            np.testing.assert_allclose(probs, archive[name + "_rate"], rtol=2e-12, atol=2e-13)
            raw = archive[name + "_count"]
            if raw.shape != probs.shape or raw.dtype.kind not in "iu" or raw.min() < 0 or raw.max() > SHOTS:
                raise ArithmeticError("A stored terminal count has an invalid shape or range")
            counts[name] = raw
    cells = sum(raw.size for raw in counts.values())
    radius = np.sqrt(np.log(2 * cells / .01) / (2 * SHOTS))
    max_count_error = max(np.max(abs(counts[name] / SHOTS - probs)) for name, probs in readings.items())
    if max_count_error > radius:
        raise ArithmeticError("Some original-denominator rates fail the simultaneous frequency test")
    for j, dist in enumerate(DISTANCES):
        t = int(np.argmax(reference[:, j]))
        at_peak = {}
        for name, rows in readings.items():
            rates = rows[t, :, j]
            excess = (counts[name][:, 1, j] - counts[name][:, 0, j]) / SHOTS
            possible = np.flatnonzero(excess >= excess.max() - 4 * radius)
            at_peak[name] = {
                "rates_empty_full": rates.tolist(), "occupied_counts": counts[name][t, :, j].tolist(),
                "count_peak": int(np.argmax(excess)), "compatible_peak_ages_99_percent": possible.tolist(),
            }
        peaks.append({"separation": dist, "receiver_port": route[dist], "peak_clock": t,
                      "forecast": reference[t, j], "sources": at_peak})
    return {
        "ring_stations": len(route), "sender_port": sender, "reset_rail": spare,
        "cooling_from_inverted_filling": {key: cooled[key] for key in ("initial_upper_count", "pair_proposals", "final_upper_count")},
        "ages_inclusive": [0, 80], "shots_per_arm": SHOTS, "seed": SEED,
        "tested_arm_probabilities": 81 * len(DISTANCES) * 2 * len(sources),
        "simultaneous_single_rate_radius_99_percent": radius,
        "maximum_single_rate_count_discrepancy": max_count_error,
        "maximum_absolute_signal_error": errors, "peaks": peaks,
        "pre_reset_bond_energy_changes": {
            str(b): -np.real(excited[b, b-1] - cold[b, b-1]) for b in (51, 55)
        },
        "maximum_single_arm_change_for_rotated_source": np.max(abs(readings["held_rotation_pi_over_four"] - readings["cold"])),
        "four_site_at_time_three": {
            "open_quadratic": exact_words(0),
            "open_density_interaction": exact_words(1),
            "open_native_blank_contacts": exact_words(0, native_diagonal=True),
            "literal_spin_ring_closure": exact_words(0, close=True),
        },
    }, counts, readings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    archive_path = HERE / "counts.npz"
    if args.write:
        out, counts, rates = run()
    else:
        with np.load(archive_path, allow_pickle=False) as archive:
            out, counts, rates = run(archive)
    for name, val in out["maximum_absolute_signal_error"].items():
        out["maximum_absolute_signal_error"][name] = float(format(val, ".3g"))
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(archive_path, **{k + "_count": v for k, v in counts.items()},
                            **{k + "_rate": v for k, v in rates.items()})
    else:
        prior = json.loads(path.read_text())
        for row, old in zip(out["peaks"], prior["peaks"], strict=True):
            if (row["separation"], row["receiver_port"], row["peak_clock"]) != (old["separation"], old["receiver_port"], old["peak_clock"]):
                raise ArithmeticError("Peak appointment or receiver assignment changed")
            np.testing.assert_allclose(row["forecast"], old["forecast"], rtol=0, atol=1e-13)
            for name, bank in row["sources"].items():
                np.testing.assert_allclose(bank["rates_empty_full"], old["sources"][name]["rates_empty_full"], atol=1e-12, rtol=0)
                np.testing.assert_array_equal(bank["occupied_counts"], old["sources"][name]["occupied_counts"])
                if bank["count_peak"] != old["sources"][name]["count_peak"]:
                    raise ArithmeticError("Peak from the raw counts differs from the summary")
                np.testing.assert_array_equal(bank["compatible_peak_ages_99_percent"], old["sources"][name]["compatible_peak_ages_99_percent"])
        for name, numbers in out["four_site_at_time_three"].items():
            for arm, saved in zip(numbers, prior["four_site_at_time_three"][name], strict=True):
                np.testing.assert_allclose(list(arm.values()), list(saved.values()), atol=1e-13, rtol=0)
    print("Verified the two-arm signal forecasts and the interacting four-site control.")


if __name__ == "__main__":
    main()
