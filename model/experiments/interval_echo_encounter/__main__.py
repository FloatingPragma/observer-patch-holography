"""Run the paired readout comparison or verify its archived probabilities."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate as echo
from model.evaluation.observer_encounter import RECEIPT as TWO_READ
from model.executable_model import carrier

HERE = Path(__file__).resolve().parent
PROTOCOL_HASH = 'ac335f9b98c634f228daa5cdaf4bd8d190670d2b5394ebdba2522d0256cdce1e'
DEGREE = 44


def bound_project(hamiltonian, vectors):
    """Uniform Legendre approximation of the negative-energy projection with midgap −5."""
    x = (hamiltonian + 5 * eye(hamiltonian.shape[0], format="csr")) / 17
    a = 5 / 17
    q = (1 - a) / (1 + a)

    def argument(v):
        return ((1 + a*a) * v - 2 * (x @ (x @ v))) / (1 - a*a)

    prev = vectors
    current = argument(vectors)
    expansion = prev + q * current
    for k in range(2, DEGREE + 1):
        following = ((2*k - 1) * argument(current) - (k - 1) * prev) / k
        expansion += q**k * following
        prev, current = current, following
    projected = vectors / 2 - (x @ expansion) / (1 + a)
    return projected, q**(DEGREE + 1) / (2*a)


def bound_report(h, source, nominated, first, target):
    indicators = np.column_stack((source, nominated))
    projected, error = bound_project(h, indicators)
    before = float(np.mean(np.sum(source * projected[:, :-1], axis=0)))
    band_seam = projected[:, -1]
    click_band = float(band_seam[target])
    amplitudes = first[target, :]
    click_prob = float(np.mean(np.abs(amplitudes)**2))
    interference = float(np.real(np.mean(amplitudes.conj() * (band_seam @ first))))
    quiet_joint = before - 2 * interference + click_prob * click_band
    after = quiet_joint + click_prob * click_band
    change = after - before
    # The unconditional change compares two normalized density matrices.
    return {
        "cutoff_shifted_energy": -5, "polynomial_degree_in_squared_energy": DEGREE,
        "analytic_operator_error_bound": float(error),
        "source_and_bypass": before, "conditional_sharp_positive": click_band,
        "conditional_sharp_quiet": quiet_joint / (1 - click_prob),
        "unread_sharp": after, "half_serviced": (before + after) / 2,
        "sharp_minus_bypass": change,
        "small_separated_weight_and_large_selected_weight":
            bool(before + error < .05 and click_band - error > .95),
        "absolute_ensemble_change_exceeds_one_per_million": bool(abs(change) - 2*error > 1e-6),
    }


def local_distribution(local_amplitudes, branch_weight=1.):
    """Unconditional joint weights over the cell's internal seams and its complementary outcome."""
    yes = np.mean(np.abs(local_amplitudes)**2, axis=1)
    return np.append(yes, branch_weight - np.sum(yes))


def evaluate_arm(sector, kappa, source, target, internal, baseline):
    h = sector.hamiltonian(kappa)
    nominated = np.zeros(sector.P)
    nominated[target] = 1
    first = expm_multiply(-1j * h, source)
    direct = expm_multiply(-4j * h, first)
    propagated_seam = expm_multiply(-4j * h, nominated)
    if not np.allclose(np.sum(np.abs(direct)**2, axis=0), 1., rtol=0., atol=1e-11):
        raise ArithmeticError("Incomplete free propagation")
    target_amps = first[target]
    p = float(np.mean(np.abs(target_amps)**2))
    yes = propagated_seam[internal, None] * target_amps
    bypass = direct[internal]
    quiet = bypass - yes

    # A direct evolution after zeroing the marked amplitude checks the branch subtraction.
    state = first[:, 0].copy()
    state[target] = 0
    explicit_quiet = expm_multiply(-4j * h, state)
    quiet_error = float(np.max(np.abs(explicit_quiet[internal] - quiet[:, 0])))

    sharp_yes = local_distribution(yes, p)
    sharp_quiet = local_distribution(quiet, 1-p)
    free = local_distribution(bypass)
    sharp = sharp_yes + sharp_quiet
    half_service = (free + sharp) / 2
    sign = local_distribution(quiet - yes)
    phase_half = (free + sign) / 2
    phase_quarter = (3*free + sign) / 4

    max_phase_error = float(max(np.max(np.abs(phase_half - sharp)),
                                np.max(np.abs(phase_quarter - half_service))))
    prior_error = float(max(abs(p - baseline["first_positive_per_trial"]),
                            abs(sum(sharp_yes[:-1]) / p - baseline["repeat_positive_given_first"])))
    if max(max_phase_error, quiet_error, prior_error) > 1e-10:
        raise ArithmeticError("Phase, explicit quiet branch or original two-read check disagrees")
    for dist in (free, sharp, half_service, phase_half, phase_quarter, sharp_yes/p, sharp_quiet/(1-p)):
        if np.min(dist) < -1e-11 or abs(np.sum(dist) - 1.) > 1e-10:
            raise ArithmeticError("Invalid recorded probability")

    delta = float(np.sum(sharp[:-1] - free[:-1]))
    result = {
        "kappa": kappa,
        "late_pair_per_enrolled_start": {
            "bypass": float(sum(free[:-1])), "half_serviced": float(sum(half_service[:-1])),
            "unread_sharp": float(sum(sharp[:-1])), "sharp_minus_bypass": delta,
            "absolute_difference_exceeds_one_per_million": bool(abs(delta) > 1e-6),
        },
        "sharp_reader": {"first_positive": p, "first_quiet": 1-p,
                         "late_pair_given_positive": float(sum(sharp_yes[:-1]) / p),
                         "late_pair_given_quiet": float(sum(sharp_quiet[:-1]) / (1-p))},
        "complete_patch_distributions": {
            "bypass": free.tolist(), "sharp_yes_joint": sharp_yes.tolist(),
            "sharp_quiet_joint": sharp_quiet.tolist(), "unread_sharp": sharp.tolist(),
            "half_serviced": half_service.tolist(), "phase_half": phase_half.tolist(),
            "phase_quarter": phase_quarter.tolist(),
        },
        "numeric_checks_below_1e_minus_10": {
            "random_sign_matches_random_read": True, "direct_quiet_branch": True,
            "reproduces_deposited_two_read": True,
        }
    }
    if kappa == 8:
        result["global_dressed_band_diagnostic"] = bound_report(h, source, nominated, first, target)
    return result


def run():
    if hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest() != PROTOCOL_HASH:
        raise ValueError("The experiment protocol was changed")
    if echo.native_law is not contrast.law or echo.KAPPA != contrast.KAPPA or echo.KAPPA != 8:
        raise ValueError("The two candidate free sectors or couplings disagree with the protocol")
    car = carrier.federation(0)
    sector = PairSector(car)
    seams = sorted(e for e in car.seams if max(e) < 12)
    internal = np.array([sector.index(*e) for e in seams])
    starts = [sector.index(*p) for p in combinations(range(12), 2) if not sector.adjacent[sector.index(*p)]]
    if len(starts) != 36 or len(seams) != 30:
        raise ArithmeticError("Enrolment differs from the frozen patch")
    source = np.zeros((sector.P, len(starts)))
    source[starts, np.arange(len(starts))] = 1
    prior = json.loads(TWO_READ.read_text())["laws"]
    arms = {name: evaluate_arm(sector, kappa, source, internal[0], internal, prior[name])
            for name, kappa in (("contrast", echo.KAPPA), ("exchange_only", 0))}
    return {
        "protocol_sha256": PROTOCOL_HASH, "carrier_seam_sha256": carrier.seam_digest(car),
        "cells": car.cells, "states": sector.P, "patch": [0], "separated_starts": len(starts),
        "nominated_seam": list(map(int, seams[0])), "durations": [1, 4],
        "outcome_order": [list(map(int, e)) for e in seams] + ["quiet"],
        "arms": arms, "bypass_yield_larger_with_contrast": bool(
            arms["contrast"]["late_pair_per_enrolled_start"]["bypass"] >
            arms["exchange_only"]["late_pair_per_enrolled_start"]["bypass"]),
    }


def verify(value, expected):
    if isinstance(value, dict):
        if value.keys() != expected.keys():
            raise ValueError("Receipt keys differ")
        for key, item in value.items():
            verify(item, expected[key])
    elif isinstance(value, list):
        if len(value) != len(expected):
            raise ValueError("Array lengths differ")
        for item, other in zip(value, expected, strict=True):
            verify(item, other)
    elif isinstance(value, float):
        if not np.isclose(value, expected, rtol=0, atol=1e-11):
            raise ValueError(f"Measured {value}, expected {expected}")
    elif type(value) is not type(expected) or value != expected:
        raise ValueError(f"Measured {value}, expected {expected}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="save measurements to a fresh file")
    args = parser.parse_args()
    result = run()
    if args.record:
        with args.record.open("x") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        verify(result, json.loads((HERE / "result.json").read_text()))
    for name, data in result["arms"].items():
        print(name, data["late_pair_per_enrolled_start"])
        if "global_dressed_band_diagnostic" in data:
            print("bound band", data["global_dressed_band_diagnostic"])
