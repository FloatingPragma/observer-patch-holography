"""Certify the eight-site bright branch and spin-return at the opposite site."""

import argparse
from collections import Counter
from dataclasses import dataclass
from functools import cache
from hashlib import sha256
import json
from math import inf, log, nextafter, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from .. import ring
from . import car

HERE = Path(__file__).resolve().parent
GROUND_CUT = -21 / 256
HOLE_CUT = -174 / 1024
STAGES = (1, -1)
TRIALS = 2**30
COMPARISONS = 16


def around(value, radius):
    return [
        nextafter(float(value) - radius, -inf),
        nextafter(float(value) + radius, inf),
    ]


def reference_check(up, down):
    """Bind the independent spin-separated chart to the source's archived order."""
    phases = []
    old_rows = []
    archive_words, archive_where = ring.sector(car.L, up, down)
    for word in car.stock(up, down)[0]:
        a, b = word & car.MASK, word >> car.L
        interleaved = 0
        parity = 0
        for site in range(car.L):
            interleaved |= ((a >> site) & 1) << (2 * site)
            if b & (1 << site):
                interleaved |= 1 << (2 * site + 1)
                parity += (a >> (site + 1)).bit_count()
        phases.append((-1) ** parity)
        old_rows.append(archive_where[interleaved])
    assert len(set(old_rows)) == len(archive_words)
    old_h = ring.hamiltonian(car.L, up, down)[old_rows][:, old_rows]
    signs = np.asarray(phases)
    assert (
        old_h != car.hamiltonian(up, down).multiply(signs).multiply(signs[:, None])
    ).nnz == 0


@dataclass
class Spectrum:
    levels: np.ndarray
    error: float
    low: list
    size: int

    def count(self, threshold):
        assert np.min(abs(self.levels - threshold)) > self.error
        return int(np.count_nonzero(self.levels < threshold))


@cache
def spectrum(up, down, k):
    h, inject, rounding = car.block(up, down, k)
    n = len(h)
    if n == 0:
        return None
    assert n <= 640
    values, vectors = eigh(h, driver="evd")
    assert np.max(abs(values)) < 8
    gamma = (16 * n + 64) * car.EPS
    gamma /= 1 - gamma
    assert np.linalg.norm(vectors) < 1.00001 * np.sqrt(n)
    gram = np.linalg.norm(vectors.conj().T @ vectors - np.eye(n))
    residual = np.linalg.norm((vectors * values) @ vectors.conj().T - h)
    hnorm = np.linalg.norm(h)
    scale = float(np.max(abs(values)))
    # Bound every complex BLAS dot by a real absolute product, then use the
    # Frobenius product bound. This also covers storage and norm roundoff.
    defect = (gram + gamma * (1.001 * n + np.sqrt(n))) * (1 + gamma)
    factor = (residual + gamma * (1.001 * n * scale + hnorm)) * (1 + gamma)
    assert defect < 1 / 2
    # For the polar unitary Q, ||V-Q||_F <= ||V*V-I||_F.
    # Weyl therefore counts the whole exact Hermitian block.
    error = float(factor + (2 + defect) * scale * defect + rounding + 1e-13)
    assert error < 2e-8
    low = []
    for i in range(min(4, n)):
        v = vectors[:, i]
        v = v / np.linalg.norm(v)
        residual_col = float(np.linalg.norm(h @ v - values[i] * v))
        radius = residual_col + gamma * (hnorm + abs(values[i])) + rounding + 1e-13
        neighbors = abs(np.delete(values, i) - values[i])
        gap = float(min(neighbors) - error) if len(neighbors) else 1.0
        eigenvector_error = 2 * radius / gap + 2e-12 if gap > radius else None
        low.append(
            {
                "energy": float(values[i]),
                "radius": float(radius),
                "error": eigenvector_error,
                "vector": inject @ v,
            }
        )
    return Spectrum(values, error, low, n)


def sectors(up, down):
    while up <= car.L and down >= 0:
        car.check_symmetries(up, down)
        reference_check(up, down)
        momenta = [spectrum(up, down, k) for k in range(car.L)]
        assert sum(s.size for s in momenta if s is not None) == len(
            car.stock(up, down)[0]
        )
        for k in range(1, car.L):
            a, b = momenta[k], momenta[-k]
            assert (a is None) == (b is None)
            if a is not None:
                np.testing.assert_allclose(a.levels, b.levels, atol=2e-12, rtol=0)
        # Integer translation and real H imply exact k/-k correspondence,
        # independent of this numerical consistency check.
        yield (up - down, momenta)
        up, down = up + 1, down - 1


def count_row(momenta, threshold):
    return [s.count(threshold) if s is not None else 0 for s in momenta]


def ground():
    central, spin1 = list(sectors(4, 4))[:2]
    (m0, b0), (m1, b1) = central, spin1
    assert m0 == 0 and m1 == 2
    assert count_row(b0, GROUND_CUT) == [0, 0, 0, 0, 1, 0, 0, 0]
    assert sum(count_row(b1, GROUND_CUT)) == 0
    g = b0[4].low[0]
    assert g["error"] is not None and g["error"] < 1e-6
    next_brackets = [
        around(b.levels[i], b.error)
        for b, i in [(w, int(k == 4)) for k, w in enumerate(b0) if w is not None]
    ]
    first_gap = [
        min(x[0] for x in next_brackets) - (g["energy"] + g["radius"]),
        min(x[1] for x in next_brackets) - (g["energy"] - g["radius"]),
    ]
    words, where = car.stock(4, 4)
    neel = sum(1 << (site + (site % 2) * car.L) for site in range(car.L))
    vector = g["vector"]
    weight = abs(vector[where[neel]]) ** 2
    g["neel_probability"] = float(weight)
    g["neel_probability_error"] = 2 * g["error"] + 1e-12
    g["momentum"] = 4
    assert np.linalg.norm(car.raising(4, 4) @ vector) < 1e-10
    assert np.linalg.norm(car.translation(4, 4) @ vector + vector) < 1e-10
    return g, {
        "original_energy": around(g["energy"], g["radius"]),
        "ground_gap": first_gap,
        "preparation_weight": around(weight, g["neel_probability_error"]),
        "half_filled_full_dimension": len(words),
        "sector_ground_counts": {
            "twice_spin_z_0": count_row(b0, GROUND_CUT),
            "twice_spin_z_2": count_row(b1, GROUND_CUT),
        },
    }


def cluster(low, hi):
    selected = []
    for k, block in enumerate(low):
        if block is not None:
            count = block.count(hi)
            assert count <= 2
            for i in range(count):
                selected.append((block.levels[i], k, i))
    selected.sort()
    result = []
    for e, k, i in selected:
        if not result or e - result[-1][0][0] > 1e-6:
            result.append([])
        result[-1].append((e, k, i))
    return result


def scoped_levels(spins, cutoff):
    central = spins[0][1]
    groups = cluster(central, cutoff)
    assert [len(g) for g in groups] == [2, 2, 1, 2]
    thresholds = (
        [-1.0]
        + [(p[-1][0] + q[0][0]) / 2 for p, q in zip(groups, groups[1:])]
        + [cutoff]
    )
    assert sum(count_row(central, thresholds[0])) == 0
    assembled = []
    for i, group in enumerate(groups):
        qset = sorted(k for _, k, _ in group)
        if len(group) == 2:
            assert 0 < qset[0] < 4 < qset[1] < 8 and sum(qset) == 8
        sector_counts = []
        for sz, bank in spins:
            a, b = count_row(bank, thresholds[i]), count_row(bank, thresholds[i + 1])
            diffs = [y - x for x, y in zip(a, b, strict=True)]
            sector_counts.append((sz, diffs))
        desired_spin = 1 if i in (0, 3) else 3
        counts = np.array([x[1] for x in sector_counts])
        mult = counts - np.vstack((counts[1:], np.zeros_like(counts[0])))
        expected = Counter({k: 1 for k in qset})
        assert all(
            int(mult[j, k]) == (expected[k] if sz == desired_spin else 0)
            for j, (sz, _) in enumerate(sector_counts)
            for k in range(car.L)
        )
        # Exact equality within a conjugate pair follows from real H. Each
        # momentum in the group is known to contain one line.
        vecs = [central[k].low[j] for _, k, j in group]
        brackets = [around(x["energy"], x["radius"]) for x in vecs]
        enclosure = [max(a for a, b in brackets), min(b for a, b in brackets)]
        assert enclosure[0] < enclosure[1]
        assembled.append(
            {
                "momenta": qset,
                "energy": enclosure,
                "total_spin": str(desired_spin) + "/2",
                "spin_z_sector_counts": {str(sz): row for sz, row in sector_counts},
                "translations_by_four": (-1) ** qset[0],
                "exact_spin_selection_excludes_source": desired_spin != 1,
                "vectors": vecs,
            }
        )
    next_brackets = [
        around(b.levels[b.count(cutoff)], b.error)
        for b in central
        if b is not None and b.count(cutoff) < b.size
    ]
    next_ene = [min(a for a, b in next_brackets), min(b for a, b in next_brackets)]
    return assembled, next_ene


def plus_pump(sign, site):
    return car.field(4, 4, site, int(sign < 0), sign > 0)


def minus_pump(sign, site):
    return sign * car.field(4, 4, site, int(sign > 0), sign > 0)


def local_spin_channel(alpha, denominator, err):
    """Infer the exact rotation-covariant local map with no ground-return flag."""
    assert denominator > err
    lo, hi = (
        min(
            (alpha + a * err) / (denominator + b * err)
            for a in (-1, 1)
            for b in (-1, 1)
        ),
        max(
            (alpha + a * err) / (denominator + b * err)
            for a in (-1, 1)
            for b in (-1, 1)
        ),
    )
    lo, hi = float(lo), float(hi)
    return {
        "single_qubit_depolarizing_eta": [lo, hi],
        "Bell_overlap_interval": [(1 + 3 * lo) / 4, (1 + 3 * hi) / 4],
        "exactly_entanglement_breaking_at_this_resolution": lo >= -1 / 3
        and hi <= 1 / 3,
        "has_nonseparable_local_reference_record": bool(lo > 1 / 3),
        "unfiltered_optimal_CHSH_interval": [
            2 * np.sqrt(2) * (0 if lo * hi < 0 else min(abs(lo), abs(hi))),
            2 * np.sqrt(2) * max(abs(lo), abs(hi)),
        ],
    }


def local_spin_records(
    sign, charge, spin, single, weights, error, absolute, prep, prep_error
):
    bare = [
        local_spin_channel(float(u), float(v), error)
        for u, v in zip(spin @ weights, single @ weights, strict=True)
    ]
    gate_4 = np.asarray(charge[4] == sign)
    tests = []
    for site in (3, 5):
        joint = gate_4 * single[site]
        moment = float((joint * spin[site]) @ weights)
        weight = float(joint @ weights)
        admitted = float(joint @ absolute)
        bound = error * np.sum(absolute)
        p = prep * admitted / 2
        p_error = (prep * bound + prep_error * (admitted + bound)) / 2
        tests.append(
            {
                "charge_site": 4,
                "spin_site": site,
                "joint_native_excess_on_the_two_receiving_sites": sign,
                "conditioned_joint_gate_weight": around(weight, error),
                "joint_gate_probability_per_original_Neel_start": around(p, p_error),
                "joint_gate_refusal_per_original_Neel_start": around(1 - p, p_error),
                "reference_channel_without_global_terminal_drain": local_spin_channel(
                    moment, weight, error
                ),
            }
        )
    return {
        "local_spin_channel_at_each_single_occupation_port": bare,
        "simultaneous_signed_charge_and_neighbor_spin": tests,
    }


def sampled_pauli_bank(report, pg, norm, rng):
    interval = report["reference_channel_without_global_terminal_drain"][
        "single_qubit_depolarizing_eta"
    ]
    eta = float(np.mean(interval))
    success = float(np.mean(report["joint_gate_probability_per_original_Neel_start"]))
    stages = [1 - pg, 3 * pg / 4, pg * (1 / 4 - norm / 2), pg * norm / 2 - success]
    assert min(stages) >= 0 and success > 0
    signals = []
    ranges = 0.0
    records = []
    for axis, phase in zip("XYZ", (1, -1, 1), strict=True):
        same = (1 + phase * eta) / 4 * success
        different = (1 - phase * eta) / 4 * success
        counts = rng.multinomial(TRIALS, [*stages, same, different, different, same])
        accepted = int(counts[4:].sum())
        assert accepted > 0 and int(counts.sum()) == TRIALS
        signals.append(
            float(phase * (counts[4] + counts[7] - counts[5] - counts[6]) / accepted)
        )
        ranges += 4 / (9 * accepted)
        records.append({"axis_pair": axis + axis, "counts": counts.tolist()})
    eta_est = float(np.mean(signals))
    systematic = (interval[1] - interval[0]) / 2
    radius = sqrt(log(2 * COMPARISONS / 0.01) * ranges / 2) + systematic
    enclosure = around(eta_est, radius)
    return {
        "simulated_eta_estimate": eta_est,
        "simultaneous_99_percent_eta": enclosure,
        "simultaneous_99_percent_Bell_overlap": [(1 + 3 * t) / 4 for t in enclosure],
        "banks": records,
    }


def relative_rows(sign, g, lines, up, down, rng):
    vector = g["vector"]
    raised = car.raising(down, up)
    maps_plus = [plus_pump(sign, y) for y in range(car.L)]
    maps_minus = [minus_pump(sign, y) for y in range(car.L)]
    # The spinor intertwining holds as an integer operator, including the
    # conjugate orientation for the hole. On the full singlet the S^+ term vanishes.
    singlet_raise = car.raising(4, 4)
    for y in range(car.L):
        up_sector = car.field(5, 3, y, int(sign > 0), sign > 0)
        assert (
            raised @ maps_minus[y] != maps_plus[y] + sign * up_sector @ singlet_raise
        ).nnz == 0
    sink = [
        np.array([map_ @ vector for map_ in mats]) for mats in (maps_plus, maps_minus)
    ]
    assert np.linalg.norm(raised @ sink[1][0] - sink[0][0]) < 1e-10
    full_states = []
    errors = []
    weights = []
    for line in (lines[0], lines[-1]):
        pulse = np.zeros_like(sink[0][0])
        line_error = 0.0
        line_weights = []
        for v in line["vectors"]:
            assert v["error"] is not None
            assert v["error"] < 1e-5
            v0 = v["vector"]
            c = np.vdot(v0, sink[0][0])
            pulse += v0 * c
            # The norm of each rank-one projection differs by at most
            # twice the column error. The pump is a contraction.
            err = 2 * v["error"] + g["error"] + 1e-12
            line_error += err
            line_weights.append(around(abs(c) ** 2, (2 * abs(c) + err) * err))
        assert min(w[0] for w in line_weights) > 0.01
        translation = car.translation(up, down, 4)
        assert (
            np.linalg.norm(translation @ pulse - line["translations_by_four"] * pulse)
            < 1e-10
        )
        full_states.append(pulse)
        weights.append(line_weights)
        errors.append(line_error)
    lower, upper = full_states
    norm = float(np.linalg.norm(lower + upper) ** 2)
    assert abs(norm - (np.linalg.norm(lower) ** 2 + np.linalg.norm(upper) ** 2)) < 1e-12
    field_error = sum(errors)
    norm_error = (2 * np.sqrt(norm) + field_error) * field_error + 1e-12
    tlow = lines[0]["energy"]
    thigh = lines[-1]["energy"]
    gap = [thigh[0] - tlow[1], thigh[1] - tlow[0]]
    assert gap[0] > 0
    translation = car.translation(up, down, 4)
    assert np.linalg.norm(translation @ (lower + upper) - (lower - upper)) < 1e-10
    # Check against native sparse propagation at a fixed physical wait,
    # independent of the spectral sign-flip identity.
    evolved = expm_multiply(-149j * car.hamiltonian(up, down) / 16, lower + upper)
    predicted = (
        np.exp(-149j * np.mean(tlow)) * lower + np.exp(-149j * np.mean(thigh)) * upper
    )
    assert np.linalg.norm(evolved - predicted) < 1e-9
    occ = np.asarray(car.stock(up, down)[0])
    charge = []
    spin = []
    singleton = []
    for y in range(car.L):
        a, b = (occ >> y) & 1, (occ >> (car.L + y)) & 1
        charge.append(a + b - 1)
        spin.append(a - b)
        singleton.append((a - b) ** 2)
    rows = {}
    for label, wave, extra in (
        ("parked", lower + upper, 0.0),
        ("half_beat", lower - upper, 0.0),
        ("fixed_wait_149", predicted, 149 * (gap[1] - gap[0])),
    ):
        err = field_error + extra
        lowering = raised.T
        lower_wave = lowering @ wave
        values = [
            np.array([np.vdot(v, w) for v in basis])
            for basis, w in zip(sink, (wave, lower_wave), strict=True)
        ]
        np.testing.assert_allclose(values[0], values[1], atol=1e-10, rtol=0)
        probabilities = (
            abs(values[0]) ** 2 / 4
        )  # pump and drain each at amplitude 1/sqrt(2)
        probabilities = probabilities * g["neel_probability"]
        bound = (2 + err + g["error"]) * (err + g["error"]) + g[
            "neel_probability_error"
        ]
        amplitude_error = err + g["error"] + 1e-12
        expectation_error = (2 * np.linalg.norm(wave) + err) * err + 1e-12
        normalized_error = 2 * expectation_error / (norm - expectation_error)
        diag = (abs(wave) ** 2) / norm
        np.testing.assert_allclose(
            [sum(charge @ diag), sum(spin @ diag)], [sign, 1], atol=1e-9
        )
        rows[label] = {
            **local_spin_records(
                sign,
                np.asarray(charge),
                np.asarray(spin),
                np.asarray(singleton),
                diag,
                normalized_error,
                abs(wave) ** 2,
                g["neel_probability"],
                g["neel_probability_error"],
            ),
            "receipt_probability_per_Neel_start": [
                around(
                    x,
                    min(
                        bound,
                        g["neel_probability"]
                        * (2 * abs(a) + amplitude_error)
                        * amplitude_error
                        / 4
                        + g["neel_probability_error"]
                        * (abs(a) + amplitude_error) ** 2
                        / 4
                        + 1e-12,
                    ),
                )
                for x, a in zip(probabilities, values[0], strict=True)
            ],
            "spin_twice_z_given_window": [
                around(v, normalized_error) for v in spin @ diag
            ],
            "signed_excess_given_window": [
                around(v, normalized_error) for v in charge @ diag
            ],
            "single_occupation_given_window": [
                around(v, normalized_error) for v in singleton @ diag
            ],
        }
    incoherent = abs(lower) ** 2 + abs(upper) ** 2
    normalized_dephasing_error = 2 * norm_error / (norm - norm_error)
    fine_local = local_spin_records(
        sign,
        np.asarray(charge),
        np.asarray(spin),
        np.asarray(singleton),
        incoherent / norm,
        normalized_dephasing_error,
        incoherent,
        g["neel_probability"],
        g["neel_probability_error"],
    )
    for label, row in [*rows.items(), ("resolved_energy_then_forgotten", fine_local)]:
        for report in row["simultaneous_signed_charge_and_neighbor_spin"]:
            report["simulated_local_Pauli_comparison"] = sampled_pauli_bank(
                report, g["neel_probability"], norm, rng
            )
    dephased = [
        g["neel_probability"]
        * (abs(np.vdot(v, lower)) ** 2 + abs(np.vdot(v, upper)) ** 2)
        / 4
        for v in sink[0]
    ]
    e = field_error + g["error"] + 1e-12
    return {
        "whole_gate_two_Bright_momentum_weights": weights,
        "local_unscaled_weights_total": around(norm, norm_error),
        "relative_energy_gap": gap,
        "half_beat_native_time": [pi / gap[1], pi / gap[0]],
        "two_spin_source_and_receiver_intertwiner_checked": True,
        "fixed_wait_direct_propagation_checked": True,
        "both_bright_translation_parities": [
            lines[0]["translations_by_four"],
            lines[-1]["translations_by_four"],
        ],
        "half_beat_T4_on_reachable_support_exact_by_line_count_and_translation": True,
        "joint_window_acceptance_per_Neel_start": around(
            g["neel_probability"] * norm / 2,
            (g["neel_probability_error"] + norm_error) / 2,
        ),
        "ground_refusal_per_Neel_start": around(
            1 - g["neel_probability"], g["neel_probability_error"]
        ),
        "pump_refusal_per_Neel_start": around(
            3 * g["neel_probability"] / 4, g["neel_probability_error"]
        ),
        "window_refusal_per_Neel_start": around(
            g["neel_probability"] * (1 / 4 - norm / 2),
            (g["neel_probability_error"] + norm_error) / 2,
        ),
        "complete_drain_refusal_by_wait_and_receiver_per_Neel_start": {
            label: [
                around(
                    g["neel_probability"] * norm / 2 - sum(p) / 2,
                    (g["neel_probability_error"] + norm_error) / 2 + (p[1] - p[0]) / 2,
                )
                for p in data["receipt_probability_per_Neel_start"]
            ]
            for label, data in rows.items()
        },
        "free_readouts": rows,
        "resolved_energy_then_forgotten_receiver_local": fine_local,
        "resolved_line_dephased_receipt_by_site_per_Neel_start": [
            around(x, 2 * e) for x in dephased
        ],
    }


def run():
    ring.source()
    for x in range(car.L):
        sigma_x = np.array([[0, 1], [1, 0]])
        g = np.linalg.matrix_power(-1j * sigma_x, x)
        gp = np.linalg.matrix_power(-1j * sigma_x, (x + 1) % car.L)
        np.testing.assert_array_equal(g @ (1j * sigma_x) @ gp.conj().T, -np.eye(2))
    g, ground_record = ground()
    banks = []
    rng = np.random.default_rng(84350149)
    for sign in STAGES:
        up, down = (5, 4) if sign == 1 else (4, 3)
        cutoff = HOLE_CUT + (sign + 1) / 2
        spins = list(sectors(up, down))
        levels, next_ene = scoped_levels(spins, cutoff)
        output = relative_rows(sign, g, levels, up, down, rng)
        for row in levels:
            row.pop("vectors")
            a, b = row["energy"]
            row["reservoir_corrected_energy"] = [
                a - g["energy"] - g["radius"] - sign / 2,
                b - g["energy"] + g["radius"] - sign / 2,
            ]
        banks.append(
            {
                "sign": sign,
                "upper_window_cut": cutoff,
                "next_distinct_energy": next_ene,
                "counts_under_upper_cut": {
                    str(sz): count_row(b, cutoff) for sz, b in spins
                },
                "complete_dark_and_bright_lines": levels,
                **output,
            }
        )
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "archived_source_sha256": ring.ANCESTRY,
        "length": car.L,
        "ground": ground_record,
        "signed_branches": banks,
        "simulated_original_Starts_per_Pauli_bank": TRIALS,
        "simultaneous_local_cohort_count": COMPARISONS,
        "simulated_count_bin_order": [
            "ground_refused",
            "pump_refused",
            "energy_window_refused",
            "two_site_charge_and_single_gate_refused",
            "++",
            "+-",
            "-+",
            "--",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    a = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if a.output:
        with a.output.open("x") as out:
            out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert result == json.loads((HERE / "records.json").read_text())
        print(
            "Replayed independent signed momentum/spin counts and the eight-site half beat."
        )


if __name__ == "__main__":
    main()
