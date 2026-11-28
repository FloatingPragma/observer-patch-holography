"""Resolve signed material and the heralded guide phase at distinct physical ports."""

import argparse
from hashlib import sha256
from itertools import product
import json
from math import inf, log, nextafter, sqrt
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ....field_charge_records.asymmetric_clock import (
    STEP,
    TICKS,
    UNDERFLOW_PAD,
    interval,
    source_columns,
    step_error,
    taylor_step,
)
from ....field_charge_records.paths import enclosing_records, joint_support, local_terms
from ...fixture import FanPath, archived_fixture
from ..measure import PAULI, pauli_record
from ..word_reader import word_pairs

HERE = Path(__file__).resolve().parent
ORIGIN, TARGET = 1045, 1010
PAIR = (334, 365)
WAITS = (0, 1, 2, 4, 8, 16, 32, 64)
TRIALS = 2**32
SETTINGS = tuple(product(range(1, 4), repeat=2))
FAMILY = 2 * len(TICKS) * (len(WAITS) - 1)
SIGNS = np.array(list(product((1, -1), repeat=2)))
COEF = np.array([[8, 2, 2, 2], [2, -2, 4, 6], [-2, -4, -6, 2], [2, 6, -2, 4]])
OBSERVABLES = np.array([[np.kron(a, b) for b in PAULI] for a in PAULI])
AXES = "IXYZ"
READ_ERROR = 1e-11


def evolution(h, state, tick_list, norm_bound, input_error=0.0, scipy_ticks=None):
    """Use the archived Taylor pad, plus the incoming norm enclosure."""
    norm, width, truncation, rounding = step_error(h)
    assert norm > 0 and width > 0
    local = truncation + rounding + UNDERFLOW_PAD
    assert max(tick_list) * local < 1 / 2
    saved = {}
    state_at_tick = state
    for tick in range(max(tick_list) + 1):
        if tick in tick_list:
            extra = float(4 * tick * local) * norm_bound
            err = nextafter(input_error + extra, inf)
            if scipy_ticks is None or tick in scipy_ticks:
                independent = expm_multiply((-1j * float(tick * STEP)) * h, state)
                assert np.linalg.norm(independent - state_at_tick) < err + 1e-12
            saved[tick] = state_at_tick, err
        if tick < max(tick_list):
            state_at_tick = taylor_step(h, state_at_tick)
    return saved


def probability(amplitudes):
    return float(np.vdot(amplitudes, amplitudes).real / 2)


def pair_amplitudes(wave, block, cell, pair):
    local = np.moveaxis(wave[block.rows].reshape(*block.shape, 2), cell, -2)
    words = block.words[cell]
    return local.reshape(-1, len(words), 2)[:, [words.index(w) for w in pair], :]


def witness_numerator(moments):
    return float(np.sum(COEF * moments) / 32)


def joint_read(wave, error, stock, path, block, pair):
    amplitudes = pair_amplitudes(wave, block, path.size - 1, pair)
    p_port = probability(wave[block.rows])
    pad = (2 * np.sqrt(stock) + error) * error + READ_ERROR
    record = pauli_record(amplitudes, pad)
    p_pair = record["reference_joint_success"]
    refusals = [1 - stock, stock - p_port, p_port - p_pair]
    assert min(refusals) > -1e-12
    np.testing.assert_allclose(sum(refusals) + p_pair, 1, atol=2e-15)
    moments = np.array(record["receiver_reference_pauli_moments"])
    score = None
    score_range = None
    if p_pair > pad:
        score = witness_numerator(moments) / p_pair
        pad_ratio = (np.abs(COEF).sum() / 32 + abs(score)) * pad / (p_pair - pad)
        score_range = interval(score, pad_ratio)
    return {
        "refusal_probabilities_per_start": np.maximum(refusals, 0).tolist(),
        "target_port_probability_per_start": p_port,
        "per_joint_probability_or_pauli_moment_error": float(pad),
        "conditional_fixed_witness": score,
        "conditional_fixed_witness_enclosure": score_range,
        **record,
    }


def density(amplitudes):
    flat = amplitudes.reshape(amplitudes.shape[0], -1)
    return flat.T @ flat.conj() / 2


def local_mark_read(plain, marked, wave_error, stock, path, block):
    cell = path.size - 2
    words = block.words[cell]

    def local_matrix(wave):
        local = np.moveaxis(wave[block.rows].reshape(*block.shape, 2), cell, -2)
        return density(local.reshape(-1, len(words), 2))

    rho, mark = local_matrix(plain), local_matrix(marked)
    p, p_mark = rho.trace().real, mark.trace().real
    pad = (2 * sqrt(stock) + wave_error) * wave_error + READ_ERROR
    if min(p, p_mark) <= pad:
        assert max(p, p_mark) <= pad
        return None
    # Wave error gives a trace-norm bound for each reduced density.
    # Normalize with each preparation's measured sharp-port probability.
    conditional = rho / p
    conditional_mark = mark / p_mark
    norm_bound = 2 * pad / (p - pad) + 2 * pad / (p_mark - pad)
    diff = conditional_mark - conditional
    change = float(np.linalg.norm(diff, ord="fro"))
    sharp = np.trace(diff.reshape(len(words), 2, len(words), 2), axis1=1, axis2=3)
    laws = []
    for face, i, j, _ in word_pairs(path, words):
        selected = [2 * i, 2 * i + 1, 2 * j, 2 * j + 1]
        pair_block = diff[np.ix_(selected, selected)]
        moments = np.einsum("ijab,ba->ij", OBSERVABLES, pair_block).real
        laws.append(
            {
                "face_and_words": [face, words[i], words[j]],
                "change_in_port_conditioned_pair_pauli_moments": moments.tolist(),
            }
        )
    return {
        "port_probability_marked": float(p_mark),
        "port_probability_change_enclosure": interval(p_mark - p, 2 * pad),
        "conditional_reduced_density_difference_frobenius_enclosure": interval(
            change, norm_bound
        ),
        "any_port_conditioned_binary_local_read_response_upper": float(
            min(1, sqrt(len(rho)) * (change + norm_bound) / 2)
        ),
        "target_sharp_word_ids": list(words),
        "target_sharp_word_probability_changes_given_port": sharp.diagonal().real.tolist(),
        "target_face_reference_probability_moment_changes": laws,
        "any_recorded_conditional_moment_change_error": float(norm_bound),
    }


def capture(path, wave, pair):
    take = np.array(
        [port == ORIGIN and words[-1] in pair for port, words in path.basis]
    )
    flip = np.array(
        [port == ORIGIN and words[-1] == pair[1] for port, words in path.basis]
    )
    projected = np.where(take[:, None], wave, 0)
    marked = np.where(flip[:, None], -projected, projected)
    return take, projected, marked


def bridge_calibration(path, take, source, pair):
    bridge = path.bridge_term
    components = []
    assert path.bridges[-1] == (TARGET, ORIGIN)
    indices, addresses, heads, edges, _ = joint_support(
        path, source, path.sign, local_terms(path)
    )
    np.testing.assert_array_equal(indices, np.arange(len(path.basis)))
    enclosing_records(path, path.sign, addresses, heads, edges)
    for a in np.flatnonzero(take):
        port, words = path.basis[a]
        assert port == ORIGIN and words[-1] in pair
        b = path.ids[TARGET, words]
        for u, v in ((a, b), (b, a)):
            row = bridge.getrow(u)
            assert dict(zip(row.indices, row.data, strict=True)) == {u: 4, v: -4}
        flipped = np.flatnonzero(heads[a] != heads[b])
        np.testing.assert_array_equal(flipped, [len(edges) - 1])
        assert heads[a, -1] == (ORIGIN if path.sign == 1 else TARGET)
        assert heads[b, -1] == (TARGET if path.sign == 1 else ORIGIN)
        components.append((int(a), int(b)))
    return components


def counter(rng, record, coefficient=COEF, denominator=32, family=FAMILY):
    t = np.array(record["receiver_reference_pauli_moments"])
    p = t[0, 0]
    bound = record["per_joint_probability_or_pauli_moment_error"]
    assert p > bound
    raw_score = float(np.sum(coefficient * t) / denominator / p)
    systematic = (
        (np.abs(coefficient).sum() / denominator + abs(raw_score)) * bound / (p - bound)
    )
    x, y = SIGNS.T
    groups = []
    estimator = float(coefficient[0, 0] / denominator)
    ranges = 0.0
    for i, j in SETTINGS:
        probs = (p + x * t[i, 0] + y * t[0, j] + x * y * t[i, j]) / 4
        assert np.all(probs >= 0)
        raw = [*record["refusal_probabilities_per_start"], *probs]
        np.testing.assert_allclose(sum(raw), 1, atol=2e-15, rtol=0)
        counts = rng.multinomial(TRIALS, raw)
        accepted = int(counts[-4:].sum())
        assert accepted > 0 and int(counts.sum()) == TRIALS
        scores = coefficient[i, j] * x * y
        if i == j:
            scores = scores + coefficient[i, 0] * x + coefficient[0, j] * y
        scores = scores / denominator
        estimator += float(np.dot(counts[-4:], scores) / accepted)
        ranges += np.ptp(scores) ** 2 / accepted
        groups.append({"settings": AXES[i] + AXES[j], "counts": counts.tolist()})
    radius = float(sqrt(log(2 * family / 0.01) * ranges / 2) + systematic)
    return {
        "banks": groups,
        "fixed_witness_estimate": estimator,
        "fixed_witness_simultaneous_99pct_interval": interval(estimator, radius),
    }


def case(path, rng):
    mask = sum(path.field.bits.values())
    pair = PAIR if path.sign == 1 else tuple(mask ^ w for w in PAIR)
    source = source_columns(path)
    block = path.blocks[path.size - 2, path.exit[-2]]
    arrival = path.blocks[path.size - 1, 0]
    assert path.cells[block.cell][block.vertex] == TARGET
    assert path.cells[arrival.cell][arrival.vertex] == ORIGIN
    assert any(
        face == 0 and set(pair) == {words[i], words[j]}
        for words in (block.words[-1], arrival.words[-1])
        for face, i, j, _ in word_pairs(path, words)
    )
    assert block.words[-1] == arrival.words[-1]
    flights = []
    for tick, (wave, error) in evolution(path.h, source, TICKS, sqrt(2)).items():
        take, projected, marked = capture(path, wave, pair)
        stock = probability(projected)
        arrival_read = joint_read(projected, error, stock, path, arrival, pair)
        controls = []
        components = bridge_calibration(path, take, source, pair)
        indices = np.array(components)
        for arm, h in (
            ("native", path.h),
            ("half_bridge_rates", path.parked + path.bridge_term / 2),
            ("parked", path.parked),
            ("isolated_bridge_terms", path.bridge_term),
        ):
            results = []
            together = np.column_stack((projected, marked))
            for lag, (flown, total_error) in evolution(
                h, together, WAITS, 2 * (sqrt(stock) + error), sqrt(2) * error
            ).items():
                before, after = flown[:, :2], flown[:, 2:]
                read = joint_read(before, total_error, stock, path, block, pair)
                changed = joint_read(after, total_error, stock, path, block, pair)
                local = local_mark_read(before, after, total_error, stock, path, block)
                if arm == "parked":
                    assert read["target_port_probability_per_start"] == 0
                    assert local is None
                if arm == "isolated_bridge_terms":
                    angle = float(lag * STEP)
                    bridge_amplitude = (1 - np.exp(-8j * angle)) / 2
                    np.testing.assert_allclose(
                        before[indices[:, 1]],
                        bridge_amplitude * projected[indices[:, 0]],
                        atol=total_error + 1e-14,
                    )
                    np.testing.assert_allclose(
                        read["reference_joint_success"],
                        stock * np.sin(4 * angle) ** 2,
                        atol=1e-12,
                    )
                    if local:
                        lo, hi = local[
                            "conditional_reduced_density_difference_frobenius_enclosure"
                        ]
                        assert lo <= 0 <= hi
                results.append(
                    {
                        "wait": float(lag * STEP),
                        "old_face_read": read,
                        "old_face_read_with_mark": changed,
                        "local_target_phase_mark_comparison": local,
                        "simulated_counts": counter(rng, read)
                        if arm == "native" and lag
                        else None,
                    }
                )
            controls.append({"arm": arm, "waits": results})
        flights.append(
            {
                "first_appointment": float(tick * STEP),
                "first_gate_probability_per_start": stock,
                "first_charge_face_read": arrival_read,
                "straight_to_target_without_first_read": joint_read(
                    wave, error, 1.0, path, block, pair
                ),
                "isolated_bridge_pairs": len(components),
                "controls": controls,
            }
        )
    return {
        "charge_sign": path.sign,
        "ordered_face_words": list(pair),
        "flights": flights,
    }


def run():
    witness = np.array([2, 1, 1, -1 - 1j])
    coefficient = np.array(
        [
            [np.vdot(witness, np.kron(a, b.T) @ witness).real for b in PAULI]
            for a in PAULI
        ]
    )
    np.testing.assert_array_equal(coefficient, COEF)
    with archived_fixture() as archive:
        positive, negative = (FanPath(archive, 4, s, asymmetric=True) for s in (1, -1))
        mask = sum(positive.field.bits.values())
        perm = [
            negative.ids[p, tuple(mask ^ w for w in words)]
            for p, words in positive.basis
        ]
        for attribute in ("h", "parked", "bridge_term"):
            a, b = getattr(positive, attribute), getattr(negative, attribute)
            assert (a != b[perm][:, perm]).nnz == 0
        np.testing.assert_array_equal(negative.source()[perm], -positive.source())
        rng = np.random.default_rng(1045101098)
        cases = [case(path, rng) for path in (positive, negative)]
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "start_site_and_reference_tag": 98,
        "intermediate_material_and_face_site": ORIGIN,
        "terminal_material_site": TARGET,
        "simulated_trials_per_bank": TRIALS,
        "simultaneous_cohort_family_size": FAMILY,
        "sharp_record_bins": [
            "first_gate_refused",
            "terminal_not_at_1010",
            "at_1010_outside_old_face_pair",
            "++",
            "+-",
            "-+",
            "--",
        ],
        "exact_signed_conjugacy": True,
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as out:
            out.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(HERE.parents[2] / "replay.py"))["compare"]
        compare(
            result, json.loads((HERE / "records.json").read_text()), "separated_read"
        )
        print(
            "Replayed signed spatial separation records and all 28 witness intervals."
        )


if __name__ == "__main__":
    main()
