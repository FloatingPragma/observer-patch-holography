"""Source selection and three-mode recovery across the archived matched seam."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[1]
ARCHIVED = "experiments/interval_echo_electron_flux_binding"
TIMES = [np.pi / 8, np.pi / 4, np.pi / 2, np.pi]
NLEVELS = 3


def extract_sources(root):
    manifest = dict(line.split("  ", 1)[::-1] for line in (FIELD / "records.sha256").read_text().splitlines())
    sources = ("search.py", "bridge.py", "free_bridge.py", "free_bridge.json")
    hashes = {}
    with ZipFile(FIELD / "records.zip") as archive:
        for filename in sources:
            name = f"{ARCHIVED}/{filename}"
            content = archive.read(name)
            digest = sha256(content).hexdigest()
            if manifest[name] != digest:
                raise ValueError(f"Archived dependency mismatch for {name}")
            target = root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
            hashes[name] = digest
    (root / "reverse-engineering-reality").symlink_to(FIELD.parents[4], target_is_directory=True)
    sys.path.insert(0, str(root / ARCHIVED))
    return hashes, json.loads((root / ARCHIVED / "free_bridge.json").read_text())


def preparation(patch, states, h, levels, modes):
    from search import phase

    words = patch.field_component()
    assert states[:len(words)] == [(patch.faces[0][0], word) for word in words]
    strings = [None] + list(combinations(range(len(patch.faces)), 2))
    assert len(strings) == 191
    masks = [0 if faces is None else patch.path(*faces) for faces in strings]
    source = np.array([[phase(word, mask) for word in words] for mask in masks]) / np.sqrt(len(words))
    weights = np.square(source @ modes[:len(words)])
    upper = float(abs(h).sum(axis=1).max()) / 2
    step = np.pi / upper
    inputs = []
    for i in range(NLEVELS):
        sel = int(np.argmax(weights[:, i]))
        w = float(weights[sel, i])
        gap = float(np.min(np.abs(np.delete(levels[:NLEVELS + 1], i) - levels[i])))
        assert w > 0 and gap > 1e-8
        # A coherent uniform superposition of N delays has a geometric-series response.
        # Every other eigenvalue differs by at least the observed nearest gap and at most upper.
        den = np.sin(step * gap / 2)
        count = max(1, math.ceil(math.sqrt(99 * (1 - w) / w) / den))
        false_pass_upper = min(1.0, 1.0 / (count * den) ** 2)
        inputs.append({"level": i, "selected_face_pair": strings[sel], "ideal_herald_probability": w,
                       "ideal_no_herald_probability": 1 - w, "selected_source_three_weights": weights[sel].tolist(),
                       "nearest_gap_native": gap,
                       "finite_filter_delay_branches": count, "longest_delay_native": (count - 1) * step,
                       "finite_herald_probability_range": [w, w + (1 - w) * false_pass_upper],
                       "finite_heralded_fidelity_lower": w / (w + (1 - w) * false_pass_upper)})
    return {"method": "charged anchor, uniform accessible heads, optional dual parity string, heralded native energy",
            "seed_port": int(patch.faces[0][0]), "accessible_anchor_words": len(words),
            "source_settings": len(strings), "unmodified_source_three_weights": weights[0].tolist(),
            "finite_delay_step_native": float(step), "charged_spectrum_upper_bound_native": upper, "inputs": inputs}


def diagonal_fidelity(kraus, total):
    # Global optimum over two independent phases. Maximize the second analytically,
    # then enclose the first on a uniform grid with its analytic Lipschitz bound.
    diag = np.diagonal(kraus, axis1=1, axis2=2)
    gram = diag.conj().T @ diag
    mesh = 16384
    z = np.exp(2j * np.pi * np.arange(mesh) / mesh)
    candidates = np.trace(gram).real + 2 * ((gram[0, 1] * z).real + abs(gram[0, 2] + gram[1, 2] * z.conj()))
    best = float(max(candidates)) / (3 * total)
    err = float(2 * (abs(gram[0, 1]) + abs(gram[1, 2])) * np.pi / mesh) / (3 * total)
    return [best, min(1., best + err)]


def coherent_readout(reduced, neutral):
    # (input, output, departed neutral) -> one Kraus matrix per unresolved neutral outcome.
    kraus = np.array(reduced).transpose(2, 1, 0)
    gram = np.einsum("aji,ajk->ik", kraus.conj(), kraus)
    total = float(np.trace(gram).real)
    vectors = kraus.transpose(0, 2, 1).reshape(len(neutral), 9)
    choi = vectors.T @ vectors.conj() / 3
    assert abs(np.trace(choi).real - total / 3) < 2e-12
    partial = choi.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
    negativity = float(-np.minimum(0, np.linalg.eigvalsh(partial)).sum())
    k_ground = np.einsum("a,aji->ji", neutral, kraus)
    singular = np.linalg.svd(k_ground, compute_uv=False)
    # A receiver contraction sigma_min * K_ground^{-1} heralds the identity
    # uniformly for any source state. This bound conditions first on the departed ground click.
    match = np.square(abs(kraus)).sum(axis=0).T
    return {"three_level_acceptance_range_all_pure_inputs": np.linalg.eigvalsh(gram)[[0, -1]].tolist(),
            "equal_prior_three_level_acceptance": total / 3,
            "equal_prior_correct_label_unconditional": float(np.trace(match)) / 3,
            "correct_label_given_three_level_click_equal_prior": float(np.trace(match)) / total,
            "optimal_label_guess_given_three_level_click_equal_prior": float(match.max(axis=0).sum()) / total,
            "choi_partial_transpose_negativity_unconditional": negativity,
            "choi_partial_transpose_negativity_given_three_level_click": negativity / (total / 3),
            "choi_overlap_phase_corrected_identity_given_three_level_click_bounds": diagonal_fidelity(kraus, total),
            "neutral_ground_transfer_singular_values": singular.tolist(),
            "optimal_uniform_exact_inverse_filter_and_neutral_ground_probability": float(singular[-1] ** 2)}, match


def sector(data, archived):
    from free_bridge import free_flight
    from search import lowest

    patch, sign, g, receiving, _, h, states, _, idx, nh, prior_levels, _, neutral = data
    levels, vec, groups, residual = lowest(h)
    np.testing.assert_allclose(levels, prior_levels, atol=2e-11)
    assert all(len(group) == 1 for group in groups[:4])
    modes = vec[:, :NLEVELS].copy()
    modes *= np.sign(modes[np.argmax(abs(modes), axis=0), np.arange(NLEVELS)])
    np.testing.assert_allclose(modes.T @ modes, np.eye(NLEVELS), atol=2e-12)
    source = preparation(patch, states, h, levels, modes)
    overlaps = neutral @ modes[idx]
    derivative = -4 * np.outer(overlaps, overlaps)
    gateway_source_weights = np.square(modes[idx]).sum(axis=0)
    # The seam swaps the two gateway field words exactly at pi/8. Its arriving
    # product wavefunction puts the same neutral wavefunction at the charged gateway for every input.
    control_reduced = [np.outer(overlaps, mode[idx]) for mode in modes.T]
    seam_channel, seam_energy = coherent_readout(control_reduced, neutral)
    assert seam_channel["choi_partial_transpose_negativity_given_three_level_click"] < 2e-12
    np.testing.assert_allclose(seam_energy, np.outer(gateway_source_weights, overlaps ** 2), atol=1e-13)
    reduced = [[] for _ in TIMES]
    arrivals = [[] for _ in TIMES]
    records = [[] for _ in TIMES]
    reference_error = 0.
    max_norm_error = 0.
    for i, mode in enumerate(modes.T):
        outputs, tolerance = free_flight(h, nh, idx, np.outer(mode, neutral), TIMES)
        for ti, (left, right) in enumerate(outputs):
            right_prob = float(np.vdot(right, right).real)
            left_prob = float(np.vdot(left, left).real)
            max_norm_error = max(max_norm_error, abs(left_prob + right_prob - 1))
            projections = modes.T @ right
            energy = np.square(abs(projections)).sum(axis=1)
            joint = np.square(abs(projections @ neutral))
            other = float(right_prob - energy.sum())
            assert min(other, left_prob, right_prob) >= -2e-13 and np.all(energy >= joint - 2e-13)
            if ti < 2 and i < 2:
                old = archived["inputs"][i]["times"][ti]
                differences = [abs(right_prob - old["charge_arrival_probability"])]
                for j, out in enumerate(old["destination_targets"]):
                    differences.extend((abs(out["right_whole_charged_mode_probability_left_unread"] - energy[j]),
                                        abs(out["whole_swapped_product_probability"] - joint[j])))
                reference_error = max(reference_error, *differences)
            reduced[ti].append(projections)
            arrivals[ti].append(right_prob)
            records[ti].append({"input_level": i, "right_charge": right_prob, "left_charge": left_prob,
                                "right_three_energies": energy.tolist(), "right_other_energy": other,
                                "right_three_energies_and_left_neutral_ground": joint.tolist(),
                                "matching_level_given_right_charge": float(energy[i]) / right_prob,
                                "herald_and_arrive_matching_level": source["inputs"][i]["ideal_herald_probability"] * float(energy[i])})
        del outputs, left, right
    assert max_norm_error < 3e-12 and reference_error < 3e-12
    times = []
    for ti, t in enumerate(TIMES):
        ch, truth = coherent_readout(reduced[ti], neutral)
        np.testing.assert_allclose(truth, np.array([x["right_three_energies"] for x in records[ti]]), atol=2e-13)
        ch["correct_label_given_right_charge_equal_prior"] = float(np.trace(truth) / sum(arrivals[ti]))
        times.append({"time_native": t, "records": records[ti], "channel": ch})
    return {"sign": sign, "bridge": [int(g), int(receiving)], "joint_dimension": 2 * len(states) * len(idx),
            "charged_dimension": len(states), "neutral_dimension": len(idx), "lowest_four_native_energies": levels[:4].tolist(),
            "first_three_relative_to_ground_energy": (levels[:3] / levels[0]).tolist(),
            "preparation": source, "times": times, "max_eigenpair_residual": residual,
            "max_total_norm_error": max_norm_error, "max_archived_observable_discrepancy": reference_error,
            "series_precision": tolerance,
            "controls": {"bridge_off_arrival": 0, "seam_only_at_pi_over_8_arrival": gateway_source_weights.tolist(),
                         "seam_only_given_arrival_right_state": "fixed uniform gateway state for every input",
                         "seam_only_three_level_channel": seam_channel,
                         "seam_only_three_right_energies": seam_energy.tolist(),
                         "neutral_ground_transfer_initial_derivative": derivative.tolist(),
                         "neutral_ground_transfer_initial_derivative_rank": int(np.linalg.matrix_rank(derivative))}}


def run(root):
    hashes, reference = extract_sources(root)
    from bridge import setups

    sectors = []
    for data, old in zip(setups(), reference["charge_sectors"], strict=True):
        assert data[1] == old["sign"]
        sectors.append(sector(data, old))
    return {"construction": "Gauss-projected field extension, two fully active matched cells with one open 3–48 seam",
            "selected_modes": "first three nondegenerate charged-cell eigenvalues in each sign, without laboratory inputs",
            "verified_archived_dependencies": hashes, "sectors": sectors}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    with TemporaryDirectory() as work:
        result = run(Path(work).resolve())
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
