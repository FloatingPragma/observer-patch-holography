"""Locate the distinction between the two archived spectroscopy preparations."""

import json
from itertools import combinations
from math import sqrt

import numpy as np
from scipy.linalg import expm

from .. import raw_verify as raw
from ..transfer_clock import infer


def physical_partition(tree):
    edges, _, _ = raw.board(0)
    selected = raw.spectroscopy()
    rays = np.zeros((len(tree.states), len(selected)), dtype=int)
    for j, (_, paths) in enumerate(selected):
        indices = [tree.words.index(path) for path in paths]
        rays[indices, j] = (1, -1)
    np.testing.assert_array_equal(tree.h2 @ rays, rays)
    np.testing.assert_array_equal(rays.T @ rays, 2 * np.eye(2, dtype=int))

    # Compare written physical heads at equal terminal load, for each tip.
    difference = set()
    for port in (9, 11):
        left, right = [
            tree.states[int(np.flatnonzero((tree.sites == port) & (ray != 0)).item())]
            for ray in rays.T
        ]
        difference.update(left[1] ^ right[1])
    triangle = sorted(edges[i] for i in difference)
    if triangle != [(0, 1), (0, 5), (1, 5)]:
        raise ValueError(f"Different physical records distinguish the columns: {triangle}")

    # A material word is a tensor product of a load position and local heads.
    # The remaining unchanged collar factors are common to both preparations.
    inside = [tuple(sorted(flips & difference)) for _, flips in tree.states]
    outside = [(port, tuple(sorted(flips - difference))) for port, flips in tree.states]
    a = {label: i for i, label in enumerate(sorted(set(inside)))}
    b = {label: i for i, label in enumerate(sorted(set(outside)))}
    indexes = ([a[x] for x in inside], [b[x] for x in outside])

    def coefficients(amplitudes):
        matrix = np.zeros((len(a), len(b)), dtype=amplitudes.dtype)
        matrix[indexes] = amplitudes
        return matrix

    if len(set(zip(*indexes, strict=True))) != len(tree.states):
        raise ValueError("The partition loses a physical word")
    first, second = (coefficients(ray) for ray in rays.T)
    common = first.T @ first
    np.testing.assert_array_equal(second.T @ second, common)
    np.testing.assert_array_equal(first.T @ second, np.zeros_like(common))
    labels = []
    for matrix in (first, second):
        active = np.flatnonzero(np.any(matrix, axis=1))
        labels.append(active.item())
    if labels[0] == labels[1]:
        raise ValueError("The source triangle does not separate the ideal inputs")
    seam = edges.index((0, 1))
    old_head = np.array([int(seam in flips) for _, flips in tree.states])
    np.testing.assert_array_equal(rays.T @ (old_head[:, None] * rays), np.diag([0, 2]))
    for port in (9, 11):
        np.testing.assert_array_equal(rays.T @ ((tree.sites == port)[:, None] * rays),
                                      np.eye(2, dtype=int))

    # The antisymmetric hopping part is the signed current on a history link.
    for i, j in combinations(range(len(tree.states)), 2):
        if tree.h2[i, j]:
            pair = np.outer(rays[i], rays[j])
            np.testing.assert_array_equal(pair, pair.T)
    return triangle, [column for column, _ in selected], coefficients


def compiled_readings(tree, columns, coefficients):
    words = raw.read_json(raw.DONOR / "analyzers.json")
    clocks = raw.read_json(raw.HERE / "clock.json")["runs"]
    evidence = raw.read_json(raw.HERE / "provenance" / "terminal.json")["launches"]
    witness = raw.read_json(raw.HERE / "donor_release.json")["launches"]
    times = (0, 1, 7, 37, 150)
    measured_floor, calculated_floor, best_advantage = 1.0, 1.0, 0.0
    comparisons = 0
    unresolved = []
    for bank, bound, launch in zip(clocks, evidence, witness, strict=True):
        strength = bank["strength_in_simulator"]
        if not (strength == bound["strength_in_simulator"] == launch["strength_in_simulator"]):
            raise ValueError("Clock and spectroscopy inputs differ")
        clock = infer(bank["counts"])
        if clock is None:
            unresolved.append(strength)
            continue
        if not clock[1] <= abs(strength) <= clock[2]:
            raise ValueError("The clock interval misses its simulated strength")
        if (set(columns) != {entry["column"] for entry in bound["launches"]} or
                set(columns) != {entry["column"] for entry in launch["launches"]}):
            raise ValueError("Different preparation in the finite-error bound")
        measured_floor = min(measured_floor, *(x["sampling_all_time_floor"] for x in launch["launches"]))
        calculated_floor = min(calculated_floor, *(x["all_exposure_signed_ray_floor_per_original"]
                                                    for x in bound["launches"]))
        ratio = strength / clock[0]
        pulse = expm(-4j * ratio * words["step"] * tree.h2)
        sources = []
        for column in columns:
            psi = np.eye(len(tree.sites), dtype=complex)[0]
            for angles in words["pulses"][column]["angles"]:
                psi = pulse @ psi
                for port, angle in zip(words["ports"], angles, strict=True):
                    psi[tree.sites == port] *= np.exp(-1j * ratio * angle)
            sources.append(psi)
        for time in times:
            flight = expm(-0.5j * strength * time * tree.h2)
            reduced = []
            for psi in sources:
                x = coefficients(flight @ psi)
                reduced.append(x.T @ x.conj())
            # The common exhaustion flag conveys no information. The success-
            # conditioned advantage is an upper bound per original start.
            advantage = np.abs(np.linalg.eigvalsh(reduced[0] - reduced[1])).sum() / 4
            best_advantage = max(best_advantage, float(advantage))
            comparisons += 1
    if unresolved != [0] or comparisons != 80 or best_advantage > sqrt(1 - calculated_floor):
        raise ValueError("Invalid compiled-source comparison")
    return {"comparisons": comparisons, "inconclusive_strengths": unresolved,
            "max_equal_prior_success_at_five_waits": 0.5 + best_advantage,
            "all_time_success_upper_from_calculated_weight": 0.5 + sqrt(1 - calculated_floor),
            "all_time_success_upper_from_99_percent_witness": 0.5 + sqrt(1 - measured_floor)}


def verify():
    tree = raw.component(0)
    triangle, columns, coefficients = physical_partition(tree)
    return {"source_columns": columns, "prefix_seams": triangle,
            "ideal_success_outside_triangle": "1/2",
            "ideal_success_from_head_0_1": 1,
            "compiled": compiled_readings(tree, columns, coefficients)}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, allow_nan=False))
