"""Local spectrum, signed profiles and vacancy-transfer work in two matched cells."""

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation

HERE = Path(__file__).resolve().parent
STARTS = 2**34
ARMS = 2**16
STEP = 1 / 128
PANELS = (
    ("native", candidate.KAPPA, True),
    ("no_reference", candidate.KAPPA, False),
    ("exchange", 0, True),
)
NAMES = {"neutral": (2, 2), "positive": (2, 1), "vacancy": (0, 1)}
DIST = shortest_path(cell().adjacency(), directed=False, unweighted=True).astype(int)
EDGES = np.array(cell().seams)


def sector(plus, minus):
    words = []
    for p in combinations(range(12), plus):
        for m in combinations(sorted(set(range(12)) - set(p)), minus):
            row = np.zeros(12, dtype=np.int8)
            row[list(p)] = 1
            row[list(m)] = -1
            words.append(row)
    w = np.array(words)
    ids = {v.tobytes(): i for i, v in enumerate(w)}
    swap = []
    for a, b in EDGES:
        perm = w.copy()
        perm[:, [a, b]] = w[:, [b, a]]
        swap.append([ids[v.tobytes()] for v in perm])
    return w, np.array(swap)


def native_geometry():
    board = federation(0)
    mask = board.degree[:12] == 5
    comp = next(
        i
        for i in range(1, board.cells)
        if np.array_equal(board.degree[12 * i : 12 * i + 12], board.degree[:12])
    )
    law = candidate.native_law
    ex, co = law.local_term_parts()
    swap = np.array([law.SRC.index((b, a)) for a, b in law.SRC])
    expected = np.vstack([0.5 * (np.eye(9) - np.eye(9)[swap]), np.zeros((90, 9))])
    np.testing.assert_allclose(ex[:, :9], expected, atol=1e-13)
    np.testing.assert_allclose(
        co[:, :9],
        np.vstack([np.diag([(a - b) ** 2 for a, b in law.SRC]), np.zeros((90, 9))]),
        atol=1e-13,
    )
    return mask, comp


def operators(w, swap, mask):
    i = np.arange(len(w))
    kinetic = np.zeros((len(w), len(w)))
    for col in swap:
        kinetic[i, i] += 0.5
        kinetic[i, col] -= 0.5
    contrast = (w[:, EDGES[:, 0]] - w[:, EDGES[:, 1]]) ** 2
    reference = (w[:, mask] + 1) ** 2 - 1
    return kinetic, contrast, reference


def filtered(kinetic, contrast, reference, k, enable):
    h = kinetic.copy()
    h.flat[:: len(h) + 1] += k * (
        contrast.sum(axis=1) + int(enable) * reference.sum(axis=1)
    )
    e, v = eigh(h, driver="evd", overwrite_a=True, check_finite=False)
    target = np.rint(e[0] * 128) / 128
    delta = (e - target) * STEP / (2 * np.pi)
    gain = (np.sinc(ARMS * delta) / np.sinc(delta)) ** 4
    prob = gain.sum() / len(e)
    weight = gain / gain.sum()
    rho = (v * weight) @ v.T
    return rho, dict(
        target=float(target),
        floor=float(e[0]),
        acceptance=float(prob),
        population_within_half_bin=float(
            weight[abs(e - target) < np.pi / (ARMS * STEP)].sum()
        ),
        mean_energy=float(weight @ e),
        lowest_gap=float(e[1] - e[0]),
    )


def graph_read(w):
    occupied = (w != 0).astype(float)
    cost = occupied @ (DIST**2)
    median = cost == cost.min(axis=1)[:, None]
    median = median / median.sum(axis=1)[:, None]
    near = ~(DIST > 1)
    groups = []
    for mark in occupied:
        ports = np.flatnonzero(mark)
        seen = {ports[0]}
        for _ in ports:
            seen.update(a for a in ports if any(near[a, b] for b in seen))
        groups.append(len(seen) == len(ports))
    shells = np.stack([w @ (DIST == r).astype(int) for r in range(4)], axis=2)
    return median, shells, groups


def meter(rng, names, probs, acceptance):
    probs = np.asarray(probs, float)
    if probs.min() < -1e-10 or abs(probs.sum() - 1) > 1e-9:
        raise ArithmeticError("Invalid native read distribution")
    probs = np.maximum(0, probs)
    n = int(rng.binomial(STARTS, acceptance))
    draw = rng.multinomial(n, probs / probs.sum())
    return dict(rejected=STARTS - n, accepted=n, keys=names, counts=draw.tolist())


def word_read(rng, w, diagonal, graph, mask, acceptance):
    med, shells, groups = graph
    selections = {"intrinsic": med}
    for g in np.flatnonzero(~mask):
        selections[f"gateway_{g}"] = np.broadcast_to(np.eye(12)[g], med.shape)
    answer = {}
    for name, selection in selections.items():
        active = np.nonzero(selection)
        keys = np.column_stack([shells[active], np.asarray(groups)[active[0]]])
        labels, inv = np.unique(keys, axis=0, return_inverse=True)
        pmf = np.bincount(
            inv, weights=(diagonal[:, None] * selection)[active], minlength=len(labels)
        )
        row = meter(rng, labels.tolist(), pmf, acceptance)
        sample = np.array(row["counts"]) / row["accepted"]
        row["signed_means"] = (sample @ labels[:, :4]).tolist()
        row["connected_fraction"] = float(sample @ labels[:, 4])
        row["true_signed_means"] = (pmf @ labels[:, :4]).tolist()
        # Each shell has expectation bounded by +/-2 in the neutral preparation.
        row["shell_radius"] = float(
            4 * np.sqrt(np.log(2 * 2500 / 0.01) / (2 * row["accepted"]))
        )
        answer[name] = row
    return answer


def local_read(rng, rho, swap, contrast, reference, acceptance, k, enable):
    n = len(rho)
    diag = np.diag(rho)
    off = rho[np.arange(n)[:, None], swap.T]
    odd = (diag[:, None] - off) / 2
    terms = np.zeros((2, 3))
    for j, d in enumerate((0, 1, 4)):
        mask = contrast == d
        terms[1, j] = np.sum(np.where(mask, odd, 0)) / 30
        terms[0, j] = np.sum(np.where(mask, diag[:, None] - odd, 0)) / 30
    labels = [[bit, d] for bit in (0, 1) for d in (0, 1, 4)]
    seam = meter(rng, labels, terms.ravel(), acceptance)
    nr = reference.shape[1]
    pmf = [(diag[:, None] * (reference == v)).sum() / nr for v in (-1, 0, 3)]
    ref = meter(rng, [-1, 0, 3], pmf, acceptance)
    measured_seam = (
        np.array(seam["counts"])
        / seam["accepted"]
        @ (np.array(labels) @ np.array([1, k]))
    )
    measured_ref = np.array(ref["counts"]) / ref["accepted"] @ np.array([-1, 0, 3])
    energy = 30 * measured_seam + int(enable) * nr * k * measured_ref

    def radius(span, n):
        return span * np.sqrt(np.log(2 * 2500 / 0.01) / (2 * n))

    err = 30 * radius(1 + 4 * k, seam["accepted"]) + int(enable) * nr * k * radius(
        4, ref["accepted"]
    )
    return dict(
        seam=seam,
        reference=ref,
        energy=float(energy),
        radius=float(err),
        exchange_only_on_same_source=float(
            30 * np.array(seam["counts"])[3:].sum() / seam["accepted"]
        ),
    )


def run():
    mask, comp = native_geometry()
    rng = np.random.default_rng(14272026)
    result = dict(
        starts=STARTS, step=STEP, arms=ARMS, data_owner=0, storage_owner=comp, panels={}
    )
    for name, pop in NAMES.items():
        w, swap = sector(*pop)
        graph = graph_read(w)
        kinetic, co, ref = operators(w, swap, mask)
        result["panels"][name] = {}
        for tag, k, enable in PANELS:
            rho, spec = filtered(kinetic, co, ref, k, enable)
            support = spec["acceptance"]
            meters = local_read(rng, rho, swap, co, ref, support, k, enable)
            if abs(meters["energy"] - spec["mean_energy"]) > meters["radius"]:
                raise ArithmeticError("Independent local energy meter failed")
            profile = word_read(rng, w, np.diag(rho), graph, mask, support)
            result["panels"][name][tag] = dict(
                spectrum=spec, local_energy=meters, profile=profile
            )
        result["panels"][name]["unfiltered_profile"] = word_read(
            rng, w, np.full(len(w), 1 / len(w)), graph, mask, 1.0
        )
    result["accounting"] = {}
    for tag, _, _ in PANELS:
        rows = [result["panels"][n][tag]["local_energy"] for n in NAMES]
        n, p, v = (r["energy"] for r in rows)
        result["accounting"][tag] = dict(
            neutral_minus_positive=n - p,
            work_to_transfer_vacancy=n - p - v,
            energy_difference_radius=sum(r["radius"] for r in rows[:2]),
            work_radius=sum(r["radius"] for r in rows),
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    out = parser.parse_args().record
    result = run()
    if out:
        with out.open("x") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(
            result,
            json.loads((HERE / "records.json").read_text()),
            "Baryon partner and profile",
        )
    for k, r in result["accounting"].items():
        q = result["panels"]["neutral"][k]["profile"]["intrinsic"]
        print(k, r, q["signed_means"], q["connected_fraction"])
