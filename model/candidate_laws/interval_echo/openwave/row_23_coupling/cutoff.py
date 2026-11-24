"""Locate the first spatial step in the same effective interaction at fixed energy."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.executable_model import carrier
from . import experiment as e


def vertex_error(v, dv, g):
    value = e.dressed(v, g)
    inv = e.norm(np.linalg.inv(e.IDENTITY - v @ g))
    margin = 1 - inv * dv * e.norm(g)
    if margin <= 0:
        raise ValueError("Unresolved matched matrix")
    return value, inv * dv * (1 + e.norm(g @ value)) / margin


def scan(z, v, dv, energy, vec, ports, radius):
    local = vec[ports[:2]]
    per_mode = np.einsum("im,jm,m->mij", local, local, 1 / (z - energy))
    full = per_mode.sum(axis=0)
    cumulative = per_mode.cumsum(axis=0)
    ends = np.flatnonzero(np.r_[np.diff(energy) > 1e-8, True])
    rows = []
    for i in ends:
        high = full - cumulative[i]
        vertex, err = vertex_error(v, dv, high)
        even = abs(e.PARITY[:, 0] @ vertex @ e.PARITY[:, 0])
        if not rows:
            initial, uncertainty = even, err
        ratio = even / initial
        # Enclose the even modulus ratio simultaneously from the matrix bounds.
        lo = max(0, even - err) / (initial + uncertainty)
        hi = (even + err) / (initial - uncertainty)
        tag = (
            "above"
            if hi < 0.9 or lo > 1.1
            else "below"
            if lo > 0.9 and hi < 1.1
            else "unresolved"
        )
        k = float(np.sqrt(2 * max(0, energy[i])))
        rows.append(
            {
                "retained_modes": int(i + 1),
                "k": k,
                "k_radius": k * radius,
                "even_modulus": float(even),
                "ratio": float(ratio),
                "ratio_interval": [float(lo), float(hi)],
                "threshold": tag,
            }
        )
    above = next(j for j, r in enumerate(rows) if r["threshold"] == "above")
    below = max(j for j in range(above) if rows[j]["threshold"] == "below")
    return rows, [below, above]


def selected_reads(z, v, dv, energy, vec, ports, rows, selected, target, rng):
    local = vec[ports[:2]]
    gfull = (local / (z - energy)) @ local.T
    sources, specs = [], []
    for j in selected:
        n = rows[j]["retained_modes"]
        weights = vec[ports[2], :n]
        p = np.linalg.norm(weights)
        sources.append(vec[:, :n] @ weights / p)
        low = (local[:, :n] / (z - energy[:n])) @ local[:, :n].T
        x = (weights / (z - energy[:n])) @ local[:, :n].T / p
        specs.append((gfull - low, low, x))
    scores = []
    archive = {}
    for name, energies, vectors in [("blank", energy, vec), ("written", *target)]:
        archive[name] = e.returns(energies, vectors, np.asarray(sources), rng)
        real, imag = 2 * archive[name] / e.N - 1
        scores.append((-1j * e.WEIGHTS * np.exp(1j * z * e.TIMES)) @ (real + 1j * imag))
    blank, written = scores
    out = []
    for i, (hi, low, x) in enumerate(specs):
        answers = []
        for ext in (hi, specs[0][0]):
            t, err = vertex_error(v, dv, low + ext)
            predicted = blank[i] + x @ t @ x
            answers.append(
                {
                    "return": e.encoded(predicted),
                    "absolute_error": float(abs(predicted - written[i])),
                    "joint_bound": 2 * e.allowance(z) + np.linalg.norm(x) ** 2 * err,
                }
            )
        out.append(
            {
                "shell_index": selected[i],
                "held_raw_return": e.encoded(written[i]),
                "matched": answers[0],
                "frozen_pre_crossing": answers[1],
            }
        )
    return out, archive


def compute():
    rng = np.random.default_rng(230538)
    banks = {}
    studies = []
    previous = json.loads((e.HERE / "results.json").read_text())["studies"]
    with np.load(e.HERE / "results.npz") as raw:
        for level in (0, 1):
            board = carrier.federation(level)
            seam = carrier.EMITTER_SEAM[level]
            lap, energy, vec, ports, _, _, _ = e.geometry(board, seam)
            h = lap.copy()
            h[np.ix_(seam, seam)] -= e.candidate.KAPPA * e.IDENTITY + np.outer(
                e.PARITY[:, 1], e.PARITY[:, 1]
            )
            target = eigh(h)
            for width in (0.5, 0.25):
                base = [
                    e.measure(raw[f"level{level}_{name}"], 1j * width)[0]
                    for name in ("blank", "written")
                ]
                _, v, dv, _ = e.interaction(*base, 1j * width)
                z = -1 + 1j * width
                rows, selected = scan(
                    z,
                    v,
                    dv,
                    energy,
                    vec,
                    ports,
                    previous[level]["core"]["radius_links"],
                )
                checks, archive = selected_reads(
                    z, v, dv, energy, vec, ports, rows, selected, target, rng
                )
                banks.update(
                    {
                        f"level{level}_width{width}_{key}": value
                        for key, value in archive.items()
                    }
                )
                studies.append(
                    {
                        "cells": board.n // 12,
                        "width": width,
                        "shells": rows,
                        "crossing_indices": selected,
                        "held_checks": checks,
                    }
                )
    return {"studies": studies}, banks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result, raw = compute()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(args.output.with_suffix(".npz"), **raw)
    else:
        with np.load(e.HERE / "cutoffs.npz") as bank:
            for key in raw:
                np.testing.assert_allclose(raw[key], bank[key], atol=2, rtol=0)
        earlier = json.loads((e.HERE / "cutoffs.json").read_text())["studies"]
        for a, b in zip(earlier, result["studies"], strict=True):
            assert a["crossing_indices"] == b["crossing_indices"]
            for deposited, recomputed in zip(a["shells"], b["shells"], strict=True):
                assert (
                    deposited["threshold"] == recomputed["threshold"]
                    and deposited["retained_modes"] == recomputed["retained_modes"]
                )
                np.testing.assert_allclose(
                    deposited["ratio_interval"], recomputed["ratio_interval"], atol=2e-7
                )
            for deposited, recomputed in zip(
                a["held_checks"], b["held_checks"], strict=True
            ):
                for key in ("matched", "frozen_pre_crossing"):
                    np.testing.assert_allclose(
                        deposited[key]["absolute_error"],
                        recomputed[key]["absolute_error"],
                        atol=2e-7,
                    )
        print("Reproduced the spatial cutoff scan and independent boundary returns.")


if __name__ == "__main__":
    main()
