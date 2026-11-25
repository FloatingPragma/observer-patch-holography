"""Independent two-arm material flight at matched radar volume."""

import json

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import expm_multiply

from .. import measure as cuts
from .... import station


def physical_means(delta):
    native, reg = cuts.native, cuts.regional
    atlas = native.Atlas(1)
    board, patch = atlas.board, native.chosen(atlas)
    route = station.iso.pairing.formation.region(board)
    src = reg.bridge.beam.indices(sum(reg.bridge.beam.screens(board, route), []))
    edges = native.flux_on_edges(board, patch)
    cells = sorted(patch["cells"])
    out = np.zeros((2, 10 + len(cells)))
    for a, q in enumerate((-1, 1)):
        for odd in (False, True):
            h, native_ops = native.op_table(board, route, patch, q, odd)
            generator = -2j * h
            starts = np.zeros((h.shape[0], len(src)))
            starts[src] = np.eye(len(src))
            w = expm_multiply(generator, starts, traceA=generator.diagonal().sum())
            ops = [native_ops[name] for name in reg.NAMES]
            ops.extend(
                (sum(op for (u, _), op in edges.items() if atlas.dist[0][u] == d), 4)
                for d in range(6)
            )
            ops.extend(
                (sum(op for (u, _), op in edges.items() if u == target), 4)
                for target in cells
            )
            for j, (op, unit) in enumerate(ops):
                ends = np.unique(op.nonzero()[0])
                block = w[ends]
                acted = op[ends, :][:, ends] @ block
                observable = sum(
                    x.conj().T @ y
                    for x, y in zip(np.hsplit(block, 10), np.hsplit(acted, 10))
                ) / (10 * unit)
                out[a, j] += (np.trace(observable @ delta[int(odd)])).real
    return out[:, :4].ravel(), out[:, 4:10], out[:, 10:]


def run(regional, shells, cell_prices):
    j, s = (cuts.mid(x) for x in (regional, shells))
    raw = json.loads((cuts.survey.HERE.parent / "results.json").read_text())
    x = np.array(raw["x_real"]) + 1j * np.array(raw["x_imag"])
    rad = station.radar_menus()[6]
    ranges = station.pair(rad, x)
    ar, v = station.geometry(ranges)
    joint = station.iso.basis(True)

    def states(d):
        return np.einsum("ba,aij->bij", d.reshape(2, 9), joint)

    directions = np.stack([station.pair(rad, states(d)) for d in np.eye(18)], axis=1)
    da, dv = np.asarray(station.range_gradient(ranges)) @ directions
    summary = np.vstack((j, (s[0] + s[1])[1:]))
    summary = summary / np.linalg.norm(summary, axis=1)[:, None]
    kernel = np.linalg.svd(summary, full_matrices=True)[2][13:]
    adjust = (dv @ kernel.T) @ kernel
    adjust /= np.linalg.norm(adjust)
    all_rows = np.vstack((summary, dv / np.linalg.norm(dv)))
    fixed = np.linalg.svd(all_rows, full_matrices=True)[2][14:]
    direction = -(da @ fixed.T) @ fixed
    direction /= np.linalg.norm(direction)
    arms = []
    volumes = []
    areas = []
    step = 4e-5
    for sign in (-1, 1):
        offset = sign * step * (directions @ direction)
        shift = brentq(
            lambda y: station.geometry(ranges + offset + y * (directions @ adjust))[1]
            - v,
            -1e-6,
            1e-6,
            xtol=1e-16,
        )
        delta = states(sign * step * direction + shift * adjust)
        assert (
            np.linalg.norm(delta) < 1e-4 and np.linalg.eigvalsh(x + delta).min() > 0.02
        )
        arms.append(delta)
        area, volume = station.geometry(station.pair(rad, x + delta))
        volumes.append(float(volume))
        areas.append(float(area))
    delta = (arms[1] - arms[0]) / (2 * step)
    jr, sr, cell_flux = physical_means(delta)
    material_slope = (
        np.array(cell_prices, dtype=float) @ np.r_[jr, cell_flux.sum(axis=0)]
    )
    ref_shell = s.sum(axis=0)
    signed = float(
        (station.entropy(x + arms[1]) - station.entropy(x + arms[0])) / (2 * step)
    )
    return dict(
        arm_step=step,
        arm_volumes=volumes,
        ordinary_calibrated_area=-raw["ordinary_eta"]
        * (areas[1] - areas[0])
        / (2 * step),
        spectral_entropy=signed,
        common_cell_flux_slope=float(material_slope),
        regional_slopes=jr.tolist(),
        vacancy_and_excess_shell_slopes=sr.tolist(),
        max_relative_regional_drift=float(max(abs(jr) / np.linalg.norm(j, axis=1))),
        max_relative_equal_mix_drift=float(
            max(abs(sr.sum(axis=0)) / np.linalg.norm(ref_shell, axis=1))
        ),
    )
