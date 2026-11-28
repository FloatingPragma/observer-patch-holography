"""Measure the finite pairing dose with a separately released material beam."""

import argparse
import json
from hashlib import sha256
from os.path import relpath

import numpy as np
from matched_einstein_limit.source_implication.directional_stress.translation.counts import (
    selection_norm,
)
from matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts import (
    measure as cuts,
)
from scipy.sparse.linalg import expm_multiply

from . import finite_dose_curvature as curve

lab = curve.lab
iv = curve.iv
reg = curve.reg
N = 2**68


def active_ports(op, unit):
    """Total selection weight on each port that yields a nonzero analyzer outcome."""
    axes = curve.dose.heads.PRODUCTS.reshape(16, 4, 4)
    rows, cols = op.nonzero()
    active = np.zeros(op.shape[0] // 4)
    for i, j in sorted(
        {(a // 4, b // 4) for a, b in zip(rows, cols, strict=True) if a // 4 <= b // 4}
    ):
        block = op[4 * i : 4 * i + 4, 4 * j : 4 * j + 4].toarray() / unit
        coef = np.einsum("kij,ji->k", axes, block) / 4
        weight = (abs(coef.real) + abs(coef.imag)).sum()
        active[i] += weight
        if j != i:
            active[j] += weight
    return np.repeat(active, 4)


def flight(names=cuts.regional.NAMES, *, include_edges=True):
    atlas = cuts.native.Atlas(1)
    board, patch = atlas.board, cuts.native.chosen(atlas)
    bridge = cuts.regional.bridge
    route = bridge.native.iso.pairing.formation.region(board)
    launch = bridge.beam.indices(sum(bridge.beam.screens(board, route), []))
    edges = cuts.native.flux_on_edges(board, patch) if include_edges else {}
    effects, scales, occupancy = {}, {}, {}
    for charge in (-1, 1):
        for odd in (False, True):
            h, ops = cuts.native.op_table(board, route, patch, charge, odd)
            grid = -2j * h
            start = np.zeros((h.shape[0], len(launch)))
            start[launch] = np.eye(len(launch))
            w = expm_multiply(grid, start, traceA=grid.diagonal().sum())
            chosen = {(charge, k): ops[k] for k in names}
            if include_edges:
                chosen.update(
                    {k: (edges[k[1:]], 4) for k in cuts.CHANNELS if k[0] == charge}
                )
            for key, (op, unit) in chosen.items():
                terms = [a.conj().T @ (op @ a) / unit for a in np.hsplit(w, 10)]
                effects.setdefault(key, []).append(np.mean(terms, axis=0))
                scales[key] = max(scales.get(key, 0), selection_norm(op, unit))
                weights = active_ports(op, unit)
                hit_terms = [
                    a.conj().T @ (weights[:, None] * a) for a in np.hsplit(w, 10)
                ]
                occupancy.setdefault(key, []).append(np.mean(hit_terms, axis=0))
    return effects, scales, occupancy


def count(rng, p):
    # Bound individual binomial sizes and subtotal sums below binary64 and int64 limits.
    blocks = rng.binomial(2**50, p, size=N // 2**50).reshape(-1, 1024)
    return sum(map(int, blocks.sum(axis=1)))


def port_read(sampler, n, hit_prob, mu, scale, loglevel):
    """Read the signed terminal port, assigning zero to each missed appointment."""
    assert abs(mu) <= hit_prob + 1e-14 and 0 < hit_prob <= 1 + 1e-14
    plus_given_hit = (1 + mu / hit_prob) / 2
    plus = minus = 0
    block = 2**52
    assert n % (block * 1024) == 0
    for _ in range(n // (block * 1024)):
        arrivals = sampler.binomial(block, min(1, hit_prob), size=1024)
        chosen = sampler.binomial(arrivals, np.clip(plus_given_hit, 0, 1))
        plus += int(chosen.sum())
        minus += int((arrivals - chosen).sum())
    hit_upper = 2 * (plus + minus) / n + 4 * loglevel / n
    width = (
        scale * (np.sqrt(2 * hit_upper * loglevel / n) + 4 * loglevel / (3 * n))
        + 1e-14
        + 3e-12 * scale
    )
    value = scale * (plus - minus) / n
    return [plus, minus, n - plus - minus], reg.band(value - width, value + width)


def run():
    iv.dps = 35
    source, states, h, hband = curve.trials()
    effects, scales, occupancy = flight()
    lower = cuts.HERE / "effects.json"
    regional = cuts.survey.HERE / "translation.json"
    rank_file = cuts.HERE / "results.json"
    site, aggregate, rank = [
        json.loads(f.read_bytes()) for f in (lower, regional, rank_file)
    ]
    keys = [(q, k) for q in (-1, 1) for k in aggregate["order"]]
    maps = dict(zip(keys, aggregate["rows"], strict=True))
    maps.update({tuple(r["channel"]): r["row"] for r in site["cuts"]})
    assert site["scale"] == aggregate["scale"] == rank["scale"]
    tape = source["tapes"]

    def observations(name):
        return curve.dose.source.observed(tape[name]["incoming"], lab.SHOTS, False)

    x, ex = observations("reference")
    y, ey = observations("held")
    basis = curve.dose.heads.measure.joint_basis()
    coeff = curve.station.pair(basis, y - x)
    errors = np.repeat((ex + ey) / 2, 9)
    coords = [
        reg.band(c - e - 1e-14, c + e + 1e-14)
        for c, e in zip(coeff, errors, strict=True)
    ]
    rng = np.random.default_rng(28112083)
    precision_rng = np.random.default_rng(28112085)
    force_rng = np.random.default_rng(28112086)
    ordinary_energy = reg.band(
        *json.loads((lab.HERE / "ordinary_stationarity.json").read_bytes())[
            "parked_energy_slope_from_counts"
        ]
    )
    rows = []
    state_keys = ("reference", "held", "matched_ordinary", "ordinary")
    for key in [*keys, *cuts.CHANNELS]:
        scale = scales[key]
        effect = effects[key]
        assert scale > 0
        means = {k: float(curve.station.pair(effect, states[k])) for k in state_keys}
        hits = {k: count(rng, 0.5 + v / (2 * scale)) for k, v in means.items()}
        err = (
            scale * np.sqrt(2 * np.log(2 * len(state_keys) * len(effects) / 0.01) / N)
            + 3e-12 * scale
            + 1e-14
        )
        reads = {
            k: reg.band((2 * n - N) / N * scale - err, (2 * n - N) / N * scale + err)
            for k, n in hits.items()
        }
        prices = [iv.mpf(pair) / aggregate["scale"] for pair in maps[key]]
        predicted = sum(p * c for p, c in zip(prices, coords, strict=True))
        measured = reads["held"] - reads["reference"]
        population = means["held"] - means["reference"]
        assert (
            measured.a <= population <= measured.b
            and predicted.a <= population <= predicted.b
        )
        ordinary = hband * (reads["ordinary"] - reads["reference"])
        diff = measured - ordinary
        precision = {}
        exposures = (N, 2**72, 2**76) if key[1] == "marked_force" else (N, 2**72)
        for n in exposures:
            sampler = force_rng if n == 2**76 else precision_rng
            detailed = {}
            intervals = {}
            # Include both exposures in the common signed-score and hit-rate bound.
            loglevel = np.log(8 * len(state_keys) * len(effects) / 0.01)
            for k in state_keys:
                hit_prob = float(curve.station.pair(occupancy[key], states[k])) / scale
                mu = means[k] / scale
                detailed[k], intervals[k] = port_read(
                    sampler, n, hit_prob, mu, scale, loglevel
                )
            contrast = intervals["held"] - intervals["reference"]
            calibrated = hband * (intervals["ordinary"] - intervals["reference"])
            assert contrast.a <= population <= contrast.b
            precision[str(n)] = {
                "plus_minus_miss": detailed,
                "contrast": lab.six.bounds(contrast),
                "held_minus_equal_entropy_ordinary": lab.six.bounds(
                    contrast - calibrated
                ),
            }
            force_price = intervals["ordinary"] - intervals["reference"]
            if not force_price.a <= 0 <= force_price.b:
                precision[str(n)]["ordinary_energy_equivalent"] = lab.six.bounds(
                    ordinary_energy * contrast / force_price
                )
                precision[str(n)]["source_inferred_energy_equivalent"] = lab.six.bounds(
                    ordinary_energy * predicted / force_price
                )
        rows.append(
            {
                "channel": list(key),
                "analyzer_scale": scale,
                "positive_counts": hits,
                "resolved_port_instrument": precision,
                "predicted_from_joint_head_counts": lab.six.bounds(predicted),
                "independent_beam_contrast": lab.six.bounds(measured),
                "direct_propagation": population,
                "matched_ordinary_beam_contrast": lab.six.bounds(
                    reads["matched_ordinary"] - reads["reference"]
                ),
                "ordinary_from_counted_source_match": lab.six.bounds(ordinary),
                "held_minus_matched_ordinary": lab.six.bounds(diff),
            }
        )
    sigma = rank["min_singular_lower"] / rank["scale"]
    C = curve.station.log_blocks(states["independent"]) - curve.station.log_blocks(
        states["reference"]
    )
    K = np.linalg.norm(curve.station.pair(basis, C))
    q = source["held_connected_modular_gain"][0]
    return {
        "seed": 28112083,
        "port_resolved_seed": 28112085,
        "force_repeat_seed": 28112086,
        "shots_per_arm_instrument": N,
        "readings": rows,
        "finite_material_norm_lower_from_measured_charge": float(
            sigma * q / (K + 1e-12)
        ),
        "calibrated_material_minimum_singular_value_lower": sigma,
        "frozen_connected_score_HS_norm": float(K),
        "input_sha256": {
            relpath(p, lab.HERE): sha256(p.read_bytes()).hexdigest()
            for p in (
                lower,
                regional,
                rank_file,
                lab.HERE / "FINITE_DOSE_CURVATURE.md",
                lab.HERE / "MATERIAL_HIT_PRECISION.md",
                lab.HERE / "MATERIAL_FORCE_MASS.md",
                lab.HERE / "finite_pairing_dose.json",
                lab.HERE / "ordinary_stationarity.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dest = lab.HERE / "finite_dose_material.json"
    if cli.parse_args().write:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dest.read_bytes()), result)
    for key in (
        "predicted_from_joint_head_counts",
        "independent_beam_contrast",
        "held_minus_matched_ordinary",
    ):
        print(
            key, sum(not a <= 0 <= b for a, b in (r[key] for r in result["readings"]))
        )
