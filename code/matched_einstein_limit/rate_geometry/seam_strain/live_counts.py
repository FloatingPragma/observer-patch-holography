"""Independent synthetic clock-count banks for the primary coupled five-station survey.

Each edge trial enrolls a charge and orientation equiprobably, pumps a fresh
system and records its first return tick or one of the two timeouts. No source
derivative is treated as a probability distribution. The two acquisition banks
are independent. Both the derivatives of the geometric functions at reference
and the source-induced changes in edge means are inferred from their own counts.

The archive freezes one NumPy multinomial realization. The sequential binomial
sampler calls platform math functions, so a seed need not generate identical
integer counts on different systems. Read the archived counts with infer() to
reproduce the intervals exactly. Tests independently recompute the physical
propensities and check a fresh sample against its sampling band.
"""

import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

from mpmath import iv
import numpy as np

from . import live_bounds, live_radar, live_source, live_survey
from .source import measure as source

HERE = Path(__file__).parent
STATIONS = (0, 1, 38, 39, 5)
EDGES = tuple(combinations(range(5), 2))
ARMS = ("ordinary", "decorrelated", "hidden")
PREPS = (
    "reference",
    "ordinary",
    "decorrelated",
    "hidden_plus",
    "hidden_minus",
    "left_plus",
    "left_minus",
    "right_plus",
    "right_minus",
)
BANKS = ("volume", "boundary")
OUTCOMES = tuple(f"return_{t}" for t in range(2, 13)) + ("relay_only", "no_relay")
TICKS = tuple(range(2, 13)) + (12, 12)
N = 1 << 68
CHUNK = 1 << 51  # Integers in each multinomial are represented exactly in binary64.
DEN = 1 << 43
SEED = 867421
MEAN_ALLOWANCE = Fraction(1, 10**11)
SETTINGS = len(BANKS) * len(PREPS) * len(EDGES)


def preparations():
    """Nine normalized, independently preparable source states after the same pump."""
    board = source.formation.federation(1)
    heads, _ = source.form(board, source.formation.region(board), -1)
    identity = np.eye(4)
    result = {k: live_source.postpump(True)[k] for k in PREPS[:3]}
    # Uniform XZ/ZX setting and uniform first sign; prepare the same/opposite
    # second sign. The average of the two projectors at a setting is (I±P)/4.
    x, _, z = source.formation.PAULI
    parity = (np.kron(x, z) + np.kron(z, x)) / 2
    for sign, suffix in ((1, "plus"), (-1, "minus")):
        result[f"hidden_{suffix}"] = live_source.propagate_heads(
            (identity + sign * parity) / 4
        )
        for side, g in live_survey.rotations(heads).items():
            u = (identity - sign * 2j * g) / np.sqrt(2)
            result[f"{side}_{suffix}"] = live_source.propagate_heads(
                u @ heads @ u.conj().T
            )
    return result


def native_probabilities():
    """POVM bins after marginalizing the independent fair charge and direction enrollment."""
    pooled = np.zeros((10, len(OUTCOMES), 8, 8), complex)
    for q in (-1, 1):
        outward, returns = live_radar.directed_events(1, q, True, STATIONS)
        for e, (i, j) in enumerate(EDGES):
            for a, b in ((i, j), (j, i)):
                relay = outward[a, b].sum(axis=0)
                done = returns[a, b, 2:]
                pooled[e, :-2] += done / 4
                pooled[e, -2] += (relay - done.sum(axis=0)) / 4
                pooled[e, -1] += (np.eye(8) - relay) / 4
    return {
        name: np.einsum("etij,ji->et", pooled, rho).real
        for name, rho in preparations().items()
    }


def freeze(probabilities):
    """Dyadic propensities allow seeded replay independent of the numerical BLAS."""
    table = {}
    for name in PREPS:
        p = probabilities[name]
        if p.min() < -1e-13 or np.max(abs(p.sum(axis=1) - 1)) > 1e-12:
            raise ValueError("invalid native probability row")
        nums = np.rint(p * DEN).astype(np.int64)
        nums[:, -1] = DEN - nums[:, :-1].sum(axis=1)
        if nums.min() < 0:
            raise ValueError("invalid frozen propensities")
        table[name] = nums.tolist()
    return table


def draw_row(propensities, b, a, e, seed=SEED):
    """Draw and sum one independently seeded edge/preparation/bank lottery."""
    if (
        len(propensities) != len(OUTCOMES)
        or sum(propensities) != DEN
        or min(propensities) < 0
    ):
        raise ValueError("invalid frozen probability row")
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, b, a, e])))
    probs = np.asarray(propensities) / DEN
    # Accumulate in blocks whose totals fit int64, then in Python ints.
    total = [0] * len(OUTCOMES)
    group = 2048
    assert N % (CHUNK * group) == 0 and CHUNK * group < 2**63
    for _ in range(N // (CHUNK * group)):
        block = rng.multinomial(CHUNK, probs, size=group).sum(axis=0).tolist()
        total = [x + y for x, y in zip(total, block, strict=True)]
    return total


def draw(table):
    """Independently seeded multinomial samples, stored as integer category counts."""
    counts = {}
    for b, bank in enumerate(BANKS):
        rows = counts[bank] = {}
        for a, name in enumerate(PREPS):
            rows[name] = []
            for e, propensities in enumerate(table[name]):
                rows[name].append(draw_row(propensities, b, a, e))
    return counts


def score(row):
    if (
        len(row) != len(TICKS)
        or sum(row) != N
        or any(type(n) is not int or n < 0 for n in row)
    ):
        raise ValueError(
            "all registered trials must be assigned to exactly one clock bin"
        )
    return Fraction(sum(t * n for t, n in zip(TICKS, row, strict=True)), 4 * N)


def slope(name, means, e):
    if name in ("ordinary", "decorrelated"):
        return means[name][e] - means["reference"][e]
    divisor = iv.sqrt(2) if name == "hidden" else 2
    return (means[f"{name}_plus"][e] - means[f"{name}_minus"][e]) / divisor


def geometry(bank, radius):
    keys = live_bounds.KEYS
    means = {
        name: [
            live_bounds.rational_interval(score(row)) + radius * iv.mpf([-1, 1])
            for row in bank[name]
        ]
        for name in PREPS
    }
    n = len(EDGES)
    if any(len(rows) != n for rows in means.values()):
        raise ValueError("all ten edge settings are required")
    lengths = means["reference"]
    # Enclose the geometric gradient itself from this bank's reference lengths.
    # Each pass differentiates five independent coordinates through the exact
    # Euclidean face, filling-volume, and solid-angle functions.
    grads = [[], [], []]
    for block in (0, 1):
        axes = {
            k: [iv.mpf(int(e == 5 * block + i)) for e in range(n)]
            for i, k in enumerate(keys)
        }
        jets = live_bounds.geometry_jets(lengths, axes)
        for g, jet in zip(grads, jets, strict=True):
            g.extend(jet.g)
    d = {name: [slope(name, means, e) for e in range(n)] for name in keys}
    derivatives = {
        name: [sum(c * h for c, h in zip(g, delta, strict=True)) for g in grads]
        for name, delta in d.items()
    }
    return [jet.v for jet in jets], derivatives


def infer(counts):
    """Distribution-free simultaneous 99 percent intervals from the two count banks."""
    saved = iv.dps
    try:
        iv.dps = 35
        ep = live_bounds.endpoints
        # Range t/4 lies in [1/2,3], including both timeouts. A union bound over
        # all 180 means requires no independence between the derived slopes.
        half = iv.mpf(5) / 2 * iv.sqrt(iv.ln(2 * SETTINGS * 100) / (2 * N))
        radius = half + live_bounds.rational_interval(MEAN_ALLOWANCE)
        vb, v = geometry(counts["volume"], radius)
        ab, a = geometry(counts["boundary"], radius)
        if v["right"][1].a <= 0:
            raise ValueError("right rotation is not resolved as a volume adjuster")
        left_volume = v["left"][1]
        if left_volume.a <= 0 <= left_volume.b:
            raise ValueError("left rotation is not resolved as a volume adjuster")
        matched = {}
        areas = {}
        for arm in ARMS:
            turn = -v[arm][1] / v["right"][1]
            left_turn = -v[arm][1] / left_volume
            areas[arm] = a[arm][0] + turn * a["right"][0]
            matched[arm] = dict(
                left_control_fixed_volume_area=ep(a[arm][0] + left_turn * a["left"][0]),
                right_rotation=ep(turn),
                fixed_volume_boundary_area=ep(areas[arm]),
                fixed_volume_cone_area_per_radius_squared=ep(a[arm][2] / 3),
                rotated_cone_area_per_radius_squared=ep(
                    (a[arm][2] + turn * a["right"][2]) / 3
                ),
            )
        return dict(
            simultaneous_coverage="at least 99% for the synthetic probability table under independent sampling",
            settings=SETTINGS,
            hoeffding_half_width=ep(half),
            total_mean_half_width=ep(radius),
            geometry_order=["outer_area", "interior_volume", "solid_angle"],
            volume_bank_reference=list(map(ep, vb)),
            boundary_bank_reference=list(map(ep, ab)),
            volume_bank_slopes={name: list(map(ep, ds)) for name, ds in v.items()},
            boundary_bank_slopes={name: list(map(ep, ds)) for name, ds in a.items()},
            matched=matched,
            right_matched_area_over_ordinary={
                arm: ep(areas[arm] / areas["ordinary"]) for arm in ARMS[1:]
            },
        )
    finally:
        iv.dps = saved


def run():
    table = freeze(native_probabilities())
    counts = draw(table)
    return dict(
        protocol="synthetic independent first-return clock counts; coupled 80-cell source, four-tick vacancy pump, retained strain and original heads",
        stations=list(STATIONS),
        deadline=12,
        charge_and_direction="independent uniform per trial; four charge/direction enrollments pooled before multinomial sampling",
        edges=[[STATIONS[a], STATIONS[b]] for a, b in EDGES],
        outcomes=list(OUTCOMES),
        score_times_four=list(TICKS),
        shots_per_setting=N,
        total_shots=SETTINGS * N,
        seed=SEED,
        sampler="PCG64 and NumPy multinomial on dyadic propensities; the integer archive fixes the sample because platform math libraries can differ. Inference replays from the counts alone.",
        probability_denominator=DEN,
        probability_numerators=table,
        counts=counts,
        method="Volume bank determines the right rotation. Independently acquired boundary bank supplies outer area and apex solid angle. Each bank estimates its geometric gradient at its own reference. The simple cone changes only radius to hold cone volume fixed; the rotated cone additionally applies the outer-volume rotation.",
        hidden_scope="The measured contrast is the normalized nominal input (XZ+ZX)/(4 sqrt(2)) from positive same/opposite-parity product lotteries; the exact kernel direction has a separate deterministic enclosure.",
        mean_allowance=float(MEAN_ALLOWANCE),
        inference=infer(counts),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "live_counts.json").write_text(data)
    else:
        print(data, end="")
