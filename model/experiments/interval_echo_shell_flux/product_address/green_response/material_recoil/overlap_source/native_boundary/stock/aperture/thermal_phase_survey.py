"""Two-source-contact phase receipts and material ranges on the same thermal copies."""

import argparse
from hashlib import sha256
from itertools import product
import json

import numpy as np
from scipy import sparse

from . import phase_receipts as phase
from . import thermal_radar as radar

HERE = radar.HERE


def thermal_flags(sectors):
    # The same conserved parity sets the thermal weights and the dark flag.
    # Average over it only after forming their joint probabilities.
    arms = []
    for f, g in product(range(2), repeat=2):
        joint = {}
        for a, b in product(radar.lab.STRENGTHS, repeat=2):
            flags = np.sin(phase.dwell.DWELL * np.array([a, b]) / 2) ** 2
            p = np.prod(np.where([f, g], flags, 1 - flags))
            label = f"{-a},{-b}"
            joint[label] = p * sectors[label]
        arms.append(radar.lab.motion.binary_responses(joint))
    np.testing.assert_allclose(
        sum(arms), radar.lab.motion.binary_responses(sectors), atol=1e-15
    )
    # After summing all thermal addresses, each sharp head gives the calibrated
    # one-use gate. Marginal positions alone need not certify these joint records.
    c = phase.written_receiver()
    for h in range(4):
        a, b = np.array([h // 2, h % 2]) * c
        np.testing.assert_allclose(
            [v[h].sum() for v in arms],
            [(1 - a) * (1 - b), (1 - a) * b, a * (1 - b), a * b],
        )
    return np.concatenate(arms, axis=1)


def additional_entropy(plain, complete):
    """Eliminate a separate ordinary-calibrated area gain at each volume."""
    from mpmath import iv

    def area(arm):
        return iv.mpf(plain[arm]["area"][-1]["from_counts"])

    def slope(bank, arm):
        return iv.mpf(bank[arm]["entropy"]["from_counts"])

    ao, ah = area("ordinary"), area("product")
    so, sh = slope(plain, "ordinary"), slope(plain, "product")
    fo, fh = slope(complete, "ordinary") - so, slope(complete, "product") - sh
    base = ao * sh - ah * so
    extra = ao * fh - ah * fo
    scale = radar.standard.quotient(
        radar.lab.interval_output(-base), radar.lab.interval_output(extra)
    )
    return dict(
        needed_held_entropy_at_fixed_ordinary_gain=radar.lab.interval_output(
            -sh + so * ah / ao
        ),
        required_coefficient_on_flag_entropy=scale,
        flag_entropy_product=radar.lab.interval_output(fh),
        flag_entropy_ordinary=radar.lab.interval_output(fo),
    )


def run():
    board = radar.lab.motion.cold.federation(0)
    native = json.loads((HERE / "pair_contact.json").read_bytes())
    edges = native["physical_seams"]
    live, hops, center = radar.chart(board, edges)
    i, j, fixed, weights, error = radar.thermal(
        board, edges, reduce_sectors=thermal_flags
    )
    rng = np.random.default_rng(6413985)
    counts = {
        k: radar.binomial_hist(rng, weights[tag])
        for k, tag in (
            ("reference", "native"),
            ("native", "native"),
            ("ordinary", "ordinary"),
            ("product", "product"),
        )
    }
    area_counts = {k: radar.binomial_hist(rng, v.sum(0)) for k, v in weights.items()}
    delta = radar.lab.probe.log_factor(native["source"])
    direct = 65 * radar.ranges(live, center, i, j)
    assisted = []
    orientation = []
    for edge in edges:
        a, b = sorted(edge, key=lambda n: (hops[center, n], n))
        orientation.append([a, b])
        exit_leg = 65 * radar.ranges(live, b, i, j)
        assisted.append(direct[:, a, None] + 67 + exit_leg)
    dist = np.concatenate(
        [
            direct,
            np.minimum(direct, assisted[1]),
            np.minimum(direct, assisted[0]),
            np.minimum(np.minimum(direct, assisted[0]), assisted[1]),
        ]
    )
    assert np.all(dist[:, center] == 0)
    u, v = sparse.triu(live, 1).nonzero()
    inside = u // 12 == v // 12
    ix, jx = np.tile(i, 4), np.tile(j, 4)
    output = {}
    for reverse in (False, True):
        asc = np.argsort(dist[:, ::-1] if reverse else dist, axis=1, kind="stable")
        order = board.n - 1 - asc if reverse else asc
        rows = {}
        for r in (10, 11, 12, int(hops[center].max())):
            volume = int(np.count_nonzero(hops[center] <= r))
            surface, ball, record, seen = radar.geometry(
                order, volume, ix, jx, edges, u, v, inside
            )
            tests = radar.read_meters(
                weights, ball, record, surface, counts, area_counts, error, delta
            )
            rows[str(volume)] = dict(
                fully_contains_both=bool(seen.all()),
                distinct_balls=len(np.unique(ball)),
                source_coverage={
                    k: float(np.sum(v * seen)) for k, v in weights.items()
                },
                tests=tests,
                comparison=radar.compare(tests),
            )
            if seen.all() and volume < board.n:
                # Both complete heads and both local flag records are physically
                # visible on every trial of these volumes.
                f = np.repeat(np.arange(4), len(i))
                _, full = np.unique(4 * record + f, return_inverse=True)
                extra = radar.read_meters(
                    weights,
                    ball,
                    full.reshape(record.shape),
                    surface,
                    counts,
                    area_counts,
                    error,
                    delta,
                )
                rows[str(volume)]["with_local_receipt_flags"] = dict(
                    entropy={
                        label: reading["entropy"] for label, reading in extra.items()
                    },
                    comparison=radar.compare(extra),
                    additional_entropy=additional_entropy(tests, extra),
                )
        output["descending" if reverse else "ascending"] = rows
    return dict(
        sha256={
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "THERMAL_PHASE_SURVEY.md",
                "thermal_radar.json",
                "phase_receipts.json",
            )
        },
        center=center,
        ordinary_clamp=fixed,
        orientations=orientation,
        per_arm=radar.N,
        count_sha256={
            k: sha256(v.astype("<i8").tobytes()).hexdigest() for k, v in counts.items()
        },
        geometry_count_sha256={
            k: sha256(v.astype("<i8").tobytes()).hexdigest()
            for k, v in area_counts.items()
        },
        tests=output,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    options = cli.parse_args()
    result = run()
    path = HERE / "thermal_phase_survey.json"
    if options.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        phase.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    for name, rows in result["tests"].items():
        print(
            name,
            [
                (r, v["fully_contains_both"], v["comparison"]["common_gain_excluded"])
                for r, v in rows.items()
            ],
        )
