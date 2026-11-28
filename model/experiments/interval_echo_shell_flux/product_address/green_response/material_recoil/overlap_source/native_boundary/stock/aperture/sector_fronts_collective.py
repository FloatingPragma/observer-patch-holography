"""Test coherent addition of distant delay records on a complete native peer register."""

import argparse
import json
from fractions import Fraction as Q
from hashlib import sha256

import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import splu

from . import (
    peer_population,
    peer_register_lab,
    sector_ball_source,
    sector_fronts,
    sector_fronts_counts,
    sector_sampling,
)

HERE = sector_fronts.HERE


def report_table(subset, mats, values, divisor):
    """Hitting probabilities at each initial word with an integer residual bound."""
    step = (mats[2] @ mats[1] @ mats[0]).T.tocsr()
    den = 128**3
    survive = np.full(step.shape[0], divisor)
    survive[subset] -= 1
    system = eye(len(survive)) - step @ diags(survive / (divisor * den))
    solver = splu(system.tocsc())
    cap = np.ceil(2 * solver.solve(np.ones(len(survive)))).astype(np.int64)
    assert (
        min(cap) >= 1
        and min(divisor * den * cap - step @ (survive * cap)) >= divisor * den
    )
    targets = np.zeros((len(survive), values.shape[1]), dtype=np.int64)
    targets[subset] = values
    rhs = step @ targets
    grid = min(2**37, 2**61 // (divisor * den))
    h = np.rint(grid * solver.solve(rhs / (divisor * den))).astype(np.int64)
    assert max(int(abs(h).max()) * divisor * den, grid * den) < 2**62
    residual = h * (divisor * den) - step @ (survive[:, None] * h) - rhs * grid
    eps = Q(int(cap.max()) * int(abs(residual).max()), grid * divisor * den)
    assert np.all((h[subset] >= 0) & (h[subset] <= grid))
    return h[subset], grid, eps


def conditional_sums(values, grid, categories, counts, gamma, size):
    # Conditional on type counts, each free member is an independent
    # uniform word within its type. Protected words have exact quotas.
    protected = categories < 0
    groups = []
    for i in range(len(counts)):
        block = values[categories == i]
        groups.append([Q(int(z), grid * len(block)) for z in block.sum(0)])
    p = [[Q(int(z), grid) for z in v] for v in values[protected]]
    x, y, j = np.array(groups + p, dtype=object).T
    weights = np.r_[counts, np.full(sum(protected), size // len(values), dtype=object)]
    assert sum(weights) == size
    mx, my = weights @ x / size, weights @ y / size
    single = weights @ j / size - mx * my
    independent_given_counts = weights @ (j - x * y)
    k = len(counts)
    total = (independent_given_counts + x[:k] @ gamma @ y[:k]) / size
    return mx, my, single, total


def acquire(subset, full, mats, decoded, blank_time, mean_xy, base):
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(subset))
    tables = peer_register_lab.acquisition.returns.ticket_tables(mats)
    rng = np.random.default_rng(248631)
    result = {}

    def observe(members, ordinary):
        groups, size = members.shape
        sums = []
        for chunk in np.array_split(members, 32):
            end, _ = sector_sampling.report_from(
                tables, subset, native, rng, 8, chunk.ravel()
            )
            sums.append(decoded[end].reshape(len(chunk), size, 2).sum(1))
        xy = np.vstack(sums)
        used = size * (len(ordinary) // size)
        end, _ = sector_sampling.report_from(
            tables, subset, native, rng, 8, ordinary[:used]
        )
        control = decoded[end].reshape(-1, size, 2).sum(1)

        def measured(a):
            return np.cov(a.T, ddof=1)[0, 1] / size

        def covariance_se(a):
            return np.prod(a - a.mean(0), axis=1).std(ddof=1) / (size * np.sqrt(len(a)))

        measured_difference = measured(xy) - measured(control)
        # Center on the finite predicted means. The complete register is the
        # sampling unit. Report a diagnostic standard error across registers.
        score = np.prod(xy - size * np.array(mean_xy), axis=1) / size

        def native_times(a):
            unique, freq = np.unique(
                size * blank_time + 2 * a, axis=0, return_counts=True
            )
            return np.column_stack((unique, freq)).tolist()

        result.update(
            registers=groups,
            members=size,
            trials=groups * size,
            observed_excess_per_member=float(score.mean() - base),
            sampling_standard_error=float(score.std(ddof=1) / np.sqrt(groups)),
            source_blind_sample_covariance_excess=float(
                np.cov(xy.T, ddof=1)[0, 1] / size - base
            ),
            ordinary_registers=len(control),
            ordinary_unused_single_returns=len(ordinary) - used,
            both_arms_measured_covariance_excess=float(measured_difference),
            both_arms_sampling_standard_error=float(
                np.hypot(covariance_se(xy), covariance_se(control))
            ),
            measured_sum_delay_covariance_excess_per_member=float(
                4 * measured_difference
            ),
            peer_summed_times_and_multiplicity=native_times(xy),
            ordinary_summed_times_and_multiplicity=native_times(control),
        )

    peer_register_lab.run(observe=observe, pulse_seed=728597, register_seed=578322)
    return result


def run():
    words, *_ = sector_ball_source.thermal.apparatus()
    _, _, subset, full, mats = sector_ball_source.sources((8,))
    bank = []
    info, _ = peer_population.histogram(16, draw_bank=bank)
    size = info["size"]
    cat = peer_population.fibers.apparatus()[3]
    x, y = (words[:, [0, 3]] > 0).astype(int).T
    flags = np.column_stack((x, y, x * y))
    ordinary = Q(int(x @ y), len(x)) - Q(int(x.sum()) * int(y.sum()), len(x) ** 2)
    report, grid, eps = report_table(subset, mats, flags, 8)
    np.testing.assert_allclose(
        report.mean(0) / grid, flags.mean(0), atol=float(eps) + 1e-12, rtol=0
    )
    counts, gamma = sector_fronts_counts.moments(bank)
    outputs = {}
    for key, tab, d, err in (
        ("at_preparation", flags, 1, Q(0)),
        ("one_eighth", report, grid, eps),
    ):
        mx, my, single, total = conditional_sums(tab, d, cat, counts, gamma, size)
        sg, agg = single - ordinary, total - ordinary
        # Products for two distinct members and products of marginal means
        # incur at most 2 eps each as both tables have entries in [0,1].
        error = (4 * size + 1) * err
        cross = (agg - sg) / (size - 1)
        outputs[key] = {
            "individual_covariance_gain": float(sg),
            "whole_covariance_gain_per_member": float(agg),
            "per_distinct_pair_covariance": float(cross),
            "conditional_report_probability_error": float(err),
            "whole_per_member_error": float(error),
            "finite_table_covariance_gain_per_member_rational": str(agg),
            "mean_home_occupancies": [float(mx), float(my)],
            "all_subgroup_sizes_per_member_upper": float(max(sg, agg) + error),
        }
    bins = json.loads((HERE / "sector_passage.json").read_bytes())[
        "spatial_modular_source"
    ]["all_rates"]["connected_certificate"]
    continuous = sector_fronts_counts.cadence(
        subset, mats, flags, cat, counts, gamma, size, bins
    )
    theory = outputs["one_eighth"]
    assert continuous["minimum_collective_gain_per_positive_rate"] > 0
    assert continuous["minimum_distinct_member_gain_per_squared_rate"] > 0
    independent = {}
    for divisor in (1, 64):
        table, denom, resid = report_table(subset, mats, flags, divisor)
        gain = float(
            conditional_sums(table, denom, cat, counts, gamma, size)[3] - ordinary
        )
        allowance = float((4 * size + 1) * resid)
        for row in continuous["bins"]:
            if row["rate"][0] <= 1 / divisor <= row["rate"][1]:
                assert (
                    row["per_rate"][0] / divisor - allowance
                    <= gain
                    <= row["per_rate"][1] / divisor + allowance
                )
        independent[str(divisor)] = {"gain_per_member": gain, "error_upper": allowance}
    for row in continuous["bins"]:
        if row["rate"][0] <= 1 / 8 <= row["rate"][1]:
            assert (
                row["per_rate"][0] / 8
                <= theory["whole_covariance_gain_per_member"]
                <= row["per_rate"][1] / 8
            )

    decoded, paths = sector_fronts.lab.two_delay_channels(sector_fronts.chart)
    obs = acquire(
        subset,
        full,
        mats,
        decoded,
        paths["blank_duration"],
        theory["mean_home_occupancies"],
        float(ordinary),
    )
    return {
        "inputs_sha256": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in ("SECTOR_FRONTS.md", "peer_population.json", "sector_passage.json")
        },
        "size": size,
        "ordinary_covariance": str(ordinary),
        "accepted_families": len(bank),
        "forecasts": outputs,
        "continuous_cadence": continuous,
        "independent_rate_checks": independent,
        "arbitrary_member_read_covariance_gain_upper": float(
            Q(1, 4) - ordinary + sum(abs(v) for v in gamma.ravel()) / (4 * size)
        ),
        "native_acquisition": obs,
        "exterior_delay_channels": paths,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    out = run()
    path = HERE / "sector_fronts_collective.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        sector_fronts.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_text()), out
        )
    print(
        [
            (
                k,
                v["individual_covariance_gain"],
                v["whole_covariance_gain_per_member"],
                v["whole_per_member_error"],
            )
            for k, v in out["forecasts"].items()
        ]
    )
