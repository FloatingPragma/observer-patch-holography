"""Read a held six-marker surface after native writing on whole peer registers."""

import math
import numpy as np
from scipy.stats import beta

from . import peer_register_lab
from . import sector_sampling


def read(tables, subset, native, words, material_meter):
    # Reload the measured material-to-boundary instrument. Its inputs are
    # the actual retained home readings of these newly acquired carriers.
    rng = np.random.default_rng(519441)
    _, group = np.unique(words[:, :4], axis=0, return_inverse=True)
    surface = np.array(material_meter)[group]
    low, hi = int(surface.min()), int(surface.max())
    assert (low, hi) == (10, 14) and set(surface) == {low, hi}
    short = (surface == low).astype(int)
    n = 128
    draws = {}
    divisors = (8, 64)

    def observe(members, reference):
        b = len(members)
        for div in divisors:
            selected = members[
                np.arange(b)[:, None], rng.integers(members.shape[1], size=(b, n))
            ].ravel()
            ends, wait = sector_sampling.report_from(
                tables, subset, native, rng, div, selected
            )
            ref, refwait = sector_sampling.report_from(
                tables, subset, native, rng, div, reference[: b * n]
            )
            draws.setdefault(div, []).append(
                (
                    short[ends].reshape(b, n).sum(1),
                    int(short[ref].sum()),
                    max(wait, refwait),
                )
            )

    def assess(banks):
        rows = np.concatenate([entry[0] for entry in banks])
        b = len(rows)
        total = b * n
        x = rows / n
        ell = math.log(4 / 0.002)
        err = math.sqrt(2 * x.var(ddof=1) * ell / b) + 7 * ell / (3 * (b - 1))
        small = sum(entry[1] for entry in banks)
        ordinary = np.array(
            [
                beta.ppf(0.001, small, total - small + 1),
                beta.isf(0.001, small + 1, total - small),
            ]
        )
        prob = x.mean() + np.array([-err, err]) - ordinary[::-1]
        return dict(
            prepared_small_count=int(rows.sum()),
            ordinary_small_count=small,
            trials_per_arm=total,
            registers=b,
            sampled_members_per_register=n,
            maximum_sweeps=max(entry[2] for entry in banks),
            ordinary_probability_band=ordinary.tolist(),
            prepared_halfwidth=err,
            mean_change_in_links=(low - hi) * (x.mean() - small / total),
            interval_links=((low - hi) * prob[::-1]).tolist(),
        )

    archive = peer_register_lab.run(
        observe=observe, pulse_seed=591795, register_seed=808562
    )
    pilot = {str(div): assess(bank) for div, bank in draws.items()}
    draws = {}
    divisors = (64,)
    rng = np.random.default_rng(506475)
    members = 0
    for pulse, register in ((393411, 923841), (146701, 730405)):
        new = peer_register_lab.run(
            observe=observe, pulse_seed=pulse, register_seed=register
        )
        members += new["members"]
    return dict(
        physical_prepared_members=archive["members"],
        joint_coverage_at_least=0.992,
        rates=pilot,
        sparse_confirmation=dict(
            physical_prepared_members=members,
            coverage_at_least=0.996,
            measurement=assess(draws[64]),
        ),
    )
