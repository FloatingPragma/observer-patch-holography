"""Full entropy differences and reference selection on the protected routed regions."""

import json
from fractions import Fraction

import numpy as np
from mpmath import iv

from model.executable_model import federation
from model.experiments.interval_echo_einstein_station import shared_variance, volume
from model.experiments.interval_echo_einstein_station.collar_scaling import (
    measure as clock,
)


def table():
    *_, source = volume.radar_source_and_hand(federation(0))
    law, _ = shared_variance.source_words(source)
    integer = np.rint(2592 * law).astype(int)
    np.testing.assert_array_equal(law, integer / 2592)
    reverse, fresh, shared = integer
    assert np.all(reverse > 0) and np.all(shared > 0)
    np.testing.assert_array_equal(shared - fresh, fresh - reverse)
    return reverse, fresh, shared


def receipt():
    iv.dps = 35
    reverse, fresh, shared = table()
    delta = shared - fresh

    def ari(z, r):
        s = iv.mpf(0)
        for a, h in zip(fresh.flat, delta.flat):
            u = iv.mpf(int(a)) + z * int(h)
            v = iv.mpf(int(a)) + r * int(h)
            s += u * iv.ln(u / v) / 2592
        return s

    def slope(z):
        return -sum(
            iv.mpf(int(h)) * iv.ln((iv.mpf(int(a)) + z * int(h)) / int(a)) / 2592
            for a, h in zip(fresh.flat, delta.flat)
        )

    def fisher(z):
        return sum(
            Fraction(int(h) ** 2, 2592 * (int(a) + z * int(h)))
            for a, h in zip(fresh.flat, delta.flat)
        )

    # Convexity of J supplies its maximum at an endpoint. The per-bin
    # lower estimate is uniform even without separating a minimizing root.
    lower = sum(
        Fraction(int(h) ** 2, 2592 * int(max(a, b)))
        for h, a, b in zip(delta.flat, reverse.flat, shared.flat)
    )
    upper = max(fisher(-1), fisher(1))
    assert 0 < lower <= fisher(0) <= upper
    readings = json.loads((clock.HERE / "results.json").read_text())
    rows = []
    for panel in readings["panels"]:
        geo = panel["geometry"]
        arm = panel["settings"][0]
        n = geo["annulus_counts"][0][0]
        N = geo["cells"]
        q = Fraction(n * (n - 1), N * (N - 1))
        ell = arm["relative_radius"]
        qscale = float(q) / ell**4
        # A likelihood read gives a linear, unbiased estimate of D(P_z||P0)
        # on that actual prepared arm. All misses score zero.
        score = []
        for target, counts in zip(
            (reverse, fresh, shared), arm["arms"]["separated"]["reads"]
        ):
            pair = np.array(counts["00"], dtype=np.int64)
            scores = np.log(target.ravel() / fresh.ravel())
            probabilities = float(q) * target.ravel() / 2592
            mean = probabilities @ scores
            variance = probabilities @ scores**2 - mean**2
            measured = float(pair @ scores) / readings["shots"]
            score.append(
                dict(
                    pair_admissions=int(pair.sum()),
                    mean_nats=mean,
                    measured_nats=measured,
                    standard_error=float(np.sqrt(variance / readings["shots"])),
                )
            )
        rows.append(
            dict(
                stations=N,
                relative_radius=ell,
                joint_admission=str(q),
                squared_density_factor_per_radius_fourth=qscale,
                information_cost_reverse_fresh_shared=score,
            )
        )
    return dict(
        fisher_lower=str(lower),
        fisher_upper=str(upper),
        fisher_at_product=str(fisher(0)),
        entropy_slopes=[str(slope(z)) for z in (-1, 0, 1)],
        conditional_information_from_product=dict(
            reverse=str(ari(-1, 0)), shared=str(ari(1, 0))
        ),
        reversed_to_shared=dict(
            entropy_change=str(ari(-1, 0) - ari(1, 0)),
            linear_entropy_prediction=str(2 * slope(-1)),
            relative_entropy_remainder=str(ari(1, -1)),
        ),
        resolved_annuli=rows,
    )
