"""Entropy of the inclusive protected word-slot algebra on the seven-tick annuli.

Every occupied word is a classical lottery of sharp orthogonal writer records.
The station address, slot label and absence flag are retained. Once the external
preparation rails are traced out, the spectrum of this word memory is the complete
addressed word distribution. For a fixed local placement, observing zero or one
word carries an arm-independent entropy. Given two arrivals, only a common tag
changes the word distribution. Consequently the full regional entropy difference
is the probability of a correlated ordered placement times the conditional
two-word entropy difference. Given P_s=product+sD on a doubly admitted pair,
its second entropy derivative is -sum(D²/P_s)<0. The product has zero first
derivative by the fixed marginal laws. Mixing the shared and fresh endpoints
separately into the reversed reference gives s=2t-1 and t-1. For 0<t<=1/2,
both parameters are nonpositive and the first lies closer to the maximum. On
[1/2,1] their entropy difference decreases strictly; at t=1 it is negative.
This certifies a unique positive crossing in the recorded interval. The result
is specific to the addressed protected word registers.
"""

import argparse
from fractions import Fraction
import json
from math import inf, nextafter
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
import numpy as np

from model.executable_model import federation
from ... import shared_variance, volume


HERE = Path(__file__).resolve().parent
IV = MPIntervalContext()
IV.dps = 50


def probabilities():
    _, _, _, src = volume.radar_source_and_hand(federation(0))
    laws, _ = shared_variance.source_words(src)
    cells = np.rint(laws * 2592).astype(int)
    np.testing.assert_array_equal(laws, cells / 2592)
    reverse, product, shared = cells
    delta = product - reverse
    np.testing.assert_array_equal(delta, shared - product)
    np.testing.assert_array_equal(
        np.outer(product.sum(1), product.sum(0)), 2592 * product
    )
    if min(cells.ravel()) <= 0:
        raise ValueError("word family fails faithfulness")
    labels = np.stack([volume.label(0, 1), volume.label(2, 1), volume.label(2, -1)])
    np.testing.assert_allclose(labels @ labels.conj().T, np.eye(3), atol=1e-14)
    np.testing.assert_array_equal(
        labels[:, 0], np.zeros(3)
    )  # blank has its own eigenvector
    # Group repeated cells of the word table before evaluating logarithms.
    pairs, counts = np.unique(
        np.stack((product.ravel(), delta.ravel()), 1), axis=0, return_counts=True
    )
    return [(int(p), int(d), int(c)) for (p, d), c in zip(pairs, counts, strict=True)]


def interval(x):
    q = Fraction(x)
    return IV.mpf(q.numerator) / q.denominator


def endpoints(x):
    return [nextafter(float(x.a), -inf), nextafter(float(x.b), inf)]


def entropy_and_slope(groups, correlation):
    """H(P) and dH/ds for P_s=product+sD, with D carrying zero total weight."""
    h = slope = IV.mpf(0)
    for p, d, count in groups:
        w = (p + d * interval(correlation)) / 2592
        h -= count * w * IV.ln(w)
        slope -= count * interval(Fraction(d, 2592)) * IV.ln(w)
    return h, slope


def contrast(groups, t):
    # Mix the two independently prepared endpoints into the same reversed tag.
    return entropy_and_slope(groups, 2 * t - 1)[0] - entropy_and_slope(groups, t - 1)[0]


def bracket_crossing(groups):
    left, right = Fraction(1, 2), Fraction(1)
    # On this interval H'(-1+2t) <= 0 and H'(-1+t) >= 0, and one is strict.
    for _ in range(36):
        mid = (left + right) / 2
        change = contrast(groups, mid)
        if change.a > 0:
            left = mid
        elif change.b < 0:
            right = mid
        else:
            raise ArithmeticError("finite entropy sign unresolved at test point")
    return [str(left), str(right)]


def run():
    groups = probabilities()
    reversible, ds = entropy_and_slope(groups, -1)
    product, product_slope = entropy_and_slope(groups, 0)
    shared, _ = entropy_and_slope(groups, 1)
    g = Fraction(729, 256)
    a = Fraction(729, 1280)
    table = json.loads((HERE / "results.json").read_text())
    ancestor = json.loads((HERE / "ancestor_limit.json").read_text())
    ancestry = {x["level"]: x for x in ancestor["midpoint_surveys"]}
    finite = []
    for row in table["midpoint_surveys"]:
        level = row["level"]
        gp = Fraction(row["normalized_pair"])
        ap = Fraction(ancestry[level]["normalized_pair"])
        finite.append(
            dict(
                level=level,
                global_entropy_slope_per_fourth=endpoints(ds * interval(gp)),
                ancestor_entropy_slope_per_fourth=endpoints(ds * interval(ap)),
            )
        )
    return dict(
        domain="Entropy of all addressed protected word slots in the region including blank, single and double arrivals; other registers traced out.",
        conditional_entropies_nats={
            str(key): endpoints(value)
            for key, value in [(-1, reversible), (0, product), (1, shared)]
        },
        reversed_reference_conditional_derivative=endpoints(ds),
        product_reference_conditional_derivative=endpoints(product_slope),
        global_entropy_derivative_per_fourth_limit=endpoints(ds * interval(g)),
        ancestor_entropy_derivative_per_fourth_limit=endpoints(ds * interval(a)),
        global_finite_endpoint_entropy_difference_per_fourth_limit=endpoints(
            (shared - product) * interval(g)
        ),
        finite_endpoint_conditional_entropy_difference=endpoints(shared - product),
        zero_of_shared_minus_fresh_entropy_at_common_admixture=bracket_crossing(groups),
        finite_resolutions=finite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    text = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "word_entropy.json").write_text(text)
    else:
        print(text, end="")
