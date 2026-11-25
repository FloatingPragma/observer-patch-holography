"""Check every archived interval against direct tails and rare tails by elementary recurrence."""

import json
import math

import numpy as np
from scipy.stats import binom

from model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair import (
    measure as moving,
)
from model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair.counting import (
    binomial_limits,
)
from model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair.recoil import (
    measure as recoil,
)


def check_all(counts, n, tail):
    k = np.unique(counts)
    lo, hi = binomial_limits(k, n, tail)
    assert np.all((lo <= k / n) & (k / n <= hi))
    nz = k > 0
    interior = k < n
    assert np.all(binom.sf(k[nz] - 1, n, lo[nz]) <= tail)
    assert np.all(binom.cdf(k[interior], n, hi[interior]) <= tail)
    # The inward adjacent float brackets the forward-tail root, to numeric
    # precision. A very small absolute allowance covers the 64-step bracket.
    assert np.all(
        binom.sf(k[nz] - 1, n, np.minimum(1, np.nextafter(lo[nz], 1) + 1e-21))
        >= tail * (1 - 1e-9)
    )
    assert np.all(
        binom.cdf(k[interior], n, np.maximum(0, np.nextafter(hi[interior], 0) - 1e-21))
        >= tail * (1 - 1e-9)
    )


def test_all_used_bins():
    archive = json.loads((moving.HERE / "records.json").read_text())
    counts = [np.array(archive["isolated_bond_orders"])[:, :, :2].ravel()]
    for arm in archive["settings"].values():
        for scan in arm["scans"]:
            counts.extend(
                np.array(scan[field])[:, :, :2].ravel()
                for field in ("orders", "erased_orders")
            )
    check_all(np.concatenate(counts), moving.SHOTS, moving.BIN_TAIL)
    counts = []
    for scan in archive["independent_followup"]["scans"]:
        counts.extend(
            np.array(scan[field])[:, :, :2].ravel() for field in ("orders", "bench")
        )
    check_all(np.concatenate(counts), 2**38, moving.BIN_TAIL)
    archive = json.loads((recoil.HERE / "records.json").read_text())
    counts = []
    for row in archive["settings"].values():
        counts.extend(row["parity"] + row["species_selective_current"])
        for v in row["current"].values():
            counts.extend(v)
        for v in row["flight"].values():
            for bank in v:
                counts.extend(bank)
    for v in archive["phase_dimer_transfer_stay"].values():
        for bank in v:
            counts.extend(bank)
    check_all(counts, recoil.SHOTS, recoil.TAIL)


def elementary_tail(k, n, p, up):
    log_pmf = (
        k * math.log(n * p)
        + sum(math.log1p(-j / n) for j in range(k))
        - math.lgamma(k + 1)
        + (n - k) * math.log1p(-p)
    )
    term = math.exp(log_pmf)
    terms = [term]
    j = k
    while term > terms[0] * 1e-18:
        if up:
            term *= (n - j) / (j + 1) * p / (1 - p)
            j += 1
        else:
            if j == 0:
                break
            term *= j / (n - j + 1) * (1 - p) / p
            j -= 1
        terms.append(term)
    return math.fsum(terms)


def test_rare_endpoint_against_pmf_ratios():
    for n, k, t in ((2**28, 1000, 1e-7), (2**32, 999, 1e-15), (2**38, 5000, 1e-15)):
        lower, upper = binomial_limits([k], n, t)
        np.testing.assert_allclose(
            [elementary_tail(k, n, lower[0], True), elementary_tail(k, n, upper[0], False)],
            [t, t],
            rtol=2e-9,
            atol=0,
        )
