"""Independent local clock spectra and two-leg first-click timing on a native belt."""

import argparse
import functools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from scipy.signal import fftconvolve
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls.belt import (
    belts,
)
from model.executable_model import federation
from model.experiments.interval_echo_einstein_station.native_tidal_response.measure import (
    local_calibration,
)

HERE = Path(__file__).resolve().parent
ARMS = ("sham", "coin_even", "coin_odd", "empty_write")
TAGS = (1, 2)
WAITS = np.arange(1, 97) / 16
DEADLINES = np.array([8, 16, 64, 256, 4096])
SHOTS = 2**28
COMPARISONS = 20000  # upper bound on all clock probabilities and radar CDF points
ALPHA = 0.01 / COMPARISONS
RADAR_EPS = float(np.sqrt(np.log(2 / ALPHA) / (2 * SHOTS)))


def geometry():
    (poles, route) = next(iter(belts().items()))
    car = federation(0)
    pole = next(p for p in poles if all(route[i] in car.neighbours[p] for i in (2, 8)))
    assert all(route[(i + 1) % 10] in car.neighbours[p] for i, p in enumerate(route))
    return {
        "perimeter": list(map(int, route)),
        "parked_unit": int(pole),
        "source_contact_sites": [int(route[i]) for i in (2, 8)],
    }


def calibration():
    """Certify both laws and both material tags in the full 99-dimensional contact."""
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    _, inserts = local_calibration()
    empty = np.zeros(law.DIM)
    empty[law.SRC.index((0, 0))] = 1
    memory = (law.echo() @ empty).reshape(-1, law.N_SRC)[:, law.SRC.index((0, 0))]
    blocks = {}
    for q in TAGS:
        ends = [law.SRC.index(word) for word in ((q, 0), (0, q))]
        blank = np.eye(law.DIM)[:, ends]
        written = np.kron(memory[:, None], np.eye(law.N_SRC)[:, ends])
        for k in (candidate.KAPPA, 0):
            h = exchange + k * contrast
            block = blank.T @ h @ blank
            np.testing.assert_allclose(h @ blank, blank @ block, atol=1e-13)
            np.testing.assert_allclose(h @ written, 0, atol=1e-13)
            blocks[k, q] = block
    return blocks, inserts


def generator(k, q, arm, blocks, inserts, bond=None):
    n = 10 if bond is None else 2
    h = np.zeros((n, n))
    sites = list(range(10)) if bond is None else [bond, (bond + 1) % 10]
    for i in range(10) if bond is None else [bond]:
        if arm == "empty_write" and i == 2:
            continue
        ends = [sites.index(i), sites.index((i + 1) % 10)]
        h[np.ix_(ends, ends)] += blocks[k, q]
    source = inserts["sham" if arm == "empty_write" else arm]

    def difference(term):
        return source[term][q] - source[term][0]

    v = difference("odd_reversal") + k * difference("contrast_squared")
    h += np.diag([v if p in (2, 8) else 0.0 for p in sites])
    return h - np.eye(n) * h.trace() / n


@functools.lru_cache(maxsize=1)
def frequency_grid():
    # Each grid cell encloses possible responses throughout its frequency interval,
    # including maxima, collapsed fringes and close aliases.
    edges = np.arange(5001) * 0.01
    mid = (edges[1:] + edges[:-1]) / 2
    s = np.sin(mid[:, None] * WAITS / 2) ** 2
    pad = np.diff(edges)[:, None] * WAITS / 4  # |d sin²(w t/2)/dw| <= t/2
    return edges, mid, s, np.maximum(0.0, s - pad), np.minimum(1.0, s + pad)


def clock_fit(counts):
    if not np.any(counts):
        return {
            "frequency": None,
            "max_probability_upper": float(1 - (ALPHA / 2) ** (1 / SHOTS)),
        }
    y = counts / SHOTS
    edges, mid, s, lower, upper = frequency_grid()
    p_lo = np.zeros_like(y)
    p_hi = np.ones_like(y)
    ok = counts > 0
    p_lo[ok] = beta.ppf(ALPHA / 2, counts[ok], SHOTS - counts[ok] + 1)
    ok = counts < SHOTS
    p_hi[ok] = beta.ppf(1 - ALPHA / 2, counts[ok] + 1, SHOTS - counts[ok])
    min_amp = (p_lo / np.maximum(upper, 1e-300)).max(axis=1)
    max_amp = (p_hi / np.maximum(lower, 1e-300)).min(axis=1)
    feasible = np.flatnonzero(min_amp <= np.minimum(1.0, max_amp))
    if not len(feasible):
        raise ValueError("A single exchange frequency fails the declared count bounds")

    def score(freq):
        wave = np.sin(freq * WAITS / 2) ** 2
        amp = np.clip(np.dot(y, wave) / np.dot(wave, wave), 0, 1)
        return np.sum((y - amp * wave) ** 2)

    a = np.clip((s[feasible] @ y) / (s[feasible] ** 2).sum(axis=1), 0, 1)
    errors = ((s[feasible] * a[:, None] - y) ** 2).sum(axis=1)
    best = feasible[errors.argmin()]
    freq = minimize_scalar(
        score, bounds=(max(1e-9, mid[best] - 0.03), mid[best] + 0.03), method="bounded"
    ).x
    wave = np.sin(freq * WAITS / 2) ** 2
    return {
        "frequency": float(freq),
        "fringe": float(np.clip(np.dot(y, wave) / np.dot(wave, wave), 0, 1)),
        "frequency_enclosure": [
            float(edges[feasible[0]]),
            float(edges[feasible[-1] + 1]),
        ],
    }


def read_clocks(k, q, arm, blocks, inserts, rng, receipts):
    records = []
    for i in range(10):
        h = generator(k, q, arm, blocks, inserts, i)
        e, u = eigh(h)
        amp = np.exp(-1j * WAITS[:, None] * e) @ (u[0] * u[1])
        counts = rng.binomial(SHOTS, np.clip(abs(amp) ** 2, 0, 1).round(13))
        receipts[f"clock_{k}_{arm}_{q}_{i}"] = counts
        records.append(clock_fit(counts))
    return records


def first_click(u, start, target, steps):
    """Projectively check just the destination after each native hold."""
    a = np.eye(len(u), dtype=complex)[:, start].copy()
    p = np.zeros(steps + 1)
    for j in range(1, steps + 1):
        a = u @ a
        p[j] = abs(a[target]) ** 2
        a[target] = 0.0
    np.testing.assert_allclose(p.sum() + np.vdot(a, a).real, 1.0, atol=5e-11, rtol=0)
    return p


def radar(h, cadence):
    e, v = eigh(h)
    u = (v * np.exp(-1j * cadence * e)) @ v.T
    steps = int(DEADLINES[-1] / cadence)
    out, back = (first_click(u, x, y, steps) for x, y in ((0, 5), (5, 0)))
    total = fftconvolve(out, back)[: steps + 1].cumsum()
    return np.maximum.accumulate(
        np.clip(total[(DEADLINES / cadence).astype(int)], 0, 1)
    )


def weighted_guide(clocks, reference):
    h = np.zeros((10, 10))
    for i, (meter, blank) in enumerate(zip(clocks, reference, strict=True)):
        if meter["frequency"] is None:
            return None
        rate = meter["frequency"] / blank["frequency"]
        ends = [i, (i + 1) % 10]
        h[np.ix_(ends, ends)] += rate * np.array([[1.0, -1.0], [-1.0, 1.0]]) / 2
    return h


def reference_clock(block, receipts):
    """Independently sample the undisturbed dimer at each appointment spacing."""
    rng = np.random.default_rng(271125)
    e, v = eigh(block)
    output = []
    for spacing in (1 / 16, 0.5, 1.0):
        p = abs(np.exp(-1j * e * spacing) @ (v[0] * v[1])) ** 2
        ticks = rng.binomial(SHOTS, round(float(p), 13))
        receipts[f"reference_{spacing}"] = np.array([ticks, SHOTS - ticks])
        times = 2 * np.arcsin(
            np.sqrt(
                [
                    ticks / SHOTS,
                    beta.ppf(ALPHA / 2, ticks, SHOTS - ticks + 1),
                    beta.ppf(1 - ALPHA / 2, ticks + 1, SHOTS - ticks),
                ]
            )
        )
        output.append(
            {
                "spacing": spacing,
                "transfers": int(ticks),
                "measured_spacing": float(times[0]),
                "simultaneous_interval": times[1:].tolist(),
            }
        )
    return output


def timing_contrasts(banks):
    summary = []
    for bank in banks:
        blank = next(
            b
            for b in banks
            if b["contrast"] == bank["contrast"] and b["source"] == "sham"
        )
        for j, read in enumerate(bank["tags"]["1"]["radar"]):

            def changes(tag):
                src = bank["tags"][tag]["radar"][j]["completed_by_deadline"]
                ref = blank["tags"][tag]["radar"][j]["completed_by_deadline"]
                return (np.array(src) - ref) / SHOTS

            diff = changes("1") - changes("2")
            summary.append(
                {
                    "contrast": bank["contrast"],
                    "source": bank["source"],
                    "cadence": read["cadence"],
                    "background_subtracted_tag_difference": diff.tolist(),
                    "simultaneous_difference_error": 4 * RADAR_EPS,
                }
            )
    return summary


def run():
    blocks, inserts = calibration()
    rng = np.random.default_rng(271124)
    receipts, banks = {}, []
    for k in (candidate.KAPPA, 0):
        ref = None
        for arm in ARMS:
            tags = {}
            for q in TAGS:
                clocks = read_clocks(k, q, arm, blocks, inserts, rng, receipts)
                if q == 1:
                    if arm == "sham":
                        ref = clocks
                    forecast = weighted_guide(clocks, ref)
                h = generator(k, q, arm, blocks, inserts)
                reads = []
                for dt in (1.0, 0.5) if k else (1.0,):
                    p = radar(h, dt)
                    counts = rng.multinomial(
                        SHOTS, np.maximum(0, np.diff(np.r_[0, p, 1.0])).round(13)
                    )
                    receipts[f"radar_{k}_{arm}_{q}_{dt}"] = counts
                    observed = counts.cumsum()[:-1] / SHOTS
                    if np.max(abs(observed - p)) > RADAR_EPS:
                        raise ValueError("radar counts outside simultaneous tolerance")
                    reads.append(
                        {
                            "cadence": dt,
                            "completed_by_deadline": counts.cumsum()[:-1].tolist(),
                            "unreturned": int(counts[-1]),
                            "oracle_cdf": p.round(9).tolist(),
                        }
                    )
                tags[str(q)] = {"clock_bonds": clocks, "radar": reads}
            banks.append(
                {
                    "contrast": k,
                    "source": arm,
                    "tags": tags,
                    "unit_clock_fixed_ruler_forecast": None
                    if forecast is None
                    else radar(forecast, 1.0).round(7).tolist(),
                }
            )
    return {
        "geometry": geometry(),
        "waits": WAITS.tolist(),
        "shots": SHOTS,
        "radar_deadlines": DEADLINES.tolist(),
        "simultaneous_cdf_error": RADAR_EPS,
        "seed": 271124,
        "contact_calibration": inserts,
        "banks": banks,
        "reference_clock": reference_clock(blocks[candidate.KAPPA, 1], receipts),
        "timing_contrasts": timing_contrasts(banks),
    }, receipts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="replace the saved count bank and summary"
    )
    args = parser.parse_args()
    result, receipts = run()
    report = HERE / "results.json"
    bank = HERE / "counts.npz"
    if args.write:
        report.write_text(json.dumps(result, indent=2) + "\n")
        np.savez_compressed(bank, **receipts)
    else:
        with np.load(bank) as saved:
            assert set(saved.files) == set(receipts)
            for key, arr in receipts.items():
                np.testing.assert_array_equal(saved[key], arr)
        prior = json.loads(report.read_text())

        def compare(a, b):
            if isinstance(b, dict):
                assert set(a) == set(b)
                for key in b:
                    compare(a[key], b[key])
            elif isinstance(b, list):
                assert len(a) == len(b)
                for x, y in zip(a, b, strict=True):
                    compare(x, y)
            elif isinstance(b, float):
                np.testing.assert_allclose(a, b, rtol=1e-7, atol=2e-8)
            else:
                assert a == b

        compare(result, prior)
    print(
        f"Verified {len(receipts)} count banks, two independent clocks and round trips."
    )


if __name__ == "__main__":
    main()
