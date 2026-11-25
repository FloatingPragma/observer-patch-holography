"""Compare probe recoil with independent size records for a mobile two-load source."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from model.candidate_laws.interval_echo.openwave.row_21_boost import measure as pair
from model.candidate_laws.interval_echo.openwave.row_30.experiment import (
    native_matrix,
    select_route,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
SHOTS = 2**40
# Union bound over <= 4,000 tails (all spectra, source histograms and controls).
TAIL = 0.01 / 4000
DWELL = np.pi / 32


def bin_interval(count, n):
    return (
        float(beta.ppf(TAIL, count, n - count + 1)) if count else 0.0,
        float(beta.isf(TAIL, count + 1, n - count)) if count < n else 1.0,
    )


def counts(weights, rng):
    p = np.array(weights, dtype=float)
    if p.min() < -1e-13 or abs(p.sum() - 1) > 1e-11:
        raise ArithmeticError("incomplete probe outcomes")
    # Fix a portable binary resolution. Fill the largest bin from the complement.
    m = int(np.argmax(p))
    units = np.rint(p * 2**50).astype(np.int64)
    units[m] += 2**50 - units.sum()
    return rng.multinomial(SHOTS, units / 2**50).tolist()


def scatter(length, core, amp, visits, rng):
    # Each source configuration leaves the pointer a separate pure outgoing
    # vector. Translating both occupied ports changes only a Fourier phase.
    spec = np.zeros(length)
    spill = np.zeros(2)
    for sites, weight in core:
        occupied = np.zeros(length, dtype=bool)
        occupied[sites] = True
        pointer = np.full(length, 1 / np.sqrt(length), dtype=complex)
        for j in range(visits):
            queried = np.roll(occupied, -j)
            spill[j] += weight * (1 - abs(amp) ** 2) * sum(abs(pointer[queried]) ** 2)
            pointer[queried] *= amp
        spec += weight * abs(np.fft.fft(pointer, norm="ortho")) ** 2
    missed = 1 - sum(w for _, w in core)
    return counts([*spec, *spill, missed], rng)


def analyze(length, sep_counts, single, one, two):
    n = length // 2
    d = np.arange(1, n + 1)
    q = 2 * np.pi * np.arange(length) / length
    sizes = np.array(sep_counts[:-1], dtype=float)
    N = int(sizes.sum())
    sep = sizes / N
    # q is read from the terminal address-pointer bin. All nonforward bins
    # of the single-port standard have identical strength.
    factor = length**2 / (length - 1)
    gain = factor * sum(single[1:length]) / SHOTS
    glo, ghi = np.array(bin_interval(sum(single[1:length]), SHOTS)) * factor
    size_err = np.sqrt(np.log(1 / TAIL) / (2 * N))
    d2 = sep @ d**2
    rad = np.sqrt(d2) / 2
    rad_int = (
        np.sqrt(
            np.maximum(0, [d2 - (n * n - 1) * size_err, d2 + (n * n - 1) * size_err])
        )
        / 2
    )
    # Rounding an acquisition probability to 2^-50 has negligible effect,
    # but allocate its full possible L1 error even after admission.
    guard = length * 2**-49 * length**2 / (1 - sep_counts[-1] / sum(sep_counts))
    cosine = np.cos(np.outer(q, d)) @ sep
    # Single-port spill determines the real phase-interference term too:
    # 2 Re(amp-1) = -gain - L * P(spill | one sharp double).
    losslo, losshi = length * np.array(bin_interval(single[length], SHOTS))
    adj_int = bin_interval(sep_counts[0], N)
    matched_two = []
    reject_additive = []
    passed = []
    ratios = []
    for j in range(1, length):
        lo, hi = np.array(bin_interval(one[j], sum(one[:-1]))) * length**2
        predicted_lo = max(0.0, 2 * glo * (1 + cosine[j] - 2 * size_err) - guard)
        predicted_hi = 2 * ghi * (1 + cosine[j] + 2 * size_err) + guard
        passed.append(bool(lo <= predicted_hi and hi >= predicted_lo))
        phase_weight = 2 + 2 * np.cos(q[j])
        add_lo = max(0.0, phase_weight * predicted_lo - guard)
        add_hi = phase_weight * predicted_hi + guard
        # The correction at d=1 is gain² - A*gain*(gain+loss),
        # A=|1+exp(iq)|². Endpoint enclosure keeps the calibration errors.
        corr_lo = glo**2 - phase_weight * ghi * (ghi + losshi)
        corr_hi = ghi**2 - phase_weight * glo * (glo + losslo)
        products = [c * v for c in (corr_lo, corr_hi) for v in adj_int]
        full_lo, full_hi = (
            max(0.0, add_lo + min(products) - guard),
            add_hi + max(products) + guard,
        )
        tlo, thi = np.array(bin_interval(two[j], sum(two[:-1]))) * length**2
        matched_two.append(bool(tlo <= full_hi and thi >= full_lo))
        reject_additive.append(bool(tlo > add_hi or thi < add_lo))
        ratios.append(
            np.sqrt([max(0.0, lo - guard) / (4 * ghi), (hi + guard) / (4 * glo)])
        )
    ratios = np.array(ratios)[:n]
    below = np.flatnonzero(ratios[:, 1] < 0.9)
    i = int(below[0] + 1) if below.size else None
    # Last above-threshold bin before the first below; zero is the
    # normalization set by twice the measured single-load amplitude.
    previous = max([0] + [j for j in range(1, i or 1) if ratios[j - 1, 0] > 0.9])
    crossing = (
        None
        if i is None
        else [float(q[previous] * rad_int[0]), float(q[i] * rad_int[1])]
    )
    # In the two-contact spectrum at q=pi the entire additive response
    # vanishes, leaving only two hits on one occupied adjacent pair.
    predicted = [
        max(0.0, adj_int[0] * glo**2 - guard) / length**2,
        (adj_int[1] * ghi**2 + guard) / length**2,
    ]
    seen = bin_interval(two[n], sum(two[:-1]))
    return dict(
        radius=rad,
        radius_interval=rad_int.tolist(),
        single_strength=gain,
        single_strength_interval=[float(glo), float(ghi)],
        matched_one_contact_bins=sum(passed),
        multiple_contact_matches=sum(matched_two),
        additive_two_contact_misses=sum(reject_additive),
        nonforward_bins=length - 1,
        first_below_bin=i,
        previous_above_bin=previous,
        qR_crossing=crossing,
        first_apparent_amplitude_interval=None if i is None else ratios[i - 1].tolist(),
        adjacency=sep[0],
        two_contact_pi_clicks=two[n],
        two_contact_pi_interval=list(seen),
        two_contact_pi_prediction=predicted,
        double_visit_agrees=bool(seen[0] <= predicted[1] and seen[1] >= predicted[0]),
    )


def run():
    prior = json.loads((pair.HERE / "results.json").read_text())
    rng = np.random.default_rng(230322)
    geometries = []
    studies = []
    for level in range(3):
        board = federation(level)
        route, _ = select_route(board)
        L = len(route)
        off = [sorted(set(board.neighbours[p]) - set(route)) for p in route]
        gaps = [p for p, opts in zip(route, off) if not opts]
        geometries.append(
            dict(cells=board.cells, ports=L, missing_outside_neighbours=gaps)
        )
        if gaps:
            continue
        for k in (8, 0):
            source = next(
                s
                for s in prior["studies"]
                if s["board_cells"] == board.cells and s["contrast"] == k
            )
            if list(route) != source["route"]:
                raise ValueError("guide changed")
            amp = native_matrix(DWELL, k)[0, 0]
            single = scatter(L, [([0], 1.0)], amp, 1, rng)
            blank = scatter(L, [([], 1.0)], amp, 1, rng)
            entries = pair.transfers(L // 2, k)
            panels = []
            for winding in [0, max(p["winding"] for p in source["panels"])]:
                panel = next(p for p in source["panels"] if p["winding"] == winding)
                _, wave = pair.line(L // 2, entries, winding)
                weight = np.mean(abs(wave(16)) ** 2, axis=1)
                core = [
                    ([idx % 2, (idx % 2 + 1 + idx // 2) % L], w)
                    for idx, w in enumerate(weight)
                ]
                one = scatter(L, core, amp, 1, rng)
                two = scatter(L, core, amp, 2, rng)
                sep_counts = panel["separation"]["16"]["counts"]
                decoded = analyze(L, sep_counts, single, one, two)
                if (
                    decoded["matched_one_contact_bins"] != L - 1
                    or decoded["multiple_contact_matches"] != L - 1
                    or not decoded["double_visit_agrees"]
                ):
                    raise ValueError(
                        "independent shape prediction misses counted signal"
                    )
                panels.append(
                    dict(
                        winding=winding,
                        one_contact=one,
                        two_contact=two,
                        reading=decoded,
                    )
                )
            studies.append(
                dict(
                    contrast=k,
                    route=route,
                    blank_neighbour=[v[0] for v in off],
                    single=single,
                    empty=blank,
                    panels=panels,
                )
            )
    return dict(
        shots=SHOTS,
        ordering="Fourier bins j=0..L-1, spill 1, spill 2, source failed",
        geometry=geometries,
        studies=studies,
    )


def verify(current, saved, path="records"):
    if isinstance(current, dict):
        if current.keys() != saved.keys():
            raise ValueError(path)
        for k in current:
            verify(current[k], saved[k], path + "." + k)
    elif isinstance(current, list):
        if len(current) != len(saved):
            raise ValueError(path)
        for i, (a, b) in enumerate(zip(current, saved)):
            verify(a, b, path + f"[{i}]")
    elif isinstance(current, float):
        if not np.isclose(current, saved, atol=5e-10, rtol=1e-9):
            raise ValueError(path)
    elif current != saved:
        raise ValueError(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run()
    if args.output:
        with args.output.open("x") as fp:
            json.dump(record, fp, indent=2, allow_nan=False)
            fp.write("\n")
    else:
        verify(record, json.loads((HERE / "records.json").read_text()))
    for s in record["studies"]:
        for panel in s["panels"]:
            print(s["contrast"], panel["winding"], panel["reading"])
