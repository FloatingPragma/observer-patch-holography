"""Decode ratios, contrasts, and held timing forecasts solely from terminal counts."""

from collections import defaultdict
from math import asin, ceil, log, pi, sin, sqrt


def decode(bank):
    shots = bank["shots"]
    rings = bank["cycles"]
    bonds = [tuple(e) for e in bank["bond_labels"]]
    if (
        bank["seed"] != 271216
        or shots != 131072
        or len(bonds) != 30
        or len(rings) != 12
        or len(bank["circulation"]) != 47
        or len(bank["populations"]) != 72
        or len(bank["quadrature_and_short_flight"]) != 60
    ):
        raise ValueError("Missing acquisition")

    scans = {
        (lev, o, 0, m, t)
        for lev in (0, 1)
        for o in (0, 13)
        for m in (1, 2)
        for t in (1, 3)
    }
    expected = {(-1, 0, p, m, 1, True) for p in range(12) for m in (1, 2)}
    expected.update((-1, 0, 0, m, 1, True) for m in (-1, -2, 0))
    expected.update(
        (-1, 0, 0, m, t, flag) for m in (1, 2) for t, flag in ((0, False), (3, True))
    )
    expected.update((*scan, True) for scan in scans)
    fields = ("level", "owner", "source", "winding", "time")
    sharp = {
        (-1, 0, 0, m, t, r) for r in (1, 0.5) for m in (1, 2, -1, -2) for t in range(7)
    }
    sharp.update((*scan, 1) for scan in scans)

    def labels(rows, extra):
        return {tuple(row[k] for k in (*fields, extra)) for row in rows}

    if (
        labels(bank["circulation"], "parity") != expected
        or labels(bank["populations"], "rate") != sharp
    ):
        raise ValueError("Missing registered acquisition")
    if set((a, b, s) for a, b, s, *_ in bank["quadrature_and_short_flight"]) != {
        (*bond, s) for bond in bonds for s in (-1, 1)
    }:
        raise ValueError("Missing signed calibration")
    from .apparatus import geometry
    from model.executable_model import cell

    cycles, opposites, _ = geometry()
    if (
        bonds != list(cell().seams)
        or rings != cycles
        or bank["antipodes"] != list(opposites)
    ):
        raise ValueError("Changed detector labels")

    def freqs(counts, n):
        if (
            len(counts) != n
            or any(type(c) is not int or c < 0 for c in counts)
            or sum(counts) != shots
        ):
            raise ValueError("Incomplete terminal bank")
        return [c / shots for c in counts]

    # 564 ring means, 288 sharp fractions and 120 calibration differences.
    terms = 47 * 12 + 72 * 4 + 60 * 2

    def radius(range_width, n):
        return (
            ceil(
                1e8 * (range_width * sqrt(log(2 * terms / 0.01) / (2 * n)) + 8 / 2**36)
            )
            / 1e8
        )

    edge_width = radius(2, shots)
    width = radius(2, 5 * shots)
    fraction_width = radius(1, shots)
    calibr_error = 0.0
    for a, b, sign, balanced, drift in bank["quadrature_and_short_flight"]:
        for counts, predicted in [(balanced, sign), (drift, sign / 2)]:
            left, right, _ = freqs(counts, 3)
            calibr_error = max(calibr_error, abs(right - left - predicted))
    table = {}
    ranges = defaultdict(list)
    signed = {}
    shape_checks = []
    for row in bank["circulation"]:
        level, owner, p, m, t = (
            row[k] for k in ("level", "owner", "source", "winding", "time")
        )
        flux = {}
        for (a, b), counts in zip(bonds, row["edges"], strict=True):
            left, right, _ = freqs(counts, 3)
            flux[a, b] = right - left
            flux[b, a] = left - right
        values = [
            sum(flux[a, b] for a, b in zip(r, r[1:] + r[:1], strict=True)) / 5
            for r in rings
        ]
        pred = []
        contrasts = []
        sign = 1 if abs(m) == 1 else -1
        feed = sin(2 * pi * m / 5) / 5
        for q in range(12):
            if q == p:
                cosine = 1
            elif q == bank["antipodes"][p]:
                cosine = -1
            elif q in rings[p]:
                cosine = sign / sqrt(5)
            else:
                cosine = -sign / sqrt(5)
            pred.append(cosine * feed)
            contrasts.append(
                abs(values[q] - cosine * values[p]) > (1 + abs(cosine)) * width
            )
        resolved = abs(values[p]) > width
        if level >= 0:
            shape_checks.append(
                [
                    level,
                    owner,
                    m,
                    t,
                    values[1] / values[0] if resolved else None,
                    sum(contrasts),
                ]
            )
        if level == -1 and t == 1 and m in (1, 2):
            key = f"winding{m}"
            if resolved:
                ranges[key].extend(values[q] / values[p] for q in rings[p])
            else:
                ranges[key].append(None)
            if p == 0:
                # Ratio interval from the joint bounds on the two raw means.
                bounds = (
                    [
                        (values[1] + a * width) / (values[0] + b * width)
                        for a in (-1, 1)
                        for b in (-1, 1)
                    ]
                    if resolved
                    else [None]
                )
                signed[key] = [values[0], values[1], min(bounds), max(bounds)]
        group = (
            f"open{level}_{owner}_t{t}"
            if level >= 0
            else "unsorted"
            if not row["parity"]
            else "isolated"
        )
        n, fail, max_error = table.get(group, [0, 0, 0.0])
        errors = [
            abs(measured - expected)
            for measured, expected in zip(values, pred, strict=True)
        ]
        table[group] = [
            n + 12,
            fail + sum(e > width for e in errors),
            max(max_error, max(errors)),
        ]
    waves = defaultdict(dict)
    open_readings = []
    for row in bank["populations"]:
        v = freqs(row["counts"], 4)
        if row["level"] == -1:
            waves[row["rate"], row["winding"]][row["time"]] = v[1]
        else:
            open_readings.append(
                [row["level"], row["owner"], row["winding"], row["time"], *v]
            )
    clocks = {}
    for (rate, m), data in waves.items():
        lo, hi = [
            2 * asin(sqrt(max(0, min(1, data[1] + s * fraction_width))))
            for s in (-1, 1)
        ]
        # A global Lipschitz bound on d sin²(ωt/2)/dω supplies direct
        # simultaneous transfer from the independent t=1 bank.
        w = 2 * asin(sqrt(data[1]))
        tests = [
            abs(meas - sin(w * t / 2) ** 2)
            <= fraction_width + t * max(w - lo, hi - w) / 2
            for t, meas in data.items()
            if t >= 2
        ]
        nominal = rate * (sqrt(5) + (1 if abs(m) == 1 else -1)) / 2
        full_fit = all(
            abs(f - sin(nominal * t / 2) ** 2) <= fraction_width
            for t, f in data.items()
        )
        clocks[f"{rate}_{m}"] = [w, lo, hi, sum(tests), full_fit]
    return dict(
        calibration_max_error=calibr_error,
        calibration_radius=edge_width,
        signed_radius=width,
        population_radius=fraction_width,
        marked_pair=signed,
        all_neighbor_ratios={
            k: ([min(v), max(v)] if None not in v else None) for k, v in ranges.items()
        },
        absolute_forecasts=table,
        open_marked_ratio_and_shape_exclusions=shape_checks,
        arrival_frequencies=clocks,
        open_populations=open_readings,
    )
