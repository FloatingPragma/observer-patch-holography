"""Read flux and discrimination intervals from terminal counts alone."""

import json
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import round_ceiling, to_int

# An integer on this grid divided by GRID is exactly representable in the archive.
GRID = 1 << 40


def relative_phase(iv, reference, measured, shots):
    # Multiplying a complex reading by the conjugate reference replaces complex division.
    # Both have the same positive shot denominator, which cancels out of the angle.
    r, s = (2 * c - shots for c in reference)
    x, y = (2 * c - shots for c in measured)
    real, imag = x * r + y * s, y * r - x * s
    if not (real or imag):
        raise ValueError("phase is undefined for a zero quadrature")
    return iv.atan2(iv.mpf(imag), iv.mpf(real))


def ceiling_grid(interval):
    """Round the upper binary interval endpoint upward using its exact mantissa."""
    return to_int((interval * GRID)._mpi_[1], round_ceiling)


def flux_ratios(iv, work, bank_count):
    shots = work["shots"]
    delta = iv.sqrt(2 * iv.ln(200 * bank_count) / shots)
    sine = iv.sqrt(2) * delta
    angular_radius = iv.atan2(sine, iv.sqrt(1 - sine**2))
    records = {}
    for row in work["banks"]:
        tau = iv.mpf(str(row["dwell"]))
        total = [iv.mpf(0) for _ in work["source_order"]]
        for edge in row["edges"]:
            controls = edge["plus_counts"][0]
            for k, pair in enumerate(edge["plus_counts"]):
                pos = relative_phase(iv, controls[0], pair[0], shots)
                neg = relative_phase(iv, controls[1], pair[1], shots)
                total[k] += (neg - pos) / (2 * tau)
        err = len(row["edges"]) * 2 * angular_radius / tau
        records[tuple(row["seam"]), row["radius"], row["dwell"]] = total, err
    charge = written = covered = 0
    for (seam, radius, tau), (v, eps) in records.items():
        if radius == 1:
            continue
        inside, inner_error = records[seam, 1, tau]
        gain = (inside[1] - inside[2]) / 2
        if (gain - inner_error).a <= 0:
            raise ValueError("unresolved flux calibration")
        for z in v[1:3]:
            charge = max(charge, ceiling_grid(abs(z) / gain))
            covered = max(covered, ceiling_grid((abs(z) + eps) / (gain - inner_error)))
        for z in v[3:5]:
            written = max(written, ceiling_grid(abs(z) / gain))
    return charge / GRID, written / GRID, covered / GRID


def binomial_cdf(iv, k, n, tick):
    """Enclose P(X <= k) at exactly p=tick/GRID by a positive finite sum."""
    if k < 0:
        return iv.mpf(0)
    if k >= n or tick == 0:
        return iv.mpf(1)
    if tick == GRID:
        return iv.mpf(0)
    if k > n // 2:
        return 1 - binomial_cdf(iv, n - k - 1, n, GRID - tick)
    p = iv.mpf(tick) / GRID
    mass = (1 - p) ** n
    total = mass
    odds = p / (1 - p)
    for j in range(k):
        mass = mass * (n - j) * odds / (j + 1)
        total += mass
    return total


def root_bracket(iv, k, n, target):
    """CDF decreases in p. Locate adjacent dyadic points on opposite sides of target."""
    lo, hi = 0, GRID
    while hi - lo > 1:
        mid = (lo + hi) // 2
        cdf = binomial_cdf(iv, k, n, mid)
        if cdf.a > target.b:
            lo = mid
        elif cdf.b < target.a:
            hi = mid
        else:
            raise ArithmeticError(
                "binomial tail comparison is unresolved at this precision"
            )
    return lo, hi


def endpoint_interval(iv, k, n, tail):
    # Lower limit solves P(X <= k-1)=1-tail and upper limit solves P(X <= k)=tail.
    lo = root_bracket(iv, k - 1, n, 1 - tail)[0] if k else 0
    hi = root_bracket(iv, k, n, tail)[1] if k < n else GRID
    return [lo / GRID, hi / GRID]


def intervals(work, mobile):
    iv = MPIntervalContext()
    iv.dps = 50
    bank_count = (
        4 * len(work["source_order"]) * sum(len(row["edges"]) for row in work["banks"])
    )
    charge, written, covered = flux_ratios(iv, work, bank_count)
    comparisons = mobile["banks"]
    selected = max(
        comparisons,
        key=lambda row: abs(
            row["receiver_records"]["write02"]["endpoint_hits"]
            - row["receiver_records"]["write11"]["endpoint_hits"]
        ),
    )
    tail = 1 / iv.mpf(400 * len(comparisons))
    endpoint_intervals = {
        source: endpoint_interval(
            iv,
            selected["receiver_records"][source]["endpoint_hits"],
            mobile["shots"],
            tail,
        )
        for source in ("write02", "write11")
    }
    return dict(
        real_quadrature_banks=bank_count,
        simultaneous_flux_coverage=0.99,
        maximum_abs_charged_outer_ratio=charge,
        maximum_abs_written_extreme_outer_ratio=written,
        simultaneous_upper_on_abs_charged_outer_ratio=covered,
        endpoint_selection={k: selected[k] for k in ("seam", "sign", "time", "launch")},
        simultaneous_selected_endpoint_intervals=endpoint_intervals,
    )


if __name__ == "__main__":
    data = json.loads(Path(__file__).with_name("results.json").read_text())
    print(json.dumps(intervals(data["boundary_work"], data["mobile_probe"]), indent=2))
