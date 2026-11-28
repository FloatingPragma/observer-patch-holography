"""Exact feasibility of a local potential at selected catalogue radii."""

from collections import defaultdict
from fractions import Fraction as Q


def decimal(x):
    return Q(str(x))


def slope_events(radius, edges, weights, shape):
    """Slopes of 4π s ρ for positive length s, at one window boundary."""
    events = defaultdict(Q)
    for lo, hi, seed in zip(edges[:-1], edges[1:], weights, strict=True):
        weight = seed / (hi**2 - lo**2)
        for hop, share in shape:
            if not hop:
                inside = Q(int(lo < radius < hi)) + Q(int(radius in (lo, hi)), 2)
                events[0] += 2 * weight * share * inside / radius
                continue
            knots = sorted(
                {
                    Q(0),
                    abs(radius - lo),
                    abs(radius - hi),
                    radius,
                    radius + lo,
                    radius + hi,
                }
            )
            prefactor = weight * share / (radius * hop)
            previous = Q(0)
            for a, b in zip(knots, knots[1:]):

                def width(h):
                    return max(0, min(hi, radius + h) - max(lo, abs(radius - h)))

                slope = prefactor * hop * (width(b) - width(a)) / (b - a)
                events[a / hop] += slope - previous
                previous = slope
            events[knots[-1] / hop] -= previous
    return events


def segment_feasible(a, b, lower, upper):
    """Solve for a position in the segment and eliminate one positive amplitude."""
    left, right = Q(0), Q(1)
    for i in range(len(a)):
        for j in range(len(a)):
            x = upper[j] * a[i] - lower[i] * a[j]
            y = upper[j] * b[i] - lower[i] * b[j]
            if x < 0 and y < 0:
                return False
            if x < 0:
                left = max(left, x / (x - y))
            if y < 0:
                right = min(right, x / (x - y))
    if left > right:
        return False
    middle = (left + right) / 2
    return all((1 - middle) * x + middle * y > 0 for x, y in zip(a, b, strict=True))


def certify(rows, annuli, shape, indices):
    """Exact for the decimal map and baryonic speeds used by the held test."""
    edges, weights = [[decimal(x) for x in vector] for vector in annuli]
    hops, shares = shape
    # There are exactly 64 endpoints. Recover their integer multiplicities.
    counts = [int(round(x * 64)) for x in shares]
    assert sum(counts) == 64 and all(
        x == n / 64 for x, n in zip(shares, counts, strict=True)
    )
    shell = list(zip(map(int, hops), [Q(n, 64) for n in counts], strict=True))
    lower, upper, slopes = [], [], []
    for i in indices:
        r, obs, bar, err = [
            decimal(rows[k][i])
            for k in ("radius_kpc", "observed_kms", "baryons_kms", "error_kms")
        ]
        lower.append(max(0, obs - 2 * err) ** 2 - bar**2 - Q(1, 1000))
        upper.append((obs + 2 * err) ** 2 - bar**2 + Q(1, 1000))
        low = slope_events(Q(98, 100) * r, edges, weights, shell)
        high = slope_events(Q(102, 100) * r, edges, weights, shell)
        slopes.append({k: low[k] - high[k] for k in set(low) | set(high)})
    assert min(lower) > 0
    knots = sorted(set().union(*slopes))
    knots.append(knots[-1] + 1)  # Beyond the last corner only zero-hop stock remains.
    value = [Q(0)] * len(indices)
    slope = value.copy()
    possible = 0
    for a, b in zip(knots, knots[1:]):
        slope = [x + events.get(a, 0) for x, events in zip(slope, slopes, strict=True)]
        following = [v + m * (b - a) for v, m in zip(value, slope, strict=True)]
        possible += int(segment_feasible(value, following, lower, upper))
        value = following
    # Past the final corner the common shape is exactly proportional to s.
    assert all(v == knots[-1] * m for v, m in zip(value, slope, strict=True))
    return dict(
        ring_indices=indices,
        radius_kpc=[rows["radius_kpc"][i] for i in indices],
        intervals_including_tail=len(knots) - 1,
        feasible_intervals=possible,
    )
