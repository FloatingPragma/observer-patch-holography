"""Translation-averaged digital boundary of an independently ranged Euclidean ball."""

from fractions import Fraction


def ball_overlaps_per_pi(radius, spacing):
    """Volumes divided by pi for the ball and its translate by one grid spacing."""
    r, h = Fraction(radius), Fraction(spacing)
    if h <= 0 or 2 * r < h:
        raise ValueError("the two radius-r balls must overlap at positive spacing")
    cap = r - h / 2
    volume = 4 * r**3 / 3
    overlap = 2 * cap**2 * (r - cap / 3)
    return volume, overlap


def cut_area_per_pi(radius, spacing):
    """Mean exterior cubical face area divided by pi for a uniformly shifted grid."""
    h = Fraction(spacing)
    v, overlap = ball_overlaps_per_pi(radius, h)
    return 6 * (v - overlap) / h


def midpoint_count_per_pi(radius, spacing):
    """Mean number of actual face centers inside the smooth domain, divided by pi."""
    v, _ = ball_overlaps_per_pi(radius, spacing)
    return 3 * v / Fraction(spacing) ** 3
