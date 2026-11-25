"""Exact interior ambiguity on any tetrahedron from the certified native survey."""

from pathlib import Path
import json

import sympy as sp


def verify():
    survey = json.loads(Path(__file__).with_name("duration_geometry.json").read_text())
    volume_lower = sp.Rational(survey["common_metric_box"]["volume_lower"])
    # Edge lengths are at most three by the original capped twelve-appointment
    # radar protocol. An equilateral face maximizes area at fixed longest side.
    face_area_upper = 9 * sp.sqrt(3) / 4
    inradius_lower = 3 * volume_lower / (4 * face_area_upper)
    width = sp.Rational(1, 16)
    separation = 3 * width
    assert separation > 2 * width and separation + width < inradius_lower

    u, amplitude, source = sp.symbols("u amplitude source", real=True)
    bump = u**2 * (1 - u**2) ** 4
    assert bump.subs(u, 0) == 0
    assert all(sp.diff(bump, u, n).subs(u, 1) == 0 for n in range(4))
    integral = 4 * sp.pi * width**3 * sp.integrate(u**2 * bump, (u, 0, 1))
    base = 3 * width**2 / 8
    first, second = base * source, -base * source
    assert sp.expand((first + second) * integral) == 0
    clearance = sp.Rational(1, 20)
    assert inradius_lower - separation - width > clearance
    # Each support point lies this far from every edge line. Even with
    # maximum contraction, visiting it is slower than the straight edge.
    detour_margin = (1 - base) ** 2 * (1 + 4 * clearance**2 / 9) ** 3 - 1
    assert detour_margin > 0

    # At the first center g_ij = delta_ij, and all first derivatives vanish.
    # The only second derivatives are C delta_ij delta_kl. Compute Ricci
    # by differentiating the Christoffel formula, with an orthonormal time.
    second_radial = (
        sp.diff((1 + amplitude * bump) ** sp.Rational(2, 3), u, 2).subs(u, 0) / width**2
    )

    def d2(i, j, k, ell):
        return second_radial if i == j and k == ell else 0

    spatial_ricci = sp.Matrix(
        3,
        3,
        lambda i, j: sp.S.Half
        * sum(
            d2(k, j, k, i) + d2(i, k, j, k) - d2(i, j, k, k) - d2(k, k, i, j)
            for k in range(3)
        ),
    )
    assert spatial_ricci.subs(amplitude, first) == -source * sp.eye(3)

    # Independently expand the radial proper ruler and enclosed volume at
    # the center, using the actual line and volume elements of this metric.
    x = sp.symbols("x")
    omega = sp.series((1 + amplitude * bump) ** sp.Rational(1, 3), u, 0, 4).removeO()
    radius = sp.integrate(omega, (u, 0, x))
    volume = 3 * sp.integrate((1 + amplitude * bump) * u**2, (u, 0, x))
    coefficient = sp.limit((volume / radius**3 - 1) / (width * radius) ** 2, x, 0)
    assert sp.simplify(coefficient.subs(amplitude, first) - source / 10) == 0
    return dict(
        inradius_lower=str(inradius_lower),
        profile_radius=str(width),
        integrated_profile=str(integral),
        proper_volume_offset=str((first + second) * integral),
        detour_ratio_sixth_power_margin=str(detour_margin),
        null_ricci_at_center=str(-source),
        relative_small_ball_quadratic_coefficient=str(source / 10),
    )


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
