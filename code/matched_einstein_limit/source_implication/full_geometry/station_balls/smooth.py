"""Exact smooth chamber bounds and small-ball coefficients at an enrolled station."""

from itertools import combinations
import json

import sympy as sp

from .geometry import ALL, observed_squared, raw_bounds


def certificate():
    square = observed_squared(raw_bounds())
    c = sp.Rational(1, 10)
    shortest, longest = sp.Rational(13, 5), sp.Rational(31, 10)
    for i, j in combinations(ALL, 2):
        assert square[i, j].a > float(shortest**2)
        assert square[i, j].b < float(longest**2)
    two_tubes = 2 * (shortest - 2 * c)
    assert two_tubes > longest

    # For sectional curvature k of the reference center, geodesic polar area
    # density is f(r)^2 = r^2 -(k/3)r^4 + ... . The transverse deformation
    # multiplies this by 1+e*r^2. All first-order coefficients are exact.
    t, e, k = sp.symbols("r e k", real=True)
    f = t * (1 - k * t**2 / 6 + k**2 * t**4 / 120) * sp.sqrt(1 + e * t**2)
    area = sp.series(4 * sp.pi * f**2, t, 0, 6).removeO()
    volume = sp.integrate(area, (t, 0, t))
    dr = sp.symbols("b")
    v_fixed = sp.expand(volume.subs(t, t + dr * t**3)).coeff(t, 5)
    adjust = sp.solve(v_fixed, dr)[0]
    area_fixed = sp.expand(area.subs(t, t + adjust * t**3)).coeff(t, 4)
    ricci_spatial = -2 * sp.limit(sp.diff(f, t, 2) / f, t, 0)
    # The normal-coordinate quadratic increment of the tensor is universal
    # even when the starting metric is anisotropic. Its normal first jets
    # vanish at the center, so the change in Ricci is its linearized Ricci.
    x = sp.symbols("x0:3")
    cross = sp.Matrix.hstack(*(sp.Matrix(x).cross(sp.eye(3)[:, j]) for j in range(3)))
    a, b, c2, d, h, f2 = sp.symbols("a b c2 d h f2")
    shape = sp.Matrix([[a, b, c2], [b, d, h], [c2, h, f2]])
    delta = cross.T * shape * cross
    ricci = sp.Matrix(
        3,
        3,
        lambda i, j: sum(
            sp.diff(delta[j, m], x[i], x[m])
            + sp.diff(delta[i, m], x[j], x[m])
            - sp.diff(delta[i, j], x[m], x[m])
            - sp.diff(delta[m, m], x[i], x[j])
            for m in range(3)
        )
        / 2,
    )
    assert (delta * sp.Matrix(x)).applyfunc(sp.expand) == sp.zeros(3, 1)
    assert ricci == 3 * (shape - sp.trace(shape) * sp.eye(3))
    isotropic = ricci.subs({a: 1, d: 1, f2: 1, b: 0, c2: 0, h: 0})
    assert isotropic == -6 * sp.eye(3)
    skew = ricci.subs({a: 1, d: -1, f2: 0, b: 0, c2: 0, h: 0})
    assert skew == sp.diag(3, -3, 0)
    chi_surface = 6 * (
        2 - 5
    )  # cylinders and their attaching circles have Euler number zero
    genus = 1 - chi_surface // 2
    collar_width = sp.Integer(2)
    assert genus == 10 and 2 * collar_width > longest
    out = dict(
        smooth_three_sphere_bicollar_width=str(collar_width),
        doubled_exit_lower=str(2 * collar_width),
        station_surface_genus=genus,
        arbitrary_normal_chart_ricci_increment=[
            list(map(int, row)) for row in isotropic.tolist()
        ],
        tagged_ray_ricci_slopes_at_matching_rest_ball_volumes=[
            int(skew[0, 0]),
            -int(skew[0, 0]),
        ],
        common_station_to_outlet=str(c),
        native_length_bounds=list(map(str, (shortest, longest))),
        two_full_cylinder_lower=str(two_tubes),
        longest_direct_upper=str(longest),
        curvature_increment=str(sp.diff(3 * ricci_spatial, e)),
        ultrastatic_null_ricci_increment=str(sp.diff(ricci_spatial, e)),
        fixed_radius_area_increment=str(sp.diff(area.coeff(t, 4), e)),
        fixed_radius_volume_increment=str(sp.diff(volume.coeff(t, 5), e)),
        fixed_volume_area_increment=str(sp.diff(area_fixed, e)),
    )
    assert (
        out["ultrastatic_null_ricci_increment"] == "-6"
        and out["curvature_increment"] == "-18"
    )
    assert out["fixed_volume_area_increment"] == "12*pi/5"
    return out


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
