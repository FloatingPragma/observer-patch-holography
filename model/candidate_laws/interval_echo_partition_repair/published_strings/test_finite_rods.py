from fractions import Fraction as F
from math import comb, pi

from model.candidate_laws.interval_echo_partition_repair.published_strings import finite_rods as f
from model.candidate_laws.interval_echo_partition_repair.published_strings import geometry as g


def test_actual_adjacent_rod_polynomial_and_integer_error():
    for n in (2, 3, 8, 21):
        for z in (F(0), F(1, 7), F(1, 2)):
            a, b = F(0), F(1, n)
            for _ in range(1, n):
                a, b = b, 2*b-a-z*b/n**2
            exact = sum((-z)**k*F(comb(n+k, 2*k+1), n**(2*k+1)) for k in range(n))
            assert exact == b
            assert abs(exact-1+(1-F(1, n*n))*z/6) <= z*z/120
            assert abs(f.rod_ratio(z, n, n**8)-exact) <= F(1, n**6)


def test_positive_finite_radial_cubature_and_angle_unit():
    assert all(w > 0 for _, w in f.RADIAL)
    for degree in range(8):
        assert sum(w*t**degree for t, w in f.RADIAL) == F(1, degree+1)
    for n in (32, 128, 512):
        ell = F(1, n)
        # Certify this against the same alternating expansion at a much finer
        # scale, with both rigorously proved tail bounds added.
        finer = f.angle_gain(ell/8)
        assert abs(f.angle_gain(ell)-finer) + (ell/8)**4/16 < ell**4/16


def test_ten_finite_volume_stops_predict_the_archived_collar_transfer():
    ordinary, held = F(1, 3), F(1, 3)+F(8, 33)
    for n in (32, 128, 512):
        ell = F(1, n)
        for u in g.FRAME:
            b2 = float(g.spatial_stretch(u))
            Kunit = f.angle_gain(ell)
            rods = []
            for density in (ordinary, held):
                K = Kunit*density
                read = f.read_ball(K, ell, u)
                assert read.stop_width <= ell**5
                ideal = g.measured_ball(float(density)*2*pi, float(ell), u)
                assert abs(float(read.radius)-ideal.radius) < 12*float(ell)**5
                assert abs(4*pi*float(read.area_over_4pi)-ideal.area) < 100*float(ell)**6
                W = g.weighted_volume(float(density)*2*pi, ideal.radius, u)
                # The quadrature and the stopped-radius bias enter W at order ell^8.
                # Keep a relative numerical allowance for the smooth float comparator.
                assert abs(4*pi*float(read.weighted_over_4pi)-W) < 5*float(ell)**8+2e-14*W
                rods.append(read.area_over_4pi)
            measured = 4*pi*float((rods[1]-rods[0])/ell**4)
            forecast = -4*pi/15 * 2*pi * float(held-ordinary)*b2
            assert abs(measured-forecast) < 90*float(ell)**2


def test_root_comparison_uses_only_finite_instrument(monkeypatch):
    def forbidden(*args, **kwargs):
        raise RuntimeError("a finite rod cannot read the smooth area oracle")
    monkeypatch.setattr(g, "product_measures", forbidden)
    ell, u = F(1, 64), g.FRAME[8]
    K = f.angle_gain(ell)*F(1, 3)
    read = f.read_ball(K, ell, u)
    assert read.weighted_over_4pi > 0 and ell < read.radius
    assert f.weighted_volume(K, read.radius, u, ell=ell) == 4*pi*float(read.weighted_over_4pi)
