"""Gates of the map from four-dimensional parameters and temperature to (g3^2, x, y)."""
import math

from . import reduction as rd


def test_constants_of_the_reduction():
    assert abs(rd.C_B - 1.953808) < 2e-6
    assert abs(rd.C_F - 0.567514) < 2e-6
    assert abs(rd.C_CONST + 0.348725) < 5e-6


def test_special_values_of_F():
    m = 80.0
    assert abs(rd.F_k_m1_m2(m, m, m) - (2.0 - math.pi / math.sqrt(3.0))) < 1e-12
    r = 1.3                                       # F(m1; m2, m2) with r = m2/m1
    closed = 2.0 - 2.0 * math.sqrt(4 * r * r - 1) * math.atan(1.0 / math.sqrt(4 * r * r - 1))
    assert abs(rd.F_k_m1_m2(m, r * m, r * m) - closed) < 1e-12
    r = 1.5                                       # F(m1; m1, m2) with r = m2/m1
    closed = (1.0 - r * r * (3 - r * r) / (1 - r * r) * math.log(r)
              - 2.0 * r * math.sqrt(4 - r * r) * math.atan(math.sqrt(2 - r) / math.sqrt(2 + r)))
    assert abs(rd.F_k_m1_m2(m, m, r * m) - closed) < 1e-12


def test_linearized_form_is_independent_of_the_scale():
    a = rd.reduce(160.0, **rd.SM, linear=True, mu_over_T=math.pi)
    b = rd.reduce(160.0, **rd.SM, linear=True, mu_over_T=4 * math.pi)
    for key in ("x", "y", "g32"):
        assert abs(a[key] - b[key]) < 1e-10 * max(1.0, abs(a[key]))


def test_vacuum_relations_run_with_the_one_loop_beta_functions():
    sm = rd.SM
    g02 = 4 * math.sqrt(2.0) * sm["G_mu"] * sm["mW"] ** 2
    h2 = (sm["mH"] / sm["mW"]) ** 2
    t2 = (sm["mt"] / sm["mW"]) ** 2
    s2 = 4 * math.pi * sm["alpha_s"] / g02
    lam, gy2, gs2 = g02 * h2 / 8.0, g02 * t2 / 2.0, s2 * g02
    lo = rd.vacuum_brackets(100.0, sm["mW"], sm["mH"], sm["mt"], s2)
    hi = rd.vacuum_brackets(100.0 * math.e, sm["mW"], sm["mH"], sm["mt"], s2)
    dg, dH, dW, dS = (hi[i] - lo[i] for i in range(4))
    # mu d/dmu of g^2, nu^2, lambda, g_Y^2 in units of 1/(8 pi^2), hep-ph/9508379 section 4.1
    assert abs(dg * g02 / 2.0 - (8 * rd.NF + rd.NS - 44) / 6.0 * g02) < 1e-9
    assert abs(dH * g02 / 2.0 - (-2.25 * g02 + 6 * lam + 3 * gy2)) < 1e-9
    run_lam = (9 / 16.0 * g02 ** 2 - 4.5 * lam * g02 + 12 * lam ** 2 - 3 * gy2 ** 2 + 6 * lam * gy2) / lam
    assert abs((-dW + dg + dH) * g02 / 2.0 - run_lam) < 1e-9
    run_gy = 4.5 * gy2 - 2.25 * g02 - 8 * gs2
    assert abs((-dW + dg + 2 * dS) * g02 / 2.0 - run_gy) < 1e-9


def test_top_relation_contains_the_pole_mass_correction():
    sm = rd.SM
    g02 = 4 * math.sqrt(2.0) * sm["G_mu"] * sm["mW"] ** 2
    s2 = 4 * math.pi * sm["alpha_s"] / g02
    with_qcd = rd.vacuum_brackets(sm["mt"], sm["mW"], sm["mH"], sm["mt"], s2)[3]
    without = rd.vacuum_brackets(sm["mt"], sm["mW"], sm["mH"], sm["mt"], 0.0)[3]
    assert abs(2 * (with_qcd - without) * g02 * rd.LOOP + 8.0 / 3.0 * sm["alpha_s"] / math.pi) < 1e-12


def test_standard_model_curves_of_the_literature():
    # D'Onofrio, Rummukainen 2016, figure 1: x = 0.2899, g3^2/T = 0.39282 at 160 GeV, y = 0 at 162.1 GeV
    r = rd.reduce(160.0, **rd.SM, linear=True)
    assert abs(r["x"] - 0.2899) < 0.003
    assert abs(r["g32"] / 160.0 - 0.39282) < 0.001
    assert abs(rd.solve_T(lambda x: 0.0, rd.SM, linear=True) - 162.1) < 0.6
    # one-loop terms move the temperature scale by about 16 percent
    assert 135.0 < rd.solve_T(lambda x: 0.0, rd.SM, tree=True) < 145.0


def test_cosmic_time():
    t = rd.cosmic_time_seconds(160.0)
    assert abs(t / 9.150e-12 - 1.0) < 1e-3
    assert abs(rd.ticks(t) / 1.329e32 - 1.0) < 1e-3
