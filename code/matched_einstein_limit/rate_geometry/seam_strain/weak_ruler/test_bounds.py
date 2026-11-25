from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.integrate import dblquad, quad
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from . import bounds


def test_two_leg_simplex_integral_and_actual_matching():
    q, low = bounds.continuum(5,Fraction(5,8))
    def integrand(t,s):
        return ((5./8)**10*(s*t)**5/120**2)**2
    mass = dblquad(integrand,0,6,0,lambda s:6-s)[0]
    rebate = dblquad(lambda t,s:(6-s-t)*integrand(t,s)/2,0,6,0,lambda s:6-s)[0]
    np.testing.assert_allclose([float(q), float(3-low)], [mass,rebate], rtol=1e-11)
    assert bounds.discrete(5,Fraction(5,8),16)[1] > Fraction(29999,10000)


def test_passive_continuous_and_weak_discrete_arrivals_on_a_dimer():
    # The weak conditional generator retains the receiving-site amplitude.
    h = np.array([[.2,1.1],[1.1,-.3]])
    p = np.diag([0,1])
    g = -1j*h-p/2
    t = 2.5
    live = quad(lambda s:abs(expm(g*s)[1,0])**2,0,t)[0]
    free = quad(lambda s:abs(expm(-1j*h*s)[1,0])**2,0,t)[0]
    assert 0 < live < free
    step = 1/16
    u = expm(-1j*h*step)
    miss = np.diag([1,np.exp(-step/2)])
    wave = np.array([1.,0.])
    free_wave = wave
    possible = recorded = 0
    for _ in range(40):
        wave, free_wave = u@wave, u@free_wave
        recorded += -np.expm1(-step)*abs(wave[1])**2
        possible += step*abs(free_wave[1])**2
        wave = miss@wave
    assert 0 < recorded < possible
    assert abs(recorded-(1-np.vdot(wave,wave).real)) < 2e-14


def test_dyadic_polynomial_blocks_against_an_independent_exponential():
    rng = np.random.default_rng(834619)
    a = rng.integers(-36,37,size=(16,16))
    a = a+a.T
    weights, scalar_error = bounds.checked_coefficients()
    bands, error = bounds.arrival_grid(csr_matrix(a),np.arange(16).reshape(2,8),weights,scalar_error)
    for t in (0,1,31,85,367,384):
        actual = expm(-1j*(t/64)*a/40)[8:,:8]
        re,im = (x[t,1,0]/(1<<bounds.OUT_BITS) for x in bands)
        assert np.linalg.norm(actual-(re+1j*im),2) < float(error)


def test_occupied_source_in_both_charges_and_six_held_stations():
    fresh = bounds.run()
    saved = json.loads(Path(bounds.__file__).with_suffix(".json").read_text())
    for field in ("stations","continuous","discrete","cell_distances","intercell_operator_norm"):
        assert fresh[field] == saved[field]
    assert fresh["stations"] == [0,264,372,396,540,660]
    assert min(fresh["cell_distances"][i][j] for i,j in combinations(range(6),2)) == 5
    assert max(map(Fraction,fresh["intercell_operator_norm"].values())) == Fraction(5,8)
    assert Fraction(fresh["continuous"]["range_lower"]) > Fraction(299991,100000)
    fresh_local = fresh["original_local_six"]
    saved_local = saved["original_local_six"]
    assert fresh_local["stations"] == saved_local["stations"] == [0,39,5,7,40,42]
    for charge in ("-1","1"):
        for kind in ("continuous","spacing_1/16"):
            new = fresh_local["by_charge"][charge][kind]
            old = saved_local["by_charge"][charge][kind]
            # Independent platforms may propose slightly different Bessel
            # coefficients. Both must pass exact scalar and matrix checks.
            for entry in ("range_lower","completion_upper"):
                assert abs(Fraction(new[entry])-Fraction(old[entry])) < Fraction(1,500000)
            assert new["geometry"] == old["geometry"]
            assert Fraction(new["range_lower"]) > Fraction(29,10)
            assert Fraction(new["geometry"]["round_three_sphere_max_distance_change_strict_lower"]) > Fraction(1,5)
            assert Fraction(new["geometry"]["convex_normal_ball_abs_sectional_curvature_times_radius_squared_strict_lower"]) > Fraction(1,2)
