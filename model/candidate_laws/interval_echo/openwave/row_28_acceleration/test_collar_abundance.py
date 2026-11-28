from itertools import product

import numpy as np
from scipy.optimize import linprog

from . import collar_abundance as test_sky


def test_minimum_mass_with_two_free_background_terms():
    report = test_sky.run()
    assert report["measured_ratio"] == "8/11"
    box = report["external_box"]
    curved = report["curved_constant_vacuum"]
    assert len(curved["triples"]) == curved["total_triples"] == 20
    data = test_sky.DistanceData()
    radial = data.kind == "DH_over_rs"
    s = 1 + data.z[radial]
    dh = data.obs[radial]
    width = curved["k"] * np.sqrt(np.diag(data.cov)[radial])
    scale = (100 * box["ruler_upper"] / test_sky.C_LIGHT) ** 2
    hi, lo = (dh - width) ** -2 / scale, (dh + width) ** -2 / scale
    design = np.column_stack((s**3, s**2, np.ones(len(s)), s**4))
    fit = linprog(
        [1, 0, 0, 0],
        A_ub=np.vstack((design, -design)),
        b_ub=np.r_[hi, -lo],
        bounds=[(0, None), (None, None), (None, None), (0, box["radiation_upper"])],
    )
    assert fit.success
    assert abs(fit.fun - curved["strongest"]["physical_dust_lower"]) < 1e-10
    assert 1.23 < curved["strongest"]["extra_per_baryon_lower"] < 1.24


def test_fixed_capacity_box_endpoints():
    d = test_sky.run()
    b = d["external_box"]
    minima = [
        test_sky.required(dh, d["z"], ruler, length, baryon, nu, rad)
        for dh, ruler, length, baryon, nu, rad in product(
            [d["DH_over_rd"], b["dh_upper"]],
            [test_sky.R_DRAG, b["ruler_upper"]],
            [test_sky.L_PLANCK, b["length_lower"]],
            [test_sky.OMEGA_B, b["baryon_upper"]],
            [0, b["neutrino_upper"]],
            [0, b["radiation_upper"]],
        )
    ]
    assert 4.36 < d["minimum_required_ratio"] == min(minima) < 4.37
    assert d["minimum_DH_at_ceiling"] > b["dh_upper"]
