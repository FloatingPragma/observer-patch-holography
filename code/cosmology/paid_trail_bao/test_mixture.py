"""Independent checks of the distance likelihood and radial feasibility."""

import json

import numpy as np

from . import mixture as m
from . import rings


def test_distance_and_stocks():
    r = json.loads((m.bao.HERE / "mixture.json").read_bytes())
    d = m.bao.DistanceData()
    for fit in (r["best"], *r["profile_95"]):
        _, v = m.bao.physical_background(d, fit["omega_b"], fit["q"], fit["H0"])
        res = d.obs - v
        score = (
            res @ np.linalg.solve(d.cov, res)
            + ((fit["omega_b"] - m.BBN) / m.SIGMA) ** 2
        )
        np.testing.assert_allclose(score, fit["chi2"], atol=1e-7)
    for fit in r["profile_95"]:
        assert abs(fit["chi2"] - r["best"]["chi2"] - m.DELTA) < 5e-5
    paid, extra = m.inventory()
    assert [str(q) for q in paid] == r["native_paid"]
    assert str(extra) == "51/202"


def test_all_curve_witnesses_and_one_ring_control():
    raw = json.loads((m.bao.HERE / "rings.json").read_bytes())
    meta, groups = rings.o.catalogue()
    kept, _ = rings.o.select(meta, groups)
    pair = raw["boxes"]["two_errors"]
    assert pair["missed"]["cosmic_95"] == ["F571-8"]
    assert not pair["outer_only_missed"] and not pair["missed"]["unlimited"]
    for box, result in raw["boxes"].items():
        k, ml_min, ml_max, gas_min, gas_max = rings.o.BOXES[box]
        for name, rows in kept.items():
            fit = result["galaxies"][name]["curves"]["cosmic_95"]
            t, stars, gas = (fit[x] for x in ("t", "scaled_stellar", "scaled_gas"))
            meta_g = meta[name]
            dl = k * meta_g["e_L36_gsun"] / meta_g["L36_gsun"]
            assert (
                ml_min * max(0, 1 - dl) - 1e-8 <= stars / t <= ml_max * (1 + dl) + 1e-8
            )
            assert gas_min - 1e-8 <= gas / t <= gas_max + 1e-8
            mass = np.array(fit["scaled_enclosed_mass"])
            assert min(mass) >= 0 and np.all(np.diff(mass) >= 0)
            budget = (
                stars * meta_g["L36_gsun"]
                + gas * rings.o.sparc.GAS_HELIUM * meta_g["MHI_gsun"]
            )
            assert mass[-1] <= raw["caps"]["cosmic_95"] * budget + 1e-8
            for row, h in zip(rows, mass, strict=True):
                v, sig = row["vobs_kms"], row["e_vobs_kms"]
                model = (
                    stars * row["vdisk_kms"] * abs(row["vdisk_kms"])
                    + gas * row["vgas_kms"] * abs(row["vgas_kms"])
                    + h * 1e9 * rings.o.G / row["rad_kpc"]
                )
                tol = v * v * (fit["widening"] + 1e-8)
                assert (
                    max(0, v - k * sig) ** 2 - tol <= model <= (v + k * sig) ** 2 + tol
                )
