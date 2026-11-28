"""Independent source, aperture, nuisance and aggregation checks."""

import json
import math

import numpy as np

from .measure import BOXES, G, HERE, bound, catalogue, native_caps, run, select


def test_native_work_and_every_reported_exceedance():
    r = json.loads((HERE / "results.json").read_text())
    assert r == json.loads(json.dumps(run()))
    assert r["source_sparc_sha256"] == {
        "table1.dat": "86771a7f9f7af768ba7932dc57c6a9d81de16a287d98f4cd49bd2b4fb2a42c5c",
        "table2.dat": "7d027e515441c6b4ebbf6aadee0327e6ad81156c4e8b151af2f6d62cb44c3962",
    }
    meta, groups = catalogue()
    for galaxy in r["galaxies"]:
        assert all(p["vbul_kms"] == 0 for p in groups[galaxy["galaxy"]])
        ring = galaxy["ring"]
        m = meta[galaxy["galaxy"]]
        mass = (0.5 * m["L36_gsun"] + 1.33 * m["MHI_gsun"]) * 1e9
        vsq = 0.5 * ring["vdisk_kms"] * abs(ring["vdisk_kms"]) + ring["vgas_kms"] * abs(
            ring["vgas_kms"]
        )
        ratio = (ring["vobs_kms"] ** 2 - vsq) * ring["rad_kpc"] / (G * mass)
        np.testing.assert_allclose(
            ratio, galaxy["bounds"]["nominal"]["q_required_min"], atol=1e-12
        )
        for key, stock in r["source_caps"].items():
            assert math.isclose(stock["cap"], 17 * int(key) / 101)
            assert math.isclose(stock["cap"] - stock["net"], 51 / 202)
            for box, v in galaxy["bounds"].items():
                assert (v["q_required_min"] > stock["cap"]) == (
                    v["caps"][key]["deficit_kms"] > 0
                )
    assert [r["all"]["two_errors"]["exceed_cap"][str(n)] for n in (14, 32, 64)] == [
        12,
        1,
        0,
    ]


def test_favorable_corners_contain_interior_and_signed_gas():
    caps = native_caps()
    m = {
        "dist_mpc": 5.0,
        "e_dist_mpc": 0.3,
        "inclination_deg": 87.0,
        "e_inclination_deg": 5.0,
        "L36_gsun": 0.3,
        "MHI_gsun": 0.07,
        "e_L36_gsun": 0.02,
    }
    ring = {
        "rad_kpc": 4.0,
        "vgas_kms": -37.0,
        "vdisk_kms": 74.0,
        "vobs_kms": 115.0,
        "e_vobs_kms": 3.0,
    }
    box = BOXES["two_errors"]
    result = bound(ring, m, box, caps)
    rng = np.random.default_rng(2114)
    for _ in range(300):
        # Scan the uncollapsed physical nuisance box independently.
        mscale = rng.uniform(0.3, 0.8) * rng.uniform(
            1 - 2 * 0.02 / 0.3, 1 + 2 * 0.02 / 0.3
        )
        gas = rng.uniform(0.7, 1.5)
        dist = rng.uniform(4.4, 5.6) / 5
        incl = rng.uniform(77, 90)
        vc = (
            rng.uniform(109, 121)
            * math.sin(math.radians(87))
            / math.sin(math.radians(incl))
        )
        rs = dist * 4
        mb = dist**2 * (mscale * 0.3e9 + gas * 0.07e9 * 1.33)
        baryon2 = dist * (mscale * 74**2 - gas * 37**2)
        q = (vc * vc - baryon2) * rs / (G * mb)
        assert result["q_required_min"] <= q
        for key, stock in caps.items():
            v = math.sqrt(baryon2 + stock["cap"] * G * mb / rs)
            assert v <= result["caps"][key]["maximum_speed_kms"]


def test_example_thresholds_and_pressure_control():
    a = json.loads((HERE / "results.json").read_text())
    obj = next(x for x in a["galaxies"] if x["galaxy"] == "NGC3741")
    m = obj["meta"]
    ring = obj["ring"]
    ex = a["NGC3741_controls"]
    q = a["source_caps"]["32"]["cap"]
    box = obj["bounds"]["two_errors"]
    s = box["stars_multiplier_at_min"]
    f = ex["thresholds"]["two_errors"]["gas_multiplier_for_32_seams"]
    d = box["distance_mpc"] / m["dist_mpc"]
    budget = d * (s * ring["vdisk_kms"] ** 2 + f * ring["vgas_kms"] ** 2)
    total = (s * m["L36_gsun"] + 1.33 * f * m["MHI_gsun"]) * 1e9 * d**2
    np.testing.assert_allclose(
        budget + G * q * total / box["radius_kpc"], box["v_floor_kms"] ** 2
    )
    assert all(
        v["q_required"]["two_errors"] < q < v["q_required"]["one_error"]
        for v in ex["removed_pressure_correction"]
    )


def test_reject_galaxy_with_any_bulge_even_on_rejected_ring():
    meta, groups = catalogue()
    name = "NGC3741"
    candidates = [dict(r) for r in groups[name]]
    candidates[0].update(vbul_kms=10, e_vobs_kms=999)
    included, tally = select({name: meta[name]}, {name: candidates})
    assert not included and tally == {"resolved_bulge": 1}
