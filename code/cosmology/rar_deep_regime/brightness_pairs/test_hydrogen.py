"""Check observed gas normalizations without fitting the orbital velocities."""

import csv
from fractions import Fraction as F
from hashlib import sha256
import json

from . import hydrogen as h
from . import sensitivity as s


def test_hydrogen_source_snapshots_and_complete_selection():
    masses = s.HERE / "hydrogen_masses.tsv"
    dish = s.HERE / "ugc06667_flux.tsv"
    assert (
        sha256(masses.read_bytes()).hexdigest()
        == "8711b79057bb1a1aa33f2d9154927f8d7b0c079566f38944baaab0a1d4658c1e"
    )
    assert (
        sha256(dish.read_bytes()).hexdigest()
        == "17f0ccfd348e51c69231e1a3b8b3da605d58fad8e852aaaf37a56d643afba453"
    )
    with masses.open() as f:
        published = list(csv.DictReader(f, delimiter="\t"))
    assert len({r["agc"] for r in published}) == len(published) == 20
    original = json.loads((s.HERE / "results.json").read_text())["cuts"]["deep"][
        "pairs"
    ]
    original_names = {p["galaxy"] for p in original}
    record = h.run()
    assert record == json.loads((s.HERE / "hydrogen.json").read_text())
    assert set(record["published_matches"]) == original_names & {
        r["galaxy"] for r in published
    }
    assert len(set(record["pairs_without_published_match"])) == 60
    assert set(record["pairs_without_published_match"]) == original_names - set(
        record["published_matches"]
    )
    assert record["counts"]["rounding_and_mass_errors"][
        "varied_requiring_original_cuts"
    ] == {"overlap": 8, "smaller_faint": 1, "no_eligible_force": 1}


def test_flux_limit_from_direct_annulus_algebra():
    pairs = json.loads((s.HERE / "results.json").read_text())["cuts"]["deep"]["pairs"]
    target = next(p for p in pairs if p["galaxy"] == "UGC06667")
    table = (s.HERE.parent / "data/table2.dat").read_text().splitlines()
    f, hg, _, q = h.records(target, table)
    # At fixed d the lower shape inequality reads c*(Qf-q*Qh) <= d*(q*Dh-Df).
    assert f[0] > q * hg[0] and q * hg[1] > f[1]
    largest = F(".6") * (q * hg[1] - f[1]) / (f[0] - q * hg[0])
    data = h.run()
    observed = data["single_dish_target"]
    assert data["lower_shape_gas_domains"]["UGC06667"]["retaining_cuts"][
        "upper"
    ] == float(largest)
    catalogue_flux = F(809000000, 235600 * 18**2)
    assert observed["largest_shape_allowed_flux_Jy_km_s"] == float(
        catalogue_flux * largest
    )
    # Even if the flux AND its errors were rounded in the most helpful direction,
    # the dish/synthesis ratio cannot approach the maximum admitted gas fraction.
    assert (F("11.55") - 2 * F("2.25")) / (F("11.05") + 2 * F("0.45")) > largest
    assert (F("11.6") - 2 * F("2.2")) / catalogue_flux > largest


def test_projection_tracks_open_endpoints_and_signs():
    # d>=0.4, d<=0.6; c>=d and c+d<1.2. The result is [0.4,0.8).
    got = h.project(
        [(-F(1), F(1), 0, False), (F(1), F(1), F("1.2"), True)], s.DISK, h.Interval()
    )
    assert got == h.Interval(F(".4"), F(".8"), False, True)
    for expected, data in [
        ("larger_faint", ((F(1), F(1)), (F(1), F(1)), F(2), F(3))),
        ("smaller_faint", ((F(1), F(1)), (F(1), F(1)), F(0), F(".5"))),
    ]:
        assert h.shape(data, s.DISK, h.Interval(0, 2)) == expected
