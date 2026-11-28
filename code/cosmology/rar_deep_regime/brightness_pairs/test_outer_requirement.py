"""Check the reported astronomical requirement and reduced rational contradictions."""

from fractions import Fraction as F
import json
import math

from . import shared_outer_order as order


def test_ngc3741_lower_bound_from_catalogue():
    sample = order.sample
    reads = json.loads((sample.HERE / "outer_requirement.json").read_bytes())
    rows = reads["galaxies"]
    data = sample.measure.DATA
    strongest = max(
        (r["total_over_baryonic"][0], n, r["line"])
        for n, z in rows.items()
        for r in z["rings"]
    )
    minimum, name, line = strongest
    assert name == "NGC3741" and line == 1129
    ring = (data / "table2.dat").read_text().splitlines()[line - 1]
    galaxy = next(
        s
        for s in (data / "table1.dat").read_text().splitlines()
        if s[:11].strip() == name
    )
    assert ring[:11].strip() == name and float(ring[19:25]) == 7
    d, de, inc, ie = [
        float(galaxy[i:j]) for i, j in ((15, 21), (22, 27), (30, 34), (35, 39))
    ]
    v, e, gas, disk, bulge = [
        float(ring[i:j]) for i, j in ((26, 32), (33, 38), (39, 45), (46, 52), (53, 59))
    ]
    gmax = max(
        gas * abs(gas) + x * disk * abs(disk) + y * bulge * abs(bulge)
        for x in (0.1, 1.0)
        for y in (0.1, 1.4)
    )
    hmax = (1 + 2 * de / d) * (
        math.sin(math.radians(min(90, inc + 2 * ie))) / math.sin(math.radians(inc))
    ) ** 2
    assert math.isclose((v - 2 * e) ** 2 / (hmax * gmax), minimum, rel_tol=2e-14)
    assert 8.7082 < minimum < 8.7083
    assert (
        sum(r["total_over_baryonic"][0] > 1 for z in rows.values() for r in z["rings"])
        == 816
    )


def test_exact_reduced_certificates(monkeypatch):
    monkeypatch.setattr(order.sample.iv, "dps", 36)
    data = json.loads((order.HERE / "outer_requirement.json").read_bytes())["galaxies"]
    archived = json.loads((order.HERE / "shared_outer_order.json").read_bytes())
    for group in ("calibration", "held"):
        certificate = archived["groups"][group]["factor_two"]
        support = {tuple(pair) for pair in certificate["support_catalogue_rings"]}
        reduced = {
            n: {"rings": [r for r in galaxy["rings"] if (n, r["line"]) in support]}
            for n, galaxy in data.items()
            if any(name == n for name, _ in support)
        }
        rows, upper, labels = order.constraints([group], reduced, "factor_two")
        assert set(labels.values()) == support
        saved = certificate["isolated_witness"]
        residual = [F(0) for _ in upper]
        bound = F(0)
        for index, num, den in saved["certificate_rows"]:
            assert num > 0 and den > 0
            coefficients, rhs = rows[index]
            bound += F(num, den) * rhs
            for col, coefficient in coefficients.items():
                residual[col] += F(num, den) * coefficient
        minimum = sum(min(F(0), v * u) for v, u in zip(residual, upper, strict=True))
        assert minimum > bound and minimum - bound == F(*saved["exact_gap"])
        if group == "calibration":
            assert len(labels) == 6 and len(reduced) == 3
