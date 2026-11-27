from fractions import Fraction as F

from .collar_family import interval, verify


def test_exact_source_transfer_and_area_enclosures():
    margins, signed, area = verify()
    assert len(margins) == 6
    measured = {
        "ordinary_energy": (".29756288", ".29756345"),
        "ordinary_entropy": (".40942660", ".40948019"),
        "pairing_entropy": (".23107941", ".23113299"),
    }
    for key, (lo, hi) in measured.items():
        lower, upper = interval(signed[key])
        assert F(lo) < lower <= upper < F(hi)
    for kind, lo, hi in (
        ("ordinary", "-.01772317", "-.01772084"),
        ("pairing", "-.01000393", "-.01000160"),
    ):
        lower, upper = area[kind]
        assert F(lo) < lower <= upper < F(hi)
