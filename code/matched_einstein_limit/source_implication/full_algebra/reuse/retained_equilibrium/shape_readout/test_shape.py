import json

from . import measure as m


def test_all_edges_and_flat_realizations():
    result = json.loads((m.HERE / "results.json").read_text())
    assert m.measure() == m.measure(60) == result
    edges = result["edge_log_secants"]
    assert len(edges) == 6 and min(x[1] for x in edges.values()) < 0 < max(
        x[0] for x in edges.values()
    )
    assert result["isotropic_max_error_lower"] > 0
    assert all(lo > 0 for arm in result["gram_principal_minors"] for lo, _ in arm)


def test_equal_clock_scale_gives_isotropic_readout():
    iv = m.iv
    ranges = [iv.mpf(i) for i in range(2, 14)]
    before = m.edge_logs(ranges)
    after = m.edge_logs([3 * x for x in ranges])
    assert all(0 in (b - a - iv.ln(3)) for a, b in zip(before, after, strict=True))


def test_orientation_filter_every_matrix_unit():
    from fractions import Fraction
    from .orientations import orientation_average

    for a in range(3):
        for b in range(3):
            unit = [[int((i, j) == (a, b)) for j in range(3)] for i in range(3)]
            assert orientation_average(unit) == tuple(
                tuple(Fraction(int(i == j and a == b), 3) for j in range(3))
                for i in range(3)
            )
