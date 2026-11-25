from fractions import Fraction as F

from .coefficients import path_gap_at_readout, regional_pair_series, run
from .transport_cmi import integrated_deficit, native_path


def test_native_and_changed_law_at_depth_two():
    for strength, third_order in [(8, [F(145,96), F(145,96), F(17,96)]),
                                  (4, [F(41,96), F(41,96), F(3,32)])]:
        expected = [[-F(1,4), -F(1,4), -F(1+2*strength,4)], [0,0,0], third_order]
        for series in run(strength).values():
            assert series == expected


def test_disconnected_and_single_crossing_controls():
    bare = regional_pair_series([(1,8)], 1, (1,8), 8)
    assert bare[1:] == [[0,0,0], [0,0,0]]
    assert regional_pair_series([(0,5), (1,8)], 1, (1,8), 8) == bare
    assert regional_pair_series([(0,1), (1,8), (8,3)], 1, (8,3), 8) == bare


def test_two_contact_gap_at_fixed_temperature():
    assert path_gap_at_readout(8) > F(12, 10000)
    assert F(3, 10000) < path_gap_at_readout(4) < F(4, 10000)


def test_transport_bound_and_positive_resolvent_identity():
    for strength, beta in ((8, 1/16), (4, 1/16), (8, -1/16)):
        record = native_path(strength, beta, abs(beta))
        cmi, transfer, lower = record[4:]
        integral, error = integrated_deficit(*record[:4])
        assert 0 < transfer < 1
        assert 0 < lower < cmi
        assert abs(cmi - integral) < 2e-12
        assert error < 1e-12
