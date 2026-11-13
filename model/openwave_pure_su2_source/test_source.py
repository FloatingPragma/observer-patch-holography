from fractions import Fraction as F
from .source import certificate, cube

def test_full_cube_gap_and_noncollapse():
    c = certificate()
    assert c['vacuum_gap_interval_exact'] == ['21/2', '1024/47']
    assert c['plaquette_variance_lower_exact'] == '47/256'
    assert F(c['ground_orthogonal_weight_upper_exact']) < F(1, 1024)

def test_complete_graph_and_volume_obstruction():
    assert len(cube(3)['plaquettes']) == 36
    assert certificate(3)['strict_positive_gap_certificate']
    assert not certificate(4)['strict_positive_gap_certificate']
    assert certificate(spacing=F(1, 8))['vacuum_gap_interval_exact'] == ['84', '8192/47']
