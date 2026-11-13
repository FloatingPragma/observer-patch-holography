from fractions import Fraction as F
import sympy as sp
from . import source as s


def test_complete_transport_budget_including_empty_halo_records():
    for seed in (1, 3):
        tensor = s.staggered.fixture(seed)
        bound, actual, _ = s.exact_transport_check(tensor)
        assert 0 < actual <= bound


def test_constant_metric_component_is_retained_and_null_invisible():
    tensor = {(a,b): {} for a in range(4) for b in range(4)}
    for a in range(4):
        tensor[a,a][(0,0,0,0)] = F(7*s.staggered.ETA[a])
    C = s.colocate(tensor, (0,0,0,0))
    assert C == 7*s.decoder.ETA
    assert s.decoder.charges(C) == sp.zeros(9,1)


def test_checkerboard_cannot_be_discarded_as_zero_tensor():
    control = s.kernel_control()
    assert control['original_maximum_absolute_value_exact'] == '1'
    assert set(control['all_interpolated_values_exact']) == {'0'}


def test_component_budget_is_sharp_for_affine_data():
    tensor = {(a,b): {} for a in range(4) for b in range(4)}
    field = {(2*i+1,2*j+1,0,0): F(2*i+1+2*j+1)
             for i in range(-2,3) for j in range(-2,3)}
    tensor[0,1] = tensor[1,0] = field
    # The four corner mean is exactly zero; the selected corner is two.
    assert s.colocate(tensor, (0,0,0,0))[0,1] == 0
    assert field[1,1,0,0] == 2
