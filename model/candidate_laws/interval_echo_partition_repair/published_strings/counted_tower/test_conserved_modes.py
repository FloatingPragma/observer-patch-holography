"""Independent Fourier contractions, linearized curvature and held-out controls."""
from fractions import Fraction as F

import pytest
from sympy import Matrix, diag

from matched_einstein_limit.frame import PAIRS, design
from .. import experiment as survey
from .conserved_modes import longitudinal_prediction, reconstruct


def acquired(tensor):
    rays = survey._frames()
    rows = design([tuple(F(v).limit_denominator(3) for v in ray) for ray in rays.values()])
    components = Matrix([tensor[i, j] for i, j in PAIRS])
    return dict(zip(rays, Matrix(rows)*components, strict=True))


@pytest.mark.parametrize("c", [F(-1), F(1), F(1, 2), F(0)])
def test_six_shapes_recover_general_driven_amplitude_and_hold_out_z(c):
    tensor = Matrix([[3, 2, -4, 1], [2, -5, 6, 7], [-4, 6, 1, -3], [1, 7, -3, 4]])
    forcing = [c*tensor[0, i] + tensor[3, i] for i in range(4)]
    readings = acquired(tensor)
    expected = tuple(tensor[i, j] for i, j in PAIRS)
    assert reconstruct(readings, c, forcing) == (expected, 0)
    readings["z"] += 2
    assert reconstruct(readings, c, forcing) == (expected, 2)


def test_quarter_speed_has_two_invisible_conserved_modes():
    for jx, jy in [(1, 0), (0, 1)]:
        tensor = Matrix([[0, jx, jy, 0], [jx, -F(5, 2)*jx, -F(1, 4)*(jx+jy), -F(1, 4)*jx],
                         [jy, -F(1, 4)*(jx+jy), -F(5, 2)*jy, -F(1, 4)*jy],
                         [0, -F(1, 4)*jx, -F(1, 4)*jy, 0]])
        assert tensor.row(0)/4 + tensor.row(3) == Matrix.zeros(1, 4)
        readings = acquired(tensor)
        assert all(value == 0 for value in readings.values())
        with pytest.raises(ValueError, match="invisible"):
            reconstruct(readings, F(1, 4))


def test_directional_error_gain_for_unforced_luminal_mode():
    shapes = tuple(survey._frames())
    columns = [reconstruct({shape: int(shape == pulse) for shape in shapes}, 1)[0]
               for pulse in shapes]
    assert max(sum(map(abs, row)) for row in zip(*columns, strict=True)) == F(11, 6)


def test_moving_stream_and_shear_with_the_same_density_and_longitudinal_check():
    # Linearized Einstein tensor from the metric amplitude times cos(z-t/2).
    # The independent pressureless stream has T_ab = v_a v_b in these units.
    eta = diag(-1, 1, 1, 1)
    wave = Matrix([-F(1, 2), 0, 0, 1])
    raised = eta * wave
    norm = wave.dot(raised)

    def einstein(h):
        trace = (eta*h).trace()
        mixed = wave * (h*raised).T
        return (norm*h + trace*wave*wave.T - mixed - mixed.T
                + eta*(raised.dot(h*raised) - norm*trace))/2

    flow = Matrix([1, 0, 0, -F(1, 2)])
    stress = flow*flow.T
    metric = diag(F(5, 3), 1, 1, F(5, 3))
    metric[0, 3] = metric[3, 0] = -F(4, 3)
    material = acquired(stress)
    assert einstein(metric) == stress
    assert material["z"] == longitudinal_prediction(material["rest"], F(1, 2)) == 1

    metric[1, 2] = metric[2, 1] = 1
    geometry = acquired(einstein(metric))
    assert geometry["rest"] == material["rest"] and geometry["z"] == material["z"]
    assert geometry["xy"] - material["xy"] == 3
    # A nonzero two-direction minor excludes every common scalar coupling.
    assert geometry["xy"]*material["rest"] - geometry["rest"]*material["xy"] == 3
    coordinates, held = reconstruct(geometry, F(1, 2))
    assert held == 0 and coordinates[PAIRS.index((1, 2))] == F(3, 8)
