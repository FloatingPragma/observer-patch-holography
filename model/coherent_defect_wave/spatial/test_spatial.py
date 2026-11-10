"""Executable local constraints, complete receiver records and literal bands."""
import pytest
import sympy as s

from ..model import JointObserver
from .certificate import pair_control,band_control,cone_control
from .family import PairState,CubicFamily,dressed_swap,exchange,gamma,frame_sign
from .family import face_measurement,charge_measurement
from .dirac_route import certificate as dirac_certificate
from .dirac_route import route_step


def test_connected_gauss_preserving_receiver_protocol():
    result = pair_control()
    assert result['full_receiver_face_distributions'] == [
        {'0':'301/600','2':'1/1800','5':'1/1800','16':'0','17':'179/360'},
        {'0':'1','2':'0','5':'0','16':'0','17':'0'}]


def test_frame_intertwining_for_every_native_coordinate():
    observer = JointObserver()
    star = observer.law.vertex_stars[5]
    for g in observer.link_words:
        for port in range(12):
            for internal in range(5):
                key = (g,port,internal,0,0)
                partner = (g^star,port,internal,0,0)
                character = -1 if port == 5 and internal >= 3 else 1
                assert frame_sign(key)*frame_sign(partner)*character == 1


def test_exchange_exact_on_complex_superpositions_and_all_gauss_characters():
    state = PairState({((1,5,3,0,0),(3,6,0,1,8)):1,
                       ((6553635,6,4,1,59),(1,5,3,0,4)):s.I,
                       ((3,0,0,0,7),(6553633,5,4,1,18)):2},6)
    result = exchange(state)
    assert state.norm_squared() == result.norm_squared() == 1
    assert exchange(result,inverse=True) == state
    assert exchange(result) == dressed_swap(state)
    for cell in (0,1):
        assert exchange(gamma(state,cell)) == gamma(result,cell)
    before = {}
    after = {}
    for target,source in ((before,state),(after,result)):
        for keys,value in source.amplitudes:
            charge = sum(2 if key[2] < 3 else -3 for key in keys)
            target[charge] = s.expand(target.get(charge,0)+s.conjugate(value)*value/6)
    assert before == after


def test_literal_torus_bands_and_scoped_nonrelativistic_diagnostic():
    result = band_control()
    assert result['literal_circuit_bloch_column_checks'] == 16
    assert not result['relativistic_mass_shell_certified']


def test_separate_literal_Clifford_route_has_controlled_free_limit():
    result = dirac_certificate()
    assert result['exact_Clifford_pair_checks'] == 16
    assert result['literal_position_space_Bloch_column_checks'] == 4
    assert not result['native_U_invariant_subspace_implements_this_routing']
    assert not result['input_mass_selected_by_repair']


def test_exact_volume_independent_circuit_cone():
    result = cone_control()
    assert result['exact_backwards_cone_sizes'] == [1,64,512]
    assert not result['continuum_spacing_limit_established']


def test_empty_branches_and_rejected_uncontrolled_inputs():
    assert PairState({}).norm_squared() == s.Integer(0)
    assert not PairState({}).norm_squared().has(s.Float)
    with pytest.raises(ValueError):
        PairState({((1,0,0,0,0),(1,0,0,0,0)):0.5})
    with pytest.raises(ValueError):
        PairState({((0,0,0,0,0),(1,0,0,0,0)):1})
    with pytest.raises(ValueError):
        CubicFamily(3,periodic=True)
    with pytest.raises(ValueError):
        CubicFamily(4).cone({(4,0,0)},1)
    assert CubicFamily(4).cone(iter([(1,1,1)]),0) == {(1,1,1)}
    with pytest.raises(ValueError):
        CubicFamily(4).one_charge_step({(0,0,0):0.5})
    for bad in (s.oo,s.Symbol('amplitude'),s.nan):
        with pytest.raises(ValueError):
            PairState({((1,0,0,0,0),(1,0,0,0,0)):bad})
    for readout in (face_measurement,charge_measurement):
        for bad in (-1,2,True):
            with pytest.raises(ValueError):
                readout(PairState({}),bad)
    for angle in (s.I,s.oo,s.Symbol('angle'),0.25):
        with pytest.raises(ValueError):
            route_step({((0,0,0),0):1},4,angle)
