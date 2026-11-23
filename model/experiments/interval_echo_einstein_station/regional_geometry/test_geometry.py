import json

import numpy as np

from . import measure


def test_known_tetrahedron_with_independent_held_stations():
    xyz = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.1,.2,.3],[.2,.3,.15]])
    edges = np.array([np.linalg.norm(xyz[i]-xyz[j]) for i,j in measure.PAIRS])
    # All mixtures have the same event probability. The first component
    # scales all distance readouts at unit speed through mixture weight.
    word = np.zeros((3,3,len(edges),2))
    word[:,:,:,0] = .5
    word[:,:,:,1] = .5*edges
    word[0,0,:,1] = edges
    ref, tangent = np.zeros((3,3)), np.zeros((3,3))
    ref[0,1] = 1
    tangent[0,0], tangent[0,1] = 1,-1
    result = measure.geometry(word,ref,tangent)
    for key, expected in (("volume",1/6),("surface_area",(3+np.sqrt(3))/2),
                          ("dlog_volume",3),("dlog_area",2),("held_relative_error",0),
                          ("anchor_radius_relative_error",0)):
        np.testing.assert_allclose(result[key], expected, atol=2e-14, rtol=0)


def test_occupation_controls_and_transferred_geometry():
    data = json.loads((measure.HERE/"results.json").read_text())
    ref = np.asarray(data["reference"])
    coeff = np.asarray(data["transferred_coefficients"])
    tangent = np.einsum("a,aij->ij", coeff,measure.DIRECTIONS)
    for sign in (-1.,1.):
        state = ref+sign*measure.EPS*tangent
        assert np.all(state>0)
        np.testing.assert_allclose(state.sum(axis=0),1/3)
        np.testing.assert_allclose(state.sum(axis=1),1/3)
    for board in data["screens"]:
        assert board["blind_control_maximum"]<1e-12
        assert abs(board["local"]["first_return"]["dlog_volume"])<1e-12
        if board["projector_renewal_max_error"] is not None:
            assert board["projector_renewal_max_error"] < 1e-11
            assert board["independent_area_volume_slope_disagreement"] < 1e-12
        for result in board["transferred"].values():
            assert result["held_relative_error"] > .05
            for derivative in ("dlog_volume","dlog_area"):
                assert abs(result[derivative]-result["central_"+derivative])<1e-11
        assert abs(board["direct_entropy_slope"]-board["measured_source_charge_slope"]) < 5*board["source_slope_standard_error"]
