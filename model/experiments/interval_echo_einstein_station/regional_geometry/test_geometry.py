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


def test_both_rulers_at_fixed_volume_and_one_fixed_area():
    from itertools import combinations
    from math import factorial

    import mpmath as mp

    from . import both_volumes

    bank = json.loads((measure.HERE / "results.json").read_text())
    analysis = both_volumes.run()
    ref = np.array(bank["reference"])
    for delta in np.array(analysis["twenty_cell_transfers"]):
        np.testing.assert_allclose(delta.sum(axis=0), 0, atol=1e-14)
        np.testing.assert_allclose(delta.sum(axis=1), 0, atol=1e-14)
        assert np.min(ref - abs(delta) / 32) > 0
    # Independent Cayley-Menger finite differences on the archived word laws.
    def content(distances, sites):
        n = len(sites) - 1
        mat = mp.matrix(n + 2)
        for a, i in enumerate(sites, 1):
            mat[a, 0] = mat[0, a] = 1
            for b, j in enumerate(sites[:a-1], 1):
                mat[a, b] = mat[b, a] = distances[(j, i)] ** 2
        return mp.sqrt((-1)**(n+1) * mp.det(mat) / (2**n * factorial(n)**2))

    with mp.workdps(45):
        step = mp.mpf('1e-10')
        for raw, screen in zip(bank["screens"], analysis["boards"], strict=True):
            laws = np.array(raw["word_moments"])
            for chosen, row in enumerate(screen["records"]):
                delta = np.array(analysis["twenty_cell_transfers"][chosen])
                independent = []
                for menu in range(2):
                    shapes = []
                    for sign in (-1, 1):
                        weights = [mp.mpf(float(x)) + sign*step*mp.mpf(float(d))
                                   for x, d in zip(ref.ravel(), delta.ravel(), strict=True)]
                        distances = {}
                        for pair, events in zip(measure.PAIRS, laws[:,:,:,menu,:].transpose(2,0,1,3), strict=True):
                            chances, hand = events.reshape(9, 2).T
                            distances[pair] = mp.fdot(weights, hand) / mp.fdot(weights, chances)
                        shapes.append([content(distances, range(4)),
                                       sum(content(distances, f) for f in combinations(range(4), 3))])
                    independent.append([float(mp.log(b/a)/(2*step))
                                        for a, b in zip(*shapes, strict=True)])
                slopes = np.array(independent).T.ravel()
                np.testing.assert_allclose(slopes, row["transferred_derivatives"], rtol=0, atol=3e-14)
                assert max(abs(slopes[[0, 1, 2 + chosen]])) < 2e-13
                assert row["transferred_entropy_slope"] > 0.14
                assert abs(slopes[3 - chosen]) > 8e-9


def test_restricted_entropy_maximum_does_not_fix_either_rim():
    from . import both_volumes

    ref = np.full((3, 3), 1/9)
    for board in both_volumes.run()["boards"]:
        for chosen, row in enumerate(board["uniform_reference"]):
            delta = np.array(row["transfer"])
            np.testing.assert_allclose(delta.sum(axis=0), 0, atol=2e-15)
            np.testing.assert_allclose(delta.sum(axis=1), 0, atol=2e-15)
            assert np.min(ref - abs(delta) / 32) > 0
            assert abs(np.sum(delta * np.log(ref))) < 2e-15
            slopes = row["derivatives"]
            assert max(abs(x) for x in slopes[:2]) < 2e-14
            assert slopes[2 + chosen] > 1e-7
