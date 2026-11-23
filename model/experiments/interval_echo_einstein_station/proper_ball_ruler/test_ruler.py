import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from . import coarse, free_shapes, measure, route_audit

HERE = Path(__file__).resolve().parent


def test_range_decoder_controls_and_material_route_replay():
    pts = np.array([[0,0,0],[1,0,0],[.5,np.sqrt(3)/2,0],[.5,np.sqrt(3)/6,np.sqrt(2/3)]])
    quads = np.array(list(combinations(range(4),4)))
    for scale in (1., .5):
        d=cdist(pts*scale,pts*scale)
        v,a,err=measure.measured_hull(d)
        np.testing.assert_allclose([v,a,err], [np.sqrt(2)/12*scale**3,np.sqrt(3)*scale**2,0.],atol=1e-12)
        assert measure.quartet_bound(d,quads)[0] == 0.
    cycle=np.array([[0,1,2,1],[1,0,1,2],[2,1,0,1],[1,2,1,0]])
    assert measure.quartet_bound(cycle,quads)[0]>.17
    assert route_audit.run() == json.loads((HERE/'route_audit.json').read_text())


def test_heldball_growth_and_free_probabilities():
    large = json.loads((HERE/'coarse.json').read_text())['panels']
    for a,b in zip(large,large[1:]):
        assert a['relative_ball_radius'] == 2*b['relative_ball_radius']
        assert a['microscopic_ball_radius']*2 == b['microscopic_ball_radius']
    for p in large:
        n=np.array(p['launch_classes'])
        assert n.sum()==coarse.N_TRIALS
        assert np.isclose(np.log2(n[:2].sum()/n[0]),p['sharp_sample_dimension'])
        assert p['doubling_dimension'] < 2.1
    native=json.loads((HERE/'results.json').read_text())['cases']
    for r in native:
        m=r['fixed_volume_experiments']['both_geometry_fixed']
        assert max(abs(x) for x in m['population'][:2]) < 1e-11
        assert abs(m['population'][2]) > .56
        assert r['all_resident_words_distortion_lower_bound_from_intact_cell'] > .14
    free=json.loads((HERE/'free_shapes.json').read_text())['rows']
    assert len(free)==12
    repeated=free_shapes.observe(native[0],0)
    for field in ('exact_quantum_interior_entropy_slope','exact_sharp_interior_entropy_slope'):
        np.testing.assert_allclose(free[0][field],repeated[field],atol=1e-12)
    for r in free:
        assert max(abs(x) for x in r['exact_volume_area_slopes'])<1e-11
        assert r['max_cut_word_probability_slope']<1e-13
        assert r['exact_quantum_interior_entropy_slope']>1.14
        assert r['exact_sharp_interior_entropy_slope']>.08
        assert abs(r['exact_quantum_interior_entropy_slope']-r['coherent_interior_observations']['signed_slopes'][0])<.001
