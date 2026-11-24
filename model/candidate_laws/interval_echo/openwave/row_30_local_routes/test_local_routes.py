import json

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_local_routes.experiment import (
    HERE, gaps, grow, native_routing, shortest_internal,
)
from model.executable_model import federation
from model.executable_model.carrier import cell, icosahedron_from_coordinates, isomorphism


def test_graph_turns_match_coordinate_angles():
    car = cell()
    xy, seams = icosahedron_from_coordinates()
    xy = xy[isomorphism(car.seams,seams)]
    xy /= np.linalg.norm(xy,axis=1)[:,None]
    cosines = (1.,1/np.sqrt(5),-1/np.sqrt(5),-1.)
    for a in range(12):
        for b in range(12):
            distance = len(shortest_internal(car,a,b))-1
            np.testing.assert_allclose(xy[a]@xy[b],cosines[distance],atol=1e-12)


def test_straight_and_farthest_are_different_physical_preparations():
    car=federation(0)
    choices,across=native_routing(car)
    for entry in choices:
        assert grow(choices,across,entry,'antipodal',False)[0]==[]
        for high in (True,False):
            route,_=grow(choices,across,entry,'farthest',high)
            assert 8<=len(route)<=14
            assert all(tuple(car.degree[p])==(6,5,6) for p in route)


def test_junction_screening_and_actual_receivers():
    r=json.loads((HERE/'junction.json').read_text())
    filtered=json.loads((HERE/'filter.json').read_text())
    np.testing.assert_allclose(gaps(),[(-12,-8),(-4,-6+2*np.sqrt(17))])
    for size in (20,80,320):
        for row in r[str(size)].values():
            assert sum(row['length_histogram'].values())==3*size
            pure=row['longest_whole_trimer']
            spec=pure['native']['standing_spectrum']
            assert spec['distinct_matches']==3*(len(pure['cells'])-1)
            assert spec['max_nearest_line_error']<1e-9
            assert max(e['inside'] for e in pure['native']['interior_gaps'])==0
    for line in filtered['rows']:
        assert sum(line['all_accepted_cells'])+line['rejected']==filtered['shots']
        if line['selector']=='longest_whole_trimer':
            if line['frequency'] in (-10,-1):
                assert line['accepted']==0
            else:
                assert line['receiver_clicks']>0
