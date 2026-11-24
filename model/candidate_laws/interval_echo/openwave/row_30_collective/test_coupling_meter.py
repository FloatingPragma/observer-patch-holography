import json

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_collective.coupling_meter import (
    HERE, HAND_TRIALS, local_generators, measure,
)
from model.executable_model import federation


def test_global_probe_geometry_and_local_column_restrictions():
    r=json.loads((HERE/'coupling_meter.json').read_text())
    edges=r['addresses']
    probe,storage=edges['probe_dimer']
    car=federation(0)
    assert car.degree[probe]==car.degree[storage]
    for tag in ('condition_link','probe_phase_service'):
        u,v=edges[tag]
        assert u==probe and v!=storage and v in car.neighbours[u]
    lab=set(edges['probe_dimer']+edges['condition_link']+edges['probe_phase_service'])
    hand=set(sum(edges['hand_seams'],[]))
    assert len(hand)==4 and not lab.intersection(hand)
    h,_,residuals=local_generators()
    assert max(residuals.values())<1e-12
    assert np.allclose(np.diff(np.linalg.eigvalsh(h['exchange'])),1)
    assert np.allclose(np.diag(h['drawn_link'])-h['drawn_link'][1,1],[32,0])
    for control in ('empty_link','vacant_neighbour'):
        assert np.allclose(h[control],np.eye(2)*h[control][0,0])


def test_microscopic_forecast_from_counts_only():
    r=json.loads((HERE/'coupling_meter.json').read_text())
    counts=np.asarray(r['all_counts'])
    assert np.all(counts.sum(axis=-1)==r['trials'])
    v,step=measure(counts,np.array(r['hand_calibration'])/HAND_TRIALS)
    assert step<np.pi*.9
    estimates=r['results']
    for i,key in enumerate(['h','J','predicted_m','predicted_c_squared','exchange',
                            'drawn_link','empty_link','vacant_neighbour']):
        assert abs(v[i]-estimates[key]['value'])<1e-9
    for ref,key in [(15,'predicted_m'),(16,'predicted_c_squared'),(0,'empty_link'),
                    (0,'vacant_neighbour')]:
        lo,hi=estimates[key]['bootstrap95']
        assert lo<ref<hi
