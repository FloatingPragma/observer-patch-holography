import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import screen_ward
from . import measure


def direct_currents():
    reg = measure.regional
    board = reg.bridge.native.iso.pairing.formation.federation(1)
    route = reg.bridge.native.iso.pairing.formation.region(board)
    rings = reg.bridge.beam.screens(board, route)
    positions = reg.bridge.beam.indices(sum(rings, []))
    states, _, _ = reg.source()
    amounts = []
    for charge in (-1, 1):
        total = 0
        for odd, weight in ((False, .8657), (True, .1343)):
            h = reg.bridge.native.generator(board, route, charge, odd)
            cells = reg.bridge.operators(board, route, charge, odd)
            cur = 1j*reg.bridge.current(board, cells, route[1:3])/8
            start = np.zeros((h.shape[0], 40))
            start[positions] = np.eye(40)
            op = -2j*h
            evolved = expm_multiply(op, start, traceA=op.diagonal().sum())
            effect = (evolved.conj().T @ (cur@evolved)).reshape(10,4,10,4)
            avg = np.einsum('ijik->jk', effect)/10
            total += weight*np.trace(avg@np.array(states[int(odd)], dtype=complex)).real
        amounts.append(total)
    return np.array([*amounts, sum(amounts)/2])


def test_sharp_optics_shared_clock_and_area():
    report = measure.run()
    held = json.loads((measure.HERE/'certificate.json').read_text())
    assert report == held
    ref = json.loads((measure.regional.HERE/'energy_odd.json').read_text())
    currents = direct_currents()
    for name, cur in zip(('vacancy','excess','fair'),currents,strict=True):
        assert cur in measure.parse_interval(report['rows'][name]['energy_current'])
    e, j = ([measure.parse_interval(r[key]) for r in report['rows'].values()]
            for key in ('expansion_odd', 'energy_current'))
    normal = [e[(i+1)%3]*j[(i+2)%3] - e[(i+2)%3]*j[(i+1)%3] for i in range(3)]
    floor = measure.parse_interval(report['absolute_odd_determinant'])/sum(map(abs, normal))
    assert floor.a > 0.1229
    for name in ('equal_work_control','fixed_work_control'):
        direct = screen_ward.run(ref[name])
        tables = [direct['charges'][q]['directional_images'] for q in ('-1','1')]
        tables.append(direct['fair_untagged_images'])
        jets = np.array([[x[k]['extrapolated_arrivals'] for k in ('forward', 'reverse')] for x in tables])
        odds = jets[:,0]-jets[:,1]
        det = np.linalg.det(np.column_stack((odds[:,0],currents,odds[:,1])))
        assert det in measure.parse_interval(report['absolute_odd_determinant'])
        cols = np.vstack([odds,jets[2,0]-(jets[0,0]+jets[1,0])/2])
        matrix = np.column_stack((cols[:,2],cols[:,1],cols[:,0],np.r_[currents,0]))
        u = np.linalg.solve(matrix, -cols[:,3])[0]
        assert u in measure.parse_interval(report['differential_mixture_required_clock'][name])
