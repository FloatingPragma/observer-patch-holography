"""Massive source mean-current continuity on the same retained field."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import joint_program as joined
import verified

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def positions(state):
    answer=np.zeros(96)
    for (x,field),v in state.items():answer[x]+=float(np.vdot(v,v).real)
    return answer


def run():
    pins=verified.verify_sources()
    path=Path(__file__);pins[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    graph,B,A,labels=joined.massive.af.data()
    neighbours,edges=joined.massive.af.mb.common.geometry_tables()
    rows=[];worst=0.
    for mass in (F(0),F(1)):
        state=joined.massive.af.initial();mu=float(mass/8)
        for step in range(3):
            current=np.zeros(288)
            for (x,field),v in state.items():
                mean=np.sum(v)/12
                kicked=np.exp(1j*mu)*v+(np.exp(-1j*mu)-np.exp(1j*mu))*mean
                reflected=np.sum(kicked)/6-kicked
                for direction in range(6):
                    edge=int(edges[x,direction]);sign=1 if x==graph.seams[edge][0] else -1
                    current[edge]+=sign*abs(reflected[2*direction+1])**2
            following=joined.massive.step(state,F(1,8),mass,.17,.41)
            discrepancy=positions(following)-positions(state)+B@current
            error=float(np.max(abs(discrepancy)));worst=max(worst,error)
            assert error<2e-12
            rows.append({'mass':str(mass),'step':step,'current288':current.tolist(),'continuity_residual':error})
            state=following
    # A wrong-current law gives a positive discrepancy on a fully admitted
    # initial source. This is an operator/record identity with falsifying data.
    before=joined.massive.af.initial()
    after=joined.massive.step(before,F(1,8),F(1),.17,.41)
    broken=float(np.max(abs(positions(after)-positions(before))))
    assert broken>.49
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in pins.items())
    result={'schema':'oph.common.massive_current.opposing_review.v1',
            'checks':{'same_mass_field_source':True,'all288_signed_current_coefficients':True,
                      'ordinary_mean_continuity_not_modular_average':True,
                      'omitted_current_law_rejected':True,'source_pins_unchanged':True},
            'records':rows,'maximum_continuity_residual':worst,
            'omitted_current_residual':broken,
            'instrument':'Native signed moving-slot current after the same on-site mass/reflection; the joint event/current record is an explicit additional intervention. It is not inferred by differencing centered modular field means across wrap boundaries.',
            'scope':'One-particle finite common source. A mean continuity identity is not a two-charge force, homogeneous field limit or species identification.',
            'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'CURRENT_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'maximum_residual':worst,'omitted_current_residual':broken}))


if __name__=='__main__':run()
