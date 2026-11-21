"""Directly measure the frequency-zero convention with an alternative local dimer."""
import argparse
import gzip
import json

import numpy as np

from .channel import FreeChannel, apply, matrices, preparation, standing
from .experiment import HERE, SHOTS, PatchClock, canonical, same_receipt
from .midpoint import direct_line


def run():
    with gzip.open(HERE/'midpoint.json.gz','rt') as fp:
        prior=json.load(fp)
    rng=np.random.default_rng(202611213033)
    clock=PatchClock(rng)
    car,cells,sites,_=preparation()
    d,o=matrices(car,sites)
    main=FreeChannel(d,o)
    refsites=[2,11]
    dr,off=matrices(car,refsites)
    dimer=FreeChannel(dr,off)
    ref=np.ones(2)/np.sqrt(2)
    assert np.linalg.norm(apply(dr,off,ref))<1e-12
    lines=[]
    for old in prior['lines']:
        j=round(old['q']*len(cells)/np.pi)
        wave,_=standing(len(cells),j*np.pi/len(cells),old['band'],d,o)
        line=direct_line(rng,clock,main,dimer,wave,ref)
        lines.append(dict(q=old['q'],band=old['band'],target=old['target'],
                          shift=line['signed_frequency']-old['direct'],**line))
    maxerr=max(abs(x['shift']+4) for x in lines)
    return dict(shots=SHOTS,hand_seams=prior['hand_seams'],clock_calibration=clock.calibration,
                reference_ports=refsites,parking=[4,11],bridge=[3,2],lines=lines,
                max_shift_error=maxerr,rest_signed_frequency=lines[0]['signed_frequency'],
                passes=bool(maxerr<.03))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    r=canonical(run())
    output=HERE/'reference_control.json.gz'
    if args.write:
        with output.open('xb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as f:
            f.write(json.dumps(r,separators=(',',':'),allow_nan=False).encode())
    else:
        with gzip.open(output,'rt') as f:
            assert same_receipt(r,json.load(f))
    print(r['max_shift_error'],r['rest_signed_frequency'],r['passes'])
