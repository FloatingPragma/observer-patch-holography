"""A one-station comparison with the free scalar on the same grating and counter.

Registered after the paired limits, before this independent acquisition.
Both apertures count every slot.  At one quarter and one half of a chart
period the unboosted harmonics are in the first half-turn.  Reject a common
cone only if the simultaneous intervals exclude ratio one at both returns.
"""
import argparse
import hashlib
import json

import numpy as np

from . import fringe
from . import fringe_resonance as resonance
from . import fringe_walk as w

N=4096*w.PERIOD
TIMES=(.25,.5)
MENUS=("free","resonant")
SHOTS=N*2**22
ALPHA=.01
SEED=921634
ARCHIVE=fringe.HERE/"fringe_reference_counts.npz"
SUMMARY=fringe.HERE/"fringe_reference_result.json"


def panels():
    k=2*np.pi/N
    bound=w.BrightWalk(np.array([[0,0,0],[0,0,k]]),mu=resonance.RES)
    for t in TIMES:
        ticks=int(t*N)
        states=[np.linalg.matrix_power(w.law.walk_matrix([0,0,z,0]),ticks)@w.law.EQUAL for z in (0,k)]
        for menu in MENUS:
            pair=states if menu=="free" else bound.wave(ticks)
            cross=np.vdot(*pair)
            p=np.clip((1+(fringe.PHASES*cross).real)/N,0,1)
            yield f"{menu}_{t}",p


def run(record=False):
    fringe.plan()
    w.optical_reference()
    forecasts=list(panels())
    if record:
        rng=np.random.default_rng(SEED)
        raw={}
        for key,p in forecasts:
            enrolled=rng.binomial(SHOTS,.5,size=4)
            hits=rng.binomial(enrolled,p)
            raw[key]=np.column_stack((hits,enrolled-hits,SHOTS-enrolled)).astype(np.int64)
        with ARCHIVE.open("xb") as stream:
            np.savez_compressed(stream,**raw)
    L=np.log(2*12*len(forecasts)/ALPHA)
    phase={}
    with np.load(ARCHIVE,allow_pickle=False) as raw:
        w.require(sorted(raw.files)==sorted(k for k,_ in forecasts),"incomplete scalar comparison")
        for key,p in forecasts:
            counts=raw[key]
            interval=fringe.read(counts,SHOTS,L)
            lo,hi=fringe.count_box(counts/SHOTS,SHOTS,L)
            expected=np.column_stack((p/2,(1-p)/2,np.full(4,.5)))
            w.require(np.all((lo<=expected)&(expected<=hi)) and interval is not None,
                      "unresolved or inconsistent scalar comparison")
            phase[key]=np.arccos(interval[0,::-1])
    ratios={}
    for t in TIMES:
        scalar=phase[f"free_{t}"]
        w.require(scalar[0]>0,"unresolved scalar phase")
        ratios[str(t)]=(phase[f"resonant_{t}"]/scalar[::-1]).round(7).tolist()
    result=dict(phase_ratio=ratios,common_rate_rejected=all(x[1]<1 for x in ratios.values()),
                panels=4*len(forecasts),draws=4*len(forecasts)*SHOTS,
                sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest())
    if record:
        with SUMMARY.open("x",encoding="utf8",newline="\n") as stream:
            stream.write(json.dumps(result,indent=1,sort_keys=True)+"\n")
    w.require(result==json.loads(SUMMARY.read_text(encoding="utf8")),"changed scalar comparison")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true")
    print(run(parser.parse_args().record)["phase_ratio"])
