"""Time-doubling test of a single-frequency fringe after ordinary marked capture."""
import argparse
import hashlib
import json

import numpy as np

from . import fringe
from . import fringe_walk as w

PROTO=fringe.HERE/"fringe_flight_protocol.json"
ARCHIVE=fringe.HERE/"fringe_flight_counts.npz"
SUMMARY=fringe.HERE/"fringe_flight_result.json"


def compose(left,right):
    a,da=left
    b,db=right
    return a@b,da@b+a@db


def velocities(k):
    """Eigenphase derivative from the special-unitary trace, for axial 0 < k < pi."""
    vectors=np.zeros((len(k),3))
    vectors[:,2]=k
    walk=w.BrightWalk(vectors)
    a,b=walk.a,walk.b
    da=-np.sin(k)/4
    db=-a*da/b
    df=np.empty_like(walk.routing)
    df[:,0,0],df[:,1,1]=da,-da
    df[:,0,1]=df[:,1,0]=db
    derivative=df*walk.coin
    cycle=walk.identity,np.zeros_like(derivative)
    for age in range(w.OPEN):
        theta=w.law.hop_angle(age,w.MEMBER.maturation)
        stay=np.exp(-1j*theta)*np.cos(theta)
        hop=1j*np.exp(-1j*theta)*np.sin(theta)
        cycle=compose((stay*walk.identity*walk.coin+hop*walk.mature,hop*derivative),cycle)
    count=w.PERIOD-w.OPEN
    block=walk.mature,derivative
    while count:
        if count%2:
            cycle=compose(block,cycle)
        block=compose(block,block)
        count//=2
    # det(cycle)=-i. Removing its common phase gives an SU(2) matrix.
    T,dT=[np.exp(1j*np.pi/4)*part for part in cycle]
    sin=np.hypot(T[:,0,0].imag,abs(T[:,1,0]))
    w.require(np.all(sin>1e-12),"quadrature node at an unresolved crossing")
    v=(dT[:,0,0]+dT[:,1,1]).real/(2*w.PERIOD*sin)
    w.require(np.all(np.isfinite(v)) and max(abs(v))<.500001,"invalid axial velocity")
    return v


def weak_curve(count, times):
    accum=np.zeros(len(times))
    for first in range(0,count,4096):
        k=np.pi*(np.arange(first,min(first+4096,count))+.5)/count
        v=velocities(k)
        accum+=np.cos(2*np.pi*times[:,None]*v).sum(axis=1)
    return accum/count


def predictions(config):
    for n in config["sizes"]:
        k=w.angles((1,1,n))
        walk=w.BrightWalk(k)
        analyzer=np.exp(-2j*np.pi*np.arange(n)/n)
        for t in config["time_ratios"]:
            ticks=int(n*t)
            w.require(ticks==n*t and ticks%2==0,"not a local optical return")
            amplitudes=np.fft.ifft(walk.wave(ticks)[0,0],axis=0)
            mass=(abs(amplitudes)**2).sum(axis=-1)
            w.require(abs(mass.sum()-1)<3e-9,"incomplete marked kernel")
            cross=mass@analyzer
            hit=(1+(fringe.PHASES*cross).real)/(2*n)
            p=np.column_stack((hit,.5-hit,np.full(4,.5)))
            yield f"{n}_{t}",n,p


def run(record=False):
    fringe.plan()
    pins=json.loads((fringe.HERE/"fringe_flight_sources.json").read_text(encoding="utf8"))
    names=("fringe_flight_protocol.json","fringe_flight.py","fringe_sources.json")
    w.require(set(pins)==set(names) and all(
        pins[n]==hashlib.sha256((fringe.HERE/n).read_bytes().replace(b"\r\n",b"\n")).hexdigest()
        for n in names),"changed flight source")
    cfg=json.loads(PROTO.read_text(encoding="utf8"))
    w.optical_reference()
    cases=list(predictions(cfg))
    if record:
        rng=np.random.default_rng(cfg["seed"])
        counts={k:fringe.sample(p,n*cfg["shots_per_site"],rng) for k,n,p in cases}
        with ARCHIVE.open("xb") as stream:
            np.savez_compressed(stream,**counts)
    times=np.array(cfg["time_ratios"])
    coarse,limit=[weak_curve(n,times) for n in cfg["quadrature_nodes"]]
    L=np.log(2*12*len(cases)/cfg["alpha"])
    rows={}
    with np.load(ARCHIVE,allow_pickle=False) as archive:
        w.require(sorted(archive.files)==sorted(k for k,_,_ in cases),"incomplete time-doubling menu")
        for key,n,p in cases:
            shots=n*cfg["shots_per_site"]
            data=archive[key]
            interval=fringe.read(data,shots,L)
            lo,hi=fringe.count_box(data/shots,shots,L)
            w.require(np.all((lo<=p)&(p<=hi)),"finite marked flight disagrees with counts")
            rows[key]=None if interval is None else interval[0]
    tests={}
    for n in cfg["sizes"]:
        for t,dual in zip(cfg["time_ratios"][:-1],cfg["time_ratios"][1:],strict=True):
            first,second=[rows[f"{n}_{v}"] for v in (t,dual)]
            w.require(dual==2*t,"intervals are not doubled")
            tests[f"{n}_{t}"]=(None if first is None or second is None else
                                round(float(second[0]+1-2*max(abs(first))**2),7))
    last=[rows[f"{max(cfg['sizes'])}_{v}"] for v in cfg["time_ratios"]]
    offset=(None if any(x is None for x in last) else
            max(float(max(abs(bounds-value))) for bounds,value in zip(last,limit,strict=True)))
    change=float(max(abs(coarse-limit)))
    result=dict(contrast_intervals={k:None if x is None else x.round(7).tolist() for k,x in rows.items()},
                excess_lower=tests,weak_curve=limit.round(7).tolist(),grid_change=round(change,8),
                finest_interval_distance=None if offset is None else round(offset,7),
                weak_screen=bool(offset is not None and offset<.04 and change<.003),
                panels=4*len(cases),draws=sum(4*n*cfg["shots_per_site"] for _,n,_ in cases),
                sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest())
    if record:
        with SUMMARY.open("x",encoding="utf8",newline="\n") as stream:
            stream.write(json.dumps(result,indent=1,sort_keys=True,allow_nan=False)+"\n")
    w.require(result==json.loads(SUMMARY.read_text(encoding="utf8")),"changed flight receipt")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true")
    summary=run(parser.parse_args().record)
    print("Lower bounds on the time-doubling excess",summary["excess_lower"])
