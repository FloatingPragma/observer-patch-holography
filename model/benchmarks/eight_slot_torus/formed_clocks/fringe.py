"""One-station fringe visibility after native capture and after coherent loading."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import fringe_walk as walk
from .formed_clock_readout import ordinary_birth
from .fringe_walk import require

HERE = Path(__file__).resolve().parent
SETTINGS = HERE/"fringe_protocol.json"
COUNTS = HERE/"fringe_counts.npz"
SUMMARY = HERE/"fringe_result.json"
SOURCES = ("fringe_protocol.json","fringe_resonance_protocol.json","fringe_resonance.py",
           "fringe_walk.py","fringe.py","../law.py","../formation.py",
           "formed_clock_readout.py","clock_comparison.py")
PHASES = np.array([1,1j,-1,-1j])


def plan():
    pins = json.loads((HERE/"fringe_sources.json").read_text(encoding="utf8"))
    require(set(pins)==set(SOURCES),"incomplete fringe source register")
    for name in SOURCES:
        data=(HERE/name).read_bytes().replace(b"\r\n",b"\n")
        require(hashlib.sha256(data).hexdigest()==pins[name],f"changed source: {name}")
    return json.loads(SETTINGS.read_text(encoding="utf8"))


def forecasts(cfg):
    walk.optical_reference()
    refused, source = ordinary_birth(5,0.)
    admission = float(np.vdot(source,source).real)
    require(np.allclose(refused,[0,0,.5],rtol=0,atol=3e-13) and abs(admission-.5)<3e-13,
            "incorrect heralded source")
    refused_off, off = ordinary_birth(5,0.,vertex=False)
    require(np.allclose(refused_off,[0,0,1],rtol=0,atol=3e-13) and not np.any(off),
            "vertex-off source formed a pair")
    for geometry in ("strip","cube"):
        for size in cfg[geometry+"_sizes"]:
            shape = (1,1,size) if geometry=="strip" else (size,)*3
            volume = math.prod(shape)
            for order in cfg[geometry+"_orders"]:
                prediction = walk.contrasts(shape,order,cfg["boosts"],cfg["windows"])
                for window in cfg["windows"]:
                    for boost in cfg["boosts"]:
                        marked, pure = prediction[window,boost]
                        for menu in cfg["menus"]:
                            born = 0. if menu=="vertex_off" else .5
                            cross = pure if menu=="coherent" else marked
                            local = born*(1+(PHASES*cross).real)/volume
                            require(np.min(local)>=-1e-11 and max(local)<=born+1e-11,
                                    "invalid local fringe probability")
                            local=np.clip(local,0,born)
                            p=np.column_stack((local,born-local,np.full(4,1-born)))
                            key=f"{geometry}_{size}_{''.join(map(str,order))}_{window}_{boost}_{menu}"
                            yield key, volume, p


def sample(p, shots, rng):
    # Bound each multinomial draw for platforms with 32-bit C long.
    counts = np.zeros(p.shape,dtype=np.int64)
    while shots:
        block = min(shots,2**30)
        for i,row in enumerate(p):
            counts[i]+=rng.multinomial(block,row)
        shots-=block
    return counts


def registered_counts(data,shots):
    require(isinstance(data,np.ndarray) and data.shape==(4,3) and data.dtype==np.int64
            and type(shots) is int and 0<shots<=2**60 and np.all((0<=data)&(data<=shots))
            and np.all(data.sum(axis=1)==shots),"invalid or missing fringe outcomes")


def count_box(f,shots,logtails):
    # Invert the two-sided Bernstein bound |f-p| <= sqrt(2pL/S) + 2L/(3S).
    c=2*logtails/shots
    t=f+c/3
    upper=t+c/2+np.sqrt(c*t+c*c/4)
    radius=np.sqrt(c*upper)+c/3+2e-14
    return np.maximum(0.,f-radius),np.minimum(1.,f+radius)


def read(data,shots,logtails):
    registered_counts(data,shots)
    f=data[:,0]/shots
    lo,hi=count_box(f,shots,logtails)
    parts=[]
    for i,j in ((0,2),(3,1)):
        if lo[i]+lo[j]<=0:
            return None
        parts.append(((lo[i]-hi[j])/(lo[i]+hi[j]),(hi[i]-lo[j])/(hi[i]+lo[j])))
    return np.array(parts)


def run(record=False):
    cfg=plan()
    cases=list(forecasts(cfg))
    if record:
        rng=np.random.default_rng(cfg["seed"])
        fresh={k:sample(p,v*cfg["shots_per_volume"],rng) for k,v,p in cases}
        with COUNTS.open("xb") as output:
            np.savez_compressed(output,**fresh)
    logtails=np.log(2*12*len(cases)/cfg["alpha"])
    observed={}
    unresolved=0
    with np.load(COUNTS,allow_pickle=False) as raw:
        require(sorted(raw.files)==sorted(k for k,_,_ in cases),"incomplete detector menu")
        for k,v,p in cases:
            shots=v*cfg["shots_per_volume"]
            data=raw[k]
            out=read(data,shots,logtails)
            lo,hi=count_box(data/shots,shots,logtails)
            require(np.all((lo<=p)&(p<=hi)) and np.all(data[p==0]==0),
                    f"failed simultaneous finite-record forecast: {k}")
            if k.endswith("vertex_off"):
                require(out is None,"resolved absent source")
            else:
                if out is None:
                    unresolved+=1
                observed[k]=out

    D=walk.envelope_coefficient()
    sizes=cfg["strip_sizes"][-2:]
    rows={}
    for n in sizes:
        for window in cfg["windows"]:
            vals={}
            for menu in ("written","coherent"):
                first,second=[observed[f"strip_{n}_001_{window}_{b}_{menu}"] for b in (0,1)]
                if first is None or second is None:
                    vals[menu]=None
                else:
                    # Rectangular confidence interval for boosted minus original contrast.
                    diff=np.column_stack((second[:,0]-first[:,1],second[:,1]-first[:,0]))
                    vals[menu]=dict(first=first.round(6).tolist(),boosted=second.round(6).tolist(),
                                    difference=diff.round(6).tolist(),
                                    separated=bool(np.any((diff[:,0]>0)|(diff[:,1]<0))))
            ticks=walk.appointment(n,window)
            k=np.array([0.,0.,2*np.pi/n])
            envelope=[]
            for b in cfg["boosts"]:
                freq=ticks*(2*b+1)*np.dot(k,k)
                parabolic=np.exp(-1j*D*freq)
                always_mature=np.exp(-1j*freq/(8*np.tan(walk.MU)))
                bound=sum(walk.finite_envelope(j*k,ticks,D) for j in (b,b+1))
                small=walk.BrightWalk(np.array([b*k,(b+1)*k])).wave(ticks)
                exact=np.vdot(*small)
                require(abs(exact-parabolic)<=bound,"finite envelope bound failed")
                envelope.append(dict(periodic_error=round(float(abs(exact-parabolic)),8),
                                     finite_bound=round(bound,8),
                                     one_opening_limit_error=round(float(abs(exact-always_mature)),8)))
            vals["envelope"]=envelope
            rows[f"{n}_{window}"]=vals

    require(all(v["written"] is not None and not v["written"]["separated"] for v in rows.values()),
            "marked source records contradict boost invariance")
    require(any(v["coherent"] is not None and v["coherent"]["separated"] for k,v in rows.items()
                if "diffusive" in k),"coherent boost separation unresolved")
    result=dict(sha256=hashlib.sha256(COUNTS.read_bytes()).hexdigest(),panels=4*len(cases),
                draws=sum(4*v*cfg["shots_per_volume"] for _,v,_ in cases),unresolved=unresolved,
                periodic_curvature=round(D,9),mature_curvature=round(1/(8*np.tan(walk.MU)),9),
                windows=rows)
    if record:
        with SUMMARY.open("x",encoding="utf8",newline="\n") as output:
            output.write(json.dumps(result,indent=1,sort_keys=True,allow_nan=False)+"\n")
    require(json.loads(SUMMARY.read_text(encoding="utf8"))==result,"changed counts or inferred result")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true",help="write previously absent acquisition and summary")
    r=run(parser.parse_args().record)
    print(f"{r['panels']} complete panels; age-period curvature {r['periodic_curvature']:.9f}")
