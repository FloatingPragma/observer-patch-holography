"""Observer comparison of a fixed contact phase and a resonant age cycle."""
import argparse
import hashlib
import json

import numpy as np

from . import fringe
from . import fringe_walk as w

SPEC=w.PERIOD
RES=np.pi/(4*SPEC)
PROTOCOL=fringe.HERE/"fringe_resonance_protocol.json"
ARCHIVE=fringe.HERE/"fringe_resonance_counts.npz"
SUMMARY=fringe.HERE/"fringe_resonance_result.json"


def speed():
    cycle=w.cycle_jet(RES)
    rest=np.exp(-1j*SPEC*RES)
    w.require(np.max(abs(cycle[0]-rest*np.eye(2)))<2e-10,"incorrect resonance")
    value=abs(cycle[1,1,0])/SPEC
    # Reverse triangle inequality for the geometric sum of the mature hops.
    lower=(np.sin((SPEC-w.OPEN)*RES)/np.sin(RES)-w.OPEN)/(2*SPEC)
    upper=(SPEC-w.OPEN+sum(abs(np.sin(w.law.hop_angle(j,w.MEMBER.maturation)))
                           for j in range(w.OPEN)))/(2*SPEC)
    w.require(0<lower<=value<=upper<.5,"resonant velocity outside analytic bounds")
    return float(value),float(lower),float(upper)


def kg_fringe(boost, time, mass, velocity):
    orders=2*np.pi*np.array([boost,boost+1])
    energy=np.hypot(mass,velocity*orders)
    sine=time*np.sinc(energy*time/np.pi)
    s=np.cos(energy*time)-1j*mass*sine
    tangent=velocity*orders*sine
    return np.conj(s[0])*s[1]+tangent[0]*tangent[1]


def family(menu,n):
    if menu=="native":
        return w.MU
    if menu=="zero_limit":
        return 1/n
    if menu=="resonant":
        return RES
    if menu=="detuned_resonant":
        return RES+1/n
    raise ValueError("unknown phase family")


def generate(cfg,v):
    for m in cfg["periods_per_side"]:
        n=m*SPEC
        directions=np.zeros((3,3))
        directions[:,2]=2*np.pi*np.arange(3)/n
        for menu in cfg["laws"]:
            walk=w.BrightWalk(directions,mu=family(menu,n))
            for t in cfg["proper_clock_fractions"]:
                ticks=w.triggered_service(int(t*n/2))
                w.require(ticks==t*n and ticks%SPEC==0,"not an optical cycle appointment")
                wave=walk.wave(ticks)
                for b in cfg["boosts"]:
                    cross=np.vdot(wave[b],wave[b+1])
                    local=(1+(fringe.PHASES*cross).real)/(2*n)
                    w.require(min(local)>-1e-13 and max(local)<=1/n+1e-13,"invalid local law")
                    local=np.clip(local,0,1/n)
                    p=np.column_stack((local,.5-local,np.full(4,.5)))
                    limit=(kg_fringe(b,t,float(menu=="detuned_resonant"),v)
                           if menu in ("resonant","detuned_resonant") else 1.+0j)
                    yield f"{m}_{menu}_{t}_{b}",n,p,limit,cross


def run(record=False):
    fringe.plan()
    cfg=json.loads(PROTOCOL.read_text(encoding="utf8"))
    w.optical_reference()
    v,lower,upper=speed()
    cases=list(generate(cfg,v))
    if record:
        rng=np.random.default_rng(cfg["seed"])
        raw={}
        for k,n,p,_,_ in cases:
            shots=n*cfg["shots_per_site"]
            # One fair enrollment and one conditional position read per original start.
            raw[k]=fringe.fair_enrollment(shots,2*p[:,0],rng)
        with ARCHIVE.open("xb") as stream:
            np.savez_compressed(stream,**raw)
    L=np.log(2*12*len(cases)/cfg["alpha"])
    peak={}
    final={}
    with np.load(ARCHIVE,allow_pickle=False) as raw:
        w.require(sorted(raw.files)==sorted(row[0] for row in cases),"incomplete resonance menus")
        for key,n,p,limit,cross in cases:
            shots=n*cfg["shots_per_site"]
            counts=raw[key]
            bounds=fringe.check_panel(counts,shots,L,p,key)
            prefix=str(n//SPEC)
            peak[prefix]=max(peak.get(prefix,0.),float(abs(cross-limit)))
            if n==max(cfg["periods_per_side"])*SPEC:
                target=np.array([limit.real,limit.imag])
                w.require(bounds is not None and np.all(bounds[:,0]-.02<=target)
                          and np.all(target<=bounds[:,1]+.02),f"finest limit gate failed: {key}")
                final[key.removeprefix(prefix+"_")]=bounds
    separated=False
    for t in cfg["proper_clock_fractions"]:
        for b in cfg["boosts"]:
            a,c=[final[f"{menu}_{t}_{b}"] for menu in ("native","resonant")]
            separated |= bool(np.any((a[:,0]>c[:,1])|(c[:,0]>a[:,1])))
    w.require(separated,"native and resonant records unresolved")
    result=dict(speed=round(v,10),speed_bounds=[round(lower,10),round(upper,10)],
                maximum_limit_residual={k:round(x,6) for k,x in peak.items()},
                finest_contrast={k:a.round(6).tolist() for k,a in final.items()},
                panels=4*len(cases),draws=sum(4*n*cfg["shots_per_site"] for _,n,_,_,_ in cases),
                sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest())
    if record:
        with SUMMARY.open("x",encoding="utf8",newline="\n") as output:
            output.write(json.dumps(result,indent=1,sort_keys=True,allow_nan=False)+"\n")
    w.require(json.loads(SUMMARY.read_text(encoding="utf8"))==result,"changed resonant records")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true")
    r=run(parser.parse_args().record)
    print(f"{r['panels']} panels; resonant speed {r['speed']}, bounds {r['speed_bounds']}")
