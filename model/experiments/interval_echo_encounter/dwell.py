"""Return spectroscopy of a binary pair and the native write that produced it."""
import argparse
from hashlib import sha256
import json
from math import ceil, lgamma, log, pi
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.special import jv

from model.experiments.seam_program.gate import archive, receipt
from .dwell_source import Station, local_blocks, require

HERE = Path(__file__).resolve().parent

FILES = ("dwell_protocol.json","dwell.py","dwell_source.py",
         "../../candidate_laws/contrast_memory_spherical_v1/law.py",
         "../../candidate_laws/contrast_memory_spherical_v1/sectors.py",
         "../../executable_model/carrier.py","../seam_program/gate.py",
         *[f"../../executable_model/carriers/FEDERATION_LEVEL{n}.npz" for n in range(3)])


def plan():
    pins = json.loads((HERE/"dwell_sources.json").read_text(encoding="utf8"))
    require(set(pins)==set(FILES),"incomplete dwell source map")
    for name in FILES:
        raw = (HERE/name).read_bytes()
        digest = sha256(raw if name.endswith(".npz") else raw.replace(b"\r\n",b"\n")).hexdigest()
        require(digest==pins[name],f"changed dwell source: {name}")
    return json.loads((HERE/"dwell_protocol.json").read_text(encoding="utf8"))


def returns(H, port, times):
    """Complete spectral measure of the sharp port word through Chebyshev moments."""
    times=np.asarray(times)
    require(times.ndim==1 and times.size>0 and np.isfinite(times).all() and min(times)>=0,
            "invalid return exposures")
    if max(times)==0:
        return np.ones(len(times),complex)
    diag = H.diagonal()
    off = np.asarray(abs(H).sum(axis=1)).ravel()-abs(diag)
    lo,hi = min(diag-off),max(diag+off)
    center,radius = (lo+hi)/2,(hi-lo)/2
    require(radius>0 and min(times)>=0,"invalid return spectrum")
    A = ((H-sparse.eye(H.shape[0])*center)/radius).tocsr()
    x = radius*max(times)
    order = max(2,ceil(x))
    # |J_n(x)| <= (x/2)^n/n!; at this order the tail ratios are <= 1/2.
    while order*log(x/2)-lgamma(order+1)>log(1e-15):
        order+=1
    if order%2:
        order+=1
    v = np.zeros(H.shape[0])
    v[port]=1.
    u = A@v
    moments = np.empty(order+1)
    moments[0],moments[1]=1.,u[port]
    first = u[port]
    for j in range(1,order//2+1):
        moments[2*j]=2*np.dot(u,u)-1
        if 2*j<order:
            nxt = 2*(A@u)-v
            moments[2*j+1]=2*np.dot(u,nxt)-first
            v,u=u,nxt
    k=np.arange(order+1)
    coeff=2*jv(k,radius*times[:,None])*((-1j)**(k%4))
    coeff[:,0]/=2
    return np.exp(-1j*times*center)*(coeff@moments)


def spectrum(H,port,kappa,gap):
    seed=np.zeros(H.shape[0])
    seed[port]=1
    values,vectors=eigsh(H,k=1,which="SA" if kappa>0 else "LA",v0=seed,tol=1e-13)
    value,v=float(values[0]),vectors[:,0]
    require(np.linalg.norm(H@v-value*v)<1e-8 and abs(value-H[port,port])<gap/2,
            "isolated spectral line not enclosed")
    return value,float(v[port]**2)


def valid_shots(shots):
    return (isinstance(shots,(int,np.integer)) and not isinstance(shots,(bool,np.bool_))
            and 0<shots<=np.iinfo(np.int32).max)


def counts_for(a,b,clock,shots,rng):
    a,b,clock=map(np.asarray,(a,b,clock))
    require(valid_shots(shots) and a.ndim==1 and a.size>0 and b.shape==a.shape
            and np.isfinite(a).all() and np.isfinite(b).all()
            and np.all(abs(a)<=1+1e-10) and np.all(abs(b)<=1+1e-10),
            "invalid return amplitudes or trial count")
    require(clock.shape==(a.size,4) and clock.dtype.kind in "fiu" and np.isfinite(clock).all()
            and np.all(clock>=0) and np.all(abs(clock.sum(axis=1)-1)<1e-12),
            "invalid reference clock probabilities")
    hit=(abs(a)**2+abs(b)**2)/2
    z=a*b.conj()
    cohorts=[]
    for quadrature in (z.real,z.imag):
        signal=np.stack(((hit+quadrature)/2,(hit-quadrature)/2,1-hit),axis=-1)
        require(signal.min()>-2e-10 and np.max(abs(signal.sum(axis=-1)-1))<1e-12,
                "incomplete signal measurement")
        joint=np.maximum(signal,0)[:,:,None]*clock[:,None,:]
        rows=[]
        for panel in joint:
            cdf=np.rint(2**32*np.cumsum(panel.ravel()/panel.sum())).astype(np.int64)
            rows.append(rng.multinomial(shots,np.diff(cdf,prepend=0)/2**32).reshape(3,4))
        cohorts.append(rows)
    return np.asarray(cohorts,dtype=np.int32)


def read(counts,shots,component_error,phase_distortion):
    require(valid_shots(shots) and np.isfinite(component_error) and component_error>0
            and (phase_distortion is None or (np.isfinite(phase_distortion) and phase_distortion>=0))
            and counts.ndim==4 and counts.shape[0]==2 and counts.shape[2:]==(3,4)
            and counts.shape[1]>=2 and counts.dtype.kind in "iu" and np.all(counts>=0)
            and np.all(counts<=shots) and np.all(counts.sum(axis=(2,3),dtype=np.int64)==shots),
            "invalid joint pair and timebase counts")
    sig=counts.sum(axis=-1,dtype=np.int64)/shots
    parts=sig[:,:,0]-sig[:,:,1]
    z=parts[0]+1j*parts[1]
    c=counts.sum(axis=(0,2),dtype=np.int64)/(2*shots)
    hand=(c[:,0]+c[:,1]-c[:,2]-c[:,3])+1j*(c[:,0]-c[:,1]+c[:,2]-c[:,3])
    result=dict(minimum_tag_probability=float(np.min(1-sig[:,:,2])),quiet_at_last=float(sig[:,-1,2].mean()),
                minimum_visibility=float(min(abs(z))),hand_elapsed=float(np.diff(np.unwrap(np.angle(hand))[[0,-1]])[0]),
                rate_interval=None)
    radius=np.sqrt(2)*component_error
    if phase_distortion is None or min(abs(z))<=radius or min(abs(hand))<=radius:
        return result
    error=np.arcsin(radius/min(abs(z)))+phase_distortion
    clock_error=np.arcsin(radius/min(abs(hand)))
    # The isolated eigenvalue increases by at most one per unit of parity coupling.
    # One native tick separates cohorts. These margins establish the winding.
    if 1+2*max(error,clock_error)>=pi or result["hand_elapsed"]<=2*clock_error:
        return result
    phase=np.unwrap(np.angle(z))
    slope=(phase[-1]-phase[0]-2*error,phase[-1]-phase[0]+2*error)
    bounds=[s/t for s in slope for t in (result["hand_elapsed"]-2*clock_error,result["hand_elapsed"]+2*clock_error)]
    result["rate_interval"]=[min(bounds),max(bounds)]
    return result


def run(spec):
    boards={}
    for level in spec["geometry_levels"]:
        b=Station(level,spec["tag"])
        isolated=[b.isolation(8,j) for j in (0,1)]
        boards[str(level)]=dict(clock_cell=b.clock_cell,clock_dimers=b.dimers,distance=b.clock_distance,
            boundary_writes=len(b.cuts),removed_mobile_ports=12,cells=b.host.cells,
            sharp_exchange=b.d,complement_saving=b.saving,positive_parities=isolated)
    b=Station(spec["flight_level"],spec["tag"])
    b.flight()
    times=np.arange(spec["last_tick"]+1,dtype=float)
    x,y=np.cos(times),np.sin(times)
    clock=np.stack([(1+i*x)*(1+j*y)/4 for i in (1,-1) for j in (1,-1)],axis=-1)
    rng=np.random.default_rng(spec["seed"])
    panels=len(times)*len(spec["strengths"])*len(spec["arms"])
    # Four signed marginals per time and setting, including the pooled clock.
    error=np.sqrt(2*np.log(8*panels/spec["alpha"])/spec["shots"])+1e-8
    result=dict(geometry=boards,pair_words=b.sector.P,component_error=error,arms={})
    raw={}
    for kappa in spec["strengths"]:
        local_blocks(kappa)
        for member in spec["arms"]:
            tag=f"{member}_{kappa}"
            H=b.generator(kappa,0,member)
            a=returns(H,b.p,times)
            cert,lines,phase_bound=None,[],None
            if member=="always_ticking":
                odd=a*np.exp(-1j*times)
                phase_bound=0.
            else:
                H1=b.generator(kappa,1,member)
                odd=returns(H1,b.p,times)
                cert=[b.isolation(kappa,j) for j in (0,1)]
                if all(cert):
                    lines=[spectrum(M,b.p,kappa,record["gap"]) for M,record in zip((H,H1),cert,strict=True)]
                    phase_bound=sum(np.arcsin((1-u["weight"])/u["weight"]) for u in cert)
                    for (energy,weight),u,amp in zip(lines,cert,(a,odd),strict=True):
                        require(weight>=u["weight"]-1e-10 and np.min(abs(amp)**2)>=u["survival"]-1e-9,
                                "time record violates isolated-line bound")
            rows=counts_for(a,odd,clock,spec["shots"],rng)
            raw[tag]=rows
            observed=read(rows,spec["shots"],error,phase_bound)
            delta=lines[1][0]-lines[0][0] if lines else (1. if member=="always_ticking" else None)
            result["arms"][tag]=dict(observed=observed,isolation=cert,isolated_energies=[u[0] for u in lines],
                spectral_weights=[u[1] for u in lines],branch_difference=delta,
                positive_gate=bool(kappa==8 and member=="native" and all(c["survival"]>.9 for c in cert)
                                   and observed["rate_interval"] is not None and observed["rate_interval"][1]<1.))
    return raw,result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true",help="create the two previously absent receipts")
    args=parser.parse_args()
    records,result=run(plan())
    archive(HERE/"dwell_counts.npz",records,args.record)
    receipt(HERE/"dwell_result.json",result,args.record)
    print(f"{sum(int(x.sum()) for x in records.values())} joint trials; native clocked pair:",
          result["arms"]["native_8"]["observed"]["rate_interval"])


if __name__=="__main__":
    main()
