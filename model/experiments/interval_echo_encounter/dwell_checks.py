"""Single-term interventions and independent physical reconstruction of the pair receiver."""
import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm

from model.experiments.seam_program.gate import archive, receipt
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import Carrier
from . import dwell
from .dwell_source import Station, require
from .memory import BITS, lift, supports
from .check_transport import checked_steps


class PhysicalFlight:
    """Start with the full blank graph including the tag, then replace that single contact."""
    def __init__(self,station):
        host=station.host
        ports=[p for p in range(host.n) if p//12!=station.clock_cell]
        local={p:i for i,p in enumerate(ports)}
        edges=[(local[a],local[b]) for a,b in host.seams if a in local and b in local]
        graph=Carrier("exterior including the tag",len(ports),edges)
        self.sector=PairSector(graph)
        ends=tuple(local[p] for p in station.tag)
        self.p=self.sector.index(*ends)
        self.positions=supports(self.sector,ends)
        deficits=host.degree[ports]-graph.degree
        self.boundary=deficits[self.sector.aa]+deficits[self.sector.bb]

    def generator(self,kappa,parity,member):
        require(member in ("native","always_ticking","binding_only","pause_only"),"unknown intervention")
        U=law.echo()
        g0,g1=law.generator_parts()
        h=U@(g0+kappa*g1)@U.T
        sign=np.zeros((11,2))
        for i,s in enumerate((2,-2)):
            for c in (0,1):
                sign[1+law.LABELS.index((c,s)),i]=1/np.sqrt(2)
        balanced=np.kron(sign,np.eye(9)[:,law.SRC.index((1,1)),None])
        restore=np.zeros((2,2))
        if member in ("always_ticking","pause_only"):
            restore+=4*kappa*np.eye(2)
        if member in ("always_ticking","binding_only"):
            restore+=np.array([[.5,-.5],[-.5,.5]])
        h+=balanced@restore@balanced.T
        ket=(sign[:,0]+(-1)**parity*sign[:,1])/np.sqrt(2)
        binary=np.eye(9)[:,[law.SRC.index(pair) for pair in BITS]]
        embed=np.kron(ket[:,None],binary)
        block=embed.T@h@embed
        rounded=np.rint(2*block)/2
        require(np.linalg.norm(h@embed-embed@rounded)<2e-12,"a physical control leaks from the parity sector")
        blank=binary.T@law.seam_term(kappa)@binary
        sector=self.sector
        return (sector.hamiltonian(kappa)-sparse.diags(kappa*self.boundary,dtype=float)
                +lift(rounded-blank,self.positions,sector.P)).tocsr()


def controls(spec,base):
    station=Station(base["flight_level"],base["tag"])
    station.flight()
    direct=PhysicalFlight(station)
    times=np.arange(base["last_tick"]+1)
    rng=np.random.default_rng(spec["seed"])
    # Actual dimer pulses, including the preparation of the advanced quadrature.
    single=[law.SRC.index((1,0)),law.SRC.index((0,1))]
    matrix=law.seam_term(spec["coupling"])[np.ix_(single,single)]+spec["coupling"]*np.eye(2)
    left=np.array([1.,0.])
    start=np.column_stack((left,expm(-1j*1.5*np.pi*matrix)@left))
    clock=[]
    for t in times:
        d=abs(expm(-1j*t*matrix)@start)**2
        clock.append(np.outer(d[:,0],d[:,1]).ravel())
    clock=np.array(clock)
    # Allocate the confidence budget across these two added settings.
    error=np.sqrt(2*np.log(8*len(times)*len(spec["controls"])/base["alpha"])/base["shots"])+1e-8
    arrays,out={},{}
    for member in spec["controls"]:
        pair=[dwell.returns(direct.generator(spec["coupling"],j,member),direct.p,times) for j in (0,1)]
        rows=dwell.counts_for(*pair,clock,base["shots"],rng)
        arrays[member]=rows
        info=dwell.read(rows,base["shots"],error,0. if member=="binding_only" else None)
        p=rows.sum(axis=-1,dtype=np.int64)/base["shots"]
        quadratures=p[:,-1,0]-p[:,-1,1]
        info["last_principal_phase"]=float(np.angle(quadratures[0]+1j*quadratures[1]))
        out[member]=info
    return arrays,out


def extended(base):
    """Check all amplitudes against an enclosed state evolution, retaining its common phase."""
    station=Station(base["flight_level"],base["tag"])
    station.flight()
    direct=PhysicalFlight(station)
    times=np.arange(base["last_tick"]+1)
    start=np.zeros(direct.sector.P,complex)
    start[direct.p]=1
    worst=0.
    worst_signal=worst_clock=worst_bin=0.
    # Identical parity-0 surfaces and exact scalar shifts are checked only once.
    cache={}
    with np.load(dwell.HERE/"dwell_counts.npz") as original, np.load(dwell.HERE/"dwell_controls_counts.npz") as extra:
        for kappa in base["strengths"]:
            arms=base["arms"]+(["binding_only","pause_only"] if kappa==8 else [])
            for member in arms:
                waves=[]
                for j in (0,1):
                    H=direct.generator(kappa,j,member)
                    if member in base["arms"]:
                        diff=H-station.generator(kappa,j,member)
                        require(not diff.count_nonzero(),"projected and directly assembled generators differ")
                    shift=H[direct.p,direct.p]
                    reduced=H-sparse.eye(H.shape[0],format="csr")*shift
                    reduced.eliminate_zeros()
                    reduced.sort_indices()
                    key=sha256(reduced.data.tobytes()+reduced.indices.tobytes()+reduced.indptr.tobytes()).digest()
                    if key not in cache:
                        half=reduced/2
                        off=np.asarray(abs(half).sum(axis=1)).ravel()-abs(half.diagonal())
                        center=(min(half.diagonal()-off)+max(half.diagonal()+off))/2
                        amps=np.ones(len(times),complex)
                        budget=0.
                        for t,v,err in checked_steps(half,start,set(times[1:]),degree=36):
                            amps[t]=np.dot(v,v)*np.exp(-2j*t*center)
                            budget=max(budget,3*err+2e-12)
                        cache[key]=amps,budget
                    amps,budget=cache[key]
                    direct_amp=amps*np.exp(-1j*times*shift)
                    worst=max(worst,float(np.max(abs(dwell.returns(H,direct.p,times)-direct_amp))+budget))
                    waves.append(direct_amp)
                records=original[f"{member}_{kappa}"] if member in base["arms"] else extra[member]
                same=(abs(waves[0])**2+abs(waves[1])**2)/2
                z=waves[0]*waves[1].conj()
                hand=np.column_stack([(1+c*np.cos(times))*(1+s*np.sin(times))/4 for c in (1,-1) for s in (1,-1)])
                clock_counts=records.sum(axis=(0,2),dtype=np.int64)/(2*base["shots"])
                signs=np.array([[1,1],[1,-1],[-1,1],[-1,-1]])
                clock_residual=np.max(abs((clock_counts-hand)@signs))
                panels=len(times)*(len(base["arms"])*len(base["strengths"]) if member in base["arms"] else 2)
                component_budget=np.sqrt(2*np.log(8*panels/base["alpha"])/base["shots"])+1e-8
                require(clock_residual<component_budget,"reference clock exceeds its registered error budget")
                worst_clock=max(worst_clock,clock_residual)
                for i,signal in enumerate((z.real,z.imag)):
                    observed=records[i].sum(axis=-1,dtype=np.int64)/base["shots"]
                    signal_residual=np.max(abs(observed[:,0]-observed[:,1]-signal))
                    require(signal_residual<component_budget,"signal exceeds its registered error budget")
                    worst_signal=max(worst_signal,signal_residual)
                    chance=np.column_stack(((same+signal)/2,(same-signal)/2,1-same))
                    error=np.max(abs(records[i]/base["shots"]-chance[:,:,None]*hand[:,None,:]))
                    worst_bin=max(worst_bin,error)
                    bins=8*2*len(times)*12
                    require(error<np.sqrt(np.log(2*bins/.01)/(2*base["shots"]))+1e-8,
                            "joint counts outside the independent simultaneous bin check")
    require(worst<1e-8,"full-domain numerical enclosure exceeds the acquisition budget")
    print(f"All return amplitudes agree within {worst:.3g}, including polynomial and rounding budgets.")
    print(f"Maximum absolute frequency-record residuals: signal {worst_signal:.8g}, "
          f"clock {worst_clock:.8g}, joint bin {worst_bin:.8g}.")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true")
    parser.add_argument("--extended",action="store_true")
    args=parser.parse_args()
    base=dwell.plan()
    files=("dwell_controls_protocol.json","dwell_checks.py","dwell_sources.json","check_transport.py","memory.py")
    pins=json.loads((dwell.HERE/"dwell_controls_sources.json").read_text(encoding="utf8"))
    require(set(pins)==set(files) and all(pins[n]==sha256((dwell.HERE/n).read_bytes().replace(b"\r\n",b"\n")).hexdigest()
                                          for n in files),"changed intervention source")
    spec=json.loads((dwell.HERE/"dwell_controls_protocol.json").read_text(encoding="utf8"))
    data,out=controls(spec,base)
    archive(dwell.HERE/"dwell_controls_counts.npz",data,args.record)
    receipt(dwell.HERE/"dwell_controls_result.json",out,args.record)
    if args.extended:
        extended(base)
    print("Energy and ticking interventions reproduced.")


if __name__=="__main__":
    main()
