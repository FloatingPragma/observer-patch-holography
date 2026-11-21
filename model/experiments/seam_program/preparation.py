"""Prepare a complex traveling-wave packet with native exchange, pair fuel and seam reads."""
import argparse
import hashlib

import numpy as np

from tools.strict_json import load
from model.executable_model import carrier
from model.candidate_laws.contrast_memory_spherical_v1 import reader, sectors
from .gate import archive, counters, multipliers, real_step, receipt, require
from .__main__ import HERE, ideal_law, pins
from .compiler import Compiler


def preparation_pins():
    hashes = pins()
    for name in ("preparation.py","preparation_protocol.json","sources.json"):
        hashes[name] = hashlib.sha256((HERE/name).read_bytes().replace(b"\r\n",b"\n")).hexdigest()
    return hashes


def target(spec):
    shape = tuple(spec["shape"])
    coordinates = np.moveaxis(np.indices(shape),0,-1)
    result = np.zeros((*shape,8),complex)
    equal = np.ones(8)/np.sqrt(8)
    for (x,y),(r,i) in zip(spec["modes"],spec["weights"],strict=True):
        k = 2*np.pi*np.array([x,y,0,0])/shape
        a = np.cos(k).mean()
        partner = (multipliers(k)[np.arange(8)^1]*equal-a*equal)/np.sqrt(1-a*a)
        result += (r+1j*i)*np.exp(1j*(coordinates@k))[...,None]*(equal+1j*partner)/np.sqrt(2)
    norm = np.linalg.norm(result)
    require(norm>0 and abs(result.sum())<1e-11,"target must have zero uniform amplitude")
    return result/norm


def rotation(comp,a,b,angle):
    route = comp.exchange(a,b,4)
    if not route or angle==0:
        return []
    middle = len(route)//2
    return [(edge,angle if j==middle else np.pi) for j,(edge,_) in enumerate(route)]


def make_source(comp,neg,desired,kappa,member):
    """Yield alternating local control words and one-seam detector queries."""
    require(kappa==8,"packet source requires the registered contrast coupling")
    slots = comp.slots
    require(desired.shape==(len(slots),) and abs(np.linalg.norm(desired)-1)<1e-12,"invalid launch")
    require(len(neg)==len(slots) and len(set(slots+neg))==2*len(slots)
            and comp.marker not in neg and min(neg)>=0 and max(neg)<comp.carrier.n,"rails overlap")
    require(member in ("native","omit_phase","unread_negative"),"unknown preparation")
    block = [slots[0],neg[0],slots[1],neg[1]]
    for mask in (1,2,3):
        for j in range(4):
            if j<(j^mask):
                yield "events",rotation(comp,block[j],block[j^mask],np.pi/2)
    for site in block[2:]:
        yield "exclude",site
    phase = 0.
    phases = np.zeros(len(slots))
    for j in range(1,len(slots)):
        remainder = np.linalg.norm(np.r_[desired[0],desired[j+1:]])
        angle = 2*np.arctan2(abs(desired[j]),remainder)
        for bank in (slots,neg):
            yield "events",rotation(comp,bank[0],bank[j],angle)
        phase -= (kappa+.5)*angle
        phases[j] = phase+np.pi/2
    phases[0] = phase
    if member!="omit_phase":
        for j,(current,requested) in enumerate(zip(phases,np.angle(desired),strict=True)):
            angle = float((current-requested)%(2*np.pi))/(1+kappa)
            yield "events",rotation(comp,slots[j],neg[j],angle)
    if member!="unread_negative":
        for site in neg:
            yield "exclude",site


def raw_event(events,state,edge,angle):
    # SPEC section 8 omits kappa*(n_a+n_b) from the original binary-generator expression.
    # Restore its commuting phase on every affected two-load basis state.
    events.apply(state,edge,angle)
    left,right,both = events.tables[edge]
    z = np.exp(-1j*events.kappa*angle)
    state[left] *= z
    state[right] *= z
    state[both] *= z*z


def execute_source(comp,events,neg,wave,member,initial):
    state = initial.copy()
    losses = np.zeros(2+len(neg))
    clock = np.zeros(comp.carrier.n,int)
    depth = clock.copy()
    nread,nevents,angle_sum = 0,0,0.
    sector = events.sector
    ref = comp.slots[0]
    seam_state = sector.index(comp.marker,ref)
    yes,no = (reader.gate(flag,0)>>6 for flag in (1,0))
    require((yes,no)==(29,1),"source needs agreed yes and no instruments")
    def apply(sequence):
        nonlocal nevents,angle_sum
        for edge,t in sequence:
            raw_event(events,state,edge,t)
            a,b = comp.carrier.seams[edge]
            clock[[a,b]] += 1
            depth[[a,b]] = max(depth[a],depth[b])+1
            nevents += 1
            angle_sum += t
    for action,data in make_source(comp,neg,wave,events.kappa,member):
        if action=="events":
            apply(data)
        else:
            route = rotation(comp,data,ref,np.pi)
            apply(route)
            losses[nread] = abs(state[seam_state])**2
            state[seam_state] = 0  # unnormalized no branch; the yes branch terminates with its own record
            nread += 1
            # Conservative synchronization barrier while the accepted record enables continuation.
            depth[:] = depth.max()+4
            apply(reversed(route))
    return state,losses,dict(events=nevents,angle_sum=angle_sum,reader_queries=nread,
                             port_counts=clock.tolist(),scheduled_depth=int(depth.max()))


def run(stage,root,record):
    spec = load(HERE/"preparation_protocol.json")
    require(preparation_pins()==load(HERE/"preparation_sources.json"),"changed preparation sources")
    car = carrier.federation(0)
    sector = sectors.PairSector(car)
    comp = Compiler(car,0,spec["shape"])
    neg = list(range(1+len(comp.slots),1+2*len(comp.slots)))
    wave = target(spec)
    expected_losses = np.r_[.25,.25,.25*abs(wave.ravel())**2]
    laws,quantized = {},{}
    grid = wave.copy()
    for t in range(max(spec["waits"])+1):
        if t in spec["waits"]:
            for q in spec["apertures"]:
                key = f"{t}_{q}"
                laws[key] = np.r_[expected_losses,.25*ideal_law(grid,comp,sector,q)]
                quantized[key] = counters(laws[key],spec["shots"])
        grid = real_step(grid)
    if stage=="seal" and record:
        root.mkdir(exist_ok=False,parents=True)
    quantized = archive(root/"ideal.npz",quantized,stage=="seal" and record)
    if stage=="seal":
        return
    ev = sectors.SeamEvents(sector,8)
    entry = sector.index(0,1)
    emitted,_ = sectors.emit(sector.hamiltonian(8),entry)
    residual_fuel = abs(emitted[0])**2
    bound = sectors.emitter_error_bound(*car.degree[[0,1]])
    marked = np.array([sector.index(0,j) for j in comp.slots])
    program = comp.build("nominal")
    counts,summary = {},{}
    for member in spec["members"]:
        state,losses,resources = execute_source(comp,ev,neg,wave.ravel(),member,emitted[1:])
        accepted = float(np.vdot(state,state).real+residual_fuel)
        require(abs(accepted+losses.sum()-1)<1e-10,"missing preparation outcome")
        fidelity = float(abs(np.vdot(wave.ravel(),state[marked]))**2/accepted)
        worst = 0.
        for t in range(max(spec["waits"])+1):
            if t in spec["waits"]:
                for q in spec["apertures"]:
                    localized = state.copy()
                    for edge,angle in rotation(comp,comp.slots[q],comp.slots[0],np.pi):
                        raw_event(ev,localized,edge,angle)
                    seams = abs(localized[sector.seam_states])**2
                    law = np.r_[losses,max(0.,accepted-seams.sum()),seams]
                    key = f"{t}_{q}"
                    counts[f"{member}_{key}"] = counters(law,spec["shots"])
                    worst = max(worst,float(abs(law-laws[key]).sum()/2))
            if t<max(spec["waits"]):
                for edge,q_angle in program:
                    raw_event(ev,state,edge,np.pi*q_angle/4)
        if member=="native":
            require(worst<=bound and accepted>.20,"native launch outside its predeclared budget")
        summary[member] = dict(acceptance=accepted,conditional_fidelity=fidelity,
                               max_joint_tv=worst,exceeds_emitter_bound=bool(worst>bound),resources=resources)
        print(member,accepted,fidelity,worst,resources["events"],resources["reader_queries"])
    stored = archive(root/"opened.npz",counts,record)
    for member,row in summary.items():
        # The terminal alphabet comprises no pair and every seam. All preceding bins are stops.
        accepted_counts = set()
        distances = []
        for key,forecast in quantized.items():
            opened = stored[f"{member}_{key}"]
            accepted_counts.add(int(opened[-(1+len(car.seams)):].sum(dtype=np.int64)))
            distances.append(abs(opened.astype(np.int64)-forecast).sum()/(2*spec["shots"]))
        require(len(accepted_counts)==1,"source acceptance changed with the terminal aperture")
        hits = accepted_counts.pop()
        measured_tv = float(max(distances))
        row.update(acceptance_count=hits,exposures=spec["shots"],measured_acceptance=hits/spec["shots"],
                   max_counter_tv=measured_tv,counter_exceeds_emitter_bound=measured_tv>bound)
    receipt(root/"result.json",summary,record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",choices=("seal","open"),default="open")
    parser.add_argument("--record",type=type(HERE))
    args = parser.parse_args()
    run(args.stage,args.record or HERE/"preparation_records",args.record is not None)


if __name__=="__main__":
    main()
