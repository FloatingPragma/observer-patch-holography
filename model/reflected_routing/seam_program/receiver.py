"""Read the prepared packet coherently by reversing its successful native instrument."""
import argparse
import hashlib
from pathlib import Path

import numpy as np

from tools.strict_json import load
from ...polynomial_sv import carrier, reader, sectors
from ..__main__ import archive, receipt, source_pins
from ..kernel import counters, require
from .compiler import Compiler
from .preparation import HERE, execute_source, make_source, preparation_pins, raw_event, rotation, target


def pins():
    return {name:hashlib.sha256((HERE/name).read_bytes().replace(b"\r\n",b"\n")).hexdigest()
            for name in ("receiver.py","receiver_protocol.json","preparation_sources.json")}


def ideal_vector(wave,member):
    if member=="native":
        return wave
    require(member=="omit_phase","unknown receiver")
    weights = abs(wave)
    angle = np.zeros(len(wave))
    donor_phase = 0.
    for j in range(1,len(wave)):
        remainder = np.linalg.norm(np.r_[weights[0],weights[j+1:]])
        split = 2*np.arctan2(weights[j],remainder)
        donor_phase -= 8.5*split
        angle[j] = donor_phase+np.pi/2
    angle[0] = donor_phase
    return weights*np.exp(1j*angle)


def receive(comp,events,neg,wave,member,initial):
    state = initial.copy()
    stops = []
    event_count,angle_sum = 0,0.
    word = list(make_source(comp,neg,wave,events.kappa,member))
    home = comp.slots[0]
    test = events.sector.index(comp.marker,home)
    require([reader.gate(f,0)>>6 for f in (0,1)]==[1,29],"reader lost its agreed records")
    clock = np.zeros(comp.carrier.n,int)
    depth = clock.copy()
    def execute(seq):
        nonlocal event_count,angle_sum
        for edge,t in seq:
            raw_event(events,state,edge,t)
            a,b = comp.carrier.seams[edge]
            clock[[a,b]] += 1
            depth[[a,b]] = max(depth[a],depth[b])+1
            event_count += 1
            angle_sum += t
    for action,arg in reversed(word):
        if action=="events":
            execute((e,(-t)%(2*np.pi)) for e,t in reversed(arg))
        else:
            route = rotation(comp,arg,home,np.pi)
            execute(route)
            stops.append(float(abs(state[test])**2))
            state[test] = 0
            depth[:] = depth.max()+4
            execute(reversed(route))
    return state,stops,dict(events=event_count,angle_sum=angle_sum,reader_queries=len(stops),
                            port_counts=clock.tolist(),scheduled_depth=int(depth.max()))


def run(stage,root,record):
    require(pins()==load(HERE/"receiver_sources.json") and
            preparation_pins()==load(HERE/"preparation_sources.json") and
            source_pins()==load(HERE.parent/"sources.json"),"changed phase-sensitive instrument")
    spec = load(HERE/"preparation_protocol.json")
    pairs = load(HERE/"receiver_protocol.json")["members"]
    wave = target(spec).ravel()
    forecasts = {a+"_"+b:float(abs(np.vdot(ideal_vector(wave,b),ideal_vector(wave,a)))**2/16)
                 for a,b in pairs}
    if stage=="seal" and record:
        root.mkdir(parents=True,exist_ok=False)
    receipt(root/"forecast.json",forecasts,record and stage=="seal")
    if stage=="seal":
        return
    car = carrier.federation(0)
    sector = sectors.PairSector(car)
    comp = Compiler(car,0,spec["shape"])
    neg = list(range(1+len(wave),1+2*len(wave)))
    ev = sectors.SeamEvents(sector,8)
    emitted,_ = sectors.emit(sector.hamiltonian(8),sector.index(0,1))
    bound = sectors.emitter_error_bound(*car.degree[[0,1]])
    sender = {a:execute_source(comp,ev,neg,wave,a,emitted[1:]) for a in sorted({a for a,_ in pairs})}
    actual,raw = {},{}
    for a,b in pairs:
        initial,before,send_cost = sender[a]
        state,after,read_cost = receive(comp,ev,neg,wave,b,initial)
        seams = abs(state[sector.seam_states])**2
        none = abs(emitted[0])**2+np.vdot(state,state).real-seams.sum()
        law = np.r_[before,after,max(0.,none),seams]
        key = a+"_"+b
        raw[key] = counters(law,spec["shots"])
        yes = float(seams[car.seams.index((0,1))])
        actual[key] = dict(home= yes,ideal=forecasts[key],within_emitter_bound=abs(yes-forecasts[key])<=bound,
                           sender_events=send_cost["events"],receiver=read_cost)
    matched = actual["native_native"]["home"]
    for a,b in pairs:
        row = actual[a+"_"+b]
        row["separation_from_native_match"] = abs(row["home"]-matched)
        row["exceeds_emitter_bound"] = row["separation_from_native_match"]>bound
        print(a,b,row["home"],row["ideal"],row["within_emitter_bound"])
    archive(root/"opened.npz",raw,record)
    receipt(root/"result.json",actual,record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",choices=("seal","open"),default="open")
    parser.add_argument("--record",type=Path)
    args = parser.parse_args()
    run(args.stage,args.record or HERE/"receiver_records",args.record is not None)


if __name__=="__main__":
    main()
