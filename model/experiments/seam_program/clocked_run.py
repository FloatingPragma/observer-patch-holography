"""Replay the clock-stopped scalar instrument on the complete active binary-pair sector."""
import argparse
from pathlib import Path

import numpy as np
from scipy import sparse

from tools.strict_json import load
from model.candidate_laws.contrast_memory_spherical_v1 import sectors
from .clocked_apparatus import clock_laws, patches
from .clocked_field import refinement, scalar_kernel, stopping_law
from .compiler import Compiler
from .gate import archive, counters, real_step, receipt, require
from .preparation import HERE, execute_source, make_source, preparation_pins, raw_event, rotation
from .receiver import pins as receiver_pins


def finite_source(data, host, ports, fuel):
    sector = sectors.PairSector(data)
    lost = host.degree[ports]-data.degree
    # The open contacts use their original missing-seam reference terms.
    h = sector.hamiltonian(8)-sparse.diags(8*(lost[sector.aa]+lost[sector.bb]))
    pair = sector.index(0,1)
    emitted, resonance = sectors.emit(h,pair,fuel)
    ideal = np.zeros_like(emitted)
    ideal[1+pair] = -1j
    error = float(np.linalg.norm(emitted-ideal))
    bound = sectors.emitter_error_bound(*data.degree[:2],fuel)
    require(error <= bound and abs(np.linalg.norm(emitted)-1)<1e-11, "fuel outside its unitary budget")
    return sector, emitted, dict(resonance=resonance, bound=bound, norm_error=error)


def translated_word(full, local, position):
    grid = np.asarray(full.slots).reshape((*full.shape,8))
    targets = grid[tuple(position)].ravel()
    return [instruction for a,b in zip(local.slots,targets,strict=True)
            for instruction in full.exchange(a,int(b),4)]


def receive_columns(local, events, neg, wave, states):
    """Apply the same adjoint instrument to all waits, querying each column independently."""
    states = states.copy()
    seen = []
    count = 0
    duration = 0.
    ref = local.slots[0]
    pair = events.sector.index(local.marker,ref)
    def execute(word):
        nonlocal count,duration
        for edge,angle in word:
            raw_event(events,states,edge,angle)
            count += 1
            duration += angle
    for action,arg in reversed(list(make_source(local,neg,wave,events.kappa,"native"))):
        if action == "events":
            execute((edge,(-angle)%(2*np.pi)) for edge,angle in reversed(arg))
        else:
            route = rotation(local,arg,ref,np.pi)
            execute(route)
            seen.append(abs(states[pair])**2)
            states[pair] = 0
            execute(reversed(route))
    return states, np.asarray(seen), dict(events=count, pulse_angle=duration, queries=len(seen))


def step_columns(full, nominal):
    actual = full.columns(nominal)
    count = len(full.slots)
    expected = np.column_stack([real_step(v.reshape((*full.shape,8))).ravel() for v in np.eye(count)])
    embedded = np.zeros_like(actual)
    embedded[full.slots] = expected
    error = float(np.linalg.norm(actual-embedded,ord=2))
    require(error < 1e-10, "nominal seam word differs from the scalar gate")
    return error


def replay(spec):
    host,data,ports,front,bank = patches(spec["level"])
    hits, clock_info = clock_laws(host,front,bank,spec)
    waits = np.arange(spec["alarm_after_periods"]+1)
    weights = {}
    for r in spec["thresholds"]:
        pmf,timeout = stopping_law(r,int(waits[-1]))
        require(np.array_equal(hits[r],pmf), "physical clock first-passage law differs from forecast")
        weights[r] = pmf,timeout
    sector, emitted, emitter = finite_source(data,host,ports,spec["fuel"])
    full = Compiler(data,spec["marker"],spec["shape"])
    local = Compiler(data,spec["marker"],[1,1,1,full.shape[3]])
    neg = list(range(data.n-len(local.slots),data.n))
    events = sectors.SeamEvents(sector,8)
    programs = {name:full.build(name) for name in spec["members"]}
    resources = {name:{k:v for k,v in full.resources(word).items() if k!="port_counts"}
                 for name,word in programs.items()}
    aperture = [translated_word(full,local,pos) for pos in spec["receivers"]]
    metadata = dict(active_ports=ports,data_seams=len(data.seams), data_pair_states=sector.P,
        word_operator_error_below_1e_10=step_columns(full,programs["nominal"])<1e-10,
        clock=clock_info, emitter=emitter, words=resources,
        aperture_events=[len(word) for word in aperture], refinement=refinement(),instrument={})
    away = ~sector.adjacent
    home = data.seams.index((0,1))
    opened,home_laws = {},{}
    max_nominal = 0.
    for order in spec["private_orders"]:
        wave = np.repeat(np.exp(2j*np.pi*order*np.arange(full.shape[3])/full.shape[3]),8)/np.sqrt(len(local.slots))
        launch,early,send = execute_source(local,events,neg,wave,"native",emitted[1:])
        live = abs(emitted[0])**2 + np.vdot(launch,launch).real
        require(abs(live+early.sum()-1)<1e-11, "sender lost probability")
        for name,word in programs.items():
            snapshots = np.empty((sector.P,len(waits)),complex)
            state = launch.copy()
            for n in waits:
                snapshots[:,n] = state
                if n < waits[-1]:
                    for edge,q in word:
                        events.apply(state,edge,q*np.pi/4)
            for site,(pos,route) in enumerate(zip(spec["receivers"],aperture,strict=True)):
                routed = snapshots.copy()
                for edge,q in route:
                    events.apply(routed,edge,q*np.pi/4)
                state,late,read = receive_columns(local,events,neg,wave,routed)
                seams = abs(state[sector.seam_states])**2
                quiet = abs(emitted[0])**2+(abs(state[away])**2).sum(axis=0)
                end = np.vstack((late,quiet,seams))
                require(np.max(abs(end.sum(axis=0)-live))<1e-10, "receiver lost probability")
                ideal = abs(scalar_kernel(full.shape,order,pos,waits))**2/16
                for r,(pmf,timeout) in weights.items():
                    law = np.r_[early,live*timeout,end@pmf]
                    key = f"{name}_{order}_{site}_{r}"
                    bins = counters(law,spec["shots"])
                    opened[key] = bins
                    idx = len(early)+1+len(late)+1+home
                    freq = bins[idx]/spec["shots"]
                    require(abs(freq-law[idx])<=1/spec["shots"], "counter error exceeded one exposure")
                    home_laws[key] = float(freq)
                    if name == "nominal":
                        max_nominal = max(max_nominal,abs(freq-pmf@ideal))
                metadata["instrument"][str(order)] = dict(sender_events=send["events"],
                    sender_pulse_angle=send["angle_sum"], sender_queries=send["reader_queries"],receiver=read)
    require(max_nominal <= emitter["bound"]+1/spec["shots"], "field transfer exceeds emitter and read budget")
    separations = {}
    for name in spec["members"][1:]:
        separations[name] = max(abs(v-home_laws[k.replace(name,"nominal",1)]) for k,v in home_laws.items()
                                if k.startswith(name+"_"))
    metadata.update(max_nominal_home_error=float(max_nominal), control_home_separation=separations)
    return opened,metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",type=Path,help="write a fresh archive directory")
    args = parser.parse_args()
    from .clocked_pins import verify_sources
    verify_sources()
    require(receiver_pins()==load(HERE/"receiver_sources.json") and
            preparation_pins()==load(HERE/"preparation_sources.json"), "changed source or receiver")
    spec = load(HERE/"clocked_protocol.json")["fixture"]
    opened,metadata = replay(spec)
    root = args.record or HERE/"clocked_records"
    if args.record:
        root.mkdir(parents=True,exist_ok=False)
    archive(root/"counts.npz",opened,bool(args.record))
    receipt(root/"result.json",metadata,bool(args.record))
    print(f"{len(opened)} complete distributions, maximum nominal home error {metadata['max_nominal_home_error']:.6g}")


if __name__ == "__main__":
    main()
