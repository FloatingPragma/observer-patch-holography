"""Independent detector forecasts for the routed field and the full-sector free release."""
from math import comb, factorial

import numpy as np
from scipy import sparse

from tools.strict_json import load
from .clocked_apparatus import patches
from .clocked_field import stopping_law
from .compiler import Compiler
from .clocked_run import finite_source
from .preparation import HERE, execute_source, make_source, rotation
from model.candidate_laws.contrast_memory_spherical_v1 import sectors
from model.executable_model import carrier


def logical_step(state, member):
    coin = state.copy()
    for mask in range(1,8):
        if member=="shift_only" or (mask==7 and member=="omit_mask_seven"):
            continue
        alt = coin[...,np.arange(8)^mask]
        z = np.exp(-.25j*np.pi)
        coin = ((1+z)*coin+(1-z)*alt)/2
    out = np.empty_like(coin)
    for s in range(8):
        out[...,s^1] = np.roll(coin[...,s],1 if s%2==0 else -1,axis=s//2)
    return out


def sharp_receiver(comp, neg, wave, states):
    stops = []
    home = comp.slots[0]
    def run(seq):
        for edge,t in seq:
            a,b = comp.carrier.seams[edge]
            va,vb = states[[a,b]]
            z = np.exp(-1j*t)
            common = np.exp(-8j*t)
            states[a] = common*((1+z)*va+(1-z)*vb)/2
            states[b] = common*((1-z)*va+(1+z)*vb)/2
    for action,arg in reversed(list(make_source(comp,neg,wave,8,"native"))):
        if action=="events":
            run((e,(-t)%(2*np.pi)) for e,t in reversed(arg))
        else:
            route = rotation(comp,arg,home,np.pi)
            run(route)
            stops.append(abs(states[home])**2)
            states[home] = 0
            run(reversed(route))
    out = np.zeros((len(comp.carrier.seams),states.shape[1]))
    for i,(a,b) in enumerate(comp.carrier.seams):
        if a==comp.marker:
            out[i] = abs(states[b])**2
    return np.vstack((stops,(abs(states)**2).sum(axis=0)-out.sum(axis=0),out))


def sharp_comparison():
    spec = load(HERE/"clocked_protocol.json")["fixture"]
    _,data,_,_,_ = patches(0)
    full = Compiler(data,0,spec["shape"])
    small = Compiler(data,0,[1,1,1,full.shape[-1]])
    neg = list(range(data.n-len(small.slots),data.n))
    cap = spec["alarm_after_periods"]
    max_tv = 0.
    expected_keys = set()
    with np.load(HERE/"clocked_records/counts.npz",allow_pickle=False) as archive:
        for order in spec["private_orders"]:
            wave = np.repeat(np.exp(2j*np.pi*order*np.arange(full.shape[3])/full.shape[3]),8)
            wave /= np.linalg.norm(wave)
            early = np.r_[.25,.25,.25*abs(wave)**2]
            for member in spec["members"]:
                grid = np.zeros((*full.shape,8),complex)
                grid[0,0,0] = wave.reshape(-1,8)/2
                saved = []
                for _ in range(cap+1):
                    saved.append(grid)
                    grid = logical_step(grid,member)
                for site,pos in enumerate(spec["receivers"]):
                    states = np.zeros((data.n,cap+1),complex)
                    for n,grid in enumerate(saved):
                        localized = grid.copy()
                        localized[tuple(pos)],localized[0,0,0] = grid[0,0,0],grid[tuple(pos)]
                        states[full.slots,n] = localized.ravel()
                    end = sharp_receiver(small,neg,wave,states)
                    for r in spec["thresholds"]:
                        key = f"{member}_{order}_{site}_{r}"
                        w,timeout = stopping_law(r,cap)
                        law = np.r_[early,.25*timeout,end@w]
                        assert law.min()>=-1e-13 and abs(law.sum()-1)<1e-12
                        bins = archive[key]
                        assert bins.shape==law.shape and bins.dtype==np.int32
                        assert bins.min()>=0 and bins.sum()==spec["shots"]
                        tv = abs(bins/spec["shots"]-law).sum()/2
                        assert tv <= np.sqrt(sum(data.degree[:2])-2)/(2*spec["fuel"])+len(bins)/(2*spec["shots"])
                        max_tv = max(max_tv,tv)
                        expected_keys.add(key)
        assert expected_keys==set(archive.files)
    return max_tv


def test_complete_archived_distributions_against_sharp_code():
    sharp_comparison()


def test_free_release_archive_against_contact_energies_and_polynomial_evolution():
    spec = load(HERE/"clocked_protocol.json")["fixture"]
    protocol = load(HERE/"free_marker_protocol.json")
    host, data, ports, _, _ = patches(spec["level"])
    pair, emitted, _ = finite_source(data,host,ports,spec["fuel"])
    source = Compiler(data,0,[1,1,1,spec["shape"][3]])
    events = sectors.SeamEvents(pair,8)
    neg = list(range(data.n-len(source.slots),data.n))
    orders = spec["private_orders"]
    launches, stops = [], []
    for order in orders:
        wave = np.repeat(np.exp(2j*np.pi*order*np.arange(source.shape[3])/source.shape[3]),8)
        state, early, _ = execute_source(source,events,neg,wave/np.linalg.norm(wave),"native",emitted[1:])
        launches.append(state)
        stops.append(early)
    initial = np.column_stack(launches)
    cap = spec["alarm_after_periods"]
    degree = host.degree[ports]
    expected_keys = set()
    with np.load(HERE/"free_marker_records/counts.npz",allow_pickle=False) as archive:
        for member in protocol["members"]:
            active = data
            if member=="marker_contacts_closed":
                active = carrier.Carrier("isolated marker",data.n,[edge for edge in data.seams if 0 not in edge])
            sec = sectors.PairSector(active)
            # Count open edges with exactly one occupied end, plus the original host reference.
            singles = active.degree[sec.aa]+active.degree[sec.bb]-2*sec.adjacent
            reference = 12-degree[sec.aa]-degree[sec.bb]
            h = sec.exchange()+sparse.diags(8.*(singles+reference-12))
            step = -1j*h/8
            # Hermiticity bounds the order-48 exponential remainder by ||H/8||**49 / 49!.
            radius = float(np.max(abs(step).sum(axis=1)))
            assert cap*radius**49/factorial(49)<1e-13
            state = initial.copy()
            records = []
            for n in range(cap+1):
                found = abs(state[:data.n-1])**2  # lexicographic pairs (0,1),...,(0,n-1)
                none = abs(emitted[0])**2+(abs(state[data.n-1:])**2).sum(axis=0)
                records.append(np.vstack((found,none)))
                if member!="hold" and n<cap:
                    term = state.copy()
                    for j in range(1,49):
                        term = step@term/j
                        state += term
            for r in protocol["thresholds"]:
                weights = np.zeros(cap+1)
                if r:
                    for n in range(r,cap+1):
                        weights[n] = comb(n-1,r-1)/2**n
                else:
                    weights[0] = 1
                timed_out = sum(comb(cap,j) for j in range(r))/2**cap
                end = np.tensordot(weights,records,axes=1)
                for i,order in enumerate(orders):
                    early = stops[i]
                    law = np.r_[early,(1-early.sum())*timed_out,end[:,i]]
                    key = f"{order}_{member}_{r}"
                    bins = archive[key]
                    assert bins.shape==law.shape and bins.dtype==np.int32 and bins.min()>=0
                    assert bins.sum()==spec["shots"] and abs(law.sum()-1)<1e-11
                    # Test every cumulative detector boundary to half a count, allowing numerical roundoff.
                    assert np.max(abs(bins.cumsum()-spec["shots"]*law.cumsum()))<=.50001
                    expected_keys.add(key)
        assert expected_keys==set(archive.files)
