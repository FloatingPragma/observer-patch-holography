"""Independent detector forecasts in the ideal sharp-pair code, including every failure."""
import numpy as np

from tools.strict_json import load
from .clocked_apparatus import patches
from .clocked_field import stopping_law
from .compiler import Compiler
from .preparation import HERE, make_source, rotation


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
