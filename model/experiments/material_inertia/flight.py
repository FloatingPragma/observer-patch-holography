"""Forecast an entire local-cluster position distribution from the pair-binding inbox."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from tools.strict_json import load
from . import instrument
from .__main__ import HERE,clean,frozen_files
from .challenge import challenge_inputs,custody,pins as challenge_pins,prediction_stage,write_or_check
from .source import Ring,require,spectral_measure
from .prediction import native_transport_bound
from .detector import position_bins


def pins():
    files = [HERE/n for n in ('challenge_sources.json','flight.py','flight_protocol.json','detector.py')]
    law = HERE.parents[1]/'candidate_laws'/'contrast_memory_spherical_v1'
    files += [law/'reader.py',law/'sectors.py',law/'common.py']
    return {p.relative_to(HERE.parents[2]).as_posix():
            hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files}


def probabilities(returns, start):
    # The overlap with a translated compact pattern is the Fourier inverse of its
    # complex return in every character. All other microscopic patterns form the last bin.
    compact = np.roll(abs(np.fft.ifft(returns))**2, start)
    remainder = 1-compact.sum()
    require(remainder >= -1e-10,'cluster probability exceeds one')
    return np.r_[compact,max(0.,remainder)]


def launch_probability(size, start):
    """The zero-time propagator is the identity on the compact preparation at start."""
    p = np.zeros(size + 1)
    p[start] = 1.
    return p


def forecast(spec,chemical,rate):
    size = spec['size']
    angles = 1-np.cos(2*np.pi*np.arange(size)/size)
    programs = {}
    for n in spec['orders']:
        ratio = sum(chemical['inertia'][str(n)])/2
        ticks = np.rint(np.array(spec['scaled_waits'])*ratio/rate).astype(np.int64)
        predictions = [launch_probability(size,spec['starting_port']) if tick == 0 else
                       probabilities(np.exp(-1j*tick*rate*angles/ratio),spec['starting_port'])
                       for tick in ticks]
        programs[str(n)] = (ticks,predictions)
    return programs


def collect(spec,base,programs,kappa):
    raw = {}
    bank = {}
    for n in spec['orders']:
        size = spec['size']
        ring = Ring.build(size,n)
        bins = position_bins(ring)
        require(np.array_equal(bins[ring.orbits[0]],np.arange(size)) and
                np.count_nonzero(bins!=size)==size, 'Fourier orbit differs from the public detector')
        ticks,_ = programs[str(n)]
        for member in ('contrast','assisted'):
            H = ring.hamiltonian(member,kappa)
            modes = [spectral_measure(H,ring.frame(k)) for k in range(size)]
            require(sum(len(E) for E,_ in modes)==len(ring.states),'incomplete position spectrum')
            ground = modes[0][0][0]
            key = f'{member}_{n}'
            draws = []
            for tick in ticks:
                t = tick/base['tick_divisor']
                if tick == 0:
                    p = launch_probability(size,spec['starting_port'])
                else:
                    returns = np.array([np.exp(-1j*t*(E-ground))@w for E,w in modes])
                    p = probabilities(returns,spec['starting_port'])
                edges = instrument.counts(p.cumsum(),spec['shots'])[:,1]
                draws.append(np.diff(edges,prepend=np.int32(0)))
            raw[key] = np.array(draws)
            bank[key] = (np.array([E[0]-ground for E,_ in modes]),np.array([w[0] for _,w in modes]))
    return raw,bank


def analyze(spec,base,programs,raw,bank,chemical,rate):
    result = {}
    for n in spec['orders']:
        ticks,targets = programs[str(n)]
        ratio = sum(chemical['inertia'][str(n)])/2
        predicted = (1-np.cos(2*np.pi*np.arange(spec['size'])/spec['size']))/ratio
        for member in ('contrast','assisted'):
            key = f'{member}_{n}'
            a = raw[key]
            require(a.dtype==np.int32 and a.shape==(len(ticks),spec['size']+1) and np.all(a>=0)
                    and np.all(a.sum(axis=1)==spec['shots']),'incomplete position record')
            observed = a/spec['shots']
            errors = .5*abs(observed-np.array(targets)).sum(axis=1)
            travel = 1-observed[:,-1]-observed[:,spec['starting_port']]
            # Separate source audit. Orthogonality of the characters and all omitted
            # spectral components bounds the state norm, hence the detector TV distance.
            delta,weight = bank[key]
            mismatch = delta/base['tick_divisor']-predicted*rate
            limits = [min(1.,np.sqrt(np.mean(weight*np.minimum(4.,(mismatch*tick)**2)+4*(1-weight)))
                          +(spec['size']+1)/(2*spec['shots'])) for tick in ticks]
            require(np.all(errors<=np.array(limits)+2e-8),'transport exceeds spectral/digitization bound')
            result[key] = dict(ticks=ticks.tolist(),total_variation=errors.tolist(),
                               compact_away=travel.tolist(),no_compact=observed[:,-1].tolist(),
                               passes=bool(max(errors)<=spec['max_total_variation']
                               and max(travel)>=spec['nontrivial_travel_gate']))
            if member == 'contrast':
                bounds = [native_transport_bound(chemical,spec['size'],n,int(tick),rate,base,spec['shots'])
                          for tick in ticks]
                require(np.all(errors <= np.array(bounds)+2e-8), 'flight exceeds binding-only error enclosure')
                result[key]['binding_only_tv_bound'] = bounds
    return clean(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    args = parser.parse_args()
    require(frozen_files()==load(HERE/'sources.json') and challenge_pins()==load(HERE/'challenge_sources.json')
            and pins()==load(HERE/'flight_sources.json'),'changed flight source')
    spec = load(HERE/'flight_protocol.json')
    base = load(HERE/'protocol.json')
    prospect = load(HERE/'challenge_protocol.json')
    kappa,offset = challenge_inputs()
    sealed = HERE/'records'/'challenge'
    chemical = prediction_stage(sealed,prospect,base,kappa,offset,False)
    write_or_check(sealed/'prediction.json',chemical,False)
    with np.load(sealed/'binding.npz',allow_pickle=False) as inbox:
        rate = instrument.clock_rate(inbox['clock'],base)
    programs = forecast(spec,chemical,rate)
    root = args.record or HERE/'records'/'flight'
    if args.record:
        root.mkdir(exist_ok=False,parents=True)
    raw,bank = collect(spec,base,programs,kappa)
    saved = custody(root/'counts.npz',raw,bool(args.record))
    result = analyze(spec,base,programs,saved,bank,chemical,rate)
    write_or_check(root/'result.json',result,bool(args.record))
    print(json.dumps({k:[v['passes'],max(v['total_variation']),max(v['compact_away'])] for k,v in result.items()}))


if __name__=='__main__':
    main()
