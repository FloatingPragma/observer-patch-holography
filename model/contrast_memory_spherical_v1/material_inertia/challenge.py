"""Two-stage equal-load interferometry with a law-blind prediction and a competing repair."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from tools.strict_json import load
from . import instrument
from .__main__ import HERE,check_json,clean,frozen_files
from .prediction import infer
from .source import Ring,require,spectral_measure


def pins():
    names = ('sources.json','challenge.py','prediction.py','challenge_protocol.json')
    return {n:hashlib.sha256((HERE/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}


def challenge_inputs():
    digest = hashlib.sha256((HERE.parent/'carriers/FEDERATION_LEVEL1.npz').read_bytes()).digest()
    return 4+digest[0]%8, int(digest[1]%17)-8


def amplitude(spectrum,times):
    energies,weights = spectrum
    return np.exp(-1j*np.outer(times,energies))@weights


def binding_inbox(spec,base,member,kappa,offset):
    ring = Ring.build(spec['size'],2)
    H = ring.hamiltonian(member,kappa)
    pair = spectral_measure(H,ring.frame(0))
    single = Ring.build(spec['size'],1)
    free = spectral_measure(single.hamiltonian(member,kappa),single.frame(0))
    require(len(free[0])==1 and pair[1][0]>=base['dominant_floor'], 'binding source unresolved')
    # The path with two independent singles gives the product return. Both active
    # arms have two loads; include their offsets before forming the frequency difference.
    energies = 2*(free[0][0]+offset) - (pair[0]+2*offset)
    require(max(abs(energies))/base['tick_divisor'] < np.pi-.2,'binding read aliases')
    return instrument.quadratures(amplitude((energies,pair[1]),instrument.tick_menu(base)/base['tick_divisor']),
                                  base['shots'])


def motion(spec,base,member,kappa,offset):
    size = spec['size']
    times = instrument.tick_menu(base)/base['tick_divisor']
    raw, weights = {},{}
    for n in [1,*spec['orders']]:
        ring = Ring.build(size,n)
        H = ring.hamiltonian(member,kappa)
        rest = spectral_measure(H,ring.frame(0))
        for ell in spec['momenta']:
            moving = spectral_measure(H,ring.frame(ell))
            require(max(moving[0])-min(rest[0]) < (np.pi-.2)*base['tick_divisor'] and
                    max(rest[0])-min(moving[0]) < (np.pi-.2)*base['tick_divisor'], 'motion read aliases')
            # Factor out a common phase before evaluating the two-copy overlap.
            # The parked arm has the same load in both paths.
            center = rest[0][0]+n*offset
            ref = amplitude((rest[0]+n*offset-center,rest[1]),times)
            signal = amplitude((moving[0]+n*offset-center,moving[1]),times)
            key = f'{member}_{n}_{ell}'
            raw[key] = instrument.quadratures(ref.conj()*signal,base['shots'])
            weights[key] = float(rest[1][0]*moving[1][0])
    return raw,weights


def custody(path,raw,record):
    if record:
        with path.open('xb') as out:
            np.savez_compressed(out,**raw)
    with np.load(path,allow_pickle=False) as stored:
        require(sorted(stored.files)==sorted(raw) and all(stored[n].dtype==np.int32 and
                np.array_equal(stored[n],v) for n,v in raw.items()), 'changed challenge counts')
        return dict(stored)


def write_or_check(path,result,record):
    result = clean(result)
    if record:
        with path.open('x',encoding='utf8',newline='\n') as out:
            json.dump(result,out,indent=1,allow_nan=False)
            out.write('\n')
    check_json(path,result)


def prediction_stage(root,spec,base,kappa,offset,record):
    raw = {'clock':instrument.clock(base)}
    for member in ('contrast','assisted'):
        raw[member] = binding_inbox(spec,base,member,kappa,offset)
    recorded = custody(root/'binding.npz',raw,record)
    require(np.array_equal(recorded['contrast'],recorded['assisted']), 'competing binding inboxes differ')
    apparatus = {k:base[k] for k in ('shots','clock_ticks','dyadic_orders','dominant_floor','phase_gate')}
    return infer(recorded['clock'],recorded['contrast'],spec['size'],apparatus,spec['orders'])


def score(interval,predicted,precision):
    plo,phi = predicted
    instrument.require_positive_order(plo,phi)
    require(not isinstance(precision,(bool,np.bool_)) and np.isfinite(precision) and precision>=0,
            'invalid forecast precision')
    if interval is None:
        return 'unresolved'
    low,high = instrument.positive_interval(interval)
    if high/low-1 > 2*precision:
        return 'unresolved'
    return 'excluded' if high<plo or low>phi else 'compatible'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('predict','test'),required=True)
    parser.add_argument('--record',type=Path)
    args = parser.parse_args()
    require(pins()==load(HERE/'challenge_sources.json') and frozen_files()==load(HERE/'sources.json'),
            'changed binding challenge')
    spec = load(HERE/'challenge_protocol.json')
    base = load(HERE/'protocol.json')
    kappa,offset = challenge_inputs()
    root = args.record or HERE/'records'/'challenge'
    register = bool(args.record) and args.stage=='predict'
    if register:
        root.mkdir(exist_ok=False,parents=True)
    forecast = prediction_stage(root,spec,base,kappa,offset,register)
    write_or_check(root/'prediction.json',forecast,register)
    if args.stage=='predict':
        print(json.dumps(clean(forecast),indent=1))
        return
    # The immutable forecast exists and has been replayed before any holdout generation.
    raw,weights = {},{}
    for member in ('contrast','assisted'):
        r,w = motion(spec,base,member,kappa,offset)
        raw.update(r)
        weights.update(w)
    saved = custody(root/'motion.npz',raw,bool(args.record))
    with np.load(root/'binding.npz',allow_pickle=False) as original:
        rate = instrument.clock_rate(original['clock'],base)
    lines = {key:instrument.read_line(counts,rate,base) if weights[key]>=base['dominant_floor'] else None
             for key,counts in saved.items()}
    rows = {}
    for member in ('contrast','assisted'):
        for n in spec['orders']:
            for ell in spec['momenta']:
                key = f'{member}_{n}_{ell}'
                a,b = lines[f'{member}_1_{ell}'],lines[key]
                interval = instrument.quotient(a,b) if a is not None and b is not None else None
                rows[key] = dict(inertia=interval,verdict=score(interval,forecast['inertia'][str(n)],base['relative_precision']))
    result = dict(source_coupling=kappa,per_load_shift=offset,readings=rows)
    write_or_check(root/'test.json',result,bool(args.record))
    print(json.dumps(clean(result),indent=1))


if __name__=='__main__':
    main()
