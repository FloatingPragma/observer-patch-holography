"""Reproduce the source, receiver records, spectral inventory and withheld composite predictions."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path

import numpy as np

from tools.strict_json import load
from . import federation, instrument
from .source import Ring, parameters, require, spectral_measure

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def frozen_files():
    names = ['__init__.py','__main__.py','source.py','instrument.py','federation.py','protocol.json']
    paths = [HERE/name for name in names]
    paths += [ROOT/'model/candidate_laws/contrast_memory_spherical_v1'/name for name in ('__init__.py','law.py')]
    paths += [ROOT/'model/executable_model/carrier.py']
    paths += [ROOT/'model/executable_model/carriers'/f'FEDERATION_LEVEL{i}.npz' for i in range(3)]
    paths += [ROOT/'model/__init__.py',ROOT/'model/requirements.txt',ROOT/'tools/strict_json.py']
    def digest(path):
        raw = path.read_bytes()
        return hashlib.sha256(raw if path.suffix=='.npz' else raw.replace(b'\r\n',b'\n')).hexdigest()
    return {p.relative_to(ROOT).as_posix():digest(p) for p in paths}


def cases(spec):
    for size in spec['rings']:
        for number in spec['occupations']:
            if size != max(spec['rings']) or number!=4 or spec['largest_ring_four']:
                yield size, number


def name(member,size,number,ell):
    return f'{member}_{size}_{number}_{ell}'


def acquire(spec, spectra_dir=None):
    raw = {'clock':instrument.clock(spec), **federation.collect(spec)}
    audit = {}
    times = instrument.tick_menu(spec)/spec['tick_divisor']
    if spectra_dir is not None:
        spectra_dir.mkdir(exist_ok=False,parents=True)
    for size, number in cases(spec):
        ring = Ring.build(size,number)
        for member in spec['laws']:
            H = ring.hamiltonian(member,spec['kappa'])
            k,hop = parameters(member,spec['kappa'])
            # Norm bound from up to 2N active swaps and 2N active assisted hops per row.
            cap = number*(6*k+2+4*hop)
            require(cap/spec['tick_divisor'] < np.pi-.2, 'first clock interval aliases source energies')
            dimension = 0
            inventory = {}
            for ell in range(size):
                frame = ring.frame(ell)
                E,w = spectral_measure(H,frame)
                require(abs(w.sum()-1)<3e-12, 'lost source overlap')
                dimension += len(w)
                if spectra_dir is not None:
                    inventory[str(ell)] = np.stack([E,w])
                for setting in spec['phases']:
                    if setting % size != ell:
                        continue
                    key = name(member,size,number,setting)
                    # Every mode contributes to both quadratures, including lines away from the cluster.
                    amplitude = np.exp(-1j*np.outer(times,E))@w
                    raw[key] = instrument.quadratures(amplitude,spec['shots'])
                    # A single ground eigenvector suffices. Never qualify a nearby line by tolerance.
                    ground = E[0]
                    weight = float(w[0])
                    audit[key] = (float(ground),weight)
            require(dimension==len(ring.states),'incomplete momentum decomposition')
            if spectra_dir is not None:
                np.savez_compressed(spectra_dir/name(member,size,number,'full'),**inventory)
    return raw,audit


def analyze(spec, raw):
    rate = instrument.clock_rate(raw['clock'],spec)
    phases = {}
    for size, number in cases(spec):
        for member in spec['laws']:
            for ell in spec['phases']:
                key = name(member,size,number,ell)
                phases[key] = instrument.read_phase(raw[key],spec)
    exposure = instrument.tick_menu(spec)[-2]*rate
    lines = {key:None if phase is None else phase.frequency(exposure) for key,phase in phases.items()}
    ratios = {}
    for size, number in cases(spec):
        if number==1:
            continue
        for member in spec['laws']:
            four = [phases[name(member,size,n,k)] for n,k in ((1,0),(1,1),(number,0),(number,1))]
            ratios[name(member,size,number,'ratios')] = instrument.band_ratios(four)
    predictions = {}
    for size in spec['rings']:
        measured = ratios[name('contrast',size,2,'ratios')]
        if measured is None or measured['inertia'] is None:
            continue
        r = measured['inertia']
        # Chebyshev recurrence from the two-body measurement, with no higher-body input.
        for n in (3,4):
            if (size,n) in cases(spec):
                predictions[name('contrast',size,n,'inertia')] = [v*v-1 if n==3 else v**3-2*v for v in r]
    flights = {}
    for level in range(3):
        trials = [raw[f'flight_{level}_{q}'] for q in (1,2)]
        for a in trials:
            require(a.dtype==np.int32 and a.shape==(len(spec['flight_ticks']),4,4)
                    and np.all(a>=0) and np.all(a.sum(axis=(1,2))==spec['shots']), 'invalid flight record')
        flights[str(level)] = dict(same_counts=bool(np.array_equal(*trials)),
                                   left_home_at_first_read=trials[0][:,1:,:].sum(axis=(1,2)).tolist(),
                                   changed_region_between_reads=(spec['shots']-trials[0].trace(axis1=1,axis2=2)).tolist())
    return dict(clock_units_per_tick=rate,lines=lines,ratios=ratios,
                higher_body_predictions=predictions,federations=flights)


def decisions(spec, data, audit):
    # The decoder receives only clock and pointer counts. Separately determine which returned
    # lines the source actually identifies as lowest eigenvalues at the stated weight floor.
    rate = data['clock_units_per_tick']
    normalized = rate*spec['tick_divisor']
    admitted = {}
    for key,(ground,weight) in audit.items():
        read = data['lines'][key]
        valid = weight >= spec['dominant_floor'] and read is not None
        if valid:
            require(abs(read[0]*normalized-ground) <= read[1]*normalized+3e-10,
                    'resolved line disagrees with the full spectrum: '+key)
        admitted[key] = bool(valid)
    verdicts = {}
    for size, number in cases(spec):
        if number == 1:
            continue
        for member in spec['laws']:
            key = name(member,size,number,'ratios')
            row = data['ratios'][key]
            qualified = all(admitted[name(member,size,n,k)] for n,k in ((1,0),(1,1),(number,0),(number,1)))
            verdicts[key] = instrument.rest_inertia_verdict(row if qualified else None,spec)
    return dict(dominant_ground_lines=admitted,rest_inertia=verdicts)


def clean(obj):
    if isinstance(obj,dict):
        return {k:clean(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj,(float,np.floating)):
        return float(f'{obj:.10g}')
    return obj


def check_json(path, expected):
    # Python equality accepts True in place of 1. Compare canonical typed JSON.
    stored = json.dumps(load(path), sort_keys=True, allow_nan=False, indent=2)
    replay = json.dumps(expected, sort_keys=True, allow_nan=False, indent=2)
    if stored != replay:
        changes = difflib.unified_diff(stored.splitlines(), replay.splitlines(),
                                       fromfile=str(path), tofile='replayed', lineterm='')
        raise ValueError('changed decoded result\n' + '\n'.join(changes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    parser.add_argument('--spectra-out',type=Path)
    args = parser.parse_args()
    require(load(HERE/'sources.json') == frozen_files(), 'changed source or registration')
    spec = load(HERE/'protocol.json')
    raw,audit = acquire(spec,args.spectra_out)
    root = args.record or HERE/'records'
    if args.record:
        root.mkdir(exist_ok=False,parents=True)
        np.savez_compressed(root/'counts.npz',**raw)
    with np.load(root/'counts.npz',allow_pickle=False) as archive:
        require(sorted(archive.files)==sorted(raw),'missing or extra acquisition')
        for k,v in raw.items():
            require(archive[k].dtype==np.int32 and np.array_equal(archive[k],v),'changed counts: '+k)
    result = analyze(spec,raw)
    result['audit'] = decisions(spec,result,audit)
    saved = clean(result)
    if args.record:
        (root/'result.json').write_text(json.dumps(saved,indent=1,allow_nan=False)+'\n',encoding='utf8')
    else:
        check_json(root/'result.json',saved)
    print(json.dumps(dict(rest_inertia=saved['audit']['rest_inertia'],federations=saved['federations']),indent=1))


if __name__=='__main__':
    main()
