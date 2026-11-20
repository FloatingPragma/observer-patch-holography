"""Acquire, decode and independently regenerate the entire pair-conversion experiment."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from tools.strict_json import load
from ...contrast_memory_spherical_v1 import carrier
from ..carrier import require
from . import apparatus, source

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREFIX = HERE.relative_to(ROOT).as_posix()
FILES = ['model/__init__.py','model/requirements.txt','model/contrast_memory_spherical_v1/__init__.py',
         'model/contrast_memory_spherical_v1/carrier.py','model/scalar_channel/__init__.py','model/scalar_channel/carrier.py',
         'tools/strict_json.py',*[PREFIX+'/'+f for f in ('__init__.py','__main__.py','source.py','apparatus.py','protocol.json')],
         *[f'model/contrast_memory_spherical_v1/carriers/FEDERATION_LEVEL{i}.npz' for i in range(3)]]


def pins():
    def digest(f):
        raw = (ROOT/f).read_bytes()
        return hashlib.sha256(raw if f.endswith('.npz') else raw.replace(b'\r\n',b'\n')).hexdigest()
    return {f:digest(f) for f in FILES}


def acquire(spec):
    rng = np.random.Generator(np.random.PCG64(spec['seed']))
    shots = spec['shots']
    times = np.arange(spec['samples'])*spec['clock_step']
    records = {'clock':apparatus.reference(rng,times,shots)}
    for level in spec['gap_carriers']:
        car = carrier.carrier(level)
        v,_ = source.rest_witness(car)
        Q,L = source.one_body(car)
        trajectory = expm_multiply(-1j*Q,v,start=0,stop=times[-1],num=len(times),
                                   traceA=-1j*Q.diagonal().sum())
        records[f'band_{level}'] = apparatus.path_counts(rng,trajectory@v,shots)
        if level == 'cell':
            for sign,h in [('plus',Q),('minus',Q),('neutral',L)]:
                values,U = np.linalg.eigh(h.toarray())
                amp = np.exp(-1j*np.outer(times,values)) @ abs(U[0])**2
                records[sign] = apparatus.path_counts(rng,amp,shots)

    for key in spec['carriers']:
        car = carrier.carrier(key)
        labels = source.ports_to_regions(car)
        for law in ['exchange',*spec['controls']]:
            H = source.generator(car,law)
            for mode in spec['states'] if law == 'exchange' else ['ordered']:
                loc,work = [],[]
                state = source.initial(car,key,mode)
                for wave in source.evolve(H,state,spec['waits']):
                    X,Y = wave.reshape(2,car.n,car.n)
                    require(np.max(abs(Y-Y.T))<2e-12,'lost neutral Bose symmetry')
                    loc.append(apparatus.counts(rng,source.position(X,Y,labels),shots))
                    probes,scale = source.energy_instruments(car,X,Y,law)
                    audit = np.vdot(wave,H@wave).real
                    require(abs(scale @ (probes[:,1]+2*probes[:,2])-audit)<2e-9,'wrong energy instrument')
                    work.append(apparatus.counts(rng,probes,shots))
                name = f'{key}_{law}_{mode}'
                records[name+'_position'] = np.array(loc)
                records[name+'_energy'] = np.array(work)
    return records


def decode(spec, records):
    shots = spec['shots']
    ticks = apparatus.read_clock(records['clock'],shots)
    require(len(ticks) == spec['samples'], 'wrong spectrum clock menu')
    waits = ticks[np.rint(np.array(spec['waits'])/spec['clock_step']).astype(int)]
    spectra = {name:apparatus.lines(records[name],ticks,spec)
               for name in ['plus','minus','neutral',*[f'band_{k}' for k in spec['gap_carriers']]]}
    rows = {}
    for key in spec['carriers']:
        car = carrier.carrier(key)
        regions = int(source.ports_to_regions(car).max())+1
        i,j = np.triu_indices(regions)
        for law in ['exchange',*spec['controls']]:
            for mode in spec['states'] if law=='exchange' else ['ordered']:
                name = f'{key}_{law}_{mode}'
                raw = records[name+'_position']
                energy = records[name+'_energy']
                require(raw.shape==(len(waits),regions**2+len(i)) and energy.shape[0]==len(waits)
                        and energy.shape[2]==3,'bad setting dimensions')
                for data in (raw,energy):
                    require(data.dtype.kind in 'iu' and np.all(data>=0) and np.all(data.sum(axis=-1)==shots),
                            'incomplete outcome counts')
                radiation = raw[:,regions**2:]/shots
                rows[name] = dict(neutral=radiation.sum(axis=1).tolist(),split=radiation[:,i!=j].sum(axis=1).tolist(),
                                  signed=(raw[:,:regions**2].sum(axis=1)/shots).tolist(),
                                  **apparatus.energies(car,energy,shots,law))
    return dict(receiver_waits=waits.tolist(),spectra=spectra,histories=rows)


def gates(spec, result):
    g = spec['gates']
    spectrum = result['spectra']
    root = np.sqrt(5)
    expect = dict(plus=[5-root,4,5+root,10],minus=[5-root,4,5+root,10],neutral=[0,5-root,6,5+root])
    for k in spec['gap_carriers']:
        expect[f'band_{k}'] = [5-root]
    checks = {}
    for name,expected in expect.items():
        row = spectrum[name]
        weights = ([3,5,3,1] if name in ('plus','minus') else [1,3,5,3]) if len(expected)==4 else [12]
        checks['spectrum_'+name] = (len(row['frequencies'])==len(expected) and
            np.max(abs(np.array(row['frequencies'])-expected))<g['line_error'] and
            np.max(abs(np.array(row['weights'])-np.array(weights)/12))<.01 and row['heldout_rms']<g['heldout_rms'])
    rows = result['histories']
    for key in spec['carriers']:
        def read(law,mode,obs):
            return np.asarray(rows[f'{key}_{law}_{mode}'][obs])
        r = read('exchange','ordered','neutral')
        checks[f'{key}_conversion'] = r[0]==0 and max(r)>g['min_conversion']
        checks[f'{key}_inverse'] = max(read('exchange','neutral','signed'))>g['min_inverse']
        checks[f'{key}_phase'] = max(abs(2*r-read('exchange','even','neutral')))<g['max_phase_law_error']
        checks[f'{key}_dark'] = not any(read('exchange','odd','neutral'))
        checks[f'{key}_off'] = not any(read('dephased_contact','ordered','neutral'))
        checks[f'{key}_parked'] = not any(read('parked_neutral','ordered','split'))
    if 0 in spec['carriers']:
        federation = rows['0_exchange_ordered']
        checks['separate_receiver_cells'] = max(federation['split'])>g['min_split_cells']
        checks['remote_field_energy'] = max(federation['remote_neutral_energy'])>g['min_remote_energy']
    checks['energy_balance'] = all(max(abs(np.array(r['total_energy'])-r['total_energy'][0]))<g['energy_error']
                                   for r in rows.values())
    return {k:bool(v) for k,v in checks.items()}


def compare(actual, stored):
    if isinstance(actual,dict):
        require(type(stored) is dict and actual.keys()==stored.keys(),'altered fields')
        for k,v in actual.items():
            compare(v,stored[k])
    elif isinstance(actual,list):
        require(type(stored) is list and len(actual)==len(stored),'altered observation list')
        for v,w in zip(actual,stored):
            compare(v,w)
    elif type(actual) is float:
        require(type(stored) is float and abs(stored-actual)<1e-8,'altered reading')
    else:
        require(type(stored) is type(actual) and stored==actual,'altered result')


def check_counts(root, draws, *, record=False, keys=None):
    if record:
        root.mkdir(exist_ok=False)
        with (root/'counts.npz').open('xb') as f:
            np.savez_compressed(f,**draws)
    with np.load(root/'counts.npz',allow_pickle=False) as archive:
        expected = draws.keys() if keys is None else keys
        require(sorted(archive.files)==sorted(expected),'incomplete measurements')
        saved = {}
        for k,v in draws.items():
            a = archive[k]
            require(a.dtype==np.int32 and np.array_equal(v,a),'changed counts: '+k)
            saved[k] = a
    return saved


def check_result(root, results, record):
    if record:
        with (root/'result.json').open('x',encoding='utf8',newline='\n') as f:
            json.dump(results,f,indent=1,allow_nan=False)
            f.write('\n')
    else:
        compare(results,load(root/'result.json'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--record',type=Path)
    group.add_argument('--cell-only',action='store_true',help='replay one-cell histories and resting profiles')
    args = parser.parse_args()
    require(load(HERE/'sources.json')==pins(),'source or protocol differs')
    spec = load(HERE/'protocol.json')
    if args.cell_only:
        spec = dict(spec,carriers=['cell'])
    draws = acquire(spec)
    dest = args.record or HERE/'records'
    keys = None
    if args.cell_only:
        keys = set(draws) | {'0_'+k[5:] for k in draws if k.startswith('cell_')}
    saved = check_counts(dest,draws,record=bool(args.record),keys=keys)
    results = decode(spec,saved)
    results['gates'] = gates(spec,results)
    if args.cell_only:
        full = load(dest/'result.json')
        for section in ('spectra','receiver_waits','histories'):
            stored = full[section]
            if section=='histories':
                stored = {k:v for k,v in stored.items() if k.startswith('cell_')}
            compare(results[section],stored)
        for k,v in results['gates'].items():
            if k!='energy_balance':
                compare(v,full['gates'][k])
    else:
        check_result(dest,results,bool(args.record))
    print(json.dumps(results['gates'],indent=2))


if __name__=='__main__':
    main()
