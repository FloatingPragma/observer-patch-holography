"""Offline acquisition, blind frequency readout and complete-spectrum checks."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from tools.strict_json import load

from ..carrier import require
from ..cage import construction as cell
from . import flight, instruments as ins, source

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCES = ('model/__init__.py', 'model/requirements.txt', 'model/scalar_channel/__init__.py',
           'model/scalar_channel/cage/__init__.py', 'model/scalar_channel/carrier.py', 'model/scalar_channel/cage/construction.py',
           'tools/strict_json.py', *('model/scalar_channel/refinement_spectrum/'+f for f in
           ('__init__.py','__main__.py','source.py','instruments.py','flight.py','protocol.json')),
           *(f'model/contrast_memory_spherical_v1/carriers/FEDERATION_LEVEL{i}.npz' for i in range(3)))


def pins():
    def digest(name):
        data = (ROOT/name).read_bytes()
        if not name.endswith('.npz'):
            data = data.replace(b'\r\n', b'\n')
        return hashlib.sha256(data).hexdigest()
    return {path: digest(path) for path in SOURCES}


def carrier_maps(spec):
    rows = {}
    for key in spec['carriers']:
        car = cell.cell() if key == 'cell' else cell.federation(key)
        names = spec['cell_profiles'] if key == 'cell' else spec['federation_profiles']
        profiles = {name: source.outer_profile(car, name)[2]/cell.HOP for name in names}
        if key != 'cell':
            for name, value in profiles.items():
                require(abs(value-rows['cell']['outer_rates'][name]) < 2e-12, 'carrier record law differs')
        rows[str(key)] = dict(outer_modes=sum(car['C'].shape), outer_rates=profiles)
    return rows


def experiment(spec):
    rng = np.random.Generator(np.random.PCG64(spec['seed']))
    t = np.arange(spec['times'])*spec['clock_step']/cell.HOP
    maps = carrier_maps(spec)
    # On the half-occupied dimer the clock splitter addresses its negative orbital.
    clock = ins.acquire(rng, ins.fuel_law(ins.clock_amplitude(cell.HOP*t), (0,2,1,3)), spec['shots'])
    draws = {'clock': clock}
    for name, ratio in maps['cell']['outer_rates'].items():
        p = ratio*cell.HOP
        ref = np.array([[0.,p],[p,0.]])
        draws['outer_'+name] = ins.acquire(rng, ins.fuel_law(ins.return_amplitude(ref,t), (0,2)), spec['shots'])
        for law in spec['laws']:
            for n in range(1,spec['depth']+2):
                H = source.reduced(spec['depth'], n, p, law)
                full = ins.fuel_law(ins.return_amplitude(H, t), (0,2))
                draws[f'{law}_{name}_{n}'] = ins.acquire(rng, full, spec['shots'])
    return maps, draws, flight.run(spec, rng)


def read(spec, draws):
    shots = spec['shots']
    for name, counts in draws.items():
        expected = (spec['times'], 4 if name == 'clock' else 2, 2)
        require(counts.shape == expected and counts.dtype.kind in 'iu' and
                np.all((counts >= 0) & (counts <= shots)), 'bad fuel counts: '+name)
    clock = ins.clock_read(draws['clock'], shots)
    return {name: ins.read_lines(clock, ins.quadrature(c,shots), spec)
            for name,c in draws.items() if name != 'clock'}


def compare_spectrum(row, expected, spec):
    expected = sorted(expected)
    require(np.shape(row['lines']) == (len(expected),) and
            np.max(abs(np.asarray(row['lines'])-expected), initial=0)<spec['maximum_frequency_error'],
            f"public lines {row['lines']} differ from the complete visible spectrum {expected}")
    require(np.shape(row['weights']) == (len(expected),) and np.all(np.isfinite(row['weights'])),
            'each visible line needs its own weight')
    require(abs(row['dc'] - (0 if expected else 1)) < .01, 'wrong stationary record weight')


def score(spec, maps, readings):
    depth = spec['depth']
    inv = source.inventory(depth)
    residuals = {}
    for name, p in maps['cell']['outer_rates'].items():
        ref = readings['outer_'+name]
        compare_spectrum(ref, [p] if p else [], spec)
        if p:
            require(abs(ref['weights'][0]-1) < .006, 'wrong external reference weight')
        outer = ref['lines'][0] if ref['lines'] else 0.
        for law in spec['laws']:
            errors = []
            for band in inv:
                n = band['length']
                rates = [1.] if law == 'no_children' else band['positive']
                row = readings[f'{law}_{name}_{n}']
                inner = source.chain(depth,n)[:n,n:]/cell.HOP
                _, modes = np.linalg.eigh(inner @ inner.T)
                strengths = [1.] if law == 'no_children' else abs(modes[0])**2
                if law == 'unsigned' and p:
                    expected = [abs(r+s*p) for r in rates for s in (-1,1)]
                    strengths = np.repeat(strengths,2)/2
                else:
                    expected = [math.hypot(r,p) for r in rates]
                compare_spectrum(row, expected, spec)
                require(np.allclose(row['weights'],np.asarray(strengths)[np.argsort(expected)],atol=.006,rtol=0),
                        'mode weights disagree with the source projection')
                rest = np.array(readings[f'{law}_uniform_{n}']['lines'])
                errors.extend(float(min(abs(w*w-outer**2-rest**2))) for w in row['lines'])
            residuals[law+'_'+name] = max(errors)
            if law == 'unsigned' and p:
                require(max(errors) > spec['minimum_unsigned_residual'], 'unsigned control also passed dispersion')
            else:
                require(max(errors) < spec['maximum_squared_dispersion_residual'], 'measured squared dispersion disagrees')
    light = sorted(readings[f'graded_uniform_{n}']['lines'][0] for n in range(1,depth+2))
    return dict(light=light, ratios=[x/light[0] for x in light], squared_dispersion_residuals=residuals)


def verify_saved(saved, current, path='receipt'):
    if type(current) is dict:
        require(type(saved) is dict and saved.keys()==current.keys(), 'wrong fields: '+path)
        for k in current:
            verify_saved(saved[k], current[k], path+'.'+k)
    elif type(current) is list:
        require(type(saved) is list and len(saved)==len(current), 'wrong list: '+path)
        for i,(a,b) in enumerate(zip(saved,current)):
            verify_saved(a,b,f'{path}[{i}]')
    elif type(current) is float:
        relative_only = (current == 0. or '.positive[' in path or
                         path.endswith(('.low','.upper_bound','.lower_bound','.uniform_amplitude_tail')))
        require(type(saved) is float and math.isclose(saved,current,abs_tol=0 if relative_only else 2e-8,
                                                     rel_tol=1e-8), 'wrong value: '+path)
    else:
        require(type(saved) is type(current) and saved==current, 'wrong value: '+path)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--record', type=Path, help='retain a fresh acquisition in an unused directory')
    args = cli.parse_args()
    require(load(HERE/'sources.json') == pins(), 'source freeze differs')
    spec = load(HERE/'protocol.json')
    maps, generated, transit = experiment(spec)
    where = args.record or HERE/'records'
    if args.record:
        where.mkdir(exist_ok=False)
        with (where/'counts.npz').open('xb') as f:
            np.savez_compressed(f, **generated)
        (where/'flight.json').write_text(json.dumps(transit,allow_nan=False)+'\n', encoding='utf8')
    with np.load(where/'counts.npz',allow_pickle=False) as z:
        require(sorted(z.files)==sorted(generated), 'missing or repeated record family')
        counts = {k:z[k] for k in generated}
        require(all(np.array_equal(v,generated[k]) for k,v in counts.items()), 'seeded acquisition differs')
    saved_flight = load(where/'flight.json')
    verify_saved(saved_flight, transit)
    flight.assess(spec, saved_flight)
    readings = read(spec, counts)
    result = dict(carriers=maps, reads=readings, result=score(spec,maps,readings),
                  inventory={str(d):source.inventory(d) for d in spec['depth_control']})
    if args.record:
        (where/'result.json').write_text(json.dumps(result,allow_nan=False)+'\n',encoding='utf8')
    else:
        verify_saved(load(where/'result.json'), result)
    print(json.dumps(result['result'],allow_nan=False,indent=2))


if __name__ == '__main__':
    main()
