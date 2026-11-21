"""Observable limit of a fixed root detector under subdivision."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from tools.strict_json import load

from model.candidate_laws.scalar_channel.carrier import require
from model.experiments.cage.construction import HOP
from . import instruments as ins, source
from .__main__ import HERE, pins, verify_saved


def limiting_return(t):
    def density(theta):
        lam = 5-4*math.cos(theta)
        return 2/math.pi*math.sin(theta)**2/lam*math.cos(t*math.sqrt(lam))
    value, error = quad(density, 0, math.pi, epsabs=1e-12, limit=160)
    require(error < 2e-8, 'continuum quadrature unresolved')
    return .75+value


def moment_tail(n, t):
    order = 4*n
    x = 3*abs(t)
    if x == 0:
        return 0.
    if x >= order+1:
        return 2.
    return min(2., 2*math.exp(order*math.log(x)-math.lgamma(order+1))/(1-x/(order+1)))


def run(spec):
    rng = np.random.Generator(np.random.PCG64(spec['seed']))
    times = np.arange(spec['time_samples'])*spec['clock_step']
    clock = ins.acquire(rng, ins.fuel_law(ins.clock_amplitude(times),(0,2,1,3)),spec['shots'])
    ticks = ins.clock_read(clock,spec['shots'])
    window = ticks <= spec['finite_window_end']
    limit = np.array([limiting_return(t) for t in ticks])
    records = {'clock': clock}
    rows = []
    for n in spec['chain_lengths']:
        H = source.chain(n-1,n)/HOP
        generated = ins.return_amplitude(H,times)
        counts = ins.acquire(rng,ins.fuel_law(generated,(0,2)),spec['shots'])
        records[str(n)] = counts
        err = abs(ins.quadrature(counts,spec['shots'])-limit)
        low = float(source.light_rate(n))
        lo = 3/math.sqrt(4**(n+1)-3*n-4)
        hi = 2.**(1-n)
        require(lo*(1-2e-15) <= low <= hi*(1+2e-15), 'light rate outside inverse-norm bounds')
        rows.append(dict(length=n, low=float(low), lower_bound=lo, upper_bound=hi,
                         atoms_versus_limit_over_window=float(max(err[window])),
                         atoms_versus_limit_all=float(max(err)),
                         uniform_amplitude_tail=moment_tail(n,spec['finite_window_end'])))
    return records, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    require(pins()==load(HERE/'sources.json'), 'source freeze differs')
    additions = ('limit.py','limit_protocol.json','sources.json')
    expected = {name:hashlib.sha256((HERE/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                for name in additions}
    require(load(HERE/'limit_sources.json')==expected, 'limit freeze differs')
    spec = load(HERE/'limit_protocol.json')
    records, rows = run(spec)
    where = args.record or HERE/'records'/'limit'
    if args.record:
        where.mkdir(exist_ok=False)
        with (where/'counts.npz').open('xb') as f:
            np.savez_compressed(f,**records)
        (where/'result.json').write_text(json.dumps(rows,indent=1,allow_nan=False)+'\n',encoding='utf8')
    with np.load(where/'counts.npz',allow_pickle=False) as z:
        require(sorted(z.files)==sorted(records) and all(z[k].dtype.kind in 'iu' and np.array_equal(z[k],v)
                                                       for k,v in records.items()),
                'refinement records differ')
    verify_saved(load(where/'result.json'),rows)
    require(rows[-1]['atoms_versus_limit_over_window'] < spec['maximum_deep_record_difference']
            and rows[1]['atoms_versus_limit_all'] > spec['minimum_shallow_long_time_difference'],
            'fixed-root observer comparison failed')
    print(json.dumps(rows[-1],indent=2))


if __name__ == '__main__':
    main()
