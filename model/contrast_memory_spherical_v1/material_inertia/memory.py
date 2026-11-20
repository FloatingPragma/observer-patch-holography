"""An exactly uncomputed memory conversion on the official carrier and its mass-response control."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy import sparse

from tools.strict_json import load
from .. import carrier, law, reader, sectors
from . import instrument
from .__main__ import HERE, check_json, clean, frozen_files, name
from .federation import propagate
from .memory_algebra import local_conversion, transfer_algebra
from .source import Ring, require, rotate, spectral_measure


def pins():
    files = [HERE/n for n in ('sources.json','memory_algebra.py','memory.py','memory_protocol.json')]
    files += [HERE.parent/n for n in ('sectors.py','reader.py','common.py')]
    return {p.relative_to(HERE.parents[2]).as_posix():hashlib.sha256(
        p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files}


def joined(car,kappa,active):
    sector = sectors.PairSector(car)
    baseline = sparse.block_diag((sector.hamiltonian(kappa),sectors.doubled_hamiltonian(car,kappa)),format='csr')
    # Both native blocks use the same subtracted 12*kappa. A local conversion is
    # the square of one column per seam, with the true Fock/contrast normalization.
    i,j = np.asarray(car.seams).T
    e = np.arange(len(i))
    shape = (sector.P+car.n,len(e))
    single = sparse.csr_matrix((np.ones(len(e)),(sector.seam_states,e)),shape=shape)
    double = sparse.csr_matrix((np.full(2*len(e),-1/np.sqrt(2)),
                                (sector.P+np.r_[i,j],np.r_[e,e])),shape=shape)
    addition = (single+double)@(single+double).T if active else single@single.T + double@double.T
    return sector,(baseline+addition).tocsr()


def enroll(spec,kappa):
    raw = {}
    # Check the declared switched read actually maps both logical predicates to agreement.
    for flag in (0,1):
        word = reader.gate(flag,0)
        require((word>>6)&3 == 1 and (word>>8)&7 == 7*flag, 'native public record is not an agreement')
    for key in spec['carriers']:
        car = carrier.carrier(key)
        a,b = carrier.EMITTER_SEAM[key]
        for setting in ('active','blocked'):
            sector,H = joined(car,kappa,setting=='active')
            for prep in spec['preparations']:
                psi = np.zeros(H.shape[0],complex)
                psi[sector.P+a if prep=='double' else sector.index(a,b)] = 1
                last = 0
                rows = []
                for tick in spec['ticks']:
                    psi = propagate(H,psi,(tick-last)/spec['tick_divisor'])
                    last = tick
                    p = abs(psi)**2
                    # The supplementary doubled-location meter refines the native no-pair outcome.
                    probs = np.r_[p[:sector.P][~sector.adjacent].sum(),p[sector.seam_states],p[sector.P:]]
                    require(abs(probs.sum()-1)<2e-10,'incomplete joint read')
                    edges = instrument.counts(probs.cumsum(),spec['shots'])[:,1]
                    rows.append(np.diff(edges,prepend=np.int32(0)))
                raw[f'{key}_{setting}_{prep}'] = np.asarray(rows,dtype=np.int32)
    return raw


def ternary_ring(size, number):
    """Two bitplanes encode presence and a second unit. Include every multiset at fixed load."""
    states = []
    for nd in range(number//2+1):
        for doubles in combinations(range(size),nd):
            high = sum(1<<i for i in doubles)
            for singles in combinations([i for i in range(size) if not high&(1<<i)],number-2*nd):
                low = high + sum(1<<i for i in singles)
                states.append(low | (high<<size))
    mask = (1<<size)-1
    return Ring.from_states(size,number,states,
                            lambda word:rotate(word&mask,size) | (rotate(word>>size,size)<<size))


def ternary_law(ring,local,kappa):
    size = ring.size
    mask = (1<<size)-1
    rows,cols,entries = [],[],[]
    for col,word in enumerate(ring.states):
        high = word>>size
        diagonal = 4*kappa*(ring.number+2*high.bit_count())
        for i in range(size):
            a,b = 1<<i,1<<((i+1)%size)
            x = int(bool(word&a))+int(bool(high&a))
            y = int(bool(word&b))+int(bool(high&b))
            k = 3*x+y
            diagonal += local[k,k]
            for out in np.flatnonzero(local[:,k]):
                if out==k:
                    continue
                u,v = divmod(int(out),3)
                new_lo = (word&mask&~(a|b)) | (a if u else 0) | (b if v else 0)
                new_hi = (high&~(a|b)) | (a if u==2 else 0) | (b if v==2 else 0)
                rows.append(ring.index[new_lo | (new_hi<<size)])
                cols.append(col)
                entries.append(local[out,k])
        rows.append(col)
        cols.append(col)
        entries.append(diagonal)
    return sparse.csr_matrix((entries,(rows,cols)),shape=(len(ring.states),)*2)


def run(spec,base):
    kappa = base['kappa']
    V = local_conversion()
    result = {'memory_algebra':transfer_algebra(), 'conversion':{}, 'mass_response':{}}
    raw = enroll(spec,kappa)
    for key in spec['carriers']:
        car = carrier.carrier(key)
        for prep in spec['preparations']:
            for setting in ('active','blocked'):
                keyname = f'{key}_{setting}_{prep}'
                r = raw[keyname]
                require(r.shape==(len(spec['ticks']),1+len(car.seams)+car.n) and np.all(r>=0)
                        and np.all(r.sum(axis=1)==spec['shots']), 'bad memory read')
                exchanged = r[:,1:1+len(car.seams)] if prep=='double' else r[:,1+len(car.seams):]
                result['conversion'][keyname] = exchanged.sum(axis=1).tolist()
    raw['clock'] = instrument.clock(base)
    rate = instrument.clock_rate(raw['clock'],base)
    times = instrument.tick_menu(base)/base['tick_divisor']
    for size in spec['rings']:
        reads = {}
        for number in spec['occupations']:
            ring = ternary_ring(size,number)
            H = ternary_law(ring,law.seam_term(kappa)+V,kappa)
            require(np.asarray(abs(H).sum(axis=1)).max()/base['tick_divisor'] < np.pi-.2,
                    'initial clock can alias a qutrit line')
            for ell in (0,1):
                E,w = spectral_measure(H,ring.frame(ell))
                key = name('memory',size,number,ell)
                raw[key] = instrument.quadratures(np.exp(-1j*np.outer(times,E))@w,base['shots'])
                phase = instrument.read_phase(raw[key],base)
                weight = w[0]
                admitted = phase is not None and weight>=base['dominant_floor']
                if admitted:
                    read = phase.frequency(instrument.tick_menu(base)[-2]*rate)
                    require(abs(read[0]*rate*base['tick_divisor']-E[0])<=read[1]*rate*base['tick_divisor']+3e-10,
                            'qutrit line differs from source')
                reads[number,ell] = phase if admitted else None
            if number>1:
                result['mass_response'][f'{size}_{number}'] = instrument.band_ratios(
                    [reads[n,ell] for n,ell in ((1,0),(1,1),(number,0),(number,1))])
    result['conversion_gates'] = {k:(max(v)>=16 if '_active_' in k else max(v)==0)
                                  for k,v in result['conversion'].items()}
    result['rest_inertia'] = {key:instrument.rest_inertia_verdict(row,base)
                              for key,row in result['mass_response'].items()}
    return raw,clean(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    args = parser.parse_args()
    require(frozen_files()==load(HERE/'sources.json') and pins()==load(HERE/'memory_sources.json'),
            'changed memory construction')
    raw,result = run(load(HERE/'memory_protocol.json'),load(HERE/'protocol.json'))
    root = args.record or HERE/'records'/'memory'
    if args.record:
        root.mkdir(exist_ok=False,parents=True)
        np.savez_compressed(root/'counts.npz',**raw)
        (root/'result.json').write_text(json.dumps(result,indent=1,allow_nan=False)+'\n',encoding='utf8')
    with np.load(root/'counts.npz',allow_pickle=False) as saved:
        require(sorted(saved.files)==sorted(raw) and all(saved[k].dtype==np.int32 and
                np.array_equal(saved[k],v) for k,v in raw.items()),'changed memory measurements')
    check_json(root/'result.json',result)
    print(json.dumps(result,indent=1))


if __name__=='__main__':
    main()
