"""Energy released by two opposite charges prepared in the measured lowest band."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from tools.strict_json import load
from model.executable_model import carrier
from model.candidate_laws.scalar_channel.carrier import require
from . import apparatus, source
from .__main__ import HERE, check_counts, check_result, pins


def freeze():
    return {n:hashlib.sha256((HERE/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for n in ('sources.json','rest_protocol.json','rest.py')}


def layout(car,kind):
    site = source.rest_cells(car)
    if kind=='disjoint':
        other = next(c for c in range(car.cells) if c not in site and
                     any(c in (i//12,j//12) and i//12 not in site and j//12 not in site for i,j in car.external))
        far = source.rest_cells(car,home=other,excluded=site)
        require(set(site).isdisjoint(far),'control packets overlap')
        return site,far
    return site,site


def packets(car, kind):
    site,far = layout(car,kind)
    v,_ = source.rest_witness(car)
    w = source.rest_witness(car,home=far[0],excluded=site)[0] if kind=='disjoint' else v
    x = np.outer(v,w)
    return np.r_[x.ravel(),np.zeros(car.n**2)],list(dict.fromkeys(site+far))


def menu(key):
    return ['overlap','blocked'] + (['disjoint'] if key!='cell' else [])


def collect(spec):
    rng = np.random.Generator(np.random.PCG64(spec['seed']))
    shots = spec['shots']
    times = np.arange(int(max(spec['waits'])/spec['clock_step'])+1)*spec['clock_step']
    rec = {'clock':apparatus.reference(rng,times,shots)}
    for key in spec['carriers']:
        car = carrier.carrier(key)
        labels = source.ports_to_regions(car)
        for kind in menu(key):
            law = 'dephased_contact' if kind=='blocked' else 'exchange'
            H = source.generator(car,law)
            psi,_ = packets(car,kind)
            loc,energy = [],[]
            for wave in source.evolve(H,psi,spec['waits']):
                X,Y = wave.reshape(2,car.n,car.n)
                loc.append(apparatus.counts(rng,source.position(X,Y,labels),shots))
                rows,_ = source.energy_instruments(car,X,Y,law)
                energy.append(apparatus.counts(rng,rows,shots))
            rec[f'{key}_{kind}_position'] = np.array(loc)
            rec[f'{key}_{kind}_energy'] = np.array(energy)
    return rec


def halfwidth(spec,car):
    # Bernstein for independent bounded shots across local settings. Positivity gives
    # sum of variances <= b E. Include the worst cumulative-edge digitization error.
    x = np.log(4*spec['upper_number_of_energy_estimates']/spec['simultaneous_error'])
    b,N = spec['max_bin_energy'],spec['shots']
    # Each one-particle bond operator sum has norm at most twice the maximum
    # degree. In the two-particle sector the contact blocks add at most two.
    require(spec['energy_cap'] >= 4*int(car.degree.max())+2 and b >= 4,
            'confidence budget understates the source or detector energy bound')
    num_settings = 3*len(car.seams)+car.n
    bias = 2*b*num_settings/2**30
    return float(np.sqrt(2*b*(spec['energy_cap']+bias)*x/N)+2*b*x/(3*N)+bias)


def analyze(spec,rec):
    shots = spec['shots']
    ticks = apparatus.read_clock(rec['clock'],shots)
    require(len(ticks) == round(max(spec['waits'])/spec['clock_step'])+1, 'wrong rest-pair clock menu')
    result = dict(clock=ticks[np.rint(np.array(spec['waits'])/spec['clock_step']).astype(int)].tolist(),histories={})
    settings = 0
    for key in spec['carriers']:
        car = carrier.carrier(key)
        regions = len(set(source.ports_to_regions(car)))
        bound = halfwidth(spec,car)
        for kind in menu(key):
            settings += len(spec['waits'])*(3*len(car.seams)+car.n)
            home,other = layout(car,kind)
            site = list(dict.fromkeys(home+other))
            name = f'{key}_{kind}'
            raw = rec[name+'_position']
            require(raw.shape==(len(spec['waits']),regions**2+regions*(regions+1)//2) and raw.dtype.kind in 'iu'
                    and np.all(raw>=0) and np.all(raw.sum(axis=1)==shots),'bad rest-pair position counts')
            law = 'dephased_contact' if kind=='blocked' else 'exchange'
            energies = apparatus.energies(car,rec[name+'_energy'],shots,law,site)
            require(len(energies['total_energy'])==len(raw),'bad rest-pair energy menu')
            neutral = raw[:,regions**2:].sum(axis=1)/shots
            baseline = energies['contact_energy'][0]
            # Total energy was measured before and after a source-on wait, without a switch.
            excess = np.array(energies['neutral_energy'])-baseline
            result['histories'][name] = dict(**energies,preparation_cells=site,neutral_probability=neutral.tolist(),
                energy_halfwidth=bound,field_above_initial_contact=excess.tolist(),
                simultaneous_lower_bound=(excess-2*bound).tolist())
    h = result['histories']
    # Clock quadratures, final positions and local energies use separate preparations.
    settings_total = 2*len(ticks) + len(h)*len(spec['waits']) + settings
    trials = sum(int(a.sum(dtype=np.int64)) for a in rec.values())
    require(trials == shots*settings_total, 'trial total differs from the scheduled preparations')
    result['resources'] = dict(clock_settings=2*len(ticks), position_settings=len(h)*len(spec['waits']),
                               energy_settings=settings, total_settings=settings_total, total_trials=trials)
    # Two stochastic binomial splits suffice for a three-outcome setting. Rounding
    # its conditional uniform quantile perturbs each integer CDF by <= 2^-53.
    grid_error = 2*settings*(shots+1)/2**53
    require(grid_error < spec['simultaneous_error']/2 and
            5*len(h)*len(spec['waits']) <= spec['upper_number_of_energy_estimates'],
            'energy confidence budget exhausted')
    result['quantile_total_variation_budget'] = grid_error
    result['gates'] = {str(k)+'_release':max(h[f'{k}_overlap']['simultaneous_lower_bound'])>0
                       for k in spec['carriers']}
    result['gates']['blocked'] = all(not any(h[f'{k}_blocked']['neutral_probability']) for k in spec['carriers'])
    result['gates']['disjoint'] = not any(h['0_disjoint']['neutral_probability'])
    result['gates']['balance'] = all(max(abs(np.array(r['total_energy'])-r['total_energy'][0]))
                                     <spec['absolute_balance_gate'] for r in h.values())
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    args = parser.parse_args()
    require(pins()==load(HERE/'sources.json') and freeze()==load(HERE/'rest_sources.json'),'source differs')
    spec = load(HERE/'rest_protocol.json')
    rec = collect(spec)
    root = args.record or HERE/'records'/'rest'
    saved = check_counts(root,rec,record=bool(args.record))
    result = analyze(spec,saved)
    check_result(root,result,bool(args.record))
    print(json.dumps(result['gates'],indent=2))


if __name__=='__main__':
    main()
