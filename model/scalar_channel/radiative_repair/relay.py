"""Local source switch, number-resolving receiver and subsequent conversion with work accounting."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from tools.strict_json import load
from ...contrast_memory_spherical_v1 import carrier
from ..carrier import require
from . import apparatus, source
from .__main__ import HERE, check_counts, check_result, pins


def extra_pins():
    return {name:hashlib.sha256((HERE/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for name in ('sources.json','relay_protocol.json','relay.py','subspaces.py')}


def receiver_cell(car):
    ends = carrier.EMITTER_SEAM[0]
    u,v = next((u,v) for u,v in car.external if u in ends or v in ends)
    return v//12 if u in ends else u//12


def sectors(wave, labels, receiver):
    X,Y = wave.reshape(2,len(labels),len(labels))
    inside = (labels==receiver).astype(int)
    number = inside[:,None]+inside[None,:]
    branches = []
    for h in (0,1,2):
        selected = np.stack([X if h==0 else np.zeros_like(X),Y*(number==h)])
        branches.append(selected.ravel())
    return branches


def check_switch(before, after, labels, receiver):
    """Certify the complete switch operator and its diagonal blocks in the meter basis."""
    n = len(labels)
    ports = np.flatnonzero(labels == receiver)
    pairs = ports*(n+1)
    # The first block holds signed pairs. The second holds ordered neutral coordinates.
    rows = np.r_[pairs, n*n+pairs]
    delta = (after-before).tocoo()
    order = np.argsort(delta.row)
    require(len(ports) > 0 and np.array_equal(delta.row[order], rows)
            and np.array_equal(delta.col[order], np.r_[n*n+pairs, pairs])
            and np.all(delta.data == -1), 'switch changes more than the receiver contact')
    # Its only two entries per port connect a charged coincidence (meter 0)
    # and a neutral coincidence (meter 2). Thus every P_h delta P_h is zero.


def work(car, states, mask, rng, shots):
    rows = None
    for v in states:
        probability = np.vdot(v,v).real
        if probability == 0:
            continue
        X,Y = (v/np.sqrt(probability)).reshape(2,car.n,car.n)
        probe,scale = source.energy_instruments(car,X,Y,'exchange',mask)
        if rows is None:
            rows = np.zeros_like(probe)
        rows += probability*probe
    return apparatus.counts(rng,rows,shots)


def collect(spec):
    car = carrier.federation(spec['level'])
    labels = source.ports_to_regions(car)
    receiver = receiver_cell(car)
    emitting = labels==0
    H = source.generator(car,'exchange',emitting)
    rng = np.random.Generator(np.random.PCG64(spec['seed']))
    shots = spec['shots']
    dt = spec['herald_time']
    maximum = dt+max(spec['delays'])
    rec = {'clock':apparatus.reference(rng,np.arange(8*maximum+1)/8,shots)}
    prepared = source.initial(car,0,spec['emitter'])
    incoming = next(source.evolve(H,prepared,[dt]))
    parts = sectors(incoming,labels,receiver)
    require(np.linalg.norm(sum(parts)-incoming) < 1e-13,'incomplete local number instrument')
    rec['initial_work'] = work(car,[prepared],emitting,rng,shots)
    rec['incoming_work'] = work(car,[incoming],emitting,rng,shots)
    rec['measured_work'] = work(car,parts,emitting,rng,shots)
    for setting in spec['controls']:
        mask = emitting | ((labels==receiver) & (setting=='receiver_on'))
        later = source.generator(car,'exchange',mask)
        if setting == 'receiver_on':
            check_switch(H, later, labels, receiver)
        rec[setting+'_switch_work'] = work(car,parts,mask,rng,shots)
        traces = [source.evolve(later,v,spec['delays']) for v in parts]
        histories = []
        for states in zip(*traces,strict=True):
            joint = np.array([source.position(*v.reshape(2,car.n,car.n),labels) for v in states])
            histories.append(apparatus.counts(rng,joint.ravel(),shots).reshape(joint.shape))
        rec[setting+'_joint'] = np.array(histories)
        rec[setting+'_final_work'] = work(car,states,mask,rng,shots)
    return rec


def relay_intervals(on,off,shots,alpha=.001):
    """One-sided Clopper-Pearson bounds, simultaneous over both arms and all named delays."""
    require(len(on)==len(off)>0,'empty or unequal delay menus')
    tail = alpha/(4*len(on))
    result = []
    for a,b in zip(on,off,strict=True):
        require(0<=a<=shots and 0<=b<=shots,'invalid joint count')
        lo = beta.ppf(tail,a,shots-a+1) if a else 0.
        hi = beta.ppf(1-tail,b+1,shots-b) if b<shots else 1.
        result.append(float(lo-hi-2/2**30))
    return result


def analyze(spec, rec):
    car = carrier.federation(spec['level'])
    receiver = receiver_cell(car)
    regions = car.cells
    shots = spec['shots']
    ticks = apparatus.read_clock(rec['clock'],shots)
    require(len(ticks) == 8*(spec['herald_time']+max(spec['delays']))+1, 'wrong relay clock menu')
    first = spec['delays'].index(0)
    left, right = np.triu_indices(regions)
    neutral_number = (left==receiver).astype(int) + (right==receiver)
    immediate_number = np.r_[np.zeros(regions**2, dtype=int), neutral_number]
    readings = {}
    energy = {}
    for k,data in rec.items():
        if k.endswith('_work'):
            active = np.arange(car.n)//12==0
            if k.startswith('receiver_on'):
                active |= np.arange(car.n)//12==receiver
            energy[k] = apparatus.energies(car,data,shots,'exchange',active=active)['total_energy']
    for k in spec['controls']:
        a = rec[k+'_joint']
        size = regions**2 + regions*(regions+1)//2
        require(a.shape==(len(spec['delays']),3,size) and a.dtype.kind in 'iu' and np.all(a>=0)
                and np.all(a.sum(axis=(1,2))==shots),'invalid joint receiver history')
        require(not np.any(a[first][np.arange(3)[:,None] != immediate_number]),
                'simultaneous receiver records disagree')
        both = a[:,2,:]
        arrivals = both.sum(axis=1)
        local_pair = both[:,receiver*(regions+1)]
        charged = both[:,:regions**2].sum(axis=1)
        readings[k] = dict(double_neutral_counts=arrivals.tolist(),later_signed_pair_counts=charged.tolist(),
                           later_receiver_pair_counts=local_pair.tolist(),
                           incoming_fraction=(arrivals/shots).tolist(),unconditional_local_fraction=(local_pair/shots).tolist())
    on,off = (readings[k] for k in spec['controls'])
    joint_bins = 3*(regions**2+regions*(regions+1)//2)
    grid_error = 2*len(spec['delays'])*(joint_bins-1)*(shots+1)/2**53
    require(grid_error < .001/2,'relay quantile budget exhausted')
    checks = dict(arrived=min(on['double_neutral_counts'])>=spec['min_double_arrivals'],
                  reabsorbed=on['later_receiver_pair_counts'][0]==0 and
                  max(on['later_receiver_pair_counts'])>=spec['min_returned_pairs'],
                  initially_neutral=on['later_signed_pair_counts'][0]==off['later_signed_pair_counts'][0]==0,
                  propagation_budget=abs(energy['incoming_work']-energy['initial_work'])<spec['max_energy_error'])
    # With only the sender converter left on, a returning neutral packet can convert
    # there and its charged daughters can travel back. Keep this competing history.
    checks['receiver_contrast'] = max((np.array(on['later_receiver_pair_counts'])-
                                      off['later_receiver_pair_counts'])/shots)>spec['min_returned_pairs']/shots
    for k in spec['controls']:
        checks[k+'_energy'] = abs(energy[k+'_final_work']-energy[k+'_switch_work'])<spec['max_energy_error']
    lower = relay_intervals(on['later_receiver_pair_counts'],off['later_receiver_pair_counts'],shots)
    checks['resolved_receiver_effect'] = max(lower) > 0
    return dict(receiver_cell=receiver,clock_at_meter=float(ticks[8*spec['herald_time']]),
                clock_at_read=ticks[8*(spec['herald_time']+np.array(spec['delays']))].tolist(),
                readings=readings,energy=energy,quantile_total_variation_budget=grid_error,
                local_effect_lower=lower,meter_work=energy['measured_work']-energy['incoming_work'],
                switch_work={k:energy[k+'_switch_work']-energy['measured_work'] for k in spec['controls']},
                gates={k:bool(v) for k,v in checks.items()})


def main():
    from .subspaces import certificates
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path)
    args = parser.parse_args()
    require(pins()==load(HERE/'sources.json') and extra_pins()==load(HERE/'relay_sources.json'),'source differs')
    spec = load(HERE/'relay_protocol.json')
    rec = collect(spec)
    root = args.record or HERE/'records'/'relay'
    saved = check_counts(root,rec,record=bool(args.record))
    result = analyze(spec,saved)
    result['certificates'] = certificates()
    check_result(root,result,bool(args.record))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
