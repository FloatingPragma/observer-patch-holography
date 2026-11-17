#!/usr/bin/env python3
"""Compressed composition of retained native source handoffs and native effects.

No fresh native calls are claimed here. The source records, certified Q and
executed detector tables are checked inputs; each actual address is propagated.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip, hashlib, json, math, time
import numpy as np
HERE=Path(__file__).resolve().parent
ROUTE=HERE.parent/'2026-11-16_codex_signed_pair_routing'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def propagate(name,cell,address):
    data=np.load(HERE/name/'maps.npz');q=data['certified_Q_proposal']
    d=data['dark_integer_basis'][cell:cell+12];P=d@np.linalg.inv(d.T@d)@d.T
    rho=np.zeros((12,12),complex);rho[address-cell,address-cell]=1
    assert np.linalg.norm(P@rho)==0
    escape=np.eye(12)-q.conj().T@q
    for _ in range(512):rho=q@rho@q.conj().T+np.trace(escape@rho)*np.eye(12)/12
    return rho,P
def output(flags,failure,meta):
    result={};fringe=[]
    for setting in ('X','Y'):
        law=flags@np.asarray(meta['complete_per_input_flag_output_centres'][setting])
        down=float(flags@np.asarray(meta['per_input_flag_processor_failure_centres'][setting]))
        p0,p1=(float(law[k]) for k in meta['success_outputs']);accepted=p0+p1
        result[setting]={'diagnostic_probability_centres':law.tolist(),
            'upstream_sector_processor_failure_centre':failure,
            'downstream_processor_failure_centre':down,
            'energy_success_probability_centre':accepted,
            'arithmetic_normalization_residual':float(1-failure-down-law.sum())}
        assert accepted>.99
        fringe.append((p0-p1)/accepted)
    result['complex_accepted_fringe']=fringe
    return result
def intervals(records,epsilon):
    radius=2*math.sqrt(2)*epsilon/(1-2*epsilon);phase=radius/(1-radius)
    z=complex(*records['1/10']['complex_accepted_fringe']);e=-10*math.atan2(z.imag,z.real)
    short=[e-10*phase,e+10*phase]
    z=complex(*records['4']['complex_accepted_fringe']);base=-math.atan2(z.imag,z.real)/4
    admissible=[(k,base+k*math.pi/2) for k in range(32)
                if short[0]-phase/4<=base+k*math.pi/2<=short[1]+phase/4]
    assert len(admissible)==1
    k,e=admissible[0]
    return {'complete_error_used':epsilon,'acceptance_floor':1-epsilon,
        'complex_fringe_radius':radius,'phase_error_upper':phase,
        'short_interval':short,'long_winding_from_short_interval':k,
        'long_interval':[e-phase/4,e+phase/4]}
def main():
    t=time.monotonic();reg=json.loads((HERE/'PREREG_ROUTED_JOIN.json').read_text())
    for name,digest in reg['files'].items():assert sha(HERE/name)==digest
    for name,digest in reg['route_files'].items():assert sha(ROUTE/name)==digest
    receipt=json.loads((HERE/'coupled_energy/RECEIPT.json').read_text())['scientific']
    records=json.loads(gzip.decompress((ROUTE/'RECORDS.json.gz').read_bytes()))
    sector=np.load(HERE/'poisson_energy/source_sector_maps.npz')
    result=[];densities={};detector=receipt['native_flag_controlled_words']
    for trial in records['trials']:
        if trial['law']!='native':continue
        assert trial['success']
        p,n=divmod(trial['pair_code'],256);source=trial['observed_source']
        assert p==0 and n==14 and source[p]==1 and source[n]==-1
        assert sum(abs(x) for x in source)==2 and trial['positive_frame']
        a,pa=propagate('first_cell_certified',0,p);b,pb=propagate('second_cell',12,n)
        af=np.array([np.trace(e@a).real for e in sector['first_effects']])
        bf=np.array([np.trace(e@b).real for e in sector['second_effects']])
        fa=float(np.trace(sector['first_processor_failure']@a).real)
        fb=float(np.trace(sector['second_processor_failure']@b).real)
        pair=np.array([af[i]*bf[j] for j in (0,1) for i in (0,1)])
        outs={timekey:{'single':output(af,fa,table['single']),
                      'pair':output(pair,fa+sum(af)*fb,table['pair'])}
              for timekey,table in detector.items()}
        for timekey in outs:
            for name in ('single','pair'):
                old=receipt['complete_source_coupled_records'][timekey][name]
                assert max(abs(np.asarray(old['complex_accepted_fringe'])-
                               outs[timekey][name]['complex_accepted_fringe']))<1e-12
        result.append({'trial_index':trial['index'],'input_pair_code':trial['input_pair_code'],
            'actual_routed_pair_code':trial['pair_code'],'actual_source':source,
            'route_record_sha256':hashlib.sha256(json.dumps(trial,sort_keys=True).encode()).hexdigest(),
            'native_frame':trial['positive_frame'],'complete_energy_output_centres':outs,
            'selected_population_centres':[float(np.trace(pa@a).real),float(np.trace(pb@b).real)]})
        densities['trial_%02d_first'%trial['index']]=a
        densities['trial_%02d_second'%trial['index']]=b
    assert len(result)==12
    # Eight possible capped moves, including the actual oriented frame.
    route_error=1-(1-F(1,2**32))**8
    joint=F(receipt['complete_joint_preparation_and_energy_read_error_upper'])+route_error
    assert joint<F(1,1000)
    single=F(receipt['complete_source_preparation_error_upper'][0])+F(receipt['complete_coupled_sector_reader_error_upper'])+route_error
    bounds={name:intervals({k:result[0]['complete_energy_output_centres'][k][name] for k in detector},float(err))
            for name,err in [('single',single),('pair',joint)]}
    lo,hi=bounds['single']['long_interval'];pl,ph=bounds['pair']['long_interval']
    ratio=[pl/hi,ph/lo]
    out=HERE/'routed_join';out.mkdir(exist_ok=False)
    np.savez_compressed(out/'maps.npz',**densities)
    report={'candidate':'A-old240-selected-compact-source-v1',
        'handoffs':result,'route_manifest_sha256':sha(ROUTE/'MANIFEST.sha256'),
        'coupled_receipt_sha256':sha(HERE/'coupled_energy/RECEIPT.json'),
        'route_failure_symbol':'terminal_route_timeout','route_timeout_upper':str(route_error),
        'complete_join_error_upper':str(joint),'energy_intervals':bounds,
        'pair_to_single_interval':ratio,'pair_keV_at511_anchor_interval':[511*x for x in ratio],
        'preparation_law':'Each observed coordinate propagated through the retained certified source-Q pump.',
        'centre_semantics':receipt['probability_centre_semantics'],
        'literal_new_native_calls':0,'composition':'12 compressed certified-operator/detector-table compositions',
        'retained_route_alternatives':['HOLD timeouts','forced native timeout'],
        'scope':'Conditioned on the sealed actual selected pair inputs; source-selection failure law remains in the loader packet.',
        'full_row_solved':False,'other_group_review':False}
    write(out/'RECEIPT.json',{'scientific':report,'prereg_sha256':sha(HERE/'PREREG_ROUTED_JOIN.json'),
        'maps_sha256':sha(out/'maps.npz'),'elapsed_seconds':time.monotonic()-t})
    print(json.dumps({'routed_join_pass':True,'actual_handoffs':len(result),
        'complete_error':float(joint),'ratio_interval':ratio,'seconds':time.monotonic()-t}))
if __name__=='__main__':main()
